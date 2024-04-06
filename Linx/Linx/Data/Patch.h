// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXDATA_PATCH_H
#define _LINXDATA_PATCH_H

#include "Linx/Base/TypeUtils.h"
#include "Linx/Base/mixins/Arithmetic.h"
#include "Linx/Base/mixins/Math.h"
#include "Linx/Data/Box.h"
#include "Linx/Data/impl/PatchIndexing.h"
#include "Linx/Data/mixins/Region.h"

#include <algorithm> // accumulate
#include <functional> // multiplies

namespace Linx {

/**
 * @ingroup data_classes
 * @ingroup regions
 * @brief A view of a raster region.
 * @tparam T The value type
 * @tparam TParent The parent raster or extrapolator type
 * @tparam TRegion The region type
 * @tparam IsContiguous A memory contiguity flag
 * 
 * As opposed to a raster, values of a patch are generally not contiguous in memory:
 * they are piece-wise contiguous when the region is a `Box`, and sometimes not even piece-wise contiguous.
 * 
 * Whatever the region type, patches are iterable,
 * and the iterator type depends on the parent and region types in order to maximize performance.
 * Assuming the region itself is cheap to translate, patches are cheap to translate and iterate,
 * which makes them ideal to represent sliding windows, even of arbitrary shapes (e.g. when the region is a `Mask`).
 * 
 * In-place pixel-wise operations of rasters (like arithmetic operators and math functions)
 * are applicable to patches of mutable parent.
 * 
 * @see pixelwise
 */
template <typename T, typename TParent, typename TRegion, bool IsContiguous = false>
class Patch :
    public Dimensional<TParent::Dimension>,
    public ArithmeticMixin<EuclidArithmetic, T, Patch<T, TParent, TRegion, IsContiguous>>,
    public MathFunctionsMixin<T, Patch<T, TParent, TRegion, IsContiguous>>,
    public RangeMixin<T, Patch<T, TParent, TRegion, IsContiguous>> {
public:

  /**
   * @brief The value type.
   */
  using Value = T;
  // std::conditional_t<std::is_const<TParent>::value, const typename TParent::Value, typename TParent::Value>;

  /**
   * @brief The parent type (a raster or extrapolator).
   */
  using Parent = TParent;

  /**
   * @brief The region type.
   */
  using Region = std::decay_t<TRegion>;

  /**
   * @brief The indexing strategy.
   */
  using Indexing = typename PatchTraits<std::decay_t<Parent>, Region, IsContiguous>::template Indexing<Parent, Region>;

  /**
   * @brief The iterator type.
   * @tparam The value type, can be `T` or `const T`
   */
  template <typename U>
  using Iterator = typename Indexing::template Iterator<U>;

  /// @{
  /// @group_construction

  /**
   * @brief Default constructor.
   * 
   * This is a dummy constructor defined for compatibility with `std::vector`.
   */
  Patch() : m_parent(nullptr), m_region(), m_indexing() {}

  /**
   * @brief Constructor.
   */
  Patch(Parent& parent, Region region) : m_parent(&parent), m_region(std::move(region)), m_indexing(*m_parent, m_region)
  {}
  // std::move(region) not applicable to const Region& // TODO decay?

  /// @group_properties

  /**
   * @brief Get the patch bounding box.
   */
  Box<Patch::Dimension> box() const
  {
    return box(m_region);
  }

  /**
   * @brief Get the number of pixels in the patch.
   */
  std::size_t size() const
  {
    return m_region.size();
  }

  /**
   * @brief Get the region.
   */
  const TRegion& domain() const
  {
    return m_region;
  }

  /**
   * @brief Check whether the patch is known to be contiguous in memory.
   * @return 1 if contiguity is guaranteed, 0 if non-contiguity is guaranteed, -1 if unknown
   * 
   * Although there may be false negatives (undetected contiguous patches), there are no false positives.
   * When the method returns `true`, iterating over pointers is therefore safe.
   */
  int contiguity() const
  {
    return IsContiguous ? 1 : -1; // FIXME return m_indexing.is_contiguous() for box, -1 otherwise
  }

  /**
   * @brief Access the parent raster or extrapolator.
   */
  const Parent& parent() const
  {
    return *m_parent;
  }

  /**
   * @brief Access the parent raster or extrapolator.
   */
  Parent& parent()
  {
    return *m_parent;
  }

  /// @group_elements

  /**
   * @brief Access the element at given position in the region referential.
   * 
   * Only patches based on regions with suitable `operator[]()` are supported.
   * For example, box-based patches support `Position` as an argument,
   * while line-based patches support `Index` as an argument.
   * 
   * This method is provided for convenience but is slower than iterators.
   * Generally, it should be avoided in hot paths.
   */
  template <typename TPos>
  const Value& operator[](const TPos& pos) const
  {
    return (*m_parent)[m_region[pos]];
  }

  /**
   * @copybrief operator[]()const
   */
  template <typename TPos>
  Value& operator[](const TPos& pos)
  {
    return (*m_parent)[m_region[pos]];
  }

  /**
   * @brief Constant iterator to the front pixel.
   */
  Iterator<const Value> begin() const
  {
    return m_indexing.template begin<const Value>(*m_parent, m_region);
  }

  /**
   * @brief Iterator to the front pixel.
   */
  Iterator<Value> begin()
  {
    return m_indexing.template begin<Value>(*m_parent, m_region);
  }

  /**
   * @brief Constant end iterator.
   */
  Iterator<const Value> end() const
  {
    return m_indexing.template end<const Value>(*m_parent, m_region);
  }

  /**
   * @brief End iterator.
   */
  Iterator<Value> end()
  {
    return m_indexing.template end<Value>(*m_parent, m_region);
  }

  /**
   * @brief Get a pointer to the underlying array.
   * @warning This method is only valid for contiguous patches.
   */
  const Value* data() const
  {
    return &(*begin());
  }

  /**
   * @copydoc data()const
   */
  Value* data()
  {
    return &(*begin());
  }

  /// @group_modifiers

  /**
   * @brief Translate the patch by a given vector.
   * 
   * This is cheap iff the region translation itself is cheap.
   */
  Patch& operator>>=(const Position<Patch::Dimension>& vector)
  {
    m_region += vector;
    return *this;
  }

  /**
   * @brief Translate the patch by the opposite of a given vector.
   */
  Patch& operator<<=(const Position<Patch::Dimension>& vector)
  {
    m_region -= vector;
    return *this;
  }

  /// @group_operations

  /**
   * @brief Check whether two patches are equal.
   */
  bool operator==(const Patch& other) const
  {
    return m_parent == other.m_parent && m_region == other.m_region;
  }

  /**
   * @brief Check whether two patches are different.
   */
  bool operator!=(const Patch& other) const
  {
    return not(*this == other);
  }

  /**
   * @brief Create a row at given relative position.
   * @param i The row index following the data ordering
   * @param position The row position
   * @warning This only works for box-based patches.
   */
  auto row(Index i) const
  {
    return row(Position<1> {i}); // FIXME accept other dimensions
  }

  /**
   * @copydoc row(Index)const
   */
  auto row(Index i)
  {
    return row(Position<1> {i}); // FIXME accept other dimensions
  }

  /**
   * @copydoc row(Index)const
   */
  auto row(const Position<Patch::OneLessDimension> position) const
  {
    Position<Patch::Dimension> front = insert<0>(position) + m_region.front();
    Position<Patch::Dimension> back = front;
    back[0] += m_region.length(0) - 1;
    return operator()({LINX_MOVE(front), LINX_MOVE(back)});
  }

  /**
   * @copydoc row(Index)const
   */
  auto row(const Position<Patch::OneLessDimension>& position)
  {
    Position<Patch::Dimension> front = insert<0>(position) + m_region.front();
    Position<Patch::Dimension> back = front;
    back[0] += m_region.length(0) - 1;
    return operator()({LINX_MOVE(front), LINX_MOVE(back)});
  }

  /**
   * @brief Create a cropped patch.
   */
  Patch<const T, const TParent, TRegion> operator()(const Box<Patch::Dimension>& box) const // FIXME contiguity
  {
    return Patch<const T, const TParent, TRegion>(*m_parent, m_region & box);
  }

  /**
   * @brief Create a cropped patch.
   */
  Patch operator()(const Box<Patch::Dimension>& box) // FIXME contiguity
  {
    return Patch(*m_parent, m_region & box);
  }

  /**
   * @brief Create a cropped patch.
   */
  Patch<const T, const TParent, Box<Patch::Dimension>, true>
  operator()(const Position<Patch::Dimension>& position) const
  {
    return (*this)(Box<Patch::Dimension>(position, position)); // FIXME ArrSequence<Position>?
  }

  /**
   * @brief Create a cropped patch.
   */
  Patch<T, TParent, Box<Patch::Dimension>, true> operator()(const Position<Patch::Dimension>& position)
  {
    return (*this)(Box<Patch::Dimension>(position, position)); // FIXME ArrSequence<Position>?
  }

  /// @}

private:

  /**
   * @brief The parent raster.
   */
  Parent* m_parent;

  /**
   * @brief The region.
   */
  Region m_region;

  /**
   * @brief The translation-invariant indexing helper.
   */
  Indexing m_indexing;
};

/**
 * @brief Create a patch translated by a given vector.
 */
template <typename T, typename TParent, typename TRegion, bool IsContiguous>
Patch<T, TParent, TRegion, IsContiguous>
operator>>(Patch<T, TParent, TRegion, IsContiguous> lhs, const Position<TParent::Dimension>& rhs)
{
  lhs >>= rhs;
  return lhs;
}

/**
 * @brief Create a patch translated by the opposite of a given vector.
 */
template <typename T, typename TParent, typename TRegion, bool IsContiguous>
Patch<T, TParent, TRegion, IsContiguous>
operator<<(Patch<T, TParent, TRegion, IsContiguous> lhs, const Position<TParent::Dimension>& rhs)
{
  lhs <<= rhs;
  return lhs;
}

/// @cond
namespace Internal {

template <typename T>
struct IsPatchImpl : std::false_type {};

template <typename T, typename TParent, typename TRegion, bool IsContiguous>
struct IsPatchImpl<Patch<T, TParent, TRegion, IsContiguous>> : std::true_type {};

template <typename, typename = void>
struct IsPatchableImpl : std::false_type {};

template <typename T>
struct IsPatchableImpl<T, std::void_t<decltype(std::declval<const T>().domain())>> :
    Internal::IsRegionImpl<std::decay_t<decltype(std::declval<const T>().domain())>> {}; // FIXME only for boxes?

} // namespace Internal
/// @endcond

/**
 * @brief Check whether a class is a patch.
 */
template <typename T>
constexpr bool is_patch()
{
  return Internal::IsPatchImpl<std::decay_t<T>>::value;
}

/**
 * @relatesalso Patch
 * @brief Check whether a class can be used as a patch parent.
 */
template <typename T>
constexpr bool is_patchable() // FIXME rm?
{
  return Internal::IsPatchableImpl<std::decay_t<T>>::value;
}

/**
 * @relatesalso Patch
 * @brief Get the parent raster of a patch.
 * 
 * As opposed to `Patch::parent()`, if the parent is an extrapolator,
 * then the underlying decorated raster is effectively returned.
 */
template <typename T, typename TParent, typename TRegion, bool IsContiguous>
const auto& rasterize(const Patch<T, TParent, TRegion, IsContiguous>& in)
{
  return rasterize(in.parent());
}

} // namespace Linx

#endif
