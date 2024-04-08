// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXDATA_VECTOR_H
#define _LINXDATA_VECTOR_H

#include "Linx/Base/Dimension.h"
#include "Linx/Base/TypeUtils.h"
#include "Linx/Base/mixins/DataContainer.h"
#include "Linx/Base/mixins/Math.h" // abspow

#include <numeric> // accumulate, multiplies
#include <type_traits> // conditional

namespace Linx {

/**
 * @ingroup data_classes
 * @brief N-dimensional vector, mainly intended for pixel position or image shape, i.e. set of coordinates.
 * @tparam N A non-negative dimension (0 is allowed), or -1 for variable dimension.
 * 
 * The values are stored in a `std::array<T, N>` in general (`N &ge; 0`),
 * or `std::vector<T>` for variable dimension (`N = -1`).
 *
 * Memory and services are optimized when dimension is fixed at compile-time (`N &ge; 0`).
 * 
 * @tspecialization{Position}
 */
template <typename T, Index N = 2>
class Vector :
    public Dimensional<N>,
    public DataContainer<T, StdHolder<Coordinates<T, N>>, VectorArithmetic, Vector<T, N>> {
public:

  /**
   * @brief The value type.
   */
  using Value = T;

  /**
   * @brief The container type.
   */
  using Container = DataContainer<T, StdHolder<Coordinates<T, N>>, VectorArithmetic, Vector<T, N>>;

  LINX_VIRTUAL_DTOR(Vector)
  LINX_DEFAULT_COPYABLE(Vector)
  LINX_DEFAULT_MOVABLE(Vector)

  /**
   * @brief Default constructor.
   * @warning
   * The values are unspecified.
   * To create vector 0, use `zero()` instead.
   */
  Vector() : Container() {}

  /**
   * @brief Create a vector of given dimension.
   */
  template <typename TInt, typename std::enable_if_t<std::is_integral<TInt>::value>* = nullptr>
  explicit Vector(TInt dim) : Container(std::size_t(dim))
  {}

  /**
   * @brief Create a vector from a brace-enclosed list of indices.
   */
  Vector(std::initializer_list<T> indices) : Container(indices.begin(), indices.end()) {}

  /**
   * @brief Create a vector from a range.
   */
  template <typename TRange, typename std::enable_if_t<IsRange<TRange>::value>* = nullptr>
  explicit Vector(TRange&& range) : Container(range.begin(), range.end())
  {}

  /**
   * @brief Create a vector full of `Limits::zero()'s.
   */
  static Vector<T, N> zero(Index size = std::abs(N))
  {
    Vector<T, N> res(size);
    return res.fill(Limits<T>::zero());
  }

  /**
   * @brief Create a vector full of `Limits::one()'s.
   */
  static Vector<T, N> one(Index size = std::abs(N))
  {
    Vector<T, N> res(size);
    return res.fill(Limits<T>::one());
  }

  /**
   * @brief Create a vector full of `Limits::inf()'s.
   */
  static Vector<T, N> inf(Index size = std::abs(N))
  {
    Vector<T, N> res(size);
    return res.fill(Limits<T>::inf());
  }

  /**
   * @brief Check whether the vector is zero.
   */
  bool is_zero() const
  {
    return this->contains_only(Limits<T>::zero());
  }

  /**
   * @brief Check whether the vector is one.
   */
  bool is_one() const
  {
    return this->contains_only(Limits<T>::one());
  }

  /**
   * @brief Check whether the vector is minus one.
   */
  bool is_inf() const
  {
    return this->contains_only(Limits<T>::inf());
  }
};

/**
   * @brief Create a vector of lower dimension.
   * @tparam M The new dimension; cannot be -1
   * 
   * The values up to dimension `M` are copied.
   */
template <Index M, typename T, Index N>
Vector<T, M> slice(const Vector<T, N>& vector)
{
  Vector<T, M> res(M);
  std::copy_n(vector.data(), M, res.data());
  return res;
}

/**
 * @brief Create a vector of higher dimension.
 * @tparam M The new dimension; cannot be -1
 * 
 * The values up to dimension `N` are copied.
 * Those between dimensions `N` and `M` are taken from the given padding vector.
 */
template <Index M, typename T, Index N>
Vector<T, M> extend(const Vector<T, N>& vector, Vector<T, M> padding = Vector<T, M>::zero())
{
  std::copy(vector.begin(), vector.end(), padding.begin());
  return padding;
}

/**
 * @relatesalso Vector
 * @brief Erase an element from a given vector.
 * 
 * The size of the resulting vector is that of the input vector minus one.
 */
template <Index I, typename T, Index N>
Vector<T, Dimensional<N>::OneLessDimension> erase(const Vector<T, N>& in)
{
  static_assert(I >= 0);
  constexpr auto M = Dimensional<N>::OneLessDimension;
  Vector<T, M> out(in.size() - 1);
  for (std::size_t i = 0; i < I; ++i) {
    out[i] = in[i];
  }
  for (std::size_t i = I; i < out.size(); ++i) {
    out[i] = in[i + 1];
  }
  return out;
}

/**
 * @relatesalso Vector
 * @brief Insert an element into a given vector.
 * 
 * The size of the resulting vector is that of the input vector plus one.
 */
template <Index I, typename T, Index N>
Vector<T, Dimensional<N>::OneMoreDimension> insert(const Vector<T, N>& in, T&& value = T {})
{
  static_assert(I >= 0);
  constexpr auto M = Dimensional<N>::OneMoreDimension;
  Vector<T, M> out(in.size() + 1);
  for (std::size_t i = 0; i < I; ++i) {
    out[i] = in[i];
  }
  out[I] = std::forward<T>(value);
  for (std::size_t i = I + 1; i < out.size(); ++i) {
    out[i] = in[i - 1];
  }
  return out;
}

/**
 * @ingroup data_classes
 * @brief Pixel position or image shape, i.e. set of indices.
 * @tparam N A non-negative dimension (0 is allowed), or -1 for variable dimension.
 * 
 * Anonymous brace-initialization is permitted, e.g.:
 * \code
 * Raster<float> raster({1920, 1080});
 * // Is equivalent to
 * Raster<float> raster(Position<2>({1920, 1080}));
 * \endcode
 * 
 * Classical positions are instantiated with named constructors, e.g.:
 * \code
 * auto bottom_left = Position<2>::zero();
 * auto top_right = Position<2>::inf();
 * \endcode
 * 
 * @see Box
 */
template <Index N = 2>
using Position = Vector<Index, N>;

/// @cond
namespace Internal {

template <typename>
struct IsPositionImpl : std::false_type {};

template <Index N>
struct IsPositionImpl<Position<N>> : std::true_type {};

} // namespace Internal

/**
 * @relatesalso Position
 * @brief Check whether a class is a position.
 */
template <typename T>
static constexpr bool is_position()
{
  return Internal::IsPositionImpl<std::decay_t<T>>::value;
}

/**
 * @relatesalso Position
 * @brief Get the stride along a given axis.
 */
template <Index N>
Index shape_stride(const Position<N>& shape, Index axis)
{
  return std::accumulate(shape.begin(), shape.begin() + axis, 1L, [](auto s, auto l) {
    return l > 0 ? s * l : 0; // FIXME here or for shape_size only?
  });
}

/**
 * @relatesalso Position
 * @brief Get the stride along a given axis.
 */
template <Index Axis, Index N>
Index shape_stride(const Position<N>& shape)
{
  return shape_stride(shape, Axis);
}

/**
 * @relatesalso Position
 * @brief Compute the number of pixels in a given shape.
 */
template <Index N = 2>
Index shape_size(const Position<N>& shape)
{
  const auto size = shape.size();
  return size == 0 ? 0 : shape_stride(shape, size);
}

} // namespace Linx

#endif
