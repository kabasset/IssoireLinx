/// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXDATA_BOX_H
#define _LINXDATA_BOX_H

#include "Linx/Data/mixins/Region.h"

#include <boost/operators.hpp>

namespace Linx {

/// @cond
namespace Internal {
template <Index N>
class BorderedBox; // for friendness // FIXME rm?
}
/// @endcond

/**
 * @ingroup regions
 * @brief An ND bounding box, defined by its front and back positions, both inclusive.
 * 
 * Like `Position`, this class stores no pixel values, but coordinates.
 */
template <Index N = 2>
class Box : public Dimensional<N>, public boost::additive<Box<N>, Index> {
  friend class Internal::BorderedBox<N>; // FIXME move to FilterMixin, sole user

public:

  /**
   * @brief A position iterator.
   * 
   * The scanning order maximizes data locality for row-major ordered data, like rasters.
   * 
   * That is, the increment operator is such that the corresponding offset in a raster is always increasing.
   * In particular, when screening a whole raster, pixels are visited in the storage order.
   */
  class Iterator;

  /// @{
  /// @group_construction

  Box() : Box(Position<N>::one(), Position<N>::zero()) {} // FIXME {zero(), -one()}

  /**
   * @brief Constructor.
   */
  Box(Position<N> front, Position<N> back) : m_front(LINX_MOVE(front)), m_back(LINX_MOVE(back)) {}

  /**
   * @brief Create a box from a shape (front is zero).
   */
  static Box<N> from_shape(Position<N> shape)
  {
    const auto dim = shape.size();
    return from_shape(Position<N>::zero(dim), LINX_MOVE(shape));
  }

  /**
   * @brief Create a box from a front position and shape.
   */
  static Box<N> from_shape(Position<N> front, Position<N> shape)
  {
    shape += front;
    --shape; // FIXME merge the two loops
    return {LINX_MOVE(front), LINX_MOVE(shape)};
  }

  /**
   * @brief Create a box from a radius and center position.
   */
  static Box<N> from_center(Index radius = 1, const Position<N>& center = Position<N>::zero())
  {
    return {center - radius, center + radius};
  }

  /**
   * @brief Create a conventionally unlimited box.
   * 
   * Front and back bounds along each axis are respectively 0 and inf.
   */
  static Box<N> whole()
  {
    return {Position<N>::zero(), Position<N>::inf()};
  }

  /// @group_properties

  /**
   * @brief Get the number of axes.
   */
  Index dimension() const
  {
    return m_front.size();
  }

  /**
   * @brief Get the front position.
   */
  inline const Position<N>& front() const
  {
    return m_front;
  }

  /**
   * @brief Get the back position.
   */
  inline const Position<N>& back() const
  {
    return m_back;
  }

  /**
   * @brief `Position::one()`, for compatibility with `Grid`.
   */
  static Position<N> step()
  {
    return Position<N>::one();
  }

  /**
   * @brief Compute the box shape.
   */
  Position<N> shape() const
  {
    return m_back - m_front + 1; // FIXME optimize
  }

  /**
   * @brief Compute the box size, i.e. number of positions.
   */
  Index size() const
  {
    return shape_size(shape());
  }

  /**
   * @brief Get the box length along given axis.
   */
  constexpr Index length(Index i) const
  {
    return m_back[i] - m_front[i] + 1;
  }

  /// @group_elements

  /**
   * @brief Get the absolute position given a position in the box referential.
   */
  inline Position<N> operator[](const Position<N>& p) const
  {
    return p + m_front;
  }

  /**
   * @brief Get an iterator to the beginning.
   */
  Iterator begin() const
  {
    return Iterator::begin(*this);
  }

  /**
   * @brief Get an iterator to the end.
   */
  Iterator end() const
  {
    return Iterator::end(*this);
  }

  /// @group_operations

  /**
   * @brief Check whether two boxes are equal.
   */
  bool operator==(const Box<N>& other) const
  {
    return m_front == other.m_front && m_back == other.m_back;
  }

  /**
   * @brief Check whether two boxes are different.
   */
  bool operator!=(const Box<N>& other) const
  {
    return m_front != other.m_front || m_back != other.m_back;
  }

  /**
   * @brief Check whether a position lies inside the box.
   */
  bool contains(const Position<N>& position) const
  {
    for (Index i = 0; i < dimension(); ++i) {
      if (position[i] < m_front[i] || position[i] > m_back[i]) {
        return false;
      }
    }
    return true;
  }

  /// @group_modifiers

  /**
   * @brief Flatten the box along a given axis.
   * 
   * The back of the box is set to the same coordinate as the front along the axis.
   */
  Box<N>& project(Index axis = 0)
  {
    m_back[axis] = m_front[axis];
    return *this;
  }

  /**
   * @brief Check whether the box is contained within another box.
   */
  bool operator<=(const Box<N>& rhs) const
  {
    for (Index i = 0; i < dimension(); ++i) {
      if (m_front[i] < rhs.front()[i]) {
        return false;
      }
      if (m_back[i] > rhs.back()[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Check whether the box is stricly contained within another box.
   */
  bool operator<(const Box<N>& rhs) const
  {
    for (Index i = 0; i < dimension(); ++i) {
      if (m_front[i] <= rhs.front()[i]) {
        return false;
      }
      if (m_back[i] >= rhs.back()[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Check whether the box contains another box.
   */
  bool operator>=(const Box<N>& rhs) const
  {
    for (Index i = 0; i < dimension(); ++i) {
      if (m_front[i] > rhs.front()[i]) {
        return false;
      }
      if (m_back[i] < rhs.back()[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Check whether the box stricly contains another box.
   */
  bool operator>(const Box<N>& rhs) const
  {
    for (Index i = 0; i < dimension(); ++i) {
      if (m_front[i] >= rhs.front()[i]) {
        return false;
      }
      if (m_back[i] <= rhs.back()[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Shrink the box inside another box (i.e. get the intersection of both).
   */
  Box<N>& operator&=(const Box<N>& rhs)
  {
    for (Index i = 0; i < dimension(); ++i) {
      m_front[i] = std::max(m_front[i], rhs.front()[i]);
      m_back[i] = std::min(m_back[i], rhs.back()[i]);
    }
    return *this;
  }

  /**
   * @brief Minimally grow the box to include another box (i.e. get the minimum box which contains both).
   */
  Box<N>& operator|=(const Box<N>& rhs)
  {
    for (Index i = 0; i < dimension(); ++i) {
      m_front[i] = std::min(m_front[i], rhs.front()[i]);
      m_back[i] = std::max(m_back[i], rhs.back()[i]);
    }
    return *this;
  }

  /**
   * @brief Grow the box by a given margin.
   */
  template <Index M>
  Box<N>& operator+=(const Box<M>& margin)
  {
    // FIXME allow N=-1
    m_front += extend<N>(margin.front());
    m_back += extend<N>(margin.back());
    return *this;
  }

  /**
   * @brief Shrink the box by a given margin.
   */
  template <Index M>
  Box<N>& operator-=(const Box<M>& margin)
  {
    // FIXME allow N=-1
    m_front -= extend<N>(margin.front());
    m_back -= extend<N>(margin.back());
    return *this;
  }

  /**
   * @brief Translate the box by a given vector.
   */
  template <Index M>
  Box<N>& operator+=(const Position<M>& vector)
  {
    // FIXME allow N=-1
    m_front += extend<N>(vector);
    m_back += extend<N>(vector);
    return *this;
  }

  /**
   * @brief Translate the box by the opposite of a given vector.
   */
  template <Index M>
  Box<N>& operator-=(const Position<M>& vector)
  {
    // FIXME allow N=-1
    m_front -= extend<N>(vector);
    m_back -= extend<N>(vector);
    return *this;
  }

  /**
    * @brief Add a scalar to each coordinate.
    */
  Box<N>& operator+=(Index scalar)
  {
    m_front += scalar;
    m_back += scalar;
    return *this;
  }

  /**
   * @brief Subtract a scalar to each coordinate.
   */
  Box<N>& operator-=(Index scalar)
  {
    m_front -= scalar;
    m_back -= scalar;
    return *this;
  }

  /**
   * @brief Add 1 to each coordinate.
   */
  Box<N>& operator++()
  {
    return *this += 1;
  }

  /**
   * @brief Subtract 1 to each coordinate.
   */
  Box<N>& operator--()
  {
    return *this -= 1;
  }

  /**
   * @brief Copy.
   */
  Box<N> operator+()
  {
    return *this;
  }

  /**
   * @brief Invert the sign of each coordinate.
   */
  Box<N> operator-()
  {
    return {-m_front, -m_back};
  }

  /// @}

private:

  /**
   * @brief The front position in the box.
   */
  Position<N> m_front;

  /**
   * @brief The back position in the box.
   */
  Position<N> m_back;
};

/**
 * @relatesalso Box
 * @brief Identity, for compatibility with `Region`.
 */
template <Index N>
inline const Box<N>& box(const Box<N>& region)
{
  return region;
}

/**
 * @relatesalso Box
 * @brief Create a box of higher dimension.
 */
template <Index M, Index N>
Box<M> extend(const Box<N>& in, const Position<M>& padding = Position<M>::zero())
{
  return {extend<M>(in.front(), padding), extend<M>(in.back(), padding)}; // FIXME move padding once
}

/**
 * @relatesalso Box
 * @brief Flatten the box along a given axis.
 */
template <Index N>
Box<N> project(Box<N> in, Index axis = 0)
{
  in.project(axis);
  return in;
}

/**
 * @relatesalso Box
 * @brief Erase an axis.
 * @tparam The index of the axis
 */
template <Index I, Index N>
Box<Dimensional<N>::OneLessDimension> erase(const Box<N>& in)
{
  return {erase<I>(in.front()), erase<I>(in.back())};
}

/**
 * @relatesalso Box
 * @brief Insert an axis.
 * @tparam I The index of the new axis
 * @param front The front bound along axis `I`
 * @param back The back bound along axis `I`
 */
template <Index I, Index N>
Box<Dimensional<N>::OneMoreDimension> insert(const Box<N>& in, Index front, Index back)
{
  return {insert<I>(in.front(), front), insert<I>(in.back(), back)};
}

/**
 * @relatesalso Box
 * @brief Clamp a position inside a box.
 */
template <typename T, Index N = 2>
Vector<T, N> clamp(Vector<T, N> position, const Box<N>& box)
{
  position.apply(
      [](auto p, auto f, auto b) {
        return std::clamp(p, f, b);
      },
      box.front(),
      box.back());
  return position;
}

/**
 * @relatesalso Box
 * @brief Clamp a position inside a shape.
 */
template <typename T, Index N = 2>
Vector<T, N> clamp(Vector<T, N> position, const Position<N>& shape)
{
  position.apply(
      [](auto p, auto s) {
        return std::clamp(p, Index(0), s - 1);
      },
      shape);
  return position;
}

/**
 * @relatesalso Box
 */
template <Index N, Index M>
Box<N> operator+(Box<N> in, const Box<M>& margin)
{
  in += margin;
  return in;
}

/**
 * @relatesalso Box
 */
template <Index N, Index M>
Box<N> operator-(Box<N> in, const Box<M>& margin)
{
  in -= margin;
  return in;
}

/**
 * @relatesalso Box
 */
template <Index N, Index M>
Box<N> operator+(Box<N> in, const Position<M>& vector)
{
  in += vector;
  return in;
}

/**
 * @relatesalso Box
 */
template <Index N, Index M>
Box<N> operator-(Box<N> in, const Position<M>& vector)
{
  in -= vector;
  return in;
}

/**
 * @relatesalso Box
 * @brief Get the intersection of two boxes.
 */
template <Index N>
inline Box<N> operator&(Box<N> region, const Box<N>& bounds)
{
  region &= bounds;
  return region;
}

} // namespace Linx

#include "Linx/Data/impl/BoxIterator.h"

#endif
