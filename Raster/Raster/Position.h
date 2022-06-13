// Copyright (C) 2022, Antoine Basset
// This file is part of Cnes.Raster <github.com/kabasset/Raster>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _RASTER_POSITION_H
#define _RASTER_POSITION_H

#include "Raster/Arithmetic.h"
#include "Raster/DataContainer.h"
#include "Raster/DataUtils.h"

#include <numeric> // accumulate, multiplies
#include <type_traits> // conditional

namespace Cnes {

/**
 * @relates Position
 * @brief The signed integer type which represents indices in a raster.
 */
using Index = long;

/**
 * @brief A container of coordinates.
 */
template <typename T, Index N = 2>
using Coordinates = typename std::conditional<(N == -1), std::vector<T>, std::array<T, (std::size_t)N>>::type;

/**
 * @relates Position
 * @brief The index container type.
 */
template <Index N = 2>
using Indices = Coordinates<Index, N>;

/**
 * @ingroup data_classes
 * @brief N-dimensional vector, mainly intended for pixel position or image shape, i.e. set of coordinates.
 * @tparam N A non-negative dimension (0 is allowed), or -1 for variable dimension.
 * @details
 * The values are stored in a `std::array<T, N>` in general (`N &ge; 0`),
 * or `std::vector<T>` for variable dimension (`N = -1`).
 *
 * Memory and services are optimized when dimension is fixed at compile-time (`N &ge; 0`).
 * 
 * @tspecialization{Position}
 */
template <typename T, Index N = 2>
class Vector :
    public DataContainer<T, StdHolder<Coordinates<T, N>>, Vector<T, N>>,
    public ArithmeticMixin<VectorArithmetic, T, Vector<T, N>> {
public:
  /**
   * @brief The dimension template parameter.
   */
  static constexpr Index Dim = N;

  CNES_VIRTUAL_DTOR(Vector)
  CNES_DEFAULT_COPYABLE(Vector)
  CNES_DEFAULT_MOVABLE(Vector)

  /**
   * @brief Default constructor.
   * @warning
   * The indices are unspecified.
   * To create position 0, use `zero()` instead.
   */
  Vector() : DataContainer<T, StdHolder<Coordinates<T, N>>, Vector<T, N>>() {}

  /**
   * @brief Create a position of given dimension.
   */
  explicit Vector(T dim) : DataContainer<T, StdHolder<Coordinates<T, N>>, Vector<T, N>>(dim) {}

  /**
   * @brief Create a position from a brace-enclosed list of indices.
   */
  Vector(std::initializer_list<T> indices) : DataContainer<T, StdHolder<Coordinates<T, N>>, Vector<T, N>>(indices) {}

  /**
   * @brief Create a position from an iterable.
   */
  template <typename TIterable, typename std::enable_if_t<isIterable<TIterable>::value>* = nullptr>
  explicit Vector(TIterable&& iterable) : DataContainer<T, StdHolder<Coordinates<T, N>>, Vector<T, N>>(iterable) {}

  /**
   * @brief Create position 0.
   */
  static Vector<T, N> zero() {
    Vector<T, N> res(std::abs(Dim));
    return res.fill(Limits<T>::zero());
  }

  /**
   * @brief Create a position full of 1's.
   */
  static Vector<T, N> one() {
    Vector<T, N> res(std::abs(Dim));
    return res.fill(Limits<T>::one());
  }

  /**
   * @brief Create max position (full of -1's).
   */
  static Vector<T, N> max() {
    Vector<T, N> res(std::abs(N));
    return res.fill(-Limits<T>::one());
  }

  /**
   * @brief Check whether the position is zero.
   */
  bool isZero() const {
    for (auto i : *this) {
      if (i != Limits<T>::zero()) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Check whether the position is one.
   */
  bool isOne() const {
    for (auto i : *this) {
      if (i != Limits<T>::one()) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Check whether the position is max.
   */
  bool isMax() const {
    for (auto i : *this) {
      if (i != -Limits<T>::one()) {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Create a position of lower dimension.
   * @tparam M The new dimension; cannot be -1
   * @details
   * The indices up to dimension `M` are copied.
   */
  template <Index M>
  Vector<T, M> slice() const {
    Vector<T, M> res(M);
    std::copy_n(this->data(), M, res.data());
    return res;
  }

  /**
   * @brief Create a position of higher dimension.
   * @tparam M The new dimension; cannot be -1
   * @details
   * The indices up to dimension `N` are copied.
   * Those between dimensions `N` and `M` are taken from the given position.
   */
  template <Index M>
  Vector<T, M> extend(const Vector<T, M>& padding = Vector<T, M>::zero()) const {
    auto res = padding;
    std::copy_n(this->data(), this->size(), res.data());
    return res;
  }
};

/**
 * @ingroup data_classes
 * @brief Pixel position or image shape, i.e. set of indices.
 * @tparam N A non-negative dimension (0 is allowed), or -1 for variable dimension.
 * @details
 * Anonymous brace-initialization is permitted, e.g.:
 * \code
 * Raster<float> raster({1920, 1080});
 * // Is equivalent to
 * Raster<float> raster(Position<2>({1920, 1080}));
 * \endcode
 * 
 * Classical positions are instantiated with named constructors, e.g.:
 * \code
 * auto bottomLeft = Position<2>::zero();
 * auto topRight = Position<2>::max();
 * \endcode
 * 
 * @see Region
 */
template <Index N = 2>
using Position = Vector<Index, N>;

/**
 * @relates Position
 * @brief Get the stride along a given axis.
 */
template <Index N>
Index shapeStride(const Position<N>& shape, Index axis) {
  return std::accumulate(shape.begin(), shape.begin() + axis, 1L, std::multiplies<Index> {});
}

/**
 * @relates Position
 * @brief Get the stride along a given axis.
 */
template <Index Axis, Index N>
Index shapeStride(const Position<N>& shape) {
  return shapeStride(shape, Axis);
}

/**
 * @relates Position
 * @brief Compute the number of pixels in a given shape.
 */
template <Index N = 2>
Index shapeSize(const Position<N>& shape) {
  const auto size = shape.size();
  if (size == 0) {
    return 0;
  }
  return shapeStride(shape, size);
}

} // namespace Cnes

#endif // _RASTER_POSITION_H
