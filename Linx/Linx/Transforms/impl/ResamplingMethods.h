// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXTRANSFORMS_IMPL_INTERPOLATIONMETHODS_H
#define _LINXTRANSFORMS_IMPL_INTERPOLATIONMETHODS_H

#include "Linx/Data/Raster.h"

namespace Linx {

/**
 * @ingroup resampling
 * @brief Constant, a.k.a. Dirichlet boundary conditions.
 */
template <typename T>
class Constant {
public:

  /**
   * @brief Constructor.
   */
  Constant(T v = T {}) : m_value(v) {}

  /**
   * @brief Return the value if out of bounds.
   */
  template <typename TRaster> // FIXME template <bool MaybeInDomain = true, typename TRaster>
  inline const T& at(TRaster& raster, const Position<TRaster::Dimension>& position) const
  {
    return raster.contains(position) ? raster[position] : m_value;
  }

  /**
   * @brief Get the extrapolation value.
   */
  operator T() const
  {
    return m_value;
  }

private:

  /**
   * @brief The extrapolation value.
   */
  T m_value; // FIXME private + operator T
};

/**
 * @ingroup resampling
 * @brief Nearest-neighbor interpolation or extrapolation, a.k.a. zero-flux Neumann boundary conditions.
 */
struct Nearest {
  /**
   * @brief Return the value at the nearest in-bounds position.
   */
  template <typename TRaster>
  inline const typename TRaster::value_type& at(TRaster& raster, Position<TRaster::Dimension> position) const
  {
    return raster[clamp(LINX_MOVE(position), raster.shape())];
  }

  /**
   * @brief Return the value at the nearest integer position.
   */
  template <typename T, typename TRaster>
  inline T at(TRaster& raster, const Vector<double, TRaster::Dimension>& position) const
  {
    Position<TRaster::Dimension> integral(position.size());
    integral.generate(
        [](auto e) {
          return e + .5; // FIXME std::round?
        },
        position);
    return raster[integral];
  }
};

/**
 * @ingroup resampling
 * @brief Periodic, a.k.a. symmetric or wrap-around, boundary conditions.
 */
struct Periodic {
  /**
   * @brief Return the value at the modulo position.
   */
  template <typename TRaster>
  inline const typename TRaster::value_type& at(TRaster& raster, Position<TRaster::Dimension> position) const
  {
    position.apply(
        [](auto p, auto s) {
          auto q = p % s;
          return q < 0 ? q + s : q; // Positive modulo
        },
        raster.shape());
    return raster[position];
  }
};

/**
 * @ingroup resampling
 * @brief Linear interpolation.
 */
struct Linear {
  /**
   * @brief Return the interpolated value at given index.
   */
  template <typename T, typename TRaster, typename... TIndices>
  inline T at(const TRaster& raster, const Vector<double, 1>& position, TIndices... indices) const
  {
    const auto f = floor<Index>(position.front());
    const auto d = position.front() - f;
    const T p = raster[{f, indices...}];
    const T n = raster[{f + 1, indices...}];
    return d * (n - p) + p;
  }

  /**
   * @brief Return the interpolated value at given position.
   */
  template <typename T, Index N, typename TRaster, typename... TIndices>
  inline std::enable_if_t<N != 1, T> // FIXME use if constexpr
  at(const TRaster& raster, const Vector<double, N>& position, TIndices... indices) const
  {
    const auto f = floor<Index>(position.back());
    const auto d = position.back() - f;
    const auto pos = slice<N - 1>(position);
    const auto p = at<T>(raster, pos, f, indices...);
    const auto n = at<T>(raster, pos, f + 1, indices...);
    return d * (n - p) + p;
  }
};

/**
 * @ingroup resampling
 * @brief Cubic interpolation.
 */
struct Cubic {
  /**
   * @brief Return the interpolated value at given index.
   */
  template <typename T, typename TRaster, typename... TIndices>
  inline T at(const TRaster& raster, const Vector<double, 1>& position, TIndices... indices) const
  {
    const auto f = floor<Index>(position.front());
    const auto d = position.front() - f;
    const T pp = raster[{f - 1, indices...}];
    const T p = raster[{f, indices...}];
    const T n = raster[{f + 1, indices...}];
    const T nn = raster[{f + 2, indices...}];
    return p + 0.5 * (d * (-pp + n) + d * d * (2 * pp - 5 * p + 4 * n - nn) + d * d * d * (-pp + 3 * p - 3 * n + nn));
  }

  /**
   * @brief Return the interpolated value at given position.
   */
  template <typename T, Index N, typename TRaster, typename... TIndices>
  inline std::enable_if_t<N != 1, T> // FIXME use if constexpr
  at(const TRaster& raster, const Vector<double, N>& position, TIndices... indices) const
  {
    const auto f = floor<Index>(position.back());
    const auto d = position.back() - f;
    const auto pos = slice<N - 1>(position);
    const auto pp = at<T>(raster, pos, f - 1, indices...);
    const auto p = at<T>(raster, pos, f, indices...);
    const auto n = at<T>(raster, pos, f + 1, indices...);
    const auto nn = at<T>(raster, pos, f + 2, indices...);
    return p + 0.5 * (d * (-pp + n) + d * d * (2 * pp - 5 * p + 4 * n - nn) + d * d * d * (-pp + 3 * p - 3 * n + nn));
  }
};

} // namespace Linx

#endif
