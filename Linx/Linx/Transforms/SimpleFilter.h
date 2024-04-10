// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXTRANSFORMS_SIMPLEFILTER_H
#define _LINXTRANSFORMS_SIMPLEFILTER_H

#include "Linx/Data/Raster.h"
#include "Linx/Transforms/Extrapolation.h"
#include "Linx/Transforms/mixins/Filter.h"

namespace Linx {

/// @cond
namespace Internal {

/**
 * @brief Test whether a type is a range, i.e. has `begin()` and `end()` methods.
 */
template <typename T, typename = void>
struct KernelShiftsWindow : std::false_type {};

/// @cond
template <typename T>
struct KernelShiftsWindow<T, std::void_t<typename T::ShiftsWindow>> : std::true_type {};
/// @endcond

} // namespace Internal
/// @endcond

/**
 * @ingroup filtering
 * @brief A structuring element for morphological operations.
 * @tparam TOp The morphological operator
 * @tparam T The value type
 * @tparam TWindow The type of window, e.g. `Box` or `Mask`
 */
template <typename TKernel>
class SimpleFilter : public FilterMixin<typename TKernel::Value, typename TKernel::Window, SimpleFilter<TKernel>> {
  friend class FilterMixin<typename TKernel::Value, typename TKernel::Window, SimpleFilter<TKernel>>;

public:

  /**
   * @brief The value type.
   */
  using Value = typename TKernel::Value;

  /**
   * @brief The dimension parameter.
   */
  static constexpr Index Dimension = TKernel::Dimension;

  /**
   * @brief Forwarding constructor.
   */
  template <typename... TArgs>
  explicit SimpleFilter(TArgs&&... args) : m_kernel(LINX_FORWARD(args)...)
  {}

  /**
   * @brief Get the kernel.
   */
  const TKernel& kernel() const
  {
    return m_kernel;
  }

  /**
   * @brief Get the kernel.
   */
  TKernel& kernel()
  {
    return m_kernel;
  }

protected:

  /**
   * @brief Get the window.
   */
  decltype(auto) window_impl() const
  {
    return kernel().window();
  }

  /**
   * @brief Get the bounding box of the window, extended to a given dimension.
   */
  template <Index M>
  decltype(auto) window_box() const
  {
    return extend<M>(box(kernel().window()));
  }

  /**
   * @brief Filter and crop an input raster.
   * 
   * Border regions which would require extrapolation are cropped-out,
   * such that the output domain is `in.domain() - filter.window()`.
   */
  template <typename T, Index N, typename THolder, typename TOut>
  void transform_impl(const Raster<T, N, THolder>& in, TOut& out) const
  {
    const auto region = in.domain() - window_box<N>();
    transform_monolith(in(region), out);
  }

  /**
   * @brief Filter an extrapolated raster.
   * 
   * The output raster has the same shape as the input raster.
   */
  template <typename TRaster, typename TMethod, typename TOut>
  void transform_impl(const Extrapolation<TRaster, TMethod>& in, TOut& out) const
  {
    const auto& raw = dont_extrapolate(in);
    const auto bbox = Internal::BorderedBox<TRaster::Dimension>(raw.domain(), window_box<TRaster::Dimension>());
    bbox.apply_inner_border(
        [&](const auto& ib) {
          const auto insub = raw(ib);
          if (insub.size() > 0) {
            auto outsub = out(insub.domain());
            transform_monolith(insub, outsub);
          }
        },
        [&](const auto& ib) {
          const auto insub = in(ib);
          if (insub.size() > 0) {
            auto outsub = out(insub.domain());
            transform_monolith_extrapolator(insub, outsub);
          }
        });
  }

  /**
   * @brief Filter and decimate a grid-based patch.
   * 
   * Decimation is especially useful for downsampling.
   * 
   * If input values from outside the raster domain are required,
   * then `in` must be an extrapolator.
   * On the contrary, if the box of `in` is small enough so that no extrapolated values are required,
   * then `in` can be a raw patch.
   * 
   * The output raster has the same shape as the grid.
   */
  template <typename T, typename TParent, typename TRegion, typename TOut>
  void transform_impl(const Patch<T, TParent, TRegion>& in, TOut& out) const
  {
    const auto& raw = dont_extrapolate(in);
    const auto& front = in.domain().front();
    const auto& step = in.domain().step();

    const auto grid_to_box = [&](const auto& g) {
      auto f = g.front() - front;
      for (std::size_t i = 0; i < f.size(); ++i) { // FIXME simplify with Linx
        f[i] /= step[i];
      }
      return Box<TParent::Dimension>::from_shape(f, g.shape());
    };

    const auto& window = window_box<TParent::Dimension>();
    const auto& domain = rasterize(in).domain();
    const auto bbox = Internal::BorderedBox<TParent::Dimension>(domain, window);
    // FIXME accept non-Box window, and of lower dim
    bbox.apply_inner_border(
        [&](const auto& ib) {
          const auto insub = raw(ib);
          if (insub.size() > 0) { // FIXME needed?
            auto outsub = out(grid_to_box(insub.domain()));
            transform_monolith(insub, outsub);
          }
        },
        [&](const auto& ib) {
          const auto insub = in(ib);
          if (insub.size() > 0) {
            auto outsub = out(grid_to_box(insub.domain()));
            transform_monolith(insub, outsub);
          }
        });
  }

private:

  /**
   * @brief Filter a monolithic patch (no region splitting).
   */
  template <typename TIn, typename TOut>
  void transform_monolith_extrapolator(const TIn& in, TOut& out) const
  {
    const auto extrapolated = in.parent().copy(Linx::box(in.domain()) + window_box<TIn::Dimension>());
    const auto box = extrapolated.domain() - window_box<TIn::Dimension>();
    // FIXME region - window().front()?
    transform_monolith(extrapolated(box), out);
  }

  /**
   * @brief Filter an input extrapolator or patch into an output raster or patch.
   * @param in A extrapolator or patch (or both)
   * @param out A raster or patch with compatible domain
   * 
   * If the filtering bounding box requires input values from outside the input domain,
   * then `in` must be an extrapolator.
   * If the bounding box of `in` is small enough so that no extrapolated values are required,
   * then `in` can be a raw patch.
   */
  template <typename TIn, typename TOut>
  void transform_monolith(const TIn& in, TOut& out) const
  {
    // FIXME accept any region
    auto patch = in.parent()(extend<TIn::Dimension>(window_impl()));
    auto out_it = out.begin();
    for (const auto& p : in.domain()) { // FIXME loop over out for a simpler sentinel?
      if constexpr (Internal::KernelShiftsWindow<TKernel>::value) { // FIXME ugly
        *out_it = m_kernel(in, patch, p);
      } else {
        patch >>= p;
        *out_it = m_kernel(patch);
        patch <<= p;
      }
      ++out_it;
    }
  }

private:

  /**
   * @brief The operation.
   */
  TKernel m_kernel;
};

} // namespace Linx

#endif
