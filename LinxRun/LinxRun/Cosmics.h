// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXRUN_COSMICS_H
#define _LINXRUN_COSMICS_H

#include "Linx/Data/Raster.h"
#include "Linx/Transforms/Filters.h"

namespace Linx {
namespace Cosmics {

/**
 * @brief Pearson correlation coefficient kernel.
*/
template <typename T, typename TWindow>
class PearsonCorrelation : public KernelMixin<T, TWindow, PearsonCorrelation<T, TWindow>> {
public:

  using KernelMixin<T, TWindow, PearsonCorrelation>::KernelMixin;

  void init_impl() // FIXME private
  {
    const auto mean = std::accumulate(this->m_values.begin(), this->m_values.end(), T()) / this->m_values.size();
    std::transform(this->m_values.begin(), this->m_values.end(), this->m_values.begin(), [=](auto& e) {
      return e - mean;
    });
    m_sum2 = std::inner_product(this->m_values.begin(), this->m_values.end(), this->m_values.begin(), T());
  }

  template <typename TIn>
  T operator()(const TIn& neighbors) const
  {
    const auto mean = std::accumulate(neighbors.begin(), neighbors.end(), T()) / neighbors.size();
    auto centered = neighbors;
    std::transform(centered.begin(), centered.end(), centered.begin(), [=](auto& e) {
      return e - mean;
    });
    const auto sum2 = std::inner_product(centered.begin(), centered.end(), centered.begin(), T());
    return std::inner_product(this->m_values.begin(), this->m_values.end(), centered.begin(), T()) /
        std::sqrt(m_sum2 * sum2);
  }

private:

  T m_sum2;
};

/**
 * @brief Quotient filter, i.e. minimum value of the ratio between neighbors and template, normalized.
*/
template <typename T, typename TWindow>
class QuotientFilter : public KernelMixin<T, TWindow, QuotientFilter<T, TWindow>> {
public:

  using KernelMixin<T, TWindow, QuotientFilter>::KernelMixin;

  template <typename TIn>
  T operator()(const TIn& neighbors) const
  {
    T out = std::numeric_limits<T>::max();
    T norm2 = 0;
    auto tit = this->m_values.begin();
    for (auto n : neighbors) {
      const auto q = n / *tit;
      norm2 += q * q;
      if (q < out) { // FIXME more robust? median? quantile?
        out = q;
      }
      ++tit;
    }
    return out * std::sqrt(this->m_values.size() / norm2);
  }

  void init_impl() {} // FIXME private
};

template <typename TIn, typename TPsf>
Raster<typename TPsf::Value> quotient(const TIn& in, const TPsf& psf)
{
  using T = typename TPsf::Value;
  auto filter = SimpleFilter<QuotientFilter<T, Box<2>>>(psf.domain() - (psf.shape() - 1) / 2, psf.begin(), psf.end());
  return filter * extrapolation<Nearest>(in);
}

template <typename TIn, typename TPsf>
Raster<typename TPsf::Value> match(const TIn& in, const TPsf& psf)
{
  using T = typename TPsf::Value;
  auto filter =
      SimpleFilter<PearsonCorrelation<T, Box<2>>>(psf.domain() - (psf.shape() - 1) / 2, psf.begin(), psf.end());
  return filter * extrapolation<Nearest>(in);
}

template <typename TIn>
Raster<typename TIn::Value> laplacian(const TIn& in)
{
  using T = typename TIn::Value;
  const auto filter = convolution(
      Raster<T>({3, 3}, {-1. / 6., -2. / 3., -1. / 6., -2. / 3., 10. / 3., -2. / 3., -1. / 6., -2. / 3., -1. / 6.}));
  return filter * extrapolation<Nearest>(in);
}

template <typename TIn>
Raster<typename TIn::Value> dilate(const TIn& in, Index radius = 1)
{
  using T = typename TIn::Value;
  auto filter = dilation<T>(Box<2>::from_center(radius)); // FIXME L2-ball?
  return filter * extrapolation<Nearest>(in);
}

template <typename TIn>
Raster<typename TIn::Value> blur(const TIn& in, Index radius = 1)
{
  using T = typename TIn::Value;
  auto filter = mean_filter<T>(Box<2>::from_center(radius)); // FIXME L2-ball?
  return filter * extrapolation<Nearest>(in);
}

/**
 * @brief Detect cosmic rays.
 * @param in The input 2D raster
 * @param pfa The detection probability of false alarm
 * 
 * This is a simple adaptive Laplacian thresholding.
 * The input raster is convolved with some Laplacian kernel.
 * The parameters of the background noise (empirically assumed Laplace-distributed)
 * of the filtered image are estimated to deduce the detection threshold from a PFA.
 */
template <typename TIn, typename TPsf>
Raster<char> detect(const TIn& in, const TPsf& psf, float pfa, float tq)
{
  auto laplacian_map = laplacian(in);
  Fits("/tmp/cosmic.fits").write(laplacian_map, 'a');
  float n = 0;
  float s = 0;
  for (auto e : laplacian_map) {
    if (not std::isnan(e)) {
      n += std::abs(e);
      ++s;
    }
  }
  printf("Norm: %f\n", n);

  // Empirically assume Laplace distribution
  // const auto tl = -norm<1>(laplacian_map) / laplacian_map.size() * std::log(2.0 * pfa);
  const auto tl = -n / s * std::log(2.0 * pfa);
  printf("Threshold: %f\n", tl);

  const Index radius = std::sqrt(psf.size()) / 4;
  printf("Radius: %li\n", radius);
  auto quotient_map = dilate(quotient(in, psf), radius);
  Fits("/tmp/cosmic.fits").write(quotient_map, 'a');

  Raster<char> out(in.shape());
  out.generate(
      [=](auto l, auto q) {
        return l > tl && q < tq; // FIXME compute quotient only where l > tl
      },
      laplacian_map,
      quotient_map);
  return out;
}

/**
 * @brief Compute the minimum contrast between a point and its neighbors in a mask.
 * 
 * The contrast is negative when the point intensity is higher than that of the neighborhood.
 */
template <typename TIn, typename TMask>
float min_contrast(const TIn& in, const TMask& mask, const Position<2>& p)
{
  auto out = std::numeric_limits<float>::max();
  const Position<2> dx {0, 1};
  const Position<2> dy {1, 0};
  for (const auto& neighbor : {p - dx, p + dx, p - dy, p + dy}) { // FIXME optimize?
    if (mask[neighbor]) {
      auto contrast = (in[neighbor] - in[p]) / in[neighbor]; // FIXME assumes in > 0
      if (contrast < out) {
        out = contrast;
      }
    }
  }
  return out;
}

/**
 * @brief Segment detected cosmic rays.
 * @param mask The detection map
 * @param threshold The similarity threshold
 * 
 * Given a detection map, neighbors of flagged pixels are considered as candidate cosmic rays.
 * Some similarity distance is computed in the neighborhood in order to decide
 * whether the cadidate belongs to the cosmic ray or to the background, by thresholding.
 */
template <typename TIn, typename TMask>
void segment(const TIn& in, TMask& mask, float threshold)
{
  // FIXME Mask<2>::ball<1>(1)
  auto candidates = dilation<typename TMask::Value>(Box<2>::from_center(1)) * extrapolation(mask, '\0') - mask;
  for (const auto& p : candidates.domain() - Box<2>::from_center(1)) {
    if (candidates[p]) {
      if (min_contrast(in, mask, p) < threshold) {
        mask[p] = true;
      }
    }
  }
}

} // namespace Cosmics
} // namespace Linx

#endif
