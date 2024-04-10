// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#include "Linx/Run/ProgramOptions.h"
#include "Linx/Run/Timer.h"
#include "Linx/Transforms/Filters.h"

#include <map>
#include <string>

using Image = Linx::Raster<float>;
using Duration = std::chrono::milliseconds;

void filter_monolith(Image& image, const Image& values)
{
  auto kernel = Linx::convolution(values);
  const auto extrapolation = Linx::extrapolation<Linx::Nearest>(image);
  const auto extrapolated = extrapolation.copy(image.domain() + kernel.window());
  auto patch = extrapolated(kernel.window() - kernel.window().front());
  auto out_it = image.begin();
  const auto rbegin = std::reverse_iterator(values.end());
  const auto rend = std::reverse_iterator(values.begin());
  for (const auto& p : image.domain()) {
    patch >>= p;
    *out_it = std::inner_product(rbegin, rend, patch.begin(), 0.F);
    ++out_it;
    patch <<= p;
  }
}

void filter_hardcoded(Image& image, const Image& values)
{
  auto kernel = Linx::convolution(values);
  const auto extrapolation = Linx::extrapolation<Linx::Nearest>(image);
  const auto extrapolated = extrapolation.copy(image.domain() + kernel.window());
  auto it = image.end();
  Image out(image.shape());
  auto out_it = out.begin();
  const auto inner = image.domain() - kernel.window().front();
  auto patch = extrapolated(inner);
  for (const auto& q : kernel.window()) {
    --it;
    patch >>= q;
    const auto k = *it;
    for (const auto& v : patch) {
      *out_it += k * v;
      ++out_it;
    }
    patch <<= q;
    out_it = out.begin();
  }
  image = out;
}

template <typename TDuration>
TDuration filter(Image& image, const Image& kernel, char setup)
{
  Linx::Timer<TDuration> timer;
  timer.start();
  switch (setup) {
    case '0':
      image = Linx::convolution(kernel) * Linx::extrapolation(image, 0.0F);
      break;
    case 'd':
      image = Linx::convolution(kernel) * Linx::extrapolation<Linx::Nearest>(image);
      break;
    case 's':
      image = Linx::sparse_convolution(kernel) * Linx::extrapolation<Linx::Nearest>(image);
      break;
    case 'm':
      filter_monolith(image, kernel);
      break;
    case 'h':
      filter_hardcoded(image, kernel);
      break;
    default:
      throw std::runtime_error("Case not implemented"); // FIXME CaseNotImplemented
  }
  return timer.stop();
}

int main(int argc, char const* argv[])
{
  Linx::ProgramOptions options;
  options.named("case", "Test case: d (default), m (monolith), h (hardcoded)", 'd');
  options.named("image", "Raster length along each axis", 2048L);
  options.named("kernel", "Kernel length along each axis", 5L);
  options.named("sparse", "Kernel sparsity", 0.);
  options.parse(argc, argv);
  const auto setup = options.as<char>("case");
  const auto image_diameter = options.as<Linx::Index>("image");
  const auto kernel_diameter = options.as<Linx::Index>("kernel");
  const auto kernel_sparsity = options.as<double>("sparse");

  Linx::Position<2> image_shape {image_diameter, image_diameter};
  Linx::Position<2> kernel_shape {kernel_diameter, kernel_diameter};

  std::cout << "Generating raster and kernel..." << std::endl;
  auto image = Image(image_shape).range();
  const auto mask = Image(kernel_shape).generate(Linx::ImpulseNoise<bool>(true, 1 - kernel_sparsity));
  const auto kernel = Image(kernel_shape).range() * mask;
  std::cout << "  input: " << image << std::endl;
  std::cout << "  kernel: " << kernel << " (size: " << Linx::sum(mask) << ")" << std::endl;

  std::cout << "Filtering..." << std::endl;
  const auto duration = filter<Duration>(image, kernel, setup);
  std::cout << "  output: " << image << std::endl;

  std::cout << "  Done in " << duration.count() << "ms" << std::endl;

  return 0;
}
