// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#include "Linx/Data/Raster.h"
#include "Linx/Transforms/Filters.h"
#include "Linx/Transforms/SimpleFilter.h"

#include <boost/test/unit_test.hpp>

using namespace Linx;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(SimpleFilter_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(constant0_3x3_test)
{
  Raster<int, 2> in({4, 3});
  in.fill(1);
  const auto extra = extrapolation(in, 0);
  const auto box = Box<2>::from_center(1);

  auto median_out = median_filter<int>(box) * extra;
  std::vector<int> median_expected {0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0};
  BOOST_TEST(median_out.container() == median_expected);

  auto erode_out = erosion<int>(box) * extra;
  auto min_out = minimum_filter<int>(box) * extra;
  std::vector<int> erode_expected {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  BOOST_TEST(erode_out.container() == erode_expected);
  BOOST_TEST(min_out.container() == erode_expected);

  auto dilate_out = dilation<int>(box) * extra;
  auto max_out = maximum_filter<int>(box) * extra;
  std::vector<int> dilate_expected {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  BOOST_TEST(dilate_out.container() == dilate_expected);
  BOOST_TEST(max_out.container() == dilate_expected);
}

BOOST_AUTO_TEST_CASE(pixel_test)
{
  const auto in = Raster<double>({4, 3}).range();
  const auto extra = extrapolation<Nearest>(in);
  const auto k = convolution(Raster<float>({2, 2}).range());
  const auto out = k * extra;
  for (const auto& p : in.domain()) {
    auto out_p = k * extra(p);
    BOOST_TEST(out_p == out[p]);
  }
  for (const auto& p : in.domain() - k.window()) {
    auto out_p = k * in(p);
    BOOST_TEST(out_p == out[p]);
  }
}

BOOST_AUTO_TEST_CASE(pixel_sequence_test)
{
  const auto in = Raster<double>({4, 3}).range();
  const auto extra = extrapolation<Nearest>(in);
  const auto k = convolution(Raster<float>({2, 2}).range());
  Sequence<Position<2>> positions(in.size());
  std::copy(in.domain().begin(), in.domain().end(), positions.begin());
  const auto out = k * extra;
  const auto out_seq = k * extra(positions);
  for (std::size_t i = 0; i < in.size(); ++i) {
    BOOST_TEST(out_seq[i] == out[i]);
  }
}

BOOST_AUTO_TEST_CASE(pixel_list_test)
{
  const auto in = Raster<double>({4, 3}).range();
  const auto extra = extrapolation<Nearest>(in);
  const auto k = convolution(Raster<float>({2, 2}).range());
  const auto out = k * extra;
  const auto out_seq = k * extra(Position<2>::zero(), Position<2>::one(), Position<2> {2, 2});
  for (Index i = 0; i <= 2; ++i) {
    auto out_i = out[{i, i}];
    BOOST_TEST(out_seq[i] == out_i);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
