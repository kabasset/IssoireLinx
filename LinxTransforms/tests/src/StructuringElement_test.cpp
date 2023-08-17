// @copyright 2022, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Raster>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "Linx/Data/Raster.h"
#include "Linx/Transforms/StructuringElement.h"

#include <boost/test/unit_test.hpp>

using namespace Linx;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(StructuringElement_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(constant0_3x3_test)
{
  Raster<int, 2> in({4, 3});
  in.fill(1);
  const auto extra = extrapolate(in, 0);
  const auto box = Box<2>::from_center(1);

  auto median_out = median_filter<int>(box) * extra;
  std::vector<int> median_expected {0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0};
  BOOST_TEST(median_out.container() == median_expected);

  auto erode_out = erosion<int>(box) * extra;
  std::vector<int> erode_expected {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  BOOST_TEST(erode_out.container() == erode_expected);

  auto dilate_out = dilation<int>(box) * extra;
  std::vector<int> dilate_expected {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  BOOST_TEST(dilate_out.container() == dilate_expected);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
