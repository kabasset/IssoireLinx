// @copyright 2022, Antoine Basset (CNES)
// This file is part of Raster <github.com/kabasset/Raster>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "LitlRaster/Grid.h"
#include "LitlRaster/Raster.h"

#include <boost/test/unit_test.hpp>

using namespace Litl;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Grid_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(grid_is_screened_in_order_test) {
  Grid<2> region(Box<2>({1, 2}, {6, 7}), {2, 3});
  std::vector<Position<2>> expected {{1, 2}, {3, 2}, {5, 2}, {1, 5}, {3, 5}, {5, 5}};
  std::vector<Position<2>> out;
  for (const auto& p : region) {
    out.push_back(p);
  }
  BOOST_TEST(out == expected);
}

BOOST_AUTO_TEST_CASE(grid_clamp_is_shrunk_test) {

  const Grid<1> in({{1}, {9}}, {3});
  // -+--+--+--
  BOOST_TEST(in.front()[0] == 1);
  BOOST_TEST(in.back()[0] == 7);
  BOOST_TEST(in.step()[0] == 3);

  const auto out2 = clamp(in, {{2}, {8}});
  //  --+--+-
  BOOST_TEST(out2.front()[0] == 4);
  BOOST_TEST(out2.back()[0] == 7);
  BOOST_TEST(out2.step()[0] == 3);

  const auto out4 = clamp(in, {{4}, {8}});
  //    +--+-
  BOOST_TEST(out4.front()[0] == 4);
  BOOST_TEST(out4.back()[0] == 7);
  BOOST_TEST(out4.step()[0] == 3);

  const auto out6 = clamp(in, {{6}, {8}});
  //      -+-
  BOOST_TEST(out6.front()[0] == 7);
  BOOST_TEST(out6.back()[0] == 7);
  BOOST_TEST(out6.step()[0] == 3);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
