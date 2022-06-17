// Copyright (C) 2022, Antoine Basset
// This file is part of Cnes.Raster <github.com/kabasset/Raster>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "RasterGeometry/Interpolation.h"

#include <boost/test/unit_test.hpp>

using namespace Cnes;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Interpolation_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(constant_nn_test) {
  Raster<int, 3> raster({2, 2, 2});
  raster.fill(1);

  Position<3> position {-1, -1, -1};
  Vector<double, 3> vector {.5, .5, .5};

  const auto extra = extrapolate(raster, 0);
  BOOST_TEST(extra[position] == 0);

  const auto inter = interpolate<NearestNeighbor>(raster);
  BOOST_TEST(&inter[position] < raster.data()); // Out of bounds
  BOOST_TEST(inter.at(vector) == 1);

  const auto interextra = interpolate<NearestNeighbor>(extra);
  BOOST_TEST(interextra[position] == 0);
  BOOST_TEST(interextra.at(vector) == 1);
}

BOOST_AUTO_TEST_CASE(linear_test) {

  Raster<int, 3> raster({2, 2, 2});
  raster.arange(1);

  const auto interpolator = interpolate<Linear>(raster);

  const auto front = interpolator[{0, 0, 0}];
  const auto back = interpolator[{1, 1, 1}];
  const auto center = interpolator.at({.5, .5, .5});

  BOOST_TEST(front == 1);
  BOOST_TEST(back == 8);
  BOOST_TEST(center == 4.5);
}

BOOST_AUTO_TEST_CASE(cubic_test) {

  Raster<int, 3> raster({4, 4, 4});
  raster.arange(1);

  const auto interpolator = interpolate<Cubic>(raster);

  const auto front = interpolator[{0, 0, 0}];
  const auto back = interpolator[{3, 3, 3}];
  const auto pos = interpolator.at({1, 1, 1});
  const auto center = interpolator.at({1.5, 1.5, 1.5});

  BOOST_TEST(front == 1);
  BOOST_TEST(back == 64);
  BOOST_TEST(pos == 22);
  BOOST_TEST(center == 32.5);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
