// Copyright (C) 2022, CNES
// This file is part of Raster <github.com/kabasset/Raster>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "Raster/Raster.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(RasterDemoConstructors_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(any_raster_ctors_test) {

  //! [Any raster]
  std::vector<int> vec {1, 2, 3, 4, 5, 6};
  const int* data = vec.data();

  Cnes::VecRaster<int> defaultInitialized({3, 2});
  Cnes::VecRaster<int> listInitialized({3, 2}, {1, 2, 3, 4, 5, 6});
  Cnes::VecRaster<int> copiedFromPointer({3, 2}, data);
  Cnes::VecRaster<int> copiedFromIterable({3, 2}, vec);
  Cnes::VecRaster<int> movedFromContainer({3, 2}, std::move(vec));
  //! [Any raster]

  for (std::size_t i = 0; i < 6; ++i) {
    BOOST_TEST(defaultInitialized[i] == 0);
    BOOST_TEST(listInitialized[i] == i + 1);
    BOOST_TEST(copiedFromPointer[i] == i + 1);
    BOOST_TEST(copiedFromIterable[i] == i + 1);
    BOOST_TEST(movedFromContainer[i] == i + 1);
  }
  BOOST_TEST(movedFromContainer.data() == data);

  const int* listData = listInitialized.data();
  //! [Any raster copy-move]
  Cnes::VecRaster<int> copiedFromRaster(listInitialized);
  Cnes::VecRaster<int> movedFromRaster(std::move(listInitialized));
  //! [Any raster copy-move]

  for (std::size_t i = 0; i < 6; ++i) {
    BOOST_TEST(copiedFromRaster[i] == i + 1);
    BOOST_TEST(movedFromRaster[i] == i + 1);
  }
  BOOST_TEST(movedFromRaster.data() == listData);
}

BOOST_AUTO_TEST_CASE(ptrraster_ctors_test) {

  //! [PtrRaster write]
  int data[] = {1, 2, 3, 4, 5, 6};

  Cnes::PtrRaster<int> constructed({3, 2}, data);
  auto made = Cnes::makeRaster(data, 3, 2);

  constructed[0] = 42;
  made[1] = 12;
  //! [PtrRaster write]

  BOOST_TEST(constructed[0] == 42);
  BOOST_TEST(constructed[1] == 12);
  BOOST_TEST(made[0] == 42);
  BOOST_TEST(made[1] == 12);

  //! [PtrRaster read]
  const int* cData = data;

  Cnes::PtrRaster<const int> cConstructed({3, 2}, data);
  auto cMade = Cnes::makeRaster(cData, 3, 2);
  //! [PtrRaster read]

  BOOST_TEST((made == constructed));
  BOOST_TEST((cMade == cConstructed));
}

BOOST_AUTO_TEST_CASE(vecraster_ctors_test) {
  //! [VecRaster]
  Cnes::VecRaster<int> filledWithValue({3, 2}, 42);
  //! [VecRaster]

  for (auto e : filledWithValue) {
    BOOST_TEST(e == 42);
  }
}

BOOST_AUTO_TEST_CASE(alignedraster_ctors_test) {

  //! [AlignedRaster owns]
  Cnes::AlignedRaster<int> defaultAligned({3, 2}); // SIMD-aligned
  Cnes::AlignedRaster<int> longerAligned({3, 2}, nullptr, 1024); // Longer alignment
  //! [AlignedRaster owns]

  BOOST_TEST(defaultAligned.alignment() % 16 == 0);
  BOOST_TEST(longerAligned.alignment() % 1024 == 0);

  //! [AlignedRaster shares]
  int data[] = {1, 2, 3, 4, 5, 6}; // Unknown alignment

  Cnes::AlignedRaster<int> notAligned({3, 2}, data, 1); // No alignment required

  try {
    Cnes::AlignedRaster<int> maybeAligned({3, 2}, data, 64); // Alignment required
    std::cout << "Data is aligned!" << std::endl;
  } catch (Cnes::AlignmentError& e) { // Alignment requirement not met
    BOOST_TEST(notAligned.alignment() % 64 != 0);
    std::cout << e.what() << std::endl;
  }
  //! [AlignedRaster shares]

  BOOST_TEST(notAligned.data() == data);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
