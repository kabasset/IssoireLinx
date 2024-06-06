// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXTRANSFORMS_MIXINS_KERNEL_H
#define _LINXTRANSFORMS_MIXINS_KERNEL_H

#include "Linx/Base/TypeUtils.h"
#include "Linx/Transforms/mixins/StructuringElement.h"

#include <initializer_list>

namespace Linx {

/**
 * @brief Mixin for kernel-based filters.
 */
template <typename T, typename TWindow, typename TDerived>
class KernelMixin : public StructuringElementMixin<T, TWindow, TDerived> {
public:

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  KernelMixin(TWindow window, TIt begin, TIt end) :
      StructuringElementMixin<T, TWindow, TDerived>(LINX_MOVE(window)), m_values(LINX_MOVE(begin), LINX_MOVE(end))
  {
    LINX_CRTP_DERIVED.init_impl();
  }

  /**
   * @brief List-based constructor.
   */
  KernelMixin(TWindow window, std::initializer_list<T> values) :
      KernelMixin(LINX_MOVE(window), values.begin(), values.end())
  {}

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  KernelMixin(TWindow window, TRange&& values) :
      KernelMixin(LINX_MOVE(window), LINX_FORWARD(values).begin(), LINX_FORWARD(values).end())
  {}

protected:

  /**
   * @brief The kernel values.
   */
  std::vector<T> m_values;
};

} // namespace Linx

#endif
