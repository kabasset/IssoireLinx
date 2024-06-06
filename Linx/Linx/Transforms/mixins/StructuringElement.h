// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXTRANSFORMS_MIXINS_STRUCTURINGELEMENT_H
#define _LINXTRANSFORMS_MIXINS_STRUCTURINGELEMENT_H

namespace Linx {

/**
 * @brief Mixin for footprint-based filters.
 */
template <typename T, typename TWindow, typename TDerived>
class StructuringElementMixin {
public:

  /**
   * @brief The kernel type.
   */
  using Value = T;

  /**
   * @brief The kernel dimension.
   */
  static constexpr Index Dimension = TWindow::Dimension;

  /**
   * @brief The kernel window type.
  */
  using Window = TWindow;

  /**
   * @brief Constructor.
   */
  StructuringElementMixin(TWindow window) : m_window(LINX_MOVE(window)) {}

  /**
   * @brief Get the window.
   */
  const TWindow& window() const
  {
    return m_window;
  }

private:

  /**
   * @brief The window.
   */
  TWindow m_window;
};

} // namespace Linx

#endif
