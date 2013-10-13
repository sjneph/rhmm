/*
  FILE: infinity.hpp
  AUTHOR: Shane Neph
  CREATE DATE: Sat Oct 12 22:07:12 PDT 2013
*/

//    Hidden Markov Model
//    Copyright (C) 2013 Shane Neph
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

#ifndef INFINITY_HMM_R_HPP
#define INFINITY_HMM_R_HPP

#include <climits>

namespace ci {

  namespace details {
    template <bool b, typename T>
    struct infinite {
      static const T value;
    };

    template <bool b, typename T>
    const T infinite<b, T>::value = std::numeric_limits<T>::max();

    template <typename T>
    struct infinite<true, T> {
      static const T value;
    };

    template <typename T>
    const T infinite<true, T>::value = std::numeric_limits<T>::infinity();

  } // namespace details

  //================================================
  // Get infinity representation for numeric type T
  //================================================
  template <typename T>
  inline T inf() {
    return(details::infinite< std::numeric_limits<T>::has_infinity, T >::value);
  }

} // namespace ci

#endif // INFINITY_HMM_R_HPP
