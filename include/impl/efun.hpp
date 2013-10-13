/*
  FILE: efun.hpp
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

#ifndef EFUN_HMM_R_CPP
#define EFUN_HMM_R_CPP

#include <cmath>
#include <limits>

#include "infinity.hpp"

namespace ci {

namespace hmm {

  // extended-log-sum
  template <typename T>
  inline T elnsum(T x, T y) {
    static const T infinite = inf<T>();
    if ( x == infinite ) {
      if ( y == infinite )
        return(infinite);
      return(y);
    }
    else if ( y == infinite )
      return(x);

    if ( x > y )
      return(x + std::log(1 + std::exp(y - x)));
    return(y + std::log(1 + std::exp(x - y)));
  }

  // extended-log-product
  template <typename T>
  inline T elnproduct(T x, T y) {
    static const T infinite = inf<T>();
    if ( x == infinite || y == infinite )
      return(infinite);
    return(x + y);
  }

} // namespace hmm

} // namespace ci


#endif // EFUN_HMM_R_CPP
