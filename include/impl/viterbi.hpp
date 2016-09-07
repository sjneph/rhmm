/*
  FILE: viterbi.hpp
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

#ifndef VITERBI_HMM_R_HPP
#define VITERBI_HMM_R_HPP

#include <array>
#include <cmath>
#include <vector>

#include "efun.hpp"

namespace ci {

namespace hmm {

  template <typename O, typename I, typename T, typename E, typename OutIter>
  void viterbi(const O& observed,
               const I& initial,
               const T& transition,
               const E& emission,
               OutIter out) {
    typedef typename O::value_type U;
    const std::size_t nstates = initial.size();
    std::size_t nobs = observed.size();
    std::vector<std::array<U,2>> delta(nstates); // only need [2] x [n_states] 2-d array
    std::size_t index = 0;
    for ( std::size_t i = 0; i < nstates; ++i ) {
      delta[i][0] = elnproduct(initial[i], emission[i][observed[0]]);
      if ( delta[i][0] > delta[index][0] )
        index = i;
    } // for

    *out++ = index;
    index = 0;
    U gmx = 0;
    std::size_t active = 0, passive = active + 1;
    for ( std::size_t s = 1; s < nobs; ++s ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        U mx = elnproduct(delta[0][active], transition[0][j]), tmp = mx;
        for ( std::size_t k = 1; k < nstates; ++k ) {
          tmp = elnproduct(delta[k][active], transition[k][j]);
          if ( tmp > mx )
            mx = tmp;
        } // for
        delta[j][passive] = elnproduct(mx, emission[j][observed[s]]);
        if ( delta[j][passive] > gmx || 0 == j )
          gmx = delta[j][passive], index = j;
      } // for
      *out++ = index;
      std::swap(active, passive);
    } // for
  }

} // namespace hmm

} // namespace ci

#endif // VITERBI_HMM_R_HPP
