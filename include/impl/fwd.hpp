/*
  FILE: fwd.hpp
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

#ifndef FWD_HMM_R_HPP
#define FWD_HMM_R_HPP

#include <vector>

#include "efun.hpp"
#include "infinity.hpp"

namespace ci {

namespace hmm {

  //==========================
  // forward_full() algorithm
  //  - calculates & retains all calculated alpha values up to index
  //==========================
  template <typename O, typename I, typename T, typename E, typename U>
  void forward_full(const O& observed,
                    const I& initial,
                    const T& transition,
                    const E& emission,
                    std::size_t index,
                    std::vector< std::vector<U> >& alpha) {
    if ( index < 1 )
      return;

    const std::size_t nstates = initial.size();
    for ( std::size_t i = 0; i < nstates; ++i )
      alpha[i][0] = elnproduct(initial[i], emission[i][observed[0]]);

    U tmpf = inf<U>();
    for ( std::size_t s = 1; s < index; ++s ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        tmpf = inf<U>();
        for ( std::size_t k = 0; k < nstates; ++k )
          tmpf = elnsum(tmpf, elnproduct(alpha[k][(s-1)], transition[k][j]));
        alpha[j][s] = elnproduct(tmpf, emission[j][observed[s]]);
      } // for
    } // for
  }

  //===========================
  // forward_index() algorithm
  //  - requires minimal memory to calculate alpha at a single "time" index
  //===========================
  template <typename O, typename I, typename T, typename E, typename U>
  void forward_index(const O& observed,
                     const I& initial,
                     const T& transition,
                     const E& emission,
                     std::size_t index,
                     std::vector<U>& alpha) {
    if ( index < 1 )
      return;

    const std::size_t nstates = initial.size();
    std::vector< std::vector<U> > lcl(nstates);
    for ( std::size_t i = 0; i < nstates; ++i ) {
      lcl[i].resize(2, 0);
      lcl[i][0] = elnproduct(initial[i], emission[i][observed[0]]);
    } // for

    U tmpf = inf<U>();
    std::size_t active = 0, passive = active + 1;
    for ( std::size_t s = 1; s < index; ++s ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        tmpf = inf<U>();
        for ( std::size_t k = 0; k < nstates; ++k )
          tmpf = elnsum(tmpf, elnproduct(lcl[k][active], transition[k][j]));
        lcl[j][passive] = elnproduct(tmpf, emission[j][observed[s]]);
      } // for
      std::swap(active, passive);
    } // for

    for ( std::size_t idx = 0; idx < nstates; ++idx )
      alpha[idx] = lcl[idx][active];
  }

  //==========================
  // forward_next() algorithm
  //  - calculate next forward_index() given last result
  //     giving much needed memory back to the system if used correctly
  //==========================
  template <typename O, typename I, typename T, typename E, typename U>
  void forward_next(const O& observed,
                    const I& initial,
                    const T& transition,
                    const E& emission,
                    std::size_t index,
                    std::vector<U>& alpha) {
    if ( index < 1 )
      return;

    const std::size_t nstates = initial.size();
    if ( 1 == index ) {
      for ( std::size_t i = 0; i < nstates; ++i )
        alpha[i] = elnproduct(initial[i], emission[i][observed[0]]);
      return;
    }
    const std::vector<U> lcl(alpha);

    U tmpf = inf<U>();
    for ( std::size_t j = 0; j < nstates; ++j ) {
      tmpf = inf<U>();
      for ( std::size_t k = 0; k < nstates; ++k )
        tmpf = elnsum(tmpf, elnproduct(lcl[k], transition[k][j]));
      alpha[j] = elnproduct(tmpf, emission[j][observed[index-1]]);
    } // for
  }

} // namespace hmm

} // namespace ci

#endif // FWD_HMM_R_HPP
