/*
  FILE: bkd.hpp
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


#ifndef BKD_HMM_R_HPP
#define BKD_HMM_R_HPP

#include <list>
#include <vector>

#include "efun.hpp"
#include "infinity.hpp"

namespace ci {

namespace hmm {

  //=============================
  // backward_full() algorithm()
  //  - calculates & retains all calculated beta values down to index
  //=============================
  template <typename O, typename I, typename T, typename E, typename U>
  void backward_full(const O& observed,
                     const I& initial,
                     const T& transition,
                     const E& emission,
                     std::size_t index,
                     std::vector< std::vector<U> >& beta) {
    std::size_t nobs = observed.size();
    std::size_t nstates = initial.size();
    if ( nobs < 2 )
      return;
    else if ( index > observed.size() || index < 1 )
      return;

    for ( std::size_t i = 0; i < nstates; ++i ) {
      for ( std::size_t j = nobs-1; j >= index; ) {
        beta[i][j] = 0;
        if ( 0 == j-- )
          break;
      } // for
    } // for

    for ( std::size_t s = nobs-1; s >= index; ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        U tmpf = inf<U>();
        for ( std::size_t k = 0; k < nstates; ++k ) {
          tmpf = elnsum(tmpf,
                        elnproduct(transition[j][k],
                                   elnproduct(emission[k][observed[s]],
                                              beta[k][s])));
        } // for
        beta[j][s-1] = tmpf;
      } // for
      if ( 0 == s-- )
        break;
    } // for
  }

  //==============================
  // backward_index() algorithm()
  //  - requires minimal memory to calculate beta at a single "time" index
  //==============================
  template <typename O, typename I, typename T, typename E, typename U>
  void backward_index(const O& observed,
                      const I& initial,
                      const T& transition,
                      const E& emission,
                      std::size_t index,
                      std::vector<U>& beta) {
    std::size_t nobs = observed.size();
    std::size_t nstates = initial.size();
    if ( nobs < 2 )
      return;
    else if ( index > observed.size() || index < 1 )
      return;

    std::vector< std::vector<U> > lcl(nstates);
    for ( std::size_t i = 0; i < nstates; ++i )
      for ( std::size_t j = 0; j < 2; ++j )
        lcl[i].resize(2), lcl[i][j] = 0;

    std::size_t active = 0, passive = active + 1;
    for ( std::size_t s = nobs-1; s >= index; ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        U tmpf = inf<U>();
        for ( std::size_t k = 0; k < nstates; ++k ) {
          tmpf = elnsum(tmpf,
                        elnproduct(transition[j][k],
                                   elnproduct(emission[k][observed[s]],
                                              lcl[k][active])));
        } // for
        lcl[j][passive] = tmpf;
      } // for
      std::swap(active, passive);
      if ( 0 == s-- )
        break;
    } // for

    for ( std::size_t idx = 0; idx < nstates; ++idx )
      beta[idx] = lcl[idx][active];
  }

  //==========================
  // backward_next() algorithm
  //  - calculate next backward_index() given last result
  //     giving much needed memory back to the system if used properly
  //==========================
  template <typename O, typename I, typename T, typename E, typename U>
  void backward_next(const O& observed,
                     const I& initial,
                     const T& transition,
                     const E& emission,
                     std::size_t index,
                     std::vector<U>& beta) {

    std::size_t nobs = observed.size();
    std::size_t nstates = initial.size();
    if ( nobs < 2 )
      return;
    else if ( index > observed.size() || index < 1 )
      return;

    if ( nobs == index ) {
      for ( std::size_t i = 0; i < nstates; ++i )
        beta[i] = 0;
      return;
    }
    std::vector<U> lcl(beta);
    U tmpf = inf<U>();
    for ( std::size_t j = 0; j < nstates; ++j ) {
      tmpf = inf<U>();
      for ( std::size_t k = 0; k < nstates; ++k ) {
        tmpf = elnsum(tmpf,
                      elnproduct(transition[j][k],
                                 elnproduct(emission[k][observed[index]],
                                            lcl[k])));
      } // for
      beta[j] = tmpf;
    } // for
  }

  //===========================
  // backward_enext() algorithm
  //  - calculate next backward_index() given last result
  //     but keeps intermediate results for potential matrix-LU
  //     calculations later.  Works when fully-connected model, all
  //     transition probs > 0 and matrix is invertible.  When
  //     this is the case, nothing is better in memory.
  //==========================
  template <typename O, typename I, typename T, typename E, typename U>
  void backward_enext(const O& observed,
                      const I& initial,
                      const T& transition,
                      const E& emission,
                      std::size_t index,
                      std::vector< std::vector<U> >& beta_internals) {

    const std::size_t nobs = observed.size();
    const std::size_t nstates = initial.size();
    if ( nobs < 2 )
      return;
    else if ( index > observed.size() || index < 1 )
      return;

    if ( nobs == index ) {
      for ( std::size_t i = 0; i < nstates; ++i )
        beta_internals[nstates-1][i] = 0;
      return;
    }

    std::vector<U> lcl(beta_internals[nstates-1].begin(), beta_internals[nstates-1].end());
    U tmpf = inf<U>();
    for ( std::size_t j = 0; j < nstates; ++j ) {
      tmpf = inf<U>();
      for ( std::size_t k = 0; k < nstates; ++k ) {
        tmpf = elnsum(tmpf,
                      elnproduct(transition[j][k],
                                 elnproduct(emission[k][observed[index]],
                                            lcl[k])));
        beta_internals[k][j] = tmpf;
      } // for
    } // for
  }

} // namespace hmm

} // namespace ci

#endif // BKD_HMM_R_HPP
