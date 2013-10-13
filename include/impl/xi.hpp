/*
  FILE: xi.hpp
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

#ifndef XI_HMM_R_HPP
#define XI_HMM_R_HPP

#include <vector>

#include "bkd.hpp"
#include "efun.hpp"
#include "fwd.hpp"
#include "infinity.hpp"


namespace ci {

namespace hmm {

  //=====================
  // xi_full()
  //  - Evaluate the probability of q_t being in state i and q_(t+1) being
  //     in state j given observations and model.
  //  - Computes all N*N*T probabilities and stores in probs
  //=====================
  template <typename O, typename I, typename T, typename E, typename U>
  void xi_full(const O& observed,
               const I& initial,
               const T& transition,
               const E& emission,
               std::vector< std::vector< std::vector<U> > >& probs) {

    const std::size_t nstates = initial.size();
    const std::size_t nobs = observed.size();
    if ( nobs < 1 )
      return;

    std::vector< std::vector<U> > alpha(nstates), beta(nstates);
    for ( std::size_t i = 0; i < alpha.size(); ++i )
      alpha[i].resize(nobs, 0), beta[i].resize(nobs, 0);

    forward_full(observed, initial, transition, emission, nobs, alpha);
    backward_full(observed, initial, transition, emission, 1, beta);

    U normalizer = inf<U>();
    for ( std::size_t s = 0; s < nobs-1; ++s ) {
      normalizer = inf<U>();
      for ( std::size_t i = 0; i < nstates; ++i ) {
        for ( std::size_t j = 0; j < nstates; ++j ) {
          probs[i][j][s] = elnproduct(alpha[i][s],
                                      elnproduct(transition[i][j],
                                                 elnproduct(emission[j][observed[s+1]],
                                                            beta[j][s+1])));
          normalizer = elnsum(normalizer, probs[i][j][s]);
        } // for
      } // for

      for ( std::size_t k = 0; k < nstates; ++k )
        for ( std::size_t m = 0; m < nstates; ++m )
          probs[k][m][s] = elnproduct(probs[k][m][s], -normalizer);
    } // for
  }


  //=====================
  // xi()
  //  - Evaluate the probability of q_t being in state i and q_(t+1) being
  //     in state j given observations and model.
  //  - Minimizes memory requirements if used correctly: N*N
  //=====================
  template <typename O, typename I, typename T, typename E, typename U>
  void xi(const O& observed,
          const I& initial,
          const T& transition,
          const E& emission,
          int index,
          const std::vector<U>& beta,
          std::vector<U>& alpha,
          std::vector< std::vector<U> >& probs) {

    const std::size_t nstates = initial.size();
    forward_next(observed, initial, transition, emission, index, alpha);

    U normalizer = inf<U>();
    for ( std::size_t i = 0; i < nstates; ++i ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        probs[i][j] = elnproduct(alpha[i],
                                 elnproduct(transition[i][j],
                                            elnproduct(emission[j][observed[index]],
                                                       beta[j])));
        normalizer = elnsum(normalizer, probs[i][j]);
      } // for
    } // for

    for ( std::size_t k = 0; k < nstates; ++k )
      for ( std::size_t m = 0; m < nstates; ++m )
        probs[k][m] = elnproduct(probs[k][m], -normalizer);
  }

} // namespace hmm

} // namespace ci

#endif // XI_HMM_R_HPP
