/*
  FILE: gamma.hpp
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

#ifndef GAMMA_HMM_R_HPP
#define GAMMA_HMM_R_HPP

#include <stack>
#include <vector>

#include "bkd.hpp"
#include "efun.hpp"
#include "fwd.hpp"
#include "infinity.hpp"


namespace ci {

namespace hmm {

  //=========
  // gamma_t_full()
  //  : Evaluate the probability of q_t being in state i given an observation
  //    sequence and model
  //  : Inefficient in time due to backward_index() calls; memory is good
  //  : Calculates all gam values (nstates * nobservations)
  //=========
  template <typename O, typename I, typename T, typename E, typename U>
  void gamma_t_full(const O& observed,
                    const I& initial,
                    const T& transition,
                    const E& emission,
                    std::vector < std::vector<U> >& gam) {

    std::size_t nstates = initial.size(), nobs = observed.size();
    std::vector<U> alpha(nstates, 0), beta(nstates, 0);

    U normalizer = inf<U>();
    for ( std::size_t s = 0; s < nobs; ++s ) {
      forward_next(observed, initial, transition, emission, s+1, alpha);
      backward_index(observed, initial, transition, emission, s+1, beta);

      normalizer = inf<U>();
      for ( std::size_t i = 0; i < nstates; ++i ) {
        gam[i][s] = elnproduct(alpha[i], beta[i]);
        normalizer = elnsum(normalizer, gam[i][s]);
      } // for

      for ( std::size_t j = 0; j < nstates; ++j )
        gam[j][s] = elnproduct(gam[j][s], -normalizer);
    } // for
  }

  //=========
  // gamma_m_full()
  //  : Evaluate the probability of q_t being in state i given an observation
  //    sequence and model
  //  : Inefficient in memory
  //  : Calculates all gam values (nstates * nobservations)
  //=========
  template <typename O, typename I, typename T, typename E, typename U>
  void gamma_m_full(const O& observed,
                    const I& initial,
                    const T& transition,
                    const E& emission,
                    std::vector < std::vector<U> >& gam) {

    std::size_t nstates = initial.size(), nobserved = observed.size();
    std::vector< std::vector<U> > alpha(nstates), beta(nstates);
    for ( std::size_t i = 0; i < nstates; ++i )
      alpha[i].resize(nobserved, 0), beta[i].resize(nobserved, 0);

    forward_full(observed, initial, transition, emission, nobserved, alpha);
    backward_full(observed, initial, transition, emission, 1, beta);

    U normalizer = inf<U>();
    for ( std::size_t s = 0; s < nobserved; ++s ) {
      normalizer = inf<U>();
      for ( std::size_t i = 0; i < nstates; ++i ) {
        gam[i][s] = elnproduct(alpha[i][s], beta[i][s]);
        normalizer = elnsum(normalizer, gam[i][s]);
      } // for

      for ( std::size_t j = 0; j < nstates; ++j )
        gam[j][s] = elnproduct(gam[j][s], -normalizer);
    } // for
  }

  //=========
  // gamma()
  //  : Evaluate the probability of q_t being in state i given an observation
  //    sequence and model
  //  : Calculates one time ('index') slice for gam (nstates * 1)
  //=========
  template <typename O, typename I, typename T, typename E, typename U>
  void gamma(const O& observed,
             const I& initial,
             const T& transition,
             const E& emission,
             std::size_t index,
             const std::vector<U>& beta,
             std::vector<U>& alpha,
             std::vector<U>& gam) {

    std::size_t nstates = initial.size();

    forward_next(observed, initial, transition, emission, index, alpha);

    U normalizer = inf<U>();
    for ( std::size_t i = 0; i < nstates; ++i ) {
      gam[i] = elnproduct(alpha[i], beta[i]);
      normalizer = elnsum(normalizer, gam[i]);
    } // for

    for ( std::size_t j = 0; j < nstates; ++j )
      gam[j] = elnproduct(gam[j], -normalizer);
  }

} // namespace hmm

} // namespace ci


#endif // GAMMA_HMM_R_HPP
