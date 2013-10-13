/*
  FILE: evalp.hpp
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

#ifndef EVALP_HMM_R_HPP
#define EVALP_HMM_R_HPP

#include <vector>

#include "efun.hpp"
#include "fwd.hpp"
#include "infinity.hpp"

namespace ci {

namespace hmm {

  //=====================
  // evalp() algorithm
  //   - Wraps hmm::forward() to solve "Problem 1"
  //=====================
  template <typename O, typename I, typename T, typename E>
  float evalp(const O& observed,
              const I& initial,
              const T& transition,
              const E& emission) {

    std::size_t tsize = observed.size();
    if ( tsize < 2 )
      return(inf<float>());

    const std::size_t nstates = initial.size();
    std::vector<float> alpha(nstates);

    forward_index(observed, initial,transition, emission, tsize, alpha);
    float enlp = inf<float>();
    for ( std::size_t i = 0; i < alpha.size(); ++i )
      enlp = elnsum(enlp, alpha[i]);
    return(std::exp(enlp));
  }

} // namespace hmm

} // namespace ci

#endif // EVALP_HMM_R_HPP
