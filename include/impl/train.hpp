/*
  FILE: train.hpp
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

#ifndef TRAIN_HMM_R_HPP
#define TRAIN_HMM_R_HPP

#include <cmath>
#include <vector>

#include "backcache.hpp"
#include "bkd.hpp"
#include "efun.hpp"
#include "gamma.hpp"
#include "infinity.hpp"
#include "xi.hpp"

namespace ci {

namespace hmm {


  /*
    ----------
    Interfaces
    ----------
    train_full() :
      Closest implementation to Rabiner's pseudo-code
      Unforgiving in RAM
      Can be fastest if number of observations, emission symbols
        and model states is small

    train() :
      Scales the best in time with a large number of observations,
        but takes some liberties in RAM with the number of states
        and emission symbols
      Should be the best implementation for most discrete cases

    train_mem() :
      Similar to train(), but works to reduce RAM requirements
        to a minimum.  There is a trade-off with how many
        times we must sweep through the number of observations


    ---------
    All functions take the same arguments.  Here is an example:

    template <typename O, typename I, typename T, typename E>
    void train(const O& observed,
               I& initial,
               T& transition,
               E& emission);
  */




  //=============
  // train_full()
  //   : Re-estimate model parameters
  //   : Closest to Rabiner's pseudo-code
  //   : Inefficient in memory
  template <typename O, typename I, typename T, typename E>
  void train_full(const O& observed,
                  I& initial,
                  T& transition,
                  E& emission) {

    typedef float U;
    const std::size_t nstates = initial.size();
    const std::size_t nobs = observed.size();
    const std::size_t nsymbols = emission[0].size();

    std::vector< std::vector<U> > gam(initial.size());
    for ( std::size_t i = 0; i < gam.size(); ++i )
      gam[i].resize(observed.size(), 0);
    gamma_m_full(observed, initial, transition, emission, gam);

    std::vector< std::vector< std::vector<U> > > probs(nstates);
    for ( std::size_t i = 0; i < probs.size(); ++i ) {
      probs[i].resize(nstates);
      for ( std::size_t j = 0; j < probs[i].size(); ++j )
        probs[i][j].resize(nobs, 0);
    } // for
    xi_full(observed, initial, transition, emission, probs);

    // update initial
    for ( std::size_t i = 0; i < gam.size(); ++i )
      initial[i] = std::exp(gam[i][0]);

    // update emission and transition
    U numeratorE = inf<U>(), denominatorE = inf<U>();
    U numeratorT = inf<U>(), denominatorT = inf<U>();
    std::size_t sentinel = std::max(nsymbols, nstates);
    for ( std::size_t i = 0; i < sentinel; ++i ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        numeratorE = inf<U>(), denominatorE = inf<U>();
        numeratorT = inf<U>(), denominatorT = inf<U>();
        for ( std::size_t s = 0; s < nobs-1; ++s ) {
          if ( i < nsymbols ) { // emission
            if ( observed[s] == i )
              numeratorE = elnsum(numeratorE, gam[j][s]);
            denominatorE = elnsum(denominatorE, gam[j][s]);
          }

          if ( i < nstates ) { // transition
            numeratorT = elnsum(numeratorT, probs[i][j][s]);
            denominatorT = elnsum(denominatorT, gam[i][s]);
          }
        } // for

        if ( i < nsymbols ) // emission
          emission[j][i] = elnproduct(numeratorE, -denominatorE);

        if ( i < nstates ) // transition
          transition[i][j] = elnproduct(numeratorT, -denominatorT);
      } // for
    } // for
  }

  //=========
  // train()
  //   : Re-estimate model parameters
  //   : Efficient in time & memory with regards to num observations
  //   : Best implementation for most discrete models
  template <typename O, typename I, typename T, typename E>
  void train(const O& observed,
             I& initial,
             T& transition,
             E& emission) {

    // Various constants
    typedef typename O::value_type U;
    const I& init(initial);
    const std::size_t nstates = initial.size();
    const std::size_t nobs = observed.size();
    const std::size_t nsymbols = emission[0].size();
    const std::size_t sentinel = std::max(nsymbols, nstates);
    const bool done = false;

    // Various local arrays declared
    std::vector<U> gam(nstates, 0), alphaG(gam), alphaX(gam);
    std::vector< std::vector<U> > probs(nstates);
    std::vector< std::vector<U> > numeratorT(nstates), denominatorT(nstates);
    std::vector< std::vector<U> > numeratorE(nsymbols), denominatorE(nsymbols);

    // Initialize values of local arrays
    for ( std::size_t i = 0; i < sentinel; ++i ) {
      if ( i < nstates ) {
        probs[i].resize(nstates, 0);
        numeratorT[i].resize(nstates, inf<U>());
        denominatorT[i].resize(nstates, inf<U>());
      }

      if ( i < nsymbols ) {
        numeratorE[i].resize(nstates, inf<U>());
        denominatorE[i].resize(nstates, inf<U>());
      }
    } // for

    // Prepare for back propogations
    typedef details::BackCache<O, I, T, E, U> BCache;
    BCache cache(observed, init, transition, emission);
    std::vector<U> const* beta = cache.Next();

    // First gamma() call
    gamma(observed, init, transition, emission, 1, *beta, alphaG, gam);
    if ( !beta )
      return;
    delete beta;

    beta = cache.Next(); // xi's beta stays ahead of gamma's by one
    if ( !beta )
      return;

    // First xi() call
    xi(observed, init, transition, emission, 1, *beta, alphaX, probs);

    // Update new initial state probabilities
    for ( std::size_t y = 0; y < nstates; ++y )
      initial[y] = std::exp(gam[y]);

    // Calculate intermediaries for emission and transition probabilities
    std::size_t s = 0;
    while ( !done ) {
      for ( std::size_t i = 0; i < sentinel; ++i ) {
        for ( std::size_t j = 0; j < nstates; ++j ) {
          if ( i < nsymbols ) { // emission
            if ( observed[s] == i )
              numeratorE[i][j] = elnsum(numeratorE[i][j], gam[j]);
            denominatorE[i][j] = elnsum(denominatorE[i][j], gam[j]);
          }

          if ( i < nstates ) { // transition
            numeratorT[i][j] = elnsum(numeratorT[i][j], probs[i][j]);
            denominatorT[i][j] = elnsum(denominatorT[i][j], gam[i]);
          }
        } // for 'j'
      } // for 'i'

      if ( ++s == nobs-1 )
        break;
      gamma(observed, init, transition, emission, s+1, *beta, alphaG, gam);
      delete beta;

      beta = cache.Next(); // xi's beta stays ahead of gamma's by one
      xi(observed, init, transition, emission, s+1, *beta, alphaX, probs);
    } // while !done

    if ( beta )
      delete beta;

    // Update new emission and transitional probabilities
    for ( std::size_t i = 0; i < sentinel; ++i ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        if ( i < nsymbols ) // emission
          emission[j][i] = elnproduct(numeratorE[i][j], -denominatorE[i][j]);

        if ( i < nstates ) // transition
          transition[i][j] = elnproduct(numeratorT[i][j], -denominatorT[i][j]);
      } // for 'j'
    } // for 'i'
  }

  //=============
  // train_mem()
  //   : Re-estimate model parameters
  //   : Most efficient in memory
  template <typename O, typename I, typename T, typename E>
  void train_mem(const O& observed,
                 I& initial,
                 T& transition,
                 E& emission) {

    typedef typename O::value_type U;
    const std::size_t nstates = initial.size();
    const std::size_t nobs = observed.size();
    const std::size_t nsymbols = emission[0].size();

    // Involatile copies of originals necessary
    const I init(initial);
    const T trans(transition);
    const E emis(emission);

    typedef details::BackCache<O, I, T, E, U> BCache;
    const BCache cache(observed, init, trans, emis);

    std::vector<U> gam(nstates, 0), alphaG(gam);
    std::vector< std::vector<U> > probs(nstates);

    for ( std::size_t i = 0; i < nstates; ++i )
      probs[i].resize(nstates, 0);

    // update emission and transition probabilities; update initial conditions
    U numeratorE = inf<U>(), denominatorE = inf<U>();
    U numeratorT = inf<U>(), denominatorT = inf<U>();
    std::size_t sentinel = std::max(nsymbols, nstates);
    for ( std::size_t i = 0; i < sentinel; ++i ) {
      for ( std::size_t j = 0; j < nstates; ++j ) {
        numeratorE = inf<U>(), denominatorE = inf<U>();
        numeratorT = inf<U>(), denominatorT = inf<U>();
        BCache gcache(cache);
        std::vector<U> alphaG(nstates, 0);

        // if-else is highly duplicated but worth it due to extra data copies for transition
        if ( i < nstates ) {
          BCache xcache(cache);
          std::vector<U> alphaX(nstates, 0);
          std::vector<U>* beta = xcache.Next(); // xi -> must move forward one
          if ( !beta )
            return;
          delete beta;

          for ( std::size_t s = 0; s < nobs-1; ++s ) {
            gamma(observed, init, trans, emis, s+1, gcache, alphaG, gam);

            if ( 0 == i && 0 == s && 0 == j ) { // update initial; must be the case: i < nstates
              for ( std::size_t y = 0; y < nstates; ++y )
                initial[y] = std::exp(gam[y]);
            }

            if ( i < nsymbols ) { // emission
              if ( observed[s] == i )
                numeratorE = elnsum(numeratorE, gam[j]);
              denominatorE = elnsum(denominatorE, gam[j]);
            }

            // transition
            xi(observed, init, trans, emis, s+1, xcache, alphaX, probs);
            numeratorT = elnsum(numeratorT, probs[i][j]);
            denominatorT = elnsum(denominatorT, gam[i]);
          } // for

          if ( i < nsymbols ) // emission
            emission[j][i] = elnproduct(numeratorE, -denominatorE);

          // transition
          transition[i][j] = elnproduct(numeratorT, -denominatorT);
        }
        else { // no need to make all copies/checks for transition-related items
          for ( std::size_t s = 0; s < nobs-1; ++s ) {
            gamma(observed, init, trans, emis, s+1, gcache, alphaG, gam);

            if ( i < nsymbols ) { // emission
              if ( observed[s] == i )
                numeratorE = elnsum(numeratorE, gam[j]);
              denominatorE = elnsum(denominatorE, gam[j]);
            }
          } // for

          if ( i < nsymbols ) // emission
            emission[j][i] = elnproduct(numeratorE, -denominatorE);
        }

      } // for
    } // for
  }

} // namespace hmm

} // namespace ci


#endif // TRAIN_HMM_R_HPP
