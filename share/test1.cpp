/*
  FILE: test1.cpp
  AUTHOR: Shane Neph
  CREATE DATE: Sat Oct 12 22:33:04 PDT 2013
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "hmm.hpp"

namespace {

  template <typename T>
  void make_vector(std::vector<T>& vec, const std::string& s) {
    for ( std::string::size_type i = 0; i < s.size(); ++i )
      vec.push_back(static_cast<int>(s[i]) - static_cast<int>('0'));
  }


  template <typename I, typename T, typename E>
  void dolog(I& initial, T& transition, E& emission) {

    typedef typename I::value_type U;
    for ( std::size_t i = 0; i < initial.size(); ++i ) {
      if ( 0 != initial[i] )
        initial[i] = std::log(initial[i]);
      else
        initial[i] = ci::inf<U>();
    } // for

    for ( std::size_t i = 0; i < transition.size(); ++i ) {
      for ( std::size_t j = 0; j < transition[i].size(); ++j ) {
        if ( 0 != transition[i][j] )
          transition[i][j] = std::log(transition[i][j]);
        else
          transition[i][j] = ci::inf<U>();
      } // for
    } // for

    for ( std::size_t i = 0; i < emission.size(); ++i ) {
      for ( std::size_t j = 0; j < emission[i].size(); ++j ) {
        if ( 0 != emission[i][j] )
          emission[i][j] = std::log(emission[i][j]);
        else
          emission[i][j] = ci::inf<U>();
      } // for
    } // for
  }

} // empty namespace


int main()
{
  // Observed sequence 
  typedef float T;
  std::vector<T> observed, initial, tmp;
  std::vector< std::vector<T> > transition, emission, results;
  std::string obs = "010000000010000100001000000000";
  make_vector(observed, obs);
  std::copy(observed.begin(), observed.end(), std::ostream_iterator<T>(std::cout, " "));
  std::cout << std::endl << "Length = " << observed.size() << std::endl;

  // Initial state probabilities
  initial.push_back(0.5), initial.push_back(0.5);

  // Transitional probabilities
  tmp.clear();
  tmp.push_back(0.9), tmp.push_back(0.1);
  transition.push_back(tmp);

  tmp.clear();
  tmp.push_back(0.5), tmp.push_back(0.5);
  transition.push_back(tmp);

  // Emission probabilities
  tmp.clear();
  tmp.push_back(0.2), tmp.push_back(0.3), tmp.push_back(0.5);
  emission.push_back(tmp);

  tmp.clear();
  tmp.push_back(0.5), tmp.push_back(0.2), tmp.push_back(0.3);
  emission.push_back(tmp);

  dolog(initial, transition, emission);

  // Test full forward
  std::cout << "Full Forward" << std::endl;
  std::size_t index = observed.size();
  results.resize(initial.size());
  for ( std::size_t i = 0; i < results.size(); ++i )
    results[i].resize(index);

  ci::hmm::forward_full(observed, initial, transition, emission, index, results);
  for ( std::size_t i = 0; i < results.size(); ++i ) {
    for ( std::size_t j = 0; j < results[i].size(); ++j )
      std::cout << results[i][j] << "\t";
    std::cout << std::endl;
  } // for

  // Test indexed forward
  std::cout << "Indexed Forward" << std::endl;
  std::vector<T> ok(initial.size(), 0);
  ci::hmm::forward_index(observed, initial, transition, emission, 1, ok);
  std::copy(ok.begin(), ok.end(), std::ostream_iterator<T>(std::cout, "\t"));
  std::cout << std::endl;

  // Test next forward
  ok.clear();
  ok.resize(initial.size(), 0);
  for ( std::size_t i = 1; i <= observed.size(); ++i ) {
    std::cout << "Next Forward (" << i << ")" << std::endl;
    ci::hmm::forward_next(observed, initial, transition, emission, i, ok);
    std::copy(ok.begin(), ok.end(), std::ostream_iterator<T>(std::cout, "\t"));
    std::cout << std::endl;
  }

  // Test backward
  std::cout << "Full Backward" << std::endl;
  results.clear();
  results.resize(initial.size());
  for ( std::size_t i = 0; i < results.size(); ++i )
    results[i].resize(observed.size(), 0);
  index = 1;

  ci::hmm::backward_full(observed, initial, transition, emission, index, results);
  for ( std::size_t i = 0; i < results.size(); ++i ) {
    for ( std::size_t j = 0; j < results[i].size(); ++j )
      std::cout << results[i][j] << "\t";
    std::cout << std::endl;
  } // for

  // Test indexed backward
  std::cout << "Indexed Backward" << std::endl;
  ok.clear();
  ok.resize(initial.size());
  ci::hmm::backward_index(observed, initial, transition, emission, index, ok);
  std::copy(ok.begin(), ok.end(), std::ostream_iterator<T>(std::cout, "\t"));
  std::cout << std::endl;

  // Test next backward
  ok.clear();
  ok.resize(initial.size(), 0);
  for ( std::size_t i = observed.size(); i > 0; --i ) {
    std::cout << "Next Backward (" << i << ")" << std::endl;
    ci::hmm::backward_next(observed, initial, transition, emission, i, ok);
    std::copy(ok.begin(), ok.end(), std::ostream_iterator<T>(std::cout, "\t"));
    std::cout << std::endl;
  }

  typedef std::vector<T> FOO;
  FOO const* foo;
  ci::hmm::details::BackCache< FOO, FOO, std::vector<FOO>, std::vector<FOO>, T > backcheaterA(observed, initial, transition, emission);
  for ( std::size_t i = 0; i < observed.size(); ++i ) {
    std::cout << "Cheat Backward->Forward (" << i << ")" << std::endl;
    foo = backcheaterA.Next();
    std::copy(foo->begin(), foo->end(), std::ostream_iterator<T>(std::cout, "\t"));
    std::cout << std::endl;
  }

  // Test extended next backward
  std::vector< std::vector<T> > dummy(initial.size());
  for ( std::size_t i = 0; i < dummy.size(); ++i )
    dummy[i] = std::vector<T>(initial.size(), 0);

  for ( std::size_t i = observed.size(); i > 0; --i ) {
    std::cout << "(Ext) Next Backward (" << i << ")" << std::endl;
    ci::hmm::backward_enext(observed, initial, transition, emission, i, dummy);
    for ( std::size_t j = 0; j < dummy.size(); ++j ) {
      for ( std::size_t k = 0; k < dummy[j].size(); ++k )
        std::cout << dummy[k][j] << "\t";
     std::cout << std::endl;
    } // for
    std::cout << std::endl;
  } // for


  // Problem 1
  float ans1 = ci::hmm::evalp(observed, initial, transition, emission);
  std::cout << "Answer to problem 1: " << ans1 << std::endl;

  // Problem 2
  // Test viterbi
  std::cout << "Viterbi answer to problem 2" << std::endl;
  std::ostream_iterator<T> os(std::cout, "\t");
  ci::hmm::viterbi(observed, initial, transition, emission, os);
  std::cout << std::endl;

  // Test gamma
  std::cout << "Testing Gamma" << std::endl;
  results.clear();
  results.resize(initial.size());
  for ( std::size_t i = 0; i < results.size(); ++i )
    results[i].resize(observed.size(), 0);
  ci::hmm::gamma_m_full(observed, initial, transition, emission, results);
  for ( std::size_t i = 0; i < results.size(); ++i ) {
    for ( std::size_t j = 0; j < results[i].size(); ++j )
      std::cout << results[i][j] << "\t";
    std::cout << std::endl;
  } // for

  std::cout << "New Gamma" << std::endl;
  std::vector<T> alpha(initial.size(), 0), beta(alpha.size(), 0), gam(alpha.size(), 0);
  ci::hmm::details::BackCache< FOO, FOO, std::vector<FOO>, std::vector<FOO>, T > backcheaterB(observed, initial, transition, emission);
  for ( std::size_t i = 1; i <= observed.size(); ++i ) {
    std::vector<T> const* foo = backcheaterB.Next();
    ci::hmm::gamma(observed, initial, transition, emission, i, *foo, alpha, gam);
    delete foo;
    std::copy(gam.begin(), gam.end(), std::ostream_iterator<T>(std::cout, "\t"));
    std::cout << std::endl;
  } // for

  // Test xi
  std::cout << "Testing xi" << std::endl;

  std::vector< std::vector< std::vector<T> > > v(initial.size());
  for ( std::size_t i = 0; i < v.size(); ++i ) {
    v[i].resize(initial.size());
    for ( std::size_t j = 0; j < v[i].size(); ++j )
      v[i][j].resize(observed.size(), 0); // should be observed.size()-1, but need to compare to old outputs
  }
  ci::hmm::xi_full(observed, initial, transition, emission, v);
  for ( std::size_t i = 0; i < v.size(); ++i ) {
    for ( std::size_t j = 0; j < v[i].size(); ++j ) {
      for ( std::size_t k = 0; k < v[i][j].size(); ++k )
        std::cout << v[i][j][k] << "\t";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  } // for

  std::cout << "YO NEW-XI" << std::endl;
  std::vector< std::vector<T> > v2(initial.size());
  for ( std::size_t i = 0; i < v2.size(); ++i ) {
    v2[i].resize(initial.size());
  }

  ci::hmm::details::BackCache< std::vector<T>, std::vector<T>, std::vector< std::vector<T> >, std::vector< std::vector<T> >, T > backcheater(observed, initial, transition, emission);
  alpha.resize(alpha.size(), 0), beta.resize(beta.size(), 0);
  foo = backcheater.Next(); // must move up one spot in backward->forward direction
  for ( std::size_t s = 1; s < observed.size(); ++s ) {
    std::vector<T> const* foo = backcheater.Next();
    ci::hmm::xi(observed, initial, transition, emission, s, *foo, alpha, v2);
    delete foo;
    for ( std::size_t i = 0; i < v2.size(); ++i ) {
      for ( std::size_t j = 0; j < v2[i].size(); ++j )
        std::cout << v2[i][j] << "\t";
      std::cout << std::endl;
    } // for
  } // for
  std::cout << std::endl;

  // Problem 3
  // Test training
  std::cout << "Problem 3" << std::endl;
  std::cout << "Start Initial" << std::endl;
  for ( std::size_t i = 0; i < initial.size(); ++i )
    std::cout << initial[i] << std::endl;

  std::vector<T> keepobserved(observed), keepinitial(initial);
  std::vector< std::vector<T> > keeptransition(transition), keepemission(emission);

  std::size_t numiter = 2;
  for ( std::size_t i = 0; i < numiter; ++i ) {
    ci::hmm::train_full(observed, initial, transition, emission);
//    dolog(initial, transition, emission);

    std::cout << "Iteration: " << (i+1) << std::endl;

    std::cout << "New Initial" << std::endl;
    for ( std::size_t i = 0; i < initial.size(); ++i )
      std::cout << initial[i] << std::endl;

    std::cout << "New Transition" << std::endl;
    for ( std::size_t i = 0; i < transition.size(); ++i ) {
      for ( std::size_t j = 0; j < transition[i].size(); ++j )
        std::cout << transition[i][j] << "\t";
      std::cout << std::endl;
    } // for

    std::cout << "New Emission" << std::endl;
    for ( std::size_t i = 0; i < emission.size(); ++i ) {
      for ( std::size_t j = 0; j < emission[i].size(); ++j ) {
        if ( emission[i][j] == ci::inf<T>() )
          emission[i][j] = 0;
        std::cout << emission[i][j] << "\t";
      } // for
      std::cout << std::endl;
    } // for
  } // for


  std::cout << "Finale: *********************" << std::endl;

  observed = keepobserved;
  initial = keepinitial;
  transition = keeptransition;
  emission = keepemission;

  std::cout << "New Training" << std::endl;
  for ( std::size_t i = 0; i < numiter; ++i ) {
    ci::hmm::train(observed, initial, transition, emission);

    std::cout << "Iteration " << (i+1) << std::endl;

    std::cout << "New Initial" << std::endl;
    std::cout << initial.size() << std::endl;
    std::cout.flush();
    for ( std::size_t i = 0; i < initial.size(); ++i )
      std::cout << initial[i] << std::endl;
    std::cout.flush();

    std::cout << "New Transition" << std::endl;
    for ( std::size_t i = 0; i < transition.size(); ++i ) {
      for ( std::size_t j = 0; j < transition[i].size(); ++j )
        std::cout << transition[i][j] << "\t";
      std::cout << std::endl;
    } // for

    std::cout << "New Emission" << std::endl;
    for ( std::size_t i = 0; i < emission.size(); ++i ) {
      for ( std::size_t j = 0; j < emission[i].size(); ++j ) {
        if ( emission[i][j] == ci::inf<T>() )
          emission[i][j] = 0;
        std::cout << emission[i][j] << "\t";
      } // for
      std::cout << std::endl;
    } // for
  } // for

  return(0);
}
