/*
  FILE: backcache.hpp
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


#ifndef BKD_CACHE_DETAILS_HMM_R_HPP
#define BKD_CACHE_DETAILS_HMM_R_HPP

#include <algorithm>
#include <cmath>
#include <list>
#include <vector>


namespace ci {

namespace hmm {

namespace details {

  //=============
  // BackCache<>
  //   : Reverse iterate through all observations caching at intervals
  //       with main goal of allowing efficient fwd(-bkd) iteration later
  //   : Basically, BackCache<> provides an interface (via Next()) to
  //       iterate through all betas in backward algorithm in reverse
  //       order (reverse-reverse == forward)
  //   : If fully-connected model with all nonzero transitions and all
  //       appropriate matrices are non-singular, massive memory savings
  //       (and hence, time savings for large inputs) can be had
  //
  //       These conditions are too constraining so this class is an
  //       attempt to retain some of those savings within the general case
  //
  //     Try to keep some reasonable number of items in memory at any given
  //       time without needing to traverse all observations more than twice.
  template <typename O, typename I, typename T, typename E, typename U>
  struct BackCache {

    //=============
    // Constructor
    BackCache(const O& o, const I& i, const T& t, const E& e)
                           : initialize_(true),
                             sz_(std::max(static_cast<std::size_t>(10000),
                                          static_cast<std::size_t>(std::sqrt(o.size())))),
                             observed_(o), initial_(i),
                             transition_(t), emission_(e)
      {  populate(); } // do full backward traversal

    //============
    // Destructor
    ~BackCache() {
      typedef typename std::list< std::vector<U>* >::iterator CType;
      for ( CType c = passiveItems_.begin(); c != passiveItems_.end(); ++c )
        delete *c;
      for ( CType c = activeItems_.begin(); c != activeItems_.end(); ++c )
        delete *c;
    }

    //=========
    // Next()
    //  : User must manage memory obtained
    inline std::vector<U> const* Next() {
      if ( !activeItems_.empty() ) {
        std::vector<U>* rtn = activeItems_.front();
        activeItems_.pop_front();
        return(rtn);
      }

      if ( passiveItems_.empty() ) // exhausted
        return(static_cast< std::vector<U>* >(0));
      populate();
      return(Next());
    }

    //========
    // Size()
    std::size_t Size() const
      { return(activeItems_.size() + passiveItems_.size()); }

    //==================
    // Copy Construction
    BackCache(const BackCache& b)
              : markers_(b.markers_), counters_(b.counters_),
                initialize_(b.initialize_), sz_(b.sz_), observed_(b.observed_),
                initial_(b.initial_), transition_(b.transition_),
                emission_(b.emission_) {

      typedef typename std::list< std::vector<U>* >::const_iterator CType;
      for ( CType c = b.passiveItems_.begin(); c != b.passiveItems_.end(); ++c )
        passiveItems_.push_back(new std::vector<U>(**c));
      for ( CType c = b.activeItems_.begin(); c != b.activeItems_.end(); ++c )
        activeItems_.push_back(new std::vector<U>(**c));
    }

    //===============
    // No Assignment
    void operator=(const BackCache& b); // disabled purposefully for now
 
  private:
    void populate() {
      typedef std::vector<U> V;

      typedef typename std::list< std::vector<U>* >::iterator CType;
      for ( CType c = activeItems_.begin(); c != activeItems_.end(); ++c )
        delete *c;

      activeItems_.clear();
      if ( sz_ <= 1 || observed_.empty() ) {
        for ( CType c = passiveItems_.begin(); c != passiveItems_.end(); ++c )
          delete *c;
        passiveItems_.clear();
        return;
      }

      if ( initialize_ ) { // traverse all observations and cache chosen results
        initialize_ = false;
        V beta(initial_.size(), 0);
        bool lastleg = (observed_.size() <= sz_);
        if ( !lastleg ) {
          passiveItems_.push_front(new V(beta));
          markers_.push_front(observed_.size());
          counters_.push_front(sz_);
        }

        for ( std::size_t i = observed_.size()-1, j = 1; i > 0; --i, ++j ) {
          if ( lastleg )
            activeItems_.push_front(new V(beta));
          else if ( i == sz_ ) {
            activeItems_.push_front(new V(beta));
            counters_.pop_front();
            counters_.push_front(j);
            lastleg = true;
          }
          else if ( j == sz_ ) {
            passiveItems_.push_front(new V(beta));
            markers_.push_front(i);
            counters_.push_front(sz_);
            j = 0;
          }
          backward_next(observed_, initial_, transition_, emission_, i, beta);
        } // for
        activeItems_.push_front(new V(beta));
        return;
      }
      else if ( passiveItems_.empty() ) // nothing left to do
        return;

      // pop off next item from passiveItems_ and re-populate activeItems_
      V beta = *passiveItems_.front();
      activeItems_.push_front(passiveItems_.front());
      passiveItems_.pop_front();
      std::size_t mark = markers_.front();
      markers_.pop_front();
      std::size_t count = counters_.front();
      counters_.pop_front();
      for ( std::size_t s = mark, i = count; i > 1; --i, --s ) {
        backward_next(observed_, initial_, transition_, emission_, s, beta);
        activeItems_.push_front(new V(beta));
      } // for
    }

  private:
    std::list< std::vector<U>* > passiveItems_;
    std::list< std::vector<U>* > activeItems_;
    std::list< std::size_t > markers_, counters_;
    bool initialize_;
    const std::size_t sz_;
    const O& observed_;
    const I& initial_;
    const T& transition_;
    const E& emission_;
  };

} // namespace details

} // namespace hmm

} // namespace ci

#endif // BKD_CACHE_DETAILS_HMM_R_HPP
