
#ifndef _PPSC_RANGE_HPP
#define _PPSC_RANGE_HPP

// -----------------------------------------------------------------------
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include <boost/range/irange.hpp>

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

typedef std::size_t size_type;
typedef std::ptrdiff_t ptrdiff_type;
  
boost::iterator_range< boost::range_detail::integer_iterator<size_type> > range(size_type first, size_type last);

boost::iterator_range< boost::range_detail::integer_iterator_with_step<ptrdiff_type> >
range(ptrdiff_type first, ptrdiff_type last, ptrdiff_type step);
  
// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_RANGE_HPP
