
#ifndef _PPSC_RANGE_CPP
#define _PPSC_RANGE_CPP

// -----------------------------------------------------------------------
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include <boost/range/irange.hpp>

#include "ppsc/range.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------
typedef std::size_t size_type;
typedef std::ptrdiff_t ptrdiff_type;
  
boost::iterator_range< boost::range_detail::integer_iterator<size_type> > range(size_type first, size_type last) {
  return boost::irange<size_type>(first, last);
}

boost::iterator_range< boost::range_detail::integer_iterator_with_step<ptrdiff_type> >
  range(ptrdiff_type first, ptrdiff_type last, ptrdiff_type step) {
  return boost::irange<ptrdiff_type, ptrdiff_type>(first, last, step);
}  
  
// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_RANGE_CPP
