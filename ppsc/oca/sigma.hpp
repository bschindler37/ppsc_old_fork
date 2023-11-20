
#ifndef _PPSC_OCA_SIGMA_HPP
#define _PPSC_OCA_SIGMA_HPP

// -----------------------------------------------------------------------
//
// OCA pseudo particle self energy integrators
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------

#include "ppsc/diag.hpp"

extern template class ppsc::diagram_handler<ppsc::diagram_handler_type::oca_sdh>;

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;

typedef diagram_handler<diagram_handler_type::oca_sdh> sigma_diagram_handler;

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_HPP
