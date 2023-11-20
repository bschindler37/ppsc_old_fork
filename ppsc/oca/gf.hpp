
#ifndef _PPSC_OCA_GF_HPP
#define _PPSC_OCA_GF_HPP

// -----------------------------------------------------------------------
//
// OCA single-particle Green's function integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------

#include "ppsc/diag.hpp"

extern template class ppsc::diagram_handler<ppsc::diagram_handler_type::oca_gdh>;

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;

typedef diagram_handler<diagram_handler_type::oca_gdh> gf_diagram_handler;

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_HPP
