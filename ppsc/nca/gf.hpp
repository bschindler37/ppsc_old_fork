
#ifndef _PPSC_NCA_GF_HPP
#define _PPSC_NCA_GF_HPP

// -----------------------------------------------------------------------
//
// NCA single-particle Green's function integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------

#include "ppsc/diag.hpp"

extern template class ppsc::diagram_handler<ppsc::diagram_handler_type::nca_gdh>;

// -----------------------------------------------------------------------
namespace ppsc {
namespace nca {
// -----------------------------------------------------------------------

using namespace ppsc;

typedef diagram_handler<diagram_handler_type::nca_gdh> gf_diagram_handler;

// -----------------------------------------------------------------------
} // end namespace nca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_NCA_GF_HPP
