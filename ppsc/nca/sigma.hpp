
#ifndef _PPSC_NCA_SIGMA_HPP
#define _PPSC_NCA_SIGMA_HPP

// -----------------------------------------------------------------------
//
// NCA pseudo particle self energy calculator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------

#include "ppsc/diag.hpp"

extern template class ppsc::diagram_handler<ppsc::diagram_handler_type::nca_sdh>;

// -----------------------------------------------------------------------
namespace ppsc {
namespace nca {
// -----------------------------------------------------------------------

using namespace ppsc;

typedef diagram_handler<diagram_handler_type::nca_sdh> sigma_diagram_handler;

// -----------------------------------------------------------------------
} // end namespace nca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_NCA_SIGMA_HPP
