
#ifndef _PPSC_OCA_GF_DISP_HPP
#define _PPSC_OCA_GF_DISP_HPP

// -----------------------------------------------------------------------
//
// OCA single-particle Green's function integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "data_types.hpp"
#include "gf_diag.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -- Dispatch on equilibrium / real-time

void gdiagram_dispatch(int tstp, gf_tstp_type & gtstp, gdiagram_configuration<> & diagram,
		       double beta, double h, int kt, int nomp);
  
// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_DISP_HPP
