
#ifndef _PPSC_OCA_SIGMA_DISP_HPP
#define _PPSC_OCA_SIGMA_DISP_HPP

// -----------------------------------------------------------------------
//
// OCA pseudo particle self energy integrators
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/data_types.hpp"
#include "sigma_diag.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

void sdiagram_dispatch(int tstp, gf_tstp_type & ststp, sdiagram_configuration<> & diagram,
		       double beta, double h, int kt, int nomp);

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_DISP_HPP
