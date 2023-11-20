
#ifndef _PPSC_CNTR_TYPES_HPP
#define _PPSC_CNTR_TYPES_HPP

// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling expansion cntr types
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include <vector>

#define NO_IMPLEMENTATION
#include "cntr/cntr.hpp"
#undef NO_IMPLEMENTATION

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

typedef std::complex<double> value_type;
  
typedef cntr::herm_matrix<double> gf_type;
typedef cntr::herm_pseudo<double> ppgf_type;
typedef cntr::herm_matrix_timestep<double> gf_tstp_type;
typedef cntr::function<double> function_type;
  
typedef std::vector<gf_type> gfs_type;
typedef std::vector<ppgf_type> ppgfs_type;
typedef std::vector<gf_tstp_type> gf_tstps_type;
typedef std::vector<function_type> functions_type;
  
// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_CNTR_TYPES_HPP
