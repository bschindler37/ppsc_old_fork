
#ifndef _PPSC_OCA_SIGMA_CPP
#define _PPSC_OCA_SIGMA_CPP

// -----------------------------------------------------------------------
//
// OCA pseudo particle self energy integrators
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/diag_base.hpp"
#include "ppsc/diag_impl.hpp"

#include "sigma_diag.hpp"
#include "sigma_disp.hpp"

#include "sigma.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;

// ----------------------------------------------------------------------
  void get_sectors(int sigma_sector, pp_int_type & int1, pp_int_type & int2,
		   int & sa_out, int & sb_out, int & sc_out,
		   bool & valid_diagram) {

  valid_diagram = false;
  int ss = sigma_sector;
  int sc =     int2.op2.to_sector_[ss]; if(sc == -1) return; sc_out = sc;
  int sb =     int1.op2.to_sector_[sc]; if(sb == -1) return; sb_out = sb;
  int sa =     int2.op1.to_sector_[sb]; if(sa == -1) return; sa_out = sa;
  int ss_ref = int1.op1.to_sector_[sa]; if(ss_ref != ss) return;
  valid_diagram = true;
}

// ----------------------------------------------------------------------
bool valid_sectors(int sigma_sector, pp_int_type & int1, pp_int_type & int2) {
  int sa, sb, sc; bool valid_diagram;
  get_sectors(sigma_sector, int1, int2, sa, sb, sc, valid_diagram);
  return valid_diagram;
}

// -----------------------------------------------------------------------
sdiagram_configuration<> construct_sigma_diagram(
  gf_tstp_type & ststp, int sigma_sig, int sigma_sector,
  ppgfs_type & ppGfs, pp_int_type & int1, pp_int_type & int2) {

  int ss = sigma_sector;
  int sa, sb, sc;
  bool valid_diagram;
  get_sectors(sigma_sector, int1, int2, sa, sb, sc, valid_diagram);

  assert( valid_diagram );

  auto op1 = int1.op1.M_[sa];
  auto op2 = int2.op1.M_[sb];
  auto op3 = int1.op2.M_[sc];
  auto op4 = int2.op2.M_[ss];

  // -- Diagram prefactor [Eckstein, Werner PRB 82, 115155 (2010)]
  const int fermion=-1, boson=+1, fwd=+1, bwd=-1;

  double stat_sign = +1;

  // Check for bwd propagating fermionic interactions
  if(int1.sig == fermion && int1.dir == bwd) stat_sign *= -1;
  if(int2.sig == fermion && int2.dir == bwd) stat_sign *= -1;

  // crossing of two fermionic interaction lines
  if(int1.sig == fermion && int2.sig == fermion) stat_sign *= -1;

  // extra factor not accounted for by diagram rules
  double debug_prefactor = -1.0;

  // -- test for bose hubbard model
  //if(int1.sig == boson && int2.sig == boson) debug_prefactor = 1.0; // DEBUG

  // (i)^{order} factor
  int order = 2;
  value_type I = value_type(0., 1.);

  value_type prefactor = std::pow(I, order) * stat_sign * debug_prefactor;

  //value_type prefactor = -1.0;

  /*
  std::cout << "sigma diagram: ss, sa, sb, sc = "
	    << sigma_sector << " [ "
	    << sa << ", " << sb << ", " << sc << " ] "
	    << sigma_sector
	    << std::endl;
  */

  /*
  std::cout << "i1,i2,sigma sig: "
	    << int1.sig << ", " << int2.sig << ", " << sigma_sig
	    << " debug_prefactor = " << debug_prefactor
	    << " prefactor = " << prefactor
	    << std::endl;
  */

  return sdiagram_configuration<>(
    ststp, sigma_sig,
    op1, op2, op3, op4,
    ppGfs[sa], ppGfs[sb], ppGfs[sc],
    int1.lam, int2.lam,
    ss, prefactor);
}

// -----------------------------------------------------------------------
sdiagrams_type build_all_sigma_diagrams(hilbert_space_type & hilbert_space,
  ppgfs_type & ppGfs, pp_ints_type & pp_ints) {

  sdiagrams_type sdiagrams;

  for(int ss = 0; ss < hilbert_space.ns_; ss++) { // sigma sector
    int sector_dim = hilbert_space.ssdim_[ss];
    int sector_sig = hilbert_space.sig_[ss];
    for(auto int1 : pp_ints) {
      for(auto int2 : pp_ints) {

	if(!valid_sectors(ss, int1, int2)) continue;

	gf_tstp_type sigma_tstp_diagram(0, 0, sector_dim);
	sdiagram_configuration<> diagram = oca::construct_sigma_diagram(
	  sigma_tstp_diagram, sector_sig, ss, ppGfs, int1, int2);

	sdiagrams.push_back(diagram);
      } // int 2
    } // int1
  } // ss
  return sdiagrams;
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

template class ppsc::diagram_handler<ppsc::diagram_handler_type::oca_sdh>;

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<>
class diagram_handler<diagram_handler_type::oca_sdh>::impl :
    public diagram_handler_base<
      diagram_handler<diagram_handler_type::oca_sdh>::impl,
      oca::sdiagrams_type,
      diagram_class_type::pseudo_particle_self_energy > {

public:

  // ---------------------------------------------------------------------
  oca::sdiagrams_type build_all_diagrams(hilbert_space_type & hilbert_space,
    ppgfs_type & ppGfs, gf_verts_type & gf_verts, pp_ints_type & pp_ints) {
    return oca::build_all_sigma_diagrams(hilbert_space, ppGfs, pp_ints);
  }

  // ---------------------------------------------------------------------
  void diagram_dispatch(int tstp, gf_tstp_type & gtstp,
			oca::sdiagram_configuration<> & diagram,
			double beta, double h, int kt, int nomp) {
    oca::sdiagram_dispatch(tstp, gtstp, diagram, beta, h, kt, nomp);
  }
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_CPP
