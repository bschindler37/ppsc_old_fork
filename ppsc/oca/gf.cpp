
#ifndef _PPSC_OCA_GF_CPP
#define _PPSC_OCA_GF_CPP

// -----------------------------------------------------------------------
//
// OCA single-particle Green's function integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/diag_base.hpp"
#include "ppsc/diag_impl.hpp"

#include "gf_diag.hpp"
#include "gf_disp.hpp"

#include "gf.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;

// ----------------------------------------------------------------------
void get_sectors(int sa, pp_int_type & hyb, gf_vert_type & gf_vertex,
		 int & sb_out, int & sc_out, int & sd_out, bool & valid_diagram) {

  valid_diagram = false;
  int sb =           hyb.op1.to_sector_[sa]; if(sb == -1) return; sb_out = sb;
  int sc =     gf_vertex.op1.to_sector_[sb]; if(sc == -1) return; sc_out = sc;
  int sd =           hyb.op2.to_sector_[sc]; if(sd == -1) return; sd_out = sd;
  int sa_ref = gf_vertex.op2.to_sector_[sd]; if(sa_ref != sa) return;
  valid_diagram = true;
}

// ----------------------------------------------------------------------
bool valid_sectors(int ga_sector, pp_int_type & hyb, gf_vert_type & gf_vertex) {
  int sb, sc, sd; bool valid_diagram;
  get_sectors(ga_sector, hyb, gf_vertex, sb, sc, sd, valid_diagram);
  return valid_diagram;
}

// -----------------------------------------------------------------------
gdiagram_configuration<> construct_gf_diagram(
  gf_tstp_type & gtstp, int gf_sig, int vertex_idx, int sa, ppgfs_type & ppGfs,
  pp_int_type & hyb, gf_vert_type & gf_vertex) {

  int sb, sc, sd;
  bool valid_diagram;
  get_sectors(sa, hyb, gf_vertex, sb, sc, sd, valid_diagram);

  auto op1 = hyb.op2.M_[sc];
  auto op2 = gf_vertex.op1.M_[sb];
  auto op3 = hyb.op1.M_[sa];
  auto op4 = gf_vertex.op2.M_[sd];

  // -- A bwd prop fermionic hybridization line gives a fermi factor
  // -- nb! By def. a standard "fwd" hybridization line is propagating bwd
  // -- in the OCA single-particle Green's function diagram (see gf_diag.hpp)

  int diagram_sig = 1;

  // previously there was an extra sign in the bwd prop fermionic hybr line
  // which was corrected for here. No need anymore.

  /*
  int fermion = -1;
  if(hyb.sig == fermion) {
    diagram_sig = -1 * hyb.dir;
  }
  */

  // -- new sign calc?

  // -- Diagram prefactor [Eckstein, Werner PRB 82, 115155 (2010)]
  const int fermion=-1, boson=+1, fwd=+1, bwd=-1;

  // This practially means that we always get a minus sign for a fermionic hybridization
  if(hyb.sig == fermion && hyb.dir == fwd) diagram_sig = -1;
  if(hyb.sig == fermion && hyb.dir == bwd) diagram_sig = -1;

  //if(hyb.sig == boson) diagram_sig *= 1; // DEBUG

  value_type prefactor = diagram_sig;

  /*
  std::cout << "gf diagram: sa, sb, sc, sd = "
	    << sa << ", " << sb << ", "
	    << sc << ", " << sd
	    << std::endl;
  */

  /*
  std::cout << "hyb, sig, dir, pref "
	    << hyb.sig << ", "
	    << hyb.dir << ", "
	    << prefactor
	    << std::endl;
  */

  return gdiagram_configuration<>(
    gtstp, gf_sig, // nb the sign here is inactive... FIXME?
    op1, op2, op3, op4,
    ppGfs[sa], ppGfs[sb], ppGfs[sc], ppGfs[sd],
    hyb.lam, vertex_idx, prefactor);
}

// -----------------------------------------------------------------------
gdiagrams_type build_all_gf_diagrams(
  hilbert_space_type & hilbert_space, ppgfs_type & ppGfs, gf_verts_type & gf_verts,
  pp_ints_type & pp_ints, int gf_sig=-1) {

  gdiagrams_type gdiagrams;

  for(int vertex_idx = 0; vertex_idx < gf_verts.size(); vertex_idx++) {
    auto gf_vertex = gf_verts[vertex_idx];
    for(auto hyb : pp_ints) {
      for(int sa = 0; sa < hilbert_space.ns_; sa++) { // sum over Ga sector

	if(!valid_sectors(sa, hyb, gf_vertex)) continue;

	gf_tstp_type gf_tstp_diagram(0, 0, 1);
	gdiagram_configuration<> diagram = construct_gf_diagram(
	  gf_tstp_diagram, gf_sig, vertex_idx, sa, ppGfs, hyb, gf_vertex);

	gdiagrams.push_back(diagram);

	} // sa
      } // hyb
    } // vertex_idx
  return gdiagrams;
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

template class ppsc::diagram_handler<ppsc::diagram_handler_type::oca_gdh>;

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<>
class diagram_handler<diagram_handler_type::oca_gdh>::impl :
    public diagram_handler_base<diagram_handler<diagram_handler_type::oca_gdh>::impl, oca::gdiagrams_type,
    diagram_class_type::single_particle_greens_function> {

public:

  // ---------------------------------------------------------------------
  oca::gdiagrams_type build_all_diagrams(hilbert_space_type & hilbert_space,
    ppgfs_type & ppGfs, gf_verts_type & gf_verts, pp_ints_type & pp_ints) {
    return oca::build_all_gf_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);
  }

  // ---------------------------------------------------------------------
  void diagram_dispatch(int tstp, gf_tstp_type & gtstp,
			oca::gdiagram_configuration<> & diagram,
			double beta, double h, int kt, int nomp) {
    oca::gdiagram_dispatch(tstp, gtstp, diagram, beta, h, kt, nomp);
  }
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_CPP
