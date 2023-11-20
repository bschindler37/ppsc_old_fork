
#ifndef _PPSC_NCA_GF_CPP
#define _PPSC_NCA_GF_CPP

// -----------------------------------------------------------------------
//
// NCA single particle Green's function construction
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "gf.hpp"

#include "ppsc/diag_base.hpp"
#include "ppsc/diag_impl.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace nca {
// -----------------------------------------------------------------------

using namespace ppsc;

// -----------------------------------------------------------------------
template<int ADIM=-1, int BDIM=-1>
class gdiagram_configuration {

public:

  const int adim = ADIM;
  const int bdim = BDIM;

  // -- Matrix types
  using lambda_matrix_type = mam::static_matrix_type<1, 1>;
  using dynamic_matrix_type = mam::dynamic_matrix_type;
  using g_matrix_type = mam::static_matrix_type<1, 1>;
  using ga_matrix_type = mam::static_matrix_type<ADIM, ADIM>;
  using gb_matrix_type = mam::static_matrix_type<BDIM, BDIM>;
  using op1_matrix_type = mam::static_matrix_type<BDIM, ADIM>;
  using op2_matrix_type = mam::static_matrix_type<ADIM, BDIM>;

  // -- Contour compound types
  using g_type = ppsc::cntr::timeslice_array_matrix_ref<g_matrix_type>;
  using ga_type = ppsc::cntr::herm_matrix_matrix_ref<ga_matrix_type>;
  using gb_type = ppsc::cntr::herm_matrix_matrix_ref<gb_matrix_type>;

  // -- Array matrix types
  using ga_array_type = mam::array_matrix_ref<ga_matrix_type>;
  using gb_array_type = mam::array_matrix_ref<gb_matrix_type>;

  gdiagram_configuration(gf_tstp_type & g_in, int sig,
			 dynamic_matrix_type op1_in,
			 dynamic_matrix_type op2_in,
			 ppgf_type & ga_in,
			 ppgf_type & gb_in,
			 int vertex_idx,
			 value_type prefactor=1.0) :
    g(g_in), sig(sig), op1(op1_in), op2(op2_in), ga(ga_in), gb(gb_in),
    vertex_idx(vertex_idx), prefactor(prefactor) {}

  gdiagram_configuration(const gdiagram_configuration<-1, -1> & d) :
    g(d.g), sig(d.sig), op1(d.op1), op2(d.op2), ga(d.ga), gb(d.gb),
    vertex_idx(d.vertex_idx), prefactor(d.prefactor) {}

  int out_idx() { return vertex_idx; }
  int out_size() { return g.nflavour(); }

  g_type g;
  int sig;
  ga_type ga;
  gb_type gb;
  op1_matrix_type op1;
  op2_matrix_type op2;
  int vertex_idx;
  value_type prefactor;
};

typedef std::vector<gdiagram_configuration<> > gdiagrams_type;

// ----------------------------------------------------------------------
void get_sectors(int ga_sector, gf_vert_type & gf_vertex,
		 int & sb_out, bool & valid_diagram) {

  valid_diagram = false;
  int sa = ga_sector;
  int sb = gf_vertex.op1.to_sector_[sa]; if(sb == -1) return; sb_out = sb;
  int sa_ref = gf_vertex.op2.to_sector_[sb]; if(sa_ref != sa) return;
  valid_diagram = true;
}

// ----------------------------------------------------------------------
bool valid_sectors(int ga_sector, gf_vert_type & gf_vertex) {
  int sb; bool valid_diagram;
  get_sectors(ga_sector, gf_vertex, sb, valid_diagram);
  return valid_diagram;
}

// -----------------------------------------------------------------------
gdiagram_configuration<> construct_gf_diagram(
  gf_tstp_type & gtstp, int gf_sig, int vertex_idx, int ga_sector,
  ppgfs_type & ppG, gf_vert_type & gf_vertex) {

  int sa = ga_sector;
  int sb;
  bool valid_diagram;
  get_sectors(sa, gf_vertex, sb, valid_diagram);

  assert( valid_diagram );

  auto op1 = gf_vertex.op1.M_[sa];
  auto op2 = gf_vertex.op2.M_[sb];

  value_type prefactor = value_type(0., ppG[sb].sig());

  return gdiagram_configuration<-1, -1>(
    gtstp, gf_sig, op1, op2, ppG[sa], ppG[sb], vertex_idx, prefactor);
}

// -----------------------------------------------------------------------
template<class HS>
gdiagrams_type build_all_gf_diagrams(
  HS & hilbert_space, ppgfs_type & ppGfs, gf_verts_type & gf_verts, int gf_sig=-1) {

  gdiagrams_type gdiagram_list;

  for(int vertex_idx = 0; vertex_idx < gf_verts.size(); vertex_idx++) {
    auto gf_vertex = gf_verts[vertex_idx];
    for(int sa = 0; sa < ppGfs.size(); sa++) { // sum over Ga sector

      if(!valid_sectors(sa, gf_vertex)) continue;

      gf_tstp_type gf_tstp_diagram(0, 0, 1);

      gdiagram_configuration<> diagram = construct_gf_diagram(
        gf_tstp_diagram, gf_sig, vertex_idx, sa, ppGfs, gf_vertex);

      gdiagram_list.push_back(diagram);

    } // sa
  } // vertex_idx
  return gdiagram_list;
}

// -----------------------------------------------------------------------
template<class DIAGT>
void gdiagram_mat(DIAGT & d, double beta, int kt, int nomp) {

#pragma omp parallel num_threads(nomp)
  {
    int ntau = d.g.ntau();
#pragma omp for schedule(dynamic)
    for(int m = 0; m < ntau; m++) {
      d.g.mat(m)(0,0) = value_type(0., d.gb.sig())
	* ( d.op1 * d.ga.mat(m) * d.op2 * d.gb.mat(ntau-1 - m) ).trace();
    }
  }
}

// -----------------------------------------------------------------------
template<class DIAGT>
void gdiagram_tstp(int n, DIAGT & d, double beta, double h, int kt, int nomp) {

  typedef typename DIAGT::g_matrix_type g_matrix_type;

#pragma omp parallel num_threads(nomp)
  { // omp parallel

    // -- Lesser component
#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      d.g.les(j)(0,0) = (d.op1 * d.ga.les(j, n) * d.op2 * d.gb.ret(n, j)).trace();
    }

    // -- Retarded component (from greater and lesser)
#pragma omp for schedule(dynamic)
    for(int j = 0; j < n+1; j++) {
      /*
      value_type g_gtr = (d.op1 * d.ga.ret(n, j) * d.op2 * d.gb.les(j, n)).trace();
      d.g.ret(j)(0,0) = g_gtr - d.g.les(j).adjoint()(0,0); // <-- This assumes g is its own conj!
      */

      value_type g_gtr = (d.op1 * d.ga.ret(n, j) * d.op2 * d.gb.les(j, n)).trace();
      value_type g_les = (d.op1 * d.ga.les(n, j) * d.op2 * d.gb.ret(j, n)).trace();
      d.g.ret(j)(0,0) = g_gtr - g_les; // <-- This assumes g is its own conj!

    }

    // -- TV mixed component
#pragma omp for schedule(dynamic)
    for(int m = 0; m < d.g.ntau(); m++) {
      d.g.tv(m)(0,0) = (d.op1 * d.ga.tv(n, m) * d.op2 * d.gb.vt(m, n)).trace();
    }

  } // omp parallel
}

// -----------------------------------------------------------------------
template<class DIAGT>
void gdiagram_dispatch(int tstp, gf_tstp_type & gtstp, DIAGT & diagram,
		       double beta, double h, int kt, int nomp) {

  // -- Point diagram timestep to given timestep
  new (& diagram.g) typename gdiagram_configuration<>::g_type(gtstp);

  // -- Dispatch on tstp
  if(tstp == -1) {
    //{ Timer tmr("nca::sdiagram_mat");
    gdiagram_mat(diagram, beta, kt, nomp);
    //}
  } else {
    //{ Timer tmr("nca::sdiagram_tstp");
    gdiagram_tstp(tstp, diagram, beta, h, kt, nomp);
    //}
  } // tstp

  // -- Post-multiplication by diagram prefactor
  diagram.g *= diagram.prefactor;
}

// -----------------------------------------------------------------------
} // end namespace nca
} // end namespace ppsc
// -----------------------------------------------------------------------

template class ppsc::diagram_handler<ppsc::diagram_handler_type::nca_gdh>;

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<>
class diagram_handler<diagram_handler_type::nca_gdh>::impl :
    public diagram_handler_base<diagram_handler<diagram_handler_type::nca_gdh>::impl, nca::gdiagrams_type,
    diagram_class_type::single_particle_greens_function> {

public:

  // ---------------------------------------------------------------------
  nca::gdiagrams_type build_all_diagrams(hilbert_space_type & hilbert_space,
    ppgfs_type & ppGfs, gf_verts_type & gf_verts, pp_ints_type & pp_ints) {
    return nca::build_all_gf_diagrams(hilbert_space, ppGfs, gf_verts);
  }

  // ---------------------------------------------------------------------
  void diagram_dispatch(int tstp, gf_tstp_type & gtstp,
			nca::gdiagram_configuration<> & diagram,
			double beta, double h, int kt, int nomp) {
    nca::gdiagram_dispatch(tstp, gtstp, diagram, beta, h, kt, nomp);
  }
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_NCA_GF_CPP
