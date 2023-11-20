
#ifndef _PPSC_NCA_SIGMA_CPP
#define _PPSC_NCA_SIGMA_CPP

// -----------------------------------------------------------------------
//
// NCA pseudo particle self energy construction
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "sigma.hpp"

#include "ppsc/diag_base.hpp"
#include "ppsc/diag_impl.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace nca {
// -----------------------------------------------------------------------

using namespace ppsc;

// -----------------------------------------------------------------------
template<int SDIM=-1, int ADIM=-1>
class sdiagram_configuration {

public:

  const int sdim = SDIM;
  const int adim = ADIM;

  // -- Matrix types
  //using lambda_matrix_type = mam::static_matrix_type<1, 1>;
  using lambda_matrix_type = mam::dynamic_matrix_type;
  using dynamic_matrix_type = mam::dynamic_matrix_type;
  using sigma_matrix_type = mam::static_matrix_type<SDIM, SDIM>;
  using ga_matrix_type = mam::static_matrix_type<ADIM, ADIM>;
  using op1_matrix_type = mam::static_matrix_type<SDIM, ADIM>;
  using op2_matrix_type = mam::static_matrix_type<ADIM, SDIM>;

  // -- Contour compound types
  using sigma_type = ppsc::cntr::timeslice_array_matrix_ref<sigma_matrix_type>;
  using ga_type = ppsc::cntr::herm_matrix_matrix_ref<ga_matrix_type>;
  using lambda_type = ppsc::cntr::herm_matrix_matrix_ref<lambda_matrix_type>;

  // -- Array matrix types
  using ga_array_type = mam::array_matrix_ref<ga_matrix_type>;
  using lambda_array_type = mam::array_matrix_ref<lambda_matrix_type>;

  sdiagram_configuration(gf_tstp_type & s_in, int sig,
			 dynamic_matrix_type op1_in,
			 dynamic_matrix_type op2_in,
			 ppgf_type &ga_in, gf_type &lam_in, int lam_i1, int lam_i2,
			 int sigma_sector, value_type prefactor=1.0) :
    s(s_in), sig(sig), op1(op1_in), op2(op2_in), ga(ga_in),
    lam(lam_in), lam_i1(lam_i1), lam_i2(lam_i2),
    sigma_sector(sigma_sector), prefactor(prefactor) {}

  sdiagram_configuration(gf_tstp_type & s_in, int sig,
			 dynamic_matrix_type op1_in,
			 dynamic_matrix_type op2_in,
			 ppgf_type &ga_in, lambda_type &lam_in, // mam data types for hybridiz
			 int lam_i1, int lam_i2,
			 int sigma_sector, value_type prefactor=1.0) :
    s(s_in), sig(sig), op1(op1_in), op2(op2_in), ga(ga_in),
    lam(lam_in), lam_i1(lam_i1), lam_i2(lam_i2),
    sigma_sector(sigma_sector), prefactor(prefactor) {}

  sdiagram_configuration(const sdiagram_configuration<-1, -1> & d) :
    s(d.s), sig(d.sig), op1(d.op1), op2(d.op2), ga(d.ga),
    lam(d.lam), lam_i1(d.lam_i1), lam_i2(d.lam_i2),
    sigma_sector(d.sigma_sector), prefactor(d.prefactor) {}

  int out_idx() { return sigma_sector; }
  int out_size() { return s.nflavour(); }

  sigma_type s;
  int sig;
  ga_type ga;
  lambda_type lam;
  int lam_i1, lam_i2;
  op1_matrix_type op1;
  op2_matrix_type op2;
  int sigma_sector;
  value_type prefactor;
};

typedef std::vector<sdiagram_configuration<> > sdiagrams_type;

// ----------------------------------------------------------------------
void get_sectors(int sigma_sector, pp_int_type & hyb, int & sa_out, bool & valid_diagram) {

  valid_diagram = false;
  int ss = sigma_sector;
  int sa = hyb.op2.to_sector_[ss]; if(sa == -1) return; sa_out = sa;
  int ss_ref = hyb.op1.to_sector_[sa]; if(ss_ref != ss) return;
  valid_diagram = true;
}

// ----------------------------------------------------------------------
bool valid_sectors(int sigma_sector, pp_int_type & hyb) {
  int sa; bool valid_diagram;
  get_sectors(sigma_sector, hyb, sa, valid_diagram);
  return valid_diagram;
}

// ----------------------------------------------------------------------
sdiagram_configuration<> construct_sigma_diagram(gf_tstp_type & ststp,
  int sigma_sig, int sigma_sector, ppgfs_type& ppG, pp_int_type & hyb) {

  int ss = sigma_sector;
  int sa;
  bool valid_diagram;
  get_sectors(ss, hyb, sa, valid_diagram);

  assert( valid_diagram );

  auto op1 = hyb.op1.M_[sa];
  auto op2 = hyb.op2.M_[ss];

  // -- Diagram prefactor [Eckstein, Werner PRB 82, 115155 (2010)]
  const int fermion=-1, boson=+1, fwd=+1, bwd=-1;

  double stat_sign = +1; // for everything except.

  // bwd propagating fermionic hybridization
  // => one operator commutation, i.e. a fermi minus sign

  if(hyb.sig == fermion && hyb.dir == bwd) stat_sign = -1;

  int order = 1;
  value_type I = value_type(0., 1.);
  value_type prefactor = std::pow(I, order) * stat_sign;

  //value_type prefactor = value_type(0., 1.);

  /*
  std::cout << "nca sigma hyb.idx1, hyb.idx2, prefactor = "
	    << hyb.idx1 << ", " << hyb.idx2 << ", " << prefactor << std::endl;
  */

  return sdiagram_configuration<>(ststp, sigma_sig,
    op1, op2, ppG[sa], hyb.lam, hyb.idx1, hyb.idx2, sigma_sector, prefactor);

}

// ----------------------------------------------------------------------
template<class HS>
sdiagrams_type build_all_sigma_diagrams(HS & hilbert_space,
  ppgfs_type & ppGfs, pp_ints_type & pp_ints) {

  sdiagrams_type sdiagram_list;

  for(int ss = 0; ss < hilbert_space.ns_; ss++) { // sigma sector
    int sector_dim = hilbert_space.ssdim_[ss];
    int sector_sig = hilbert_space.sig_[ss];
    for(auto hyb : pp_ints) {

      if(!valid_sectors(ss, hyb)) continue;

      gf_tstp_type sigma_tstp_diagram_dummy(0, 0, sector_dim);
      sdiagram_configuration<> diagram = construct_sigma_diagram(
        sigma_tstp_diagram_dummy, sector_sig, ss, ppGfs, hyb);

      sdiagram_list.push_back(diagram);

    }
  }
  return sdiagram_list;
}

// -----------------------------------------------------------------------
template<class DIAGT>
void sdiagram_mat(DIAGT & d, double beta, int kt, int nomp) {

#pragma omp parallel num_threads(nomp)
  {
#pragma omp for schedule(dynamic)
    for(int m = 0; m < d.s.ntau(); m++) {

      d.s.mat(m).noalias() = value_type(0., 1.)
	* d.lam.mat(m)(d.lam_i1, d.lam_i2) * (d.op1 * d.ga.mat(m) * d.op2);
    }
  }
}

// -----------------------------------------------------------------------
template<class DIAGT>
void sdiagram_tstp(int n, DIAGT & d, double beta, double h, int kt, int nomp) {

  // -- Dimensions
  int mmax = d.s.ntau() - 1;

#pragma omp parallel num_threads(nomp)
  { // omp parallel

    // -- Retarded component
#pragma omp for schedule(dynamic)
    //for(int j = 0; j < n; j++) {
    for(int j = 0; j < n+1; j++) { // use when removing separate equaltime contrib

      d.s.ret(j).noalias() = d.lam.gtr(n, j)(d.lam_i1, d.lam_i2)
	* (d.op1 * d.ga.ret(n, j) * d.op2);
    }

    // -- Lesser component
#pragma omp for schedule(dynamic)
    //for(int j = 0; j < n; j++) {
    for(int j = 0; j < n+1; j++) { // use when removing separate equaltime contrib

      d.s.les(j).noalias() = d.lam.les(j, n)(d.lam_i1, d.lam_i2)
	* (d.op1 * d.ga.les(j, n) * d.op2);
    }

    // -- TV mixed component
#pragma omp for schedule(dynamic)
    for(int m = 0; m < d.s.ntau(); m++) {

      d.s.tv(m).noalias() = d.lam.tv(n,m)(d.lam_i1, d.lam_i2)
	* (d.op1 * d.ga.tv(n, m) * d.op2);
    }

  } // omp parallel

  /*
  // -- These separate calculations are only required if we want the same
  // -- Herimcity error as the original implementation

  // -- RETARDED, Equal time component
  d.s.ret(n).noalias() = ((-d.lam.ret(n,n)).adjoint() + d.lam.les(n,n))(d.lam_i1,d.lam_i2)
    * (d.op1 * d.ga.ret(n, n) * d.op2);

  // -- LESSER, Equal time component
  d.s.les(n).noalias() = (-d.lam.les(n, n)).adjoint()(d.lam_i1,d.lam_i2)
    * (d.op1 * d.ga.les(n, n) * d.op2);
  */

}

// -----------------------------------------------------------------------
template<class DIAGT>
void sdiagram_dispatch(int tstp, gf_tstp_type & ststp, DIAGT & diagram,
		       double beta, double h, int kt, int nomp) {

  // -- Point diagram timestep to given timestep
  new (& diagram.s) typename sdiagram_configuration<>::sigma_type(ststp);

  // -- Dispatch on tstp
  if(tstp == -1) {
    //{ Timer tmr("nca::sdiagram_mat");
    sdiagram_mat(diagram, beta, kt, nomp);
    //}
  } else {
    //{ Timer tmr("nca::sdiagram_tstp");
    sdiagram_tstp(tstp, diagram, beta, h, kt, nomp);
    //}
  } // tstp

  // -- Post-multiplication by diagram prefactor
  diagram.s *= diagram.prefactor;
}

// -----------------------------------------------------------------------
} // end namespace nca
} // end namespace ppsc
// -----------------------------------------------------------------------

template class ppsc::diagram_handler<ppsc::diagram_handler_type::nca_sdh>;

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<>
class diagram_handler<diagram_handler_type::nca_sdh>::impl :
    public diagram_handler_base<diagram_handler<diagram_handler_type::nca_sdh>::impl, nca::sdiagrams_type,
    diagram_class_type::pseudo_particle_self_energy> {

public:

  // ---------------------------------------------------------------------
  nca::sdiagrams_type build_all_diagrams(hilbert_space_type & hilbert_space,
    ppgfs_type & ppGfs, gf_verts_type & gf_verts, pp_ints_type & pp_ints) {
    return nca::build_all_sigma_diagrams(hilbert_space, ppGfs, pp_ints);
  }

  // ---------------------------------------------------------------------
  void diagram_dispatch(int tstp, gf_tstp_type & gtstp,
			nca::sdiagram_configuration<> & diagram,
			double beta, double h, int kt, int nomp) {
    nca::sdiagram_dispatch(tstp, gtstp, diagram, beta, h, kt, nomp);
  }
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_NCA_SIGMA_CPP
