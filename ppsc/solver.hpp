
#ifndef _PPSC_SOLVER_HPP
#define _PPSC_SOLVER_HPP

// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling expansion solver
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include <boost/lexical_cast.hpp>

// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

#include "ppsc/nca/nca.hpp"
#include "ppsc/oca/oca.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
template<class HAM>
class solver {

public:

  template<class HILB>
  solver(int nt, int ntau, double beta, double h, int kt, int nomp,
	 HILB & hilbert_space, int expansion_order=1) :
    nt(nt), ntau(ntau), beta(beta), h(h), kt(kt), nomp(nomp),
    hilbert_space(hilbert_space), expansion_order(expansion_order),
    hamiltonian(nt, hilbert_space),
    rho(get_list_of_operators(nt, hilbert_space.ssdim_)),
    H(get_list_of_functions(nt, hilbert_space.ssdim_)),
    ppGfs(get_list_of_ppgfs(nt, ntau, hilbert_space.ssdim_, hilbert_space.sig_)),
    ppSigmas(get_list_of_ppgfs(nt, ntau, hilbert_space.ssdim_, hilbert_space.sig_)),
    ppSigmas_oca(get_list_of_ppgfs(nt, ntau, hilbert_space.ssdim_, hilbert_space.sig_)), // DEBUG
    Eint_pp(nt+2) {

    // Bug fix for the static integration object
    // this must be run outside of any openmp environment to
    // avoid a race condition inside the integrator setup
    for(int k = 1; k <= 5; k++) integration::I<double>(k);
  }
  
  // ---------------------------------------------------------------------

  void update_diagrams(pp_ints_type & pp_ints, gf_verts_type & gf_verts) {
    nca_sdh.update_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);
    nca_gdh.update_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);
    if(has_oca()) {
      oca_sdh.update_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);            
      oca_gdh.update_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);
    }
  }

  void print_hamiltonian(int tstp) {
    mam::dynamic_matrix_type tmp;
    std::cout << "Hamiltoanian " << tstp << std::endl;
    for( auto idx : range(0, H.size())){
      H[idx].get_value(tstp,tmp);
      std::cout  << tstp << " " << idx << ": " << tmp  << std::endl;
    }
  }

  void update_sigma(int tstp) {
    nca_sdh.set_pp_self_energy(ppSigmas, tstp, ntau, beta, h, kt, nomp);
    if(has_oca()) oca_sdh.add_to_pp_self_energy(ppSigmas, tstp, ntau, beta, h, kt, nomp);

    // DEBUG
    if(has_oca()) oca_sdh.set_pp_self_energy(ppSigmas_oca, tstp, ntau, beta, h, kt, nomp);
  }

  void normalize_ppgf() { normalize_ppgfs_mat(ppGfs, pp_mu, beta); }
  
  void solve_dyson(int tstp) {
    solve_pp_dyson(tstp, ppGfs, ppSigmas, H, pp_mu, beta, h, kt, nomp);
  }

  void update_density_matrix(int tstp) {
    ppsc::update_density_matrix(tstp, ppGfs, rho);
  }
  
  // ---------------------------------------------------------------------
  ppsc::gf_tstps_type get_spgf(int tstp) {

    ppsc::gf_tstps_type nca_gf_tstp_vertex_list =
      nca_gdh.evaluate_to_out_idx(tstp, ntau, beta, h, kt, nomp);
    
    if (expansion_order == 2) {	
      ppsc::gf_tstps_type oca_gf_tstp_vertex_list =
	oca_gdh.evaluate_to_out_idx(tstp, ntau, beta, h, kt, nomp);
      
      for( auto idx : range(0, nca_gf_tstp_vertex_list.size()) )
	nca_gf_tstp_vertex_list[idx].incr(oca_gf_tstp_vertex_list[idx], 1.0);
    }

    return nca_gf_tstp_vertex_list;
  }
  
  // ---------------------------------------------------------------------
  void pp_step(int tstp) {
    int n1 = (tstp <= kt && tstp >= 0 ? 0 : tstp);
    int n2 = (tstp <= kt && tstp >= 0 ? kt : tstp);

    for(int n = n1; n <= n2; n++) update_sigma(n);
    solve_dyson(tstp);
    if(tstp == -1) normalize_ppgf();
    for(int n = n1; n <= n2; n++) update_density_matrix(n);
    for(int n = n1; n <= n2; n++) hamiltonian.update_exp_vals(n, rho);
  }

  // ---------------------------------------------------------------------
  void extrapolate_timestep(int tstp) { extrapolate_ppgf_timestep(tstp, ppGfs, kt); }

  // ---------------------------------------------------------------------
  void symmetry_reduction(int tstp, double tol=1e-9) {

    if(!nca_sdh.has_symmetries()) {
      std::cout << "--> solver_base::symmetry_reduction: nca_sdh" << std::endl;
      nca_sdh.analyze_diagram_symmetries(tstp, ntau, beta, h, kt, nomp, tol);
    }

    if(!nca_gdh.has_symmetries()) {
      std::cout << "--> solver_base::symmetry_reduction: nca_gdh" << std::endl;
      nca_gdh.analyze_diagram_symmetries(tstp, ntau, beta, h, kt, nomp, tol);
    }
    
    if(has_oca() && !oca_sdh.has_symmetries()) {
      std::cout << "--> solver_base::symmetry_reduction: oca_sdh" << std::endl;
      oca_sdh.analyze_diagram_symmetries(tstp, ntau, beta, h, kt, nomp, tol);
    }

    if(has_oca() && !oca_gdh.has_symmetries()) {
      std::cout << "--> solver_base::symmetry_reduction: oca_gdh" << std::endl;
      oca_gdh.analyze_diagram_symmetries(tstp, ntau, beta, h, kt, nomp, tol);
    }

    print_diagrams_summary();
  }

  // ---------------------------------------------------------------------
  void print_diagrams_summary() {

    std::cout << "--> Number of NCA sigma diagrams: "
	      << nca_sdh.total_number() << " reduced to "
	      << nca_sdh.reduced_number() << "." << std::endl;

    if(has_oca())
      std::cout << "--> Number of OCA sigma diagrams: "
		<< oca_sdh.total_number() << " reduced to "
		<< oca_sdh.reduced_number() << "." << std::endl;
    
    std::cout << "--> Number of NCA    gf diagrams: "
	      << nca_gdh.total_number() << " reduced to "
	      << nca_gdh.reduced_number() << "." << std::endl;

    if(has_oca())
      std::cout << "--> Number of OCA    gf diagrams: "
		<< oca_gdh.total_number() << " reduced to "
		<< oca_gdh.reduced_number() << "." << std::endl;
  }

  // ---------------------------------------------------------------------
  bool has_symmetries() {
    if(has_oca()) {
      return nca_sdh.has_symmetries() && nca_gdh.has_symmetries() &&
             oca_sdh.has_symmetries() && oca_gdh.has_symmetries();
    } else {
      return nca_sdh.has_symmetries() && nca_gdh.has_symmetries();
    }
  }

  // ---------------------------------------------------------------------
  void clear_symmetries() {
    nca_sdh.clear_symmetries();
    nca_gdh.clear_symmetries();
    if(has_oca()) {
      oca_sdh.clear_symmetries();
      oca_gdh.clear_symmetries();
    }
  }

  // ---------------------------------------------------------------------
  void read_symmetries(std::string & filename) {

    std::cout << "--> Reading symmetries from file: "
	      << filename << std::endl;
    
    std::ifstream fd(filename);
    input_archive_type ia(fd);

    // -- NB! order important!
    nca_gdh.read_symmetries(ia);
    nca_sdh.read_symmetries(ia);
    oca_gdh.read_symmetries(ia);
    oca_sdh.read_symmetries(ia);

    print_diagrams_summary();    
  }
  
  // ---------------------------------------------------------------------
  void write_symmetries(std::string & filename) {

    std::cout << "--> Writing symmetries to file: "
	      << filename << std::endl;
    
    std::ofstream fd(filename);
    output_archive_type oa(fd);
    
    // -- NB! order important!
    nca_gdh.write_symmetries(oa);
    nca_sdh.write_symmetries(oa);
    oca_gdh.write_symmetries(oa);
    oca_sdh.write_symmetries(oa);
  }
  
  // ---------------------------------------------------------------------
  void update_hamiltonian(int tstp) {
    operator_type Ht;
    hamiltonian.get_hamiltonian(tstp, Ht);      
    for( auto idx : range(0, Ht.M_.size()) )
      if( Ht.to_sector_[idx] != -1 ) H[idx].set_value(tstp, Ht.M_[idx]);
  }

  void update_hamiltonian(int tstp,double nexp) {
    operator_type Ht;
    hamiltonian.get_hamiltonian(tstp, Ht,nexp);   
    for( auto idx : range(0, Ht.M_.size()) )
      if( Ht.to_sector_[idx] != -1 ) H[idx].set_value(tstp, Ht.M_[idx]);
  }

  void update_hamiltonian(int tstp,double nd_exp,double nu_exp) {
    operator_type Ht;
    hamiltonian.get_hamiltonian(tstp, Ht,nd_exp,nu_exp);
    for( auto idx : range(0, Ht.M_.size()) )
      if( Ht.to_sector_[idx] != -1 ) H[idx].set_value(tstp, Ht.M_[idx]);
  }


  void update_hamiltonian() {
    for( int tstp = -1; tstp <= nt; tstp++)
      update_hamiltonian(tstp);
  }

  // ---------------------------------------------------------------------
  void solve_atomic() { zeroth_order(ppGfs, H, pp_mu, beta); }
  void init_real_time() { set_t0_from_mat(ppGfs); }
  bool has_oca() { return expansion_order >= 2; }
  
  // ---------------------------------------------------------------------
  void load(std::string & filename) {

    std::cout << "--> Loading pseudo particles from: " << filename << std::endl;

    hid_t file_id = read_hdf5_file(filename);

    { 
      hid_t group_id = open_group(file_id, "imp");
      pp_mu = read_primitive_type<double>(group_id, "pp_mu");
      close_group(group_id);
    }
    
    for( auto idx : range(0, ppGfs.size()) ) {
      std::string name = std::string("ppG_") + boost::lexical_cast<std::string>(idx);
      hid_t group_id = open_group(file_id, name);
      //ppGfs[idx].read_from_hdf5(-1, group_id); // only read t=0
      ppGfs[idx].read_from_hdf5(group_id);
      close_group(group_id);
    }    

    for( auto idx : range(0, ppSigmas.size()) ) {
      std::string name = std::string("ppSigma_") + boost::lexical_cast<std::string>(idx);
      hid_t group_id = open_group(file_id, name);
      //ppSigmas[idx].read_from_hdf5(-1, group_id); // only read t=0
      ppSigmas[idx].read_from_hdf5(group_id);
      close_group(group_id);
    }

    close_hdf5_file(file_id);
  }

  // ---------------------------------------------------------------------
  void store(hid_t file_id, bool store_pp=false) {

    { // -- hdf group imp scope
      
      hid_t group_id = create_group(file_id, "imp");

      store_int_attribute_to_hid(group_id, "nt", nt);
      store_int_attribute_to_hid(group_id, "ntau", ntau);
      store_double_attribute_to_hid(group_id, "beta", beta);
      store_double_attribute_to_hid(group_id, "h", h);
      store_int_attribute_to_hid(group_id, "kt", kt);
      store_int_attribute_to_hid(group_id, "nomp", nomp);
      
      store_int_attribute_to_hid(group_id, "expansion_order", expansion_order);
      store_double_attribute_to_hid(group_id, "pp_mu", pp_mu);
      
      store_real_data_to_hid(group_id, "Eint_pp", Eint_pp.data(), Eint_pp.size());
      
      // -- Store Hamiltonian params and local exp values
      hamiltonian.store(group_id);
      
      // -- Hamiltonian and density-matrix storage
      // -- first map to large matrix cntr-function and then store to hdf5
      
      int tot_dim = hilbert_space.nh_;
      function_type rho_func(nt, tot_dim);
      
      for(int t=-1; t <= nt; t++) {
	mam::dynamic_matrix_type rho_mat_t(tot_dim, tot_dim);
	rho_mat_t = rho[t+1].to_dense();
	rho_func.set_value(t, rho_mat_t);
      }
      rho_func.write_to_hdf5(group_id, "rho");
    
      /// -- CLOSE IMP GROUP
      close_group(group_id);
    }

    if(store_pp) { // -- Store ppGfs and ppSigmas
      for( auto idx : range(0, ppGfs.size()) ) {
	std::string name = std::string("ppG_") + boost::lexical_cast<std::string>(idx);
	hid_t group_id = create_group(file_id, name);
	store_herm_greens_function(group_id, ppGfs[idx]);
	close_group(group_id);
      }
      for( auto idx : range(0, ppSigmas.size()) ) {
	std::string name = std::string("ppSigma_") + boost::lexical_cast<std::string>(idx);
	hid_t group_id = create_group(file_id, name);
	store_herm_greens_function(group_id, ppSigmas[idx]);
	close_group(group_id);
      }
      for( auto idx : range(0, ppSigmas_oca.size()) ) {
	std::string name = std::string("ppSigma_oca_") + boost::lexical_cast<std::string>(idx);
	hid_t group_id = create_group(file_id, name);
	store_herm_greens_function(group_id, ppSigmas_oca[idx]);
	close_group(group_id);
      }
    } // store_pp

  }

  // ---------------------------------------------------------------------
  void calculate_pp_interaction_energy() {

    std::cout << "--> ppsc::solver::calculate_pp_interaction_energy"
	      << std::endl;

    std::cout << "nt, ntau = " << nt << ", " << ntau << std::endl;
    std::cout << "ppGfs.size() = " << ppGfs.size() << std::endl;
    std::cout << "Eint_pp.size() = " << Eint_pp.size() << std::endl;
    
    for( auto t : range(0, Eint_pp.size()) ) Eint_pp[t] = 0.0;

    for( auto sector : range(0, ppGfs.size()) ) {

      std::cout << "sector = " << sector << std::endl;

      int size = ppGfs[sector].size1();
      int sig = ppGfs[sector].sig();

      std::cout << "sig = " << sig << std::endl;
      std::cout << "sig sigma = " << ppSigmas[sector].sig() << std::endl;
      
      ppgf_type GS(nt, ntau, size, sig);
      //ppgf_type SG(nt, ntau, ppGfs[sector].size1());

      // -- Convert to hermitian matrix response functs
      gf_type gs(nt, ntau, size, sig);

      gf_type gf(nt, ntau, size, sig);
      gf_type sigma(nt, ntau, size, sig);

      for(int tstp=-1; tstp <= nt; tstp++) {
	ppsc::gf_tstp_type tmp(tstp, ntau, size);

	ppGfs[sector].get_timestep(tstp, tmp);
	gf.set_timestep(tstp, tmp);

	ppSigmas[sector].get_timestep(tstp, tmp);
	sigma.set_timestep(tstp, tmp);
      }
      
      if(nt > 0) {

	::cntr::convolution(gs,
			    gf, gf, sigma, sigma,
			    integration::I<double>(kt), beta, h);

	::cntr::pseudo_convolution(GS,
				   ppGfs[sector], ppGfs[sector],
				   ppSigmas[sector], ppSigmas[sector],
				   integration::I<double>(kt), beta, h);
      }
      


      /*
      cntr::pseudo_convolution(SG,
			       ppSigmas[sector], ppSigmas[sector],
			       ppGfs[sector], ppGfs[sector],
			       integration::I<double>(kt), beta, h);
      cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> SG_ref(SG);
      */

      //Eint_pp[0] += -(GS_ref.mat(nt) + SG_ref.mat(nt)).trace().real();
      //Eint_pp[t + 1] += (GS_ref.les(t, t) + SG_ref.les(t,t)).trace().imag();		    

      /*
      cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> GS_ref(GS);
      Eint_pp[0] += -(GS_ref.mat(ntau)).trace().real();
      for( auto t : range(0, nt+1) )
	//Eint_pp[t + 1] += (GS_ref.les(t, t)).trace().real();
	Eint_pp[t + 1] += (GS_ref.les(t, t)).trace().imag();
      */

      cntr::herm_matrix_matrix_ref<mam::dynamic_matrix_type> gs_ref(gs);
      Eint_pp[0] += -(gs_ref.mat(ntau)).trace().real();
      for( auto t : range(0, nt+1) )
	Eint_pp[t + 1] += (gs_ref.les(t, t)).trace().real();
      //Eint_pp[t + 1] += (gs_ref.les(t, t)).trace().imag();
      
    }

    for(int tstp=-1; tstp <= nt; tstp++) {
      std::cout << "tstp, Eint_pp = "
		<< tstp << ", " << Eint_pp[tstp+1] << std::endl;
    }
  }
  
  // ---------------------------------------------------------------------

  int nt;
  int ntau;
  double beta;
  double h;
  int kt;
  int nomp;

  hilbert_space_type hilbert_space;
  int expansion_order;

  HAM hamiltonian;
  operators_type rho;

  functions_type H;
  ppgfs_type ppGfs;
  ppgfs_type ppSigmas;
  ppgfs_type ppSigmas_oca;

  double pp_mu = 0;
  std::vector<double> Eint_pp;

  // -- Diagram handlers
  nca::sigma_diagram_handler nca_sdh;
  oca::sigma_diagram_handler oca_sdh;

  nca::gf_diagram_handler nca_gdh;
  oca::gf_diagram_handler oca_gdh;
  
};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_SOLVER_HPP
