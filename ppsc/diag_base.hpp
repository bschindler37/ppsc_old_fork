
#ifndef _PPSC_DIAG_BASE_HPP
#define _PPSC_DIAG_BASE_HPP

// -----------------------------------------------------------------------
//
// Pseudo particle strong coupling diagram handler and symmetry analyser
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
enum diagram_class_type {
  pseudo_particle_self_energy, single_particle_greens_function
};

// -----------------------------------------------------------------------

template<class Derived, class DIAGST, diagram_class_type DIAG_CLASS>
class diagram_handler_base {
  
public:

  const diagram_class_type diagram_class = DIAG_CLASS;
  
  typedef DIAGST diagrams_type;

  // ---------------------------------------------------------------------
  // -- CRTP from Derived

  template<class HS> diagrams_type build_all_diagrams(
    HS & hilbert_space, ppgfs_type & ppGfs, gf_verts_type & gf_verts,
    pp_ints_type & pp_ints) {
    
    return static_cast<Derived *>(this)->
      build_all_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);    
  }

  template<class DIAGT> void diagram_dispatch(
    int tstp, gf_tstp_type & gf_tstp, DIAGT & diagram, double beta,
    double h, int kt, int nomp) {

    static_cast<Derived *>(this)->
      diagram_dispatch(tstp, gf_tstp, diagram, beta, h, kt, nomp);
  }
  
  // ---------------------------------------------------------------------
  void update_diagrams(hilbert_space_type & hilbert_space, ppgfs_type & ppGfs,
		       gf_verts_type & gf_verts, pp_ints_type & pp_ints) {

    diagrams = build_all_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);

    switch(diagram_class) {
      case diagram_class_type::pseudo_particle_self_energy :
	out_tstp_sizes = hilbert_space.ssdim_;
	break;
      case diagram_class_type::single_particle_greens_function :
	out_tstp_sizes = std::vector<int>(gf_verts.size(), 1);
	break;
      default:
	std::cerr << "--> diagram_handler_base: update_diagrams, error diagram_class"
		  << std::endl;
	exit(0);
    }
  }

  // ------------------------------------------------------------------
  gf_tstps_type evaluate_to_out_idx(int tstp, int ntau, double beta,
				    double h, int kt, int nomp) {
    
    gf_tstps_type gf_tstp_vertex_list =
      get_list_of_timesteps(tstp, ntau, out_tstp_sizes);

    if( !has_symmetries() ) { // No symmetries

      //std::cout << "--> Full diagram evaluation (no symmetries)" << std::endl;

      for( auto diagram : diagrams ) { // Evaluate all diagrams
	gf_tstp_type gf_tstp_diagram(tstp, ntau, diagram.out_size());
	diagram_dispatch(tstp, gf_tstp_diagram, diagram, beta, h, kt, nomp);
	gf_tstp_vertex_list[diagram.out_idx()].incr(gf_tstp_diagram, 1.0);
      }

    } else { // use symmetry reduction_map
      
      //std::cout << "--> Diagram evaluation using reduction_map" << std::endl;
      
      for( auto elem : reduction_map ) {

	int diag_idx = elem.first;
	auto diagram = diagrams[diag_idx];

	//std::cout << "diag_idx = " << diag_idx << std::endl;
	
	gf_tstp_type gf_tstp_diagram(tstp, ntau, diagram.out_size());
	diagram_dispatch(tstp, gf_tstp_diagram, diagram, beta, h, kt, nomp);
	
	for( auto el : elem.second ) {
	  int out_idx = el.first;
	  int prefactor = el.second;

	  //std::cout << "out_idx, prefactor = " << out_idx << ", " << prefactor << std::endl;
	  
	  gf_tstp_vertex_list[out_idx].incr(gf_tstp_diagram, 1.0 * prefactor);
	}
      }
    }

    return gf_tstp_vertex_list;
  }
  
  // ---------------------------------------------------------------------
  void set_pp_self_energy(ppgfs_type & Sigma,
			  int tstp, int ntau, double beta,
			  double h, int kt, int nomp) {

    assert( diagram_class == diagram_class_type::pseudo_particle_self_energy );

    gf_tstps_type sigma_tstp_sectors =
      evaluate_to_out_idx(tstp, ntau, beta, h, kt, nomp);

    for(int sector = 0; sector < sigma_tstp_sectors.size(); sector++) {
      Sigma[sector].set_timestep(tstp, sigma_tstp_sectors[sector]);
    }
  }
  
  // ---------------------------------------------------------------------
  void add_to_pp_self_energy(ppgfs_type & Sigma,
			     int tstp, int ntau, double beta,
			     double h, int kt, int nomp) {

    assert( diagram_class == diagram_class_type::pseudo_particle_self_energy );

    gf_tstps_type sigma_tstp_sectors =
      evaluate_to_out_idx(tstp, ntau, beta, h, kt, nomp);

    for(int sector = 0; sector < sigma_tstp_sectors.size(); sector++) {
      Sigma[sector].incr_timestep(tstp, sigma_tstp_sectors[sector], 1.0);
    }
  }

  // ---------------------------------------------------------------------
  void analyze_diagram_symmetries(int tstp, int ntau, double beta, double h,
				  int kt, int nomp, double tol=1e-12) {

    //std::cout << "--> diagram_handler_base: Diagram symmetry analysis" << std::endl;

    gf_tstps_type tstp_diagrams =
      evaluate_all_diagrams(tstp, ntau, beta, h, kt, nomp);

    // -- Debug printing
    /*
    for( int i1 = 0; i1 < tstp_diagrams.size(); i1++) {
      std::cout << std::string(72, '=') << std::endl;
      std::cout << "--> Diagram: " << i1 << std::endl;
      std::cout << "out_idx = " << diagrams[i1].out_idx() << std::endl;
      
      cntr::timeslice_array_matrix_ref<mam::dynamic_matrix_type> tstp1(tstp_diagrams[i1]);
      std::cout << "is_zero = " << tstp1.is_zero() << std::endl;

      if(tstp1.is_zero()) continue;
      
      for( int i2 = 0; i2 < i1; i2++) {
	cntr::timeslice_array_matrix_ref<mam::dynamic_matrix_type> tstp2(tstp_diagrams[i2]);
	
	std::cout << "equal to diagram nr. " << i2 << " ? = "
	  // << i1 << ", " << i2 << ", "
		  << tstp1.is_equal_to(tstp2, tol) << std::endl;
      }
    }
    */

    construct_reduction_map(tstp_diagrams, tol);
    //print_report();
  }

  // ---------------------------------------------------------------------
  gf_tstps_type evaluate_all_diagrams(int tstp, int ntau, double beta,
				      double h, int kt, int nomp) {
    // -- setup timestep storage
    std::vector<int> diagrams_sector_sizes;
    for( auto diagram : diagrams ) diagrams_sector_sizes.push_back(diagram.out_size());
    gf_tstps_type tstp_diagrams = get_list_of_timesteps(tstp, ntau, diagrams_sector_sizes);

    // -- Eval all diagrams at tstp
    for( int sd_idx = 0; sd_idx < diagrams.size(); sd_idx++)
      diagram_dispatch(tstp, tstp_diagrams[sd_idx], diagrams[sd_idx], beta, h, kt, nomp);      

    return tstp_diagrams;
  }
  
  // ---------------------------------------------------------------------
  void construct_reduction_map(gf_tstps_type & tstp_diagrams, double tol=1e-12) {

    //std::cout << "--> construct_reduction_map" << std::endl;
    
    reduction_map.clear();
    
    std::cout << "--> diag_base::construct_reduction_map:" << std::endl;
    
    // loop over all diagrams
    for(int i1 = 0; i1 < tstp_diagrams.size(); i1++) {

      cntr::timeslice_array_matrix_ref<mam::dynamic_matrix_type> tstp1(tstp_diagrams[i1]);

      std::cout << "i1, is_zero = " << i1 << ", " << tstp1.is_zero(tol) << std::endl;
      
      if(tstp1.is_zero(tol)) continue; // Do not account for zero diagrams
      
      int i2;
      bool is_equal = false;
      
      for( auto elem : reduction_map ) { // chech equality with all preceding (unique) diagrams
	
	i2 = elem.first; // take out idx for, unique diagram
	cntr::timeslice_array_matrix_ref<mam::dynamic_matrix_type> tstp2(tstp_diagrams[i2]);
	
	double abs_diff = tstp1.abs_diff(tstp2);
	std::cout << "cf with unique diagram: i2, abs_diff, is_equal = "
		  << i2 << ", " << abs_diff << ", "
		  << tstp1.is_equal_to(tstp2, tol) << std::endl;

	// -- per component diffs for debugging...
	double abs_diff_ret = tstp1.ret.abs_diff(tstp2.ret);
	double abs_diff_les = tstp1.les.abs_diff(tstp2.les);
	double abs_diff_tv = tstp1.tv.abs_diff(tstp2.tv);
	double abs_diff_mat = tstp1.mat.abs_diff(tstp2.mat);
	std::cout << "abs_diff, ret, les, tv, mat = "
		  << abs_diff_ret << ", "
		  << abs_diff_les << ", "
		  << abs_diff_tv << ", "
		  << abs_diff_mat << std::endl;

	is_equal = tstp1.is_equal_to(tstp2, tol);
	if(is_equal) break;      
      }
      
      if(is_equal) { // found equivalent diagram

	if( reduction_map[i2].count(diagrams[i1].out_idx()) == 0 ) {
	  // new out sector, set prefactor to one
	  reduction_map[i2][diagrams[i1].out_idx()] = 1;

	} else {
	  // already existing out sector
	  // increase its prefactor by one
	  reduction_map[i2][diagrams[i1].out_idx()] += 1;
	}
      }  else {
	// new unique diagram found, add it to the reduced set of diagrams
	reduction_map[i1][diagrams[i1].out_idx()] = 1;
      } 
    }

    print_map();
  }

  // ---------------------------------------------------------------------
  void print_map() {
    std::cout << "--> diagram_reduction_base: reduction_map" << std::endl;
    for( auto elem : reduction_map ) {
      std::cout << "didx = " << elem.first << std::endl;
      for( auto e : elem.second )
	std::cout << "out_idx, prefactor = " << e.first << ", " << e.second << std::endl;
    }    
  }

  // ---------------------------------------------------------------------
  void print_report() {
    std::cout << "--> Diagram handler symmetry report" << std::endl;
    print_map();
    std::cout << "Number of full diagrams are " << total_number()
	      << " which is reduced to " << reduced_number() << "." << std::endl; 
  }
  
  // ---------------------------------------------------------------------

  diagrams_type & get_diagrams() { return diagrams; }
  int size() { return diagrams.size(); }
  bool has_symmetries() { return !reduction_map.empty(); }
  void clear_symmetries() { reduction_map.clear(); }

  // ---------------------------------------------------------------------

  int total_number() { return diagrams.size(); } 
  int reduced_number() {
    if( has_symmetries() ) return reduction_map.size();
    else                   return total_number();
  } 

  // ---------------------------------------------------------------------

  void read_symmetries(input_archive_type & ia) { ia(reduction_map); }
  void write_symmetries(output_archive_type & oa) { oa(reduction_map); }
  
  // ---------------------------------------------------------------------
  
  diagrams_type diagrams;
  std::map< int, std::map<int, int> > reduction_map;
  std::vector<int> out_tstp_sizes;
  
};
  
// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_DIAG_BASE_HPP
