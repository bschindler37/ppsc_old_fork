
#ifndef _PPSC_DIAG_HPP
#define _PPSC_DIAG_HPP

// -----------------------------------------------------------------------
//
// Diagram handler PIMPL interface
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/ppsc.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
enum diagram_handler_type {
  nca_gdh, nca_sdh, oca_gdh, oca_sdh
};

// -----------------------------------------------------------------------
template<diagram_handler_type DHT>
class diagram_handler {

public:
  
  diagram_handler();
  ~diagram_handler();

  void update_diagrams(hilbert_space_type & hilbert_space, ppgfs_type & ppGfs,
		       gf_verts_type & gf_verts, pp_ints_type & pp_ints);

  void analyze_diagram_symmetries(int tstp, int ntau, double beta, double h,
				  int kt, int nomp, double tol=1e-12);  
  
  gf_tstps_type evaluate_to_out_idx(int tstp, int ntau, double beta,
				    double h, int kt, int nomp);

 void add_to_pp_self_energy(ppgfs_type & Sigma, int tstp, int ntau,
			    double beta, double h, int kt, int nomp);

  void set_pp_self_energy(ppgfs_type & Sigma, int tstp, int ntau, double beta,
			  double h, int kt, int nomp);

  bool has_symmetries();
  void clear_symmetries();
  int total_number();
  int reduced_number();
  void print_report();

  void read_symmetries(input_archive_type & ia);
  void write_symmetries(output_archive_type & oa);
  
private:
  class impl;
  std::unique_ptr<impl> pimpl;

};

// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_DIAG_HPP
