
#ifndef _PPSC_DIAG_CPP
#define _PPSC_DIAG_CPP

// -----------------------------------------------------------------------
//
// Diagram handler PIMPL interface
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/diag.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -- PIMPL methods
  
template<diagram_handler_type DHT>
diagram_handler<DHT>::diagram_handler() : pimpl(new impl()) {};

template<diagram_handler_type DHT>
diagram_handler<DHT>::~diagram_handler() {}

// -----------------------------------------------------------------------

template<diagram_handler_type DHT>
void diagram_handler<DHT>::update_diagrams(
  hilbert_space_type & hilbert_space, ppgfs_type & ppGfs,
  gf_verts_type & gf_verts, pp_ints_type & pp_ints) {
  pimpl->update_diagrams(hilbert_space, ppGfs, gf_verts, pp_ints);
}  

template<diagram_handler_type DHT>
void diagram_handler<DHT>::analyze_diagram_symmetries(
  int tstp, int ntau, double beta, double h, int kt, int nomp, double tol) {
  pimpl->analyze_diagram_symmetries(tstp, ntau, beta, h, kt, nomp, tol);
}
  
template<diagram_handler_type DHT>
gf_tstps_type diagram_handler<DHT>::evaluate_to_out_idx(
  int tstp, int ntau, double beta, double h, int kt, int nomp) {
  return pimpl->evaluate_to_out_idx(tstp, ntau, beta, h, kt, nomp);
}

template<diagram_handler_type DHT>
void diagram_handler<DHT>::add_to_pp_self_energy(
  ppgfs_type & Sigma, int tstp, int ntau, double beta, double h, int kt, int nomp) {
  pimpl->add_to_pp_self_energy(Sigma, tstp, ntau, beta, h, kt, nomp);
}

template<diagram_handler_type DHT>
void diagram_handler<DHT>::set_pp_self_energy(
  ppgfs_type & Sigma, int tstp, int ntau, double beta, double h, int kt, int nomp) {
  pimpl->set_pp_self_energy(Sigma, tstp, ntau, beta, h, kt, nomp);
}
  
// -----------------------------------------------------------------------

template<diagram_handler_type DHT>
bool diagram_handler<DHT>::has_symmetries() { return pimpl->has_symmetries(); }
template<diagram_handler_type DHT>
void diagram_handler<DHT>::clear_symmetries() { return pimpl->clear_symmetries(); }
template<diagram_handler_type DHT>
int diagram_handler<DHT>::total_number() { return pimpl->total_number(); }
template<diagram_handler_type DHT>
int diagram_handler<DHT>::reduced_number() { return pimpl->reduced_number(); }
template<diagram_handler_type DHT>
void diagram_handler<DHT>::print_report() { return pimpl->print_report(); }

// -----------------------------------------------------------------------

template<diagram_handler_type DHT>
void diagram_handler<DHT>::read_symmetries(input_archive_type & ia) {
  pimpl->read_symmetries(ia);
}

template<diagram_handler_type DHT>
void diagram_handler<DHT>::write_symmetries(output_archive_type & oa) {
  pimpl->write_symmetries(oa);
}
  
// -----------------------------------------------------------------------
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_DIAG_CPP
