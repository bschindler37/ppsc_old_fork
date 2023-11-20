
#ifndef _PPSC_HILBERT_SPACE_SINGLE_BAND_FERMI_DIAG_HPP		// helpful macro for debugging
#define _PPSC_HILBERT_SPACE_SINGLE_BAND_FERMI_DIAG_HPP

// -----------------------------------------------------------------------
//
// Single band fermionic hilbert space with colinear spin
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "hilbert_space_base.hpp"

// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hilbert_spaces {
// -----------------------------------------------------------------------

class single_band_fermi_diag : public hilbert_space_base {	// single_band_fermi_diag is a derived class of hilbert_space_base

public:
  
  typedef hilbert_space_base base_type;

  void init() {
    int norb = 1;
    int spin_degeneracy = 2;
    bool spin_conservation = true;
    base_type::init(norb, spin_degeneracy, spin_conservation);
  }
  
};

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_SINGLE_BAND_FERMI_DIAG_HPP
// -----------------------------------------------------------------------
