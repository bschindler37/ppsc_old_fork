
#ifndef _PPSC_HILBERT_SPACE_TWO_BAND_FERMI_DIAG_HPP
#define _PPSC_HILBERT_SPACE_TWO_BAND_FERMI_DIAG_HPP

// -----------------------------------------------------------------------
//
// Two band fermionic hilbert space
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

class two_band_fermi_diag : public hilbert_space_base {

public:
  
  typedef hilbert_space_base base_type;

  void init() {
    int norb = 2;
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
#endif // _PPSC_HILBERT_SPACE_TWO_BAND_FERMI_DIAG_HPP
// -----------------------------------------------------------------------
