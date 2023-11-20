
#ifndef _PPSC_HILBERT_SPACE_THREE_BAND_FERMI_DIAG_HPP
#define _PPSC_HILBERT_SPACE_THREE_BAND_FERMI_DIAG_HPP

// -----------------------------------------------------------------------
//
// Three band fermionic hilbert space
//
// Author: P. Werner, philipp.werner@gmail.com (2016)
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

class three_band_fermi_diag : public hilbert_space_base {

public:
  
  typedef hilbert_space_base base_type;

  void init() {
    int norb = 3;
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
#endif // _PPSC_HILBERT_SPACE_THREE_BAND_FERMI_DIAG_HPP
// -----------------------------------------------------------------------
