
#ifndef _PPSC_OCA_GF_DISP_CPP
#define _PPSC_OCA_GF_DISP_CPP

// -----------------------------------------------------------------------
//
// OCA single-particle Green's function integrator
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "gf.hpp"

#include "gf_diag.hpp"
#include "gf_mat.hpp"
#include "gf_tstp.hpp"

#include "gf_disp.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// -- Dispatch on equilibrium / real-time
  
template<class DIAGT>
void gdiagram_dispatch_impl(int tstp, gf_tstp_type & gtstp, DIAGT & diagram,
		       double beta, double h, int kt, int nomp) {

  // -- Point diagram timestep to given timestep
  new (& diagram.g) typename gdiagram_configuration<>::g_type(gtstp);

  // -- Dispatch on tstp
  if(tstp == -1) {
    //{ Timer tmr("oca::gdiagram_mat");
    gdiagram_mat(diagram, beta, kt, nomp);
    //}
  } else {
    //{ Timer tmr("oca::gdiagram_tstp");
    gdiagram_tstp(tstp, diagram, beta, h, kt, nomp);
    //}
  }

  // -- Post-multiplication by diagram prefactor
  diagram.g *= diagram.prefactor;
  //std::cout << "--> OCA gf diagram.prefactor = "
  //	    << diagram.prefactor << std::endl;
  
}

// -----------------------------------------------------------------------
// -- Dynamic dispatch on dynamically sized diagrams
// -- calling statically sized diagram integrators ( dyn -> static conv)
  
void gdiagram_dispatch(
  int tstp, gf_tstp_type & gtstp, gdiagram_configuration<> & diagram,
  double beta, double h, int kt, int nomp) {

  int adim = diagram.ga.nflavour();
  int bdim = diagram.gb.nflavour();
  int cdim = diagram.gc.nflavour();
  int ddim = diagram.gd.nflavour();

  typedef std::array<int, 4> dims_type;
  dims_type dims = {adim, bdim, cdim, ddim};

  // -- Dispatch on dimensions
  // -- this can be made more beautiful but hey ... KISS
  
  // ---------------------------------------------------------------------
  if( dims == dims_type({1, 1, 1, 1}) ) {
    
    gdiagram_configuration<1, 1, 1, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({1, 2, 1, 2}) ){

    gdiagram_configuration<1, 2, 1, 2> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);
    
  // ---------------------------------------------------------------------
  } else if( dims == dims_type({2, 1, 2, 1}) ){

    gdiagram_configuration<2, 1, 2, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({2, 2, 2, 2}) ){

    gdiagram_configuration<2, 2, 2, 2> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  // -- THREE BAND DIAGRAMS
  // ---------------------------------------------------------------------

  /* // -- Inactive by default
    
  // ---------------------------------------------------------------------
  // -- 1, 3, 3, 3 all permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({1, 3, 3, 3}) ){

    gdiagram_configuration<1, 3, 3, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 1, 3, 3}) ){

    gdiagram_configuration<3, 1, 3, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 3, 1, 3}) ){

    gdiagram_configuration<3, 3, 1, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 3, 3, 1}) ){

    gdiagram_configuration<3, 3, 3, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);
    
  // ---------------------------------------------------------------------
  // -- 3, 9, 9, 9 all permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 9, 9}) ){

    gdiagram_configuration<3, 9, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 9, 9}) ){

    gdiagram_configuration<9, 3, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 3, 9}) ){

    gdiagram_configuration<9, 9, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 9, 3}) ){

    gdiagram_configuration<9, 9, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  // -- 1, 3, 9, 3 all cyclic permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({1, 3, 9, 3}) ){

    gdiagram_configuration<1, 3, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 1, 3, 9}) ){

    gdiagram_configuration<3, 1, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 1, 3}) ){

    gdiagram_configuration<9, 3, 1, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 3, 1}) ){

    gdiagram_configuration<3, 9, 3, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);
    
  // ---------------------------------------------------------------------
  // -- 3, 3, 9, 9 all cyclic permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 3, 9, 9}) ){

    gdiagram_configuration<3, 3, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 3, 9}) ){

    gdiagram_configuration<9, 3, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 3, 3}) ){

    gdiagram_configuration<9, 9, 3, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 9, 3}) ){

    gdiagram_configuration<3, 9, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);
    
  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 3, 9}) ){

    gdiagram_configuration<3, 9, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 9, 3}) ){

    gdiagram_configuration<9, 3, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);
    
  // ---------------------------------------------------------------------
  // -- 9, 9, 9, 9
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 9, 9}) ){

    gdiagram_configuration<9, 9, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    gdiagram_dispatch_impl(tstp, gtstp, static_diagram, beta, h, kt, nomp);
    
  */
    
  // ---------------------------------------------------------------------
  } else {
    std::cout << "--> gdiagram_dispatch: Diagram dims not compiled" << std::endl;
    std::cout << "adim, bdim, cdim, ddim = "
	      << adim << ", " << bdim << ", "
	      << cdim << ", " << ddim << std::endl;
    exit(0);
  }
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_GF_DISP_CPP
