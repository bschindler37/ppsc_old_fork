
#ifndef _PPSC_OCA_SIGMA_DISP_CPP
#define _PPSC_OCA_SIGMA_DISP_CPP

// -----------------------------------------------------------------------
//
// OCA pseudo particle self energy integrators
//
// Author: Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//
// -----------------------------------------------------------------------

#include "ppsc/timer.hpp"

#include "sigma.hpp"

#include "sigma_diag.hpp"
#include "sigma_mat.hpp"
#include "sigma_tstp.hpp"

#include "sigma_disp.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
namespace oca {
// -----------------------------------------------------------------------

using namespace ppsc;

// -----------------------------------------------------------------------
template<class DIAGT>
void sdiagram_dispatch_impl(int tstp, gf_tstp_type & ststp, DIAGT & diagram,
		       double beta, double h, int kt, int nomp) {

  // -- Point diagram timestep to given timestep
  new (& diagram.s) typename sdiagram_configuration<>::sigma_type(ststp);

  // -- Dispatch on tstp
  if(tstp == -1) {
    //{ Timer tmr("oca::sdiagram_mat");
    sdiagram_mat(diagram, beta, kt, nomp);
    //}
  } else {
    //{ Timer tmr("oca::sdiagram_tstp");
    sdiagram_tstp(tstp, diagram, beta, h, kt, nomp);
    //}
  } // tstp

  // -- Post-multiplication by diagram prefactor
  diagram.s *= diagram.prefactor;
  //std::cout << "--> OCA sigma diagram.prefactor = "
  //	    << diagram.prefactor << std::endl;
}

// -----------------------------------------------------------------------
void sdiagram_dispatch(
  int tstp, gf_tstp_type & ststp, sdiagram_configuration<> & diagram,
  double beta, double h, int kt, int nomp) {

  int sdim = diagram.s.nflavour();
  int adim = diagram.ga.nflavour();
  int bdim = diagram.gb.nflavour();
  int cdim = diagram.gc.nflavour();

  typedef std::array<int, 4> dims_type;
  dims_type dims = {sdim, adim, bdim, cdim};

  // -- Dispatch on dimensions
  // -- this can be made more beautiful but hey ... KISS

  // ---------------------------------------------------------------------
  if( dims == dims_type({1, 1, 1, 1}) ) {

    sdiagram_configuration<1, 1, 1, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({1, 2, 1, 2}) ){

    sdiagram_configuration<1, 2, 1, 2> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({2, 1, 2, 1}) ){

    sdiagram_configuration<2, 1, 2, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({2, 2, 2, 2}) ){

    sdiagram_configuration<2, 2, 2, 2> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  // -- THREE BAND DIAGRAMS
  // ---------------------------------------------------------------------

  /* // -- Inactive by default

  // ---------------------------------------------------------------------
  // -- 1, 3, 3, 3 all permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({1, 3, 3, 3}) ){

    sdiagram_configuration<1, 3, 3, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 1, 3, 3}) ){

    sdiagram_configuration<3, 1, 3, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 3, 1, 3}) ){

    sdiagram_configuration<3, 3, 1, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 3, 3, 1}) ){

    sdiagram_configuration<3, 3, 3, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  // -- 1, 3, 9, 3 all cyclic permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({1, 3, 9, 3}) ){

    sdiagram_configuration<1, 3, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 1, 3, 9}) ){

    sdiagram_configuration<3, 1, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);


  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 1, 3}) ){

    sdiagram_configuration<9, 3, 1, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);


    // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 3, 1}) ){

    sdiagram_configuration<3, 9, 3, 1> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);


  // ---------------------------------------------------------------------
  // -- 3, 3, 9, 9 all cyclic permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 3, 9, 9}) ){

    sdiagram_configuration<3, 3, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);


    // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 9, 3}) ){

    sdiagram_configuration<3, 9, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);


  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 3, 9}) ){

    sdiagram_configuration<9, 3, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 3, 3}) ){

    sdiagram_configuration<9, 9, 3, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  // -- 3, 9, 9, 9 all cyclic permutations
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({3, 9, 9, 9}) ){

    sdiagram_configuration<3, 9, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 3, 9, 9}) ){

    sdiagram_configuration<9, 3, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 3, 9}) ){

    sdiagram_configuration<9, 9, 3, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 9, 3}) ){

    sdiagram_configuration<9, 9, 9, 3> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  // ---------------------------------------------------------------------
  // -- 9, 9, 9, 9
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  } else if( dims == dims_type({9, 9, 9, 9}) ){

    sdiagram_configuration<9, 9, 9, 9> static_diagram(diagram); // cast from dynamic -> fixed size
    sdiagram_dispatch_impl(tstp, ststp, static_diagram, beta, h, kt, nomp);

  */

  // ---------------------------------------------------------------------
  } else {
    std::cout << "--> sdiagram_dispatch: Diagram dims not compiled" << std::endl;
    std::cout << "sdim, adim, bdim, cdim = "
	      << sdim << ", " << adim << ", "
	      << bdim << ", " << cdim << std::endl;
    exit(0);
  }
}

// -----------------------------------------------------------------------
} // end namespace oca
} // end namespace ppsc
// -----------------------------------------------------------------------

#endif // _PPSC_OCA_SIGMA_DISP_CPP
