#include <Eigen/Core>
#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <sys/stat.h>

using namespace std;

#ifndef CNTR_USE_OMP
#define CNTR_USE_OMP
#endif

#ifndef CNTR_USE_MPI
#define CNTR_USE_MPI
#endif

#define NCA_SOLVER_ASSERT_0 0
#define NCA_SOLVER_ASSERT_1 0

#include <cntr/cntr.hpp>
#include <cntr/utils/read_inputfile.hpp>

// -----------------------------------------------------------------------

#include "./ppsc/ppsc.hpp"
#include "./ppsc/solver.hpp"

// In ppsc: one needs to program hilbert spaces and hamiltonians, the solver is
// then build on that Here I include the relevant hilbert spaces and
// hamiltonian. see also documentation in these files NOTE: the hilbert space
// contains also an assumption of a block structure -- all pseudoparticle Gs and
// Sigmas will be block diagonal due to symmetry. the action ("interaction
// lines") bewlo must not mess up this blockstructure If they break the symmetry
// (e.g., if the Hamailtonian has a spin symmetry but there ius a term
// cup_dag*Delta*c_do in the action, then pseudo_Gs will not be diagonal in spin
// ... then one must choose a Hilbert space with a lower symmetry
//

#include "./ppsc/hamiltonians/single_band_hubbard.hpp"
#include "./ppsc/hilbert_spaces/single_band_fermi_diag.hpp"

// Function to display the progress bar
#include "progress_bar.hpp"

// in ppsc, we assume the action to a sum of terms like:
//
// S = - \int dt1 dt2 Vbar(t1) Delta(t1,t2) V(t2)
// where Vbar and V are operators in the hilbert space, and Delta is a
// contour-ordered functions; In multi orbtial systems with matrix Delta
// S = - \sum_{ij} Vbar_i(t1) * Delta_{ij}(t1,t2) * V_j(t2)
// each pair (ij) is a separate line .. see multi-orbitals examples ...
//
// the following function generates a list of all interaction terms
//
// -----------------------------------------------------------------------
template <class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type &Delta_up,
                               ppsc::gf_type &Delta_up_cc, HILB &hil_) {

  // spin (u)p/(d)own and (c)reation/(a)nihilation operators copied from the
  // hilbert space. ppsc::operator_type is some kind of blockmatrix, which maps
  // each sector of the hilbert space to precisely one sector

  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  int fermion = -1, fwd = +1, bwd = -1;
  //
  // note 1:
  // in ppsc, diagrams are generate only with forward pointing lines, so for
  // each term in the action line (unless it is symmetric) we generate TWO
  // interaction lines e.g., the term  cup_dag(t)*Delta(t,t')*cup(t') in the
  // action will generate lines
  //   v1bar(t)*A1(t,t')*v1(t'), v2bar(t)*A2(t,t')*v2(t') with
  //   v1bar = cup_dag, v1=cup, A1(t,t')=Delta(t,t') "forward"
  //   v2bar = cup, v2=cup_dag, A2(t,t')=Delta(t',t) "backward"
  // the reverse function A2(t,t')=Delta(t',t) will later have to be constructed
  // explicitly from Delta and stored in Delta_up_cc, but for this there are
  // some nice routines ...
  //
  // note 2:
  // ppsc::pp_int_type stores only a reference to Delta, but does not copy the
  // data
  ppsc::pp_ints_type pp_ints;
  pp_ints.push_back(
      ppsc::pp_int_type(Delta_up, cuc, cua, fermion, fwd)); // spin up fwd
  pp_ints.push_back(
      ppsc::pp_int_type(Delta_up_cc, cua, cuc, fermion, bwd)); // spin up bwd

  pp_ints.push_back(ppsc::pp_int_type(Delta_up, cdc, cda, fermion,
                                      fwd)); // spin do fwd // assuming spin sym
  pp_ints.push_back(ppsc::pp_int_type(Delta_up_cc, cda, cdc, fermion,
                                      bwd)); // spin do bwd // assuming spin sym

  return pp_ints;
}

// -----------------------------------------------------------------------
// get_gf_verts contains the information about correlation functions which
// should be computed If you later want to compute a contour-ordered correlation
// function C(t,t') = -i <TC A(t) B(t')>, then generate an entry
// ppsc::gf_vert_type(i, j, A, B) in the vector gf_verts
// i and j are just labels which can be used to write the result
// into the (i,j) component of a contour Green's function
template <class HILB> ppsc::gf_verts_type get_gf_verts(HILB &hil_) {

  // this function is already specific to the hilbert space used below
  // (single_band_fermi_diag) in the way in which operators are stored in hil_
  // single_band_fermi_diag has a list of operators with flavors
  // (orbital=0,spin,[annihilaton|creation])
  ppsc::operator_type cua = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  ppsc::operator_type cda = hil_.c_op_[hil_.flavor(0, 1, 0)];
  ppsc::operator_type cdc = hil_.c_op_[hil_.flavor(0, 1, 1)];

  ppsc::gf_verts_type gf_verts;
  // just the spin up GF will be measured, as there is spin symmetry
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, cua, cuc));
  return gf_verts;
}

/* ---------------------------------------------------------------------
Self-consisteny "Bethe with field":

in gerenal we start from a Hopping Hamiltonian
H0=\sum_{<ij>} J0 e^{I*A(t)*r_{ij}} c_i^dag c_j,  r_{ij}=r_i-r_j, J0=real

The hybridiazation function obtained by integrating out the bath:
Delta_{ii}(t,t') = J0^2 sum_{a,b} e^{I*A(t)*r_{ia}} G_{ab}^{[i]}(t,t')
e^{I*A(t')*r_{bi}}, with the cavity G^{[i]}, and the sum is over NN of i

For the Bethe lattice,
Delta_{ii}(t,t') = J0^2 sum_{a} e^{I*A(t)*r_{ia}} G_{aa}(t,t')
e^{I*A(t')*r_{ai}}, where we also assumed G_{aa}=G_{aa}^{[i]} due to inf
coordination number

We then assume that 1/2 of the NN are in the direction of the field
(r_{ia}=-r_{ai}=+1) and 1/2 are against (r_{ia}=-r_{ai}=-1), and G_{aa}==G due
to transl invariance:

Delta(t,t') = [DL(t,t')+DR(t,t')]/2
DL(t,t') = J0^2 e^{+I*A(t)} G(t,t') e^{-I*A(t')}
DR(t,t') = J0^2 e^{-I*A(t)} G(t,t') e^{+I*A(t')}
 ---------------------------------------------------------------------*/

void init_hybridization(int tstp, cntr::herm_matrix<double> &D,
                        cntr::function<double> &J, double eps, double mu,
                        double beta, double h) {
  // compute Delta for one non-inetracting bath site with energy level eps
  // -->  Delta(t,t') = J(t) g(t-t') J*(t)
  // input is J(t) = J * exp(i A(t))
  cdmatrix eps_matrix(1, 1);
  eps_matrix(0, 0) = eps;

  cntr::herm_matrix_timestep<double> g(tstp, D.ntau(), 1, -1);
  // initialize g
  cntr::green_from_H(tstp, g, mu, eps_matrix, beta, h);

  g.left_multiply(tstp, J, 1.0);
  g.right_multiply_hermconj(tstp, J, 1.0);
  D.set_timestep(tstp, g);
}

// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter, tstp, order, nomp;
  double J0_real, eps_bath;

  int store_pp, read_eq_sym, read_rt_sym, read_state_from_file, text_output,
      atomic_limit;
  double beta, h, errmax, dmfterr, dmfterr_equil, linear_mixing;
  bool matsubara_converged = false;
  string filename, out_prefix, logfile;
  try {
    // ---------------------------------------------------------------------
    // READ GENERAL INPUT (NOT YET MICROSCOPIC PARAMETERS)
    {
      if (argc < 3)
        throw("COMMAND LINE ARGUMENT MISSING");
      // scan the input file, double underscores to avoids mismatch
      // the meaning of some variables will be explained where they are used ...
      // e.g. find_param fooks for tyhe first line starting with  "__nt=" in the
      // input file argv[1] and interprets the second word in that line as the
      // value nt
      find_param(argv[1], "__nt=", nt);
      find_param(argv[1], "__ntau=", ntau);
      find_param(argv[1], "__beta=", beta);
      find_param(argv[1], "__h=", h);
      find_param(argv[1], "__itermax=", itermax);
      find_param(argv[1], "__errmax=", errmax);
      find_param(argv[1], "__iter_rtime=", iter_rtime);
      // note: second order (OCA) is very slow and can basically be used obnly
      // in few cases basically single orbital calculations or very short times
      // so order is always one ...
      find_param(argv[1], "__kt=", kt);
      find_param(argv[1], "__order=", order);
      // the following forces "full output" .. needed if output
      // should be used for restart a different simulation
      find_param(argv[1], "__store_pp=", store_pp);
      find_param(argv[1], "__linear_mixing=", linear_mixing);
      // symmetry reduction by Hugo Strand:
      // Block structure of Hilbert space is set by hand anyway.
      // in addition, however, some diagrams are identical because, e.g.,
      // spin_up=spin_do symmetry reduction tries to detect this automatically.
      // However, in particular for real time it may be that terms in H break
      // the symmetry only at later times, so if the algorithm does not detect
      // it at early times it may produce wrong results. so use woth case.
      // Anyway, symmetry reduction is mostly for OCA; for NCA main effor its
      // solution of pp Dyson, which is soved for all blocks anyway Hence, for
      // simplicity I have commented out all parts related to symmetry reduction
      // find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      // find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__read_state_from_file=", read_state_from_file);
      // omp is used for loop over pseudoparticle Dyson in different sectors
      // this is not an efficient paralellization (unless sectors are of similar
      // size) (e.g., sor SIAM, there are 4 equasl sectors, no nomp=4 could be
      // useful)
      find_param(argv[1], "__nomp=", nomp);
      find_param(argv[1], "__text_output=", text_output);
    }
    out_prefix = argv[2];
    logfile = out_prefix + "log.out";
    ofstream logf(logfile);

    // ---------------------------------------------------------------------
    // Main setup:
    // choose a hilbert space:
    typedef ppsc::hilbert_spaces::single_band_fermi_diag hilbert_space_type;
    // choose a hamiltonian which is build on that hilbert space
    typedef ppsc::hamiltonians::single_band_hubbard<hilbert_space_type>
        hamiltonian_type;
    // the solver is then build accordiung to the hamiltonian
    // (e.g. it inherits the blockstructure, and functions to initialize the
    // Hamiltonian)
    typedef ppsc::solver<hamiltonian_type> solver_type;

    // init tghe Hilbert space:
    // each hilbert space should have a function  init() which defines the
    // sectors, initializes a set of operators, etc.

    hilbert_space_type hilbert_space;
    hilbert_space.init();
    // set up the solver:
    // nt: number of real times steps [0...nt] (nt+1 points)
    // ntau:number of imag times steps [0...ntau] (ntau+1 points)
    // beta: inverse temperature
    // h: timestep
    // kt: solve order for the integration in Dyson equastion
    //     should alwyas be maximum (5)
    // nomp:omp paralellization: currentkly only OCS, should implement some
    // paralellization of multi-orbital NCA as well!
    //  hilbert_space: the hulber space ... must be previously initialized with
    //  init!
    // order: 1 for NCA (see above for comments on OCA)
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);

    {
      // -- Setup local Hamiltonian
      // the solver now contains an instance of the Hamiltonian
      // this also contains vectors that store time dependent parameters
      // (dsee explanation in the hamiltonian class)
      // in any case, these vectors have nt+2 entries e.g.:
      // imp.hamiltonian.U[0] is for the matsubara branch
      // imp.hamiltonian.U[i+1] = U(i) on timestep i on real times branch,
      // i=0...,nt
      find_param(argv[1], "__mu=", imp.hamiltonian.mu);
      // e.g. find_param_tvector fooks for the first line starting with  "__U="
      // in the input file argv[1]: if the second words starts with --, like
      // __U= --U.txt, then U.txt should be the input file cobntaining the nt+2
      // U values; otherwise, second word is interpreted as a number and U(t) is
      // set to a constant
      find_param_tvector(argv[1], "__U=", imp.hamiltonian.U, nt);
      // e.g., by overwriting the entry U[0] with a different values one could
      // implement a quench (will probably not be used ...)
      find_param(argv[1], "__U0=", imp.hamiltonian.U[0]);
      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.eps_up, nt);
      find_param_tvector(argv[1], "__eps=", imp.hamiltonian.eps_do, nt);
      find_param(argv[1], "__eps_bath=", eps_bath);
      find_param(argv[1], "__J0_real=", J0_real);
      // this function tells the solver to generate the Hamiltonian matrix
      // from the parameters which are now set.
      imp.update_hamiltonian();
    }
    // ---------------------------------------------------------------------
    // -- setup single particle greens functions and hybridizations
    // using of course NESSi

    cntr::herm_matrix<double> Gloc_up(nt, ntau, 1, -1);
    cntr::herm_matrix<double> Delta_up(nt, ntau, 1, -1);
    Delta_up.clear();
    cntr::herm_matrix<double> Delta_up_cc(nt, ntau, 1, -1);
    Delta_up_cc.clear();
    // read in vector potential and initialize the hopping with a Peierls phase:
    vector<double> vector_potential(nt + 2);
    find_param_tvector(argv[1], "__vector_potential=", vector_potential, nt);
    find_param(argv[1], "__atomic_limit=", atomic_limit);
    cntr::function<double> J(nt);
    cdouble J0(J0_real, 0.0);
    cdmatrix Jt(1, 1);
    Jt(0, 0) = J0;
    J.set_value(-1, Jt);
    for (int tstp = 0; tstp <= nt; tstp++) {
      // note: indexing in vector_potential shifted by 1, because
      // vector_potential[0] is imag branch value peierls = exp(ii*A(t))
      Jt(0, 0) = J0 * cdouble(cos(vector_potential[tstp + 1]),
                              sin(vector_potential[tstp + 1]));
      J.set_value(tstp, Jt);
    }

    // ---------------------------------------------------------------------

    if (read_state_from_file) {

      // the data can be restored form a previous simulation
      // (e.g., to have a better starting point for an iteration)
      // this works only if input files are consistent (same nt, ntau etc!)
      // also, data_ppsc.h5 must contain the full state (see description
      // at the output generation)
      std::string filename = "data_ppsc.h5";
      std::cout << "--> Reading state from file: " << filename << std::endl;

      Gloc_up.read_from_hdf5(filename.c_str(), "g");
      Delta_up.read_from_hdf5(filename.c_str(), "d");
      Delta_up_cc.read_from_hdf5(filename.c_str(), "dcc");
      imp.load(filename);
      // treatment of symmetries is described below
      // if(read_eq_sym) {
      //      std::string filename = "sym_eq.txt";
      //	      imp.read_symmetries(filename);
      //}

      imp.update_density_matrix(-1);
      imp.hamiltonian.update_exp_vals(-1, imp.rho);

      // std::cout << "pp_mu = " << imp.pp_mu << std::endl;
    }

    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      if (!read_state_from_file) {
        // generate the atomic limit solution for the paseudo Greensfunctions.
        imp.solve_atomic();
      } else {
        gtmp.set_timestep(-1, Gloc_up); // this avoids one iteration
      }
      cout << "Matsubara: maximum " << itermax << " iterations" << endl;

      for (iter = 1; iter <= itermax; iter++) {

        // -- Construct interactions and verticies
        // onbkly references to Delta_up and Delta_up_cc are stored
        ppsc::pp_ints_type pp_ints =
            get_pp_ints(Delta_up, Delta_up_cc, imp.hamiltonian.hil_);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

        // tell the solver which interaction lines to be taken, and which
        // correlation functions to be measured:
        imp.update_diagrams(pp_ints, gf_verts);

        // if(iter == 1 && read_eq_sym && !imp.has_symmetries()) {
        //   std::string filename = "sym_eq.txt";
        //  imp.read_symmetries(filename);
        // }
        //  the actual calulation:
        imp.pp_step(-1);

        // compute the correlation functions:
        // result written to a vector of cntr_herm_matrix_timeslice
        // here: gf_tstps[0] contains the impurity greensfunction at timeslice
        // -1
        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);

        // when we update G, we use a linear mixing with the parameter
        // G = linear_mixing * old + (1-linear_mixing) * new
        // -- linear mixing in Gloc
        ppsc::gf_tstp_type gloc_old(-1, ntau, 1);
        ppsc::gf_tstp_type gloc_mix(-1, ntau, 1);
        gloc_mix.clear();
        Gloc_up.get_timestep(-1, gloc_old);
        gloc_mix.incr(gloc_old, linear_mixing);
        gloc_mix.incr(gf_tstps[0], 1.0 - linear_mixing);
        Gloc_up.set_timestep(-1, gloc_mix);

        // -- Check error
        dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc_up);
        gtmp.set_timestep(-1, Gloc_up);

        // -- Update Hybridization: Delta from one nonint. bath site
        // compute at timestep tstp=-1
        if (!atomic_limit) {
          init_hybridization(-1, Delta_up, J, eps_bath, imp.hamiltonian.mu,
                             beta, h);
        }
        // now we have to construct the function with reversed time arguments:
        // Delta_up_cc(t,t')=Delta_up(t',t)
        ppsc::set_bwd_from_fwd(-1, Delta_up_cc, Delta_up);
        logf << "iter:  " << iter << " err: " << dmfterr_equil << endl;

        displayProgressBar(iter, itermax);

        if (dmfterr_equil < errmax) {

          // if(!imp.has_symmetries()) {
          //     imp.symmetry_reduction(-1);
          //     std::string filename = "sym_eq.txt";
          //     imp.write_symmetries(filename);
          //   }
          matsubara_converged = true;
          break;
        }
      }
      if (iter > itermax) {
        cerr << "WARNING: Matsubara not converged  after " << itermax
             << "steps ... abort" << endl;
        cerr << "skip real-time calculation " << endl;
      }
    }
    cout << endl; // close the progress bar

    // ---------------------------------------------------------------------

    // if(nt > 0) imp.clear_symmetries();

    // if(read_rt_sym) {
    // std::string filename = "sym_rt.txt";
    // imp.read_symmetries(filename);
    // }
    //  ---------------------------------------------------------------------
    //  	START ...
    //   the self-consistemnt equastions have to be solved bby a "global
    //   iteration" on timesteps [0...kt]
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      cout << "START: maximum " << itermax << " iterations" << endl;
      for (iter = 1; iter <= itermax; iter++) {
        // again update list of interaction lines and vcertices ...
        // this would not have to be done again, as it has not changed
        ppsc::pp_ints_type pp_ints =
            get_pp_ints(Delta_up, Delta_up_cc, imp.hamiltonian.hil_);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
        // the calculations ...
        // when pp_step(tstp) is called with timstep tstp=kt,
        // pseudo=particle Gs are determined for times 0...kt
        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(kt);
        // computation of Greens fumnctions and Delta:
        // this must be done for all timesteps 0...kt which are updated in this
        // iteration:
        for (int n = 0; n <= kt; n++) {
          ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
          Gloc_up.set_timestep(n, gf_tstps[0]);
          // update Hybridization
          if (!atomic_limit) {
            init_hybridization(n, Delta_up, J, eps_bath, imp.hamiltonian.mu,
                               beta, h);
          }

          // Delta with reversed arguments
          ppsc::set_bwd_from_fwd(n, Delta_up_cc, Delta_up);
        }
        dmfterr = cntr::distance_norm2(kt, gtmp, Gloc_up);
        gtmp.set_timestep(kt, Gloc_up);

        logf << "START: iter:  " << iter << " err: " << dmfterr << endl;
        displayProgressBar(iter, itermax);
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
      cout << endl; // close the progress bar
    }
    // ---------------------------------------------------------------------
    // -- Reduce number of diagrams using symmetries
    // if(matsubara_converged && nt > 0 && !imp.has_symmetries()) {
    //  imp.symmetry_reduction(kt);
    //  std::string filename = "sym_rt.txt";
    //  imp.write_symmetries(filename);
    //}
    // ---------------------------------------------------------------------
    // 	REALTIME: ONLY FOR NT>0
    // here a fixed number of iterations is done at each timestep
    // to update information at given timestep
    cout << "REALTIME: nt= " << nt << " steps " << endl;
    for (tstp = kt + 1; tstp <= nt; tstp++) {
      cntr::herm_matrix_timestep<double> gtmp(tstp, ntau, 1, -1);
      imp.extrapolate_timestep(tstp - 1);
      for (iter = 1; iter <= iter_rtime; iter++) {
        ppsc::pp_ints_type pp_ints =
            get_pp_ints(Delta_up, Delta_up_cc, imp.hamiltonian.hil_);
        ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

        imp.update_diagrams(pp_ints, gf_verts);
        imp.pp_step(tstp);

        ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
        Gloc_up.set_timestep(tstp, gf_tstps[0]);
        if (!atomic_limit) {
          init_hybridization(tstp, Delta_up, J, eps_bath, imp.hamiltonian.mu,
                             beta, h);
        }
        ppsc::set_bwd_from_fwd(tstp, Delta_up_cc, Delta_up);

        dmfterr = cntr::distance_norm2(tstp, gtmp, Gloc_up);
        gtmp.set_timestep(tstp, Gloc_up);
        logf << "tstp:  " << tstp << " iter: " << iter << " err: " << dmfterr
             << endl;
      }
      displayProgressBar(tstp, nt);
    }
    cout << endl; // close the progress bar
    // computing different observables:
    // ---------------------------------------------------------------------
    // other local exp values have been calculated in pp_step and are stored in
    // Hamiltonian [X] = Q=pesudoparticle partition function,nup,ndo,docc,Eint)
    // are in imp.hamiltonian.[X]_exp
    if (text_output) {
      filename = out_prefix + "Gloc.out";
      Gloc_up.print_to_file(filename.c_str());
      filename = out_prefix + "obs.out";
      ofstream ofile(filename);
      ofile.precision(14);
      ofile << "# tstp <nup> <ndo> <docc> <Eint> <Q> A(t) E(t)" << endl;
      double intEj = 0.0;
      double Et = 0;
      for (int tstp = -1; tstp <= nt; tstp++) {
        ofile << tstp << " ";
        ofile << imp.hamiltonian.nu_exp[tstp + 1] << " ";
        ofile << imp.hamiltonian.nd_exp[tstp + 1] << " ";
        ofile << imp.hamiltonian.docc_exp[tstp + 1] << " ";
        ofile << imp.hamiltonian.Eint_exp[tstp + 1] << " ";
        ofile << imp.hamiltonian.Q_exp[tstp + 1] << " ";
        ofile << vector_potential[tstp + 1] << " ";
        ofile << endl;
      }
      ofile.close();
    }

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      filename = out_prefix + "data_ppsc.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;

      imp.store(file_id, store_pp);

      group_id = create_group(file_id, "g");
      store_herm_greens_function(group_id, Gloc_up);
      close_group(group_id);

      group_id = create_group(file_id, "d");
      store_herm_greens_function(group_id, Delta_up);
      close_group(group_id);

      group_id = create_group(file_id, "dcc");
      store_herm_greens_function(group_id, Delta_up_cc);
      close_group(group_id);
      group_id = create_group(file_id, "scalars");
      store_real_data_to_hid(group_id, "A", vector_potential.data(),
                             vector_potential.size());
      store_double_attribute_to_hid(group_id, "error", dmfterr_equil);
      close_group(group_id);
      close_hdf5_file(file_id);
    }

    logf.close();

  } // try
  catch (char *message) {
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << "CDMFT input_file [ --test ]\n" << endl;
  } catch (...) {
    cerr << "unspecified exception " << endl;
    cerr << "\nCDMFT input_file [ --test ]\n" << endl;
  }
  return 0;
}
