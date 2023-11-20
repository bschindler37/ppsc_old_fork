// -----------------------------------------------------------------------

#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// -----------------------------------------------------------------------

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

#include "./ppsc/hilbert_spaces/n_channel_one_spin.hpp"
#include "./ppsc/hamiltonians/n_channel_siam.hpp"

// -----------------------------------------------------------------------

using namespace std;

// -----------------------------------------------------------------------


#define NCHANNEL 1



class single_particle_greensfunction_type {

public:

  single_particle_greensfunction_type(int nt, int ntau,int nchannel) : nt(nt), ntau(ntau),nchannel(nchannel)
  {}
  void init_size() {
    gup.resize(nchannel);
    gdo.resize(nchannel);
    for(int m=0;m<nchannel;m++){
        gup[m]=cntr::herm_matrix<double>(nt,ntau,1,-1);
        gdo[m]=cntr::herm_matrix<double>(nt,ntau,1,-1);
    }
  }
  void update(int tstp, ppsc::gf_tstps_type & gf_tstps) {
    int idx=0;
    for(int m=0;m<nchannel;m++){
      this->gup[m].set_timestep(tstp, gf_tstps[idx]);
      idx++;
      this->gdo[m].set_timestep(tstp, gf_tstps[idx]);
      idx++;
    }
  }

  void store(hid_t file_id) {
    hid_t group_id;
    char flag[20];
    for(int m=0;m<nchannel;m++){
        sprintf(flag,"gup%d",m);
        group_id = create_group(file_id, flag);
        store_herm_greens_function(group_id, this->gup[m]);
        close_group(group_id);
        sprintf(flag,"gdo%d",m);
        group_id = create_group(file_id, flag);
        store_herm_greens_function(group_id, this->gdo[m]);
        close_group(group_id);
    }
  }

  void load(std::string filename) {
    char flag[20];
    for(int m=0;m<nchannel;m++){
        sprintf(flag,"gup%d",m);
        this->gup[m].read_from_hdf5(filename.c_str(), flag);
        sprintf(flag,"gdo%d",m);
        this->gdo[m].read_from_hdf5(filename.c_str(), flag);
    }
  }
  int nt, ntau,nchannel;
  vector<cntr::herm_matrix<double> > gup, gdo;
    
};

class hybridization_function_type {

public:
  
  hybridization_function_type(int nt, int ntau,int nchannel) :
    nt(nt), ntau(ntau),nchannel(nchannel)
  {}    
  void init_size() {
    Delta.resize(nchannel);
    Delta_cc.resize(nchannel);
    J.resize(nchannel);
    Jcc.resize(nchannel);
    kk.resize(nt+1);
    for(int m=0;m<nchannel;m++){
        Delta[m]=cntr::herm_matrix<double>(nt,ntau,1,-1);
        Delta_cc[m]=cntr::herm_matrix<double>(nt,ntau,1,-1);
        J[m]=cntr::function<double>(nt,1);
        Jcc[m]=cntr::function<double>(nt,1);
    }
  }
 
    void init_box(double beta,double h,double W,double omega) {
    // compute hyb ... to be replaced by box
    cntr::herm_matrix<double> tmp(nt, ntau, 1, -1);
    cntr::herm_matrix<double> DL(nt, ntau, 1, -1);
    cntr::herm_matrix<double> DR(nt, ntau, 1, -1);
    cntr::function<double> jtmp(nt,1);
    cntr::function<double> jcctmp(nt,1);
    cntr::smooth_box dos(-W,W,W*10.0);
    cntr::green_equilibrium(tmp,dos,beta,h,100,20);
    for(int m=0;m<nchannel;m++){
        DL=tmp;
        DR=tmp;
        jtmp=J[m];
        jcctmp=Jcc[m];
        for(int tstp=-1;tstp<=nt;tstp++){
            double cx=cos(kk[tstp+1]*cos(omega*tstp*h));
            double cy=sin(kk[tstp+1]*cos(omega*tstp*h));
            jtmp[tstp]=jtmp[tstp]*cdouble(cx,cy);
            jcctmp[tstp]=jcctmp[tstp]*cdouble(cx,-cy);
            DL.left_multiply(tstp,jtmp,1.0);
            DL.right_multiply(tstp,jcctmp,0.5);
            DR.left_multiply(tstp,jcctmp,1.0);
            DR.right_multiply(tstp,jtmp,0.5);
            Delta[m].set_timestep(tstp,DL);
            Delta[m].incr_timestep(tstp,DR,1.0);
            ppsc::set_bwd_from_fwd(tstp, Delta_cc[m], Delta[m]);
        }
    }
    }

  void store(hid_t file_id) {
    hid_t group_id;
    char flag[20];
    for(int m=0;m<nchannel;m++){
        sprintf(flag,"d%d",m);
        group_id = create_group(file_id, flag);
        store_herm_greens_function(group_id, this->Delta[m]);
        close_group(group_id);
        sprintf(flag,"dcc%d",m);
        group_id = create_group(file_id, flag);
        store_herm_greens_function(group_id, this->Delta_cc[m]);
        close_group(group_id);
    }
  }
  void load(std::string filename) {
    char flag[20];
    for(int m=0;m<nchannel;m++){
        sprintf(flag,"d%d",m);
        this->Delta[m].read_from_hdf5(filename.c_str(), flag);
        sprintf(flag,"dcc%d",m);
        this->Delta_cc[m].read_from_hdf5(filename.c_str(), flag);
    }
  }
  int nt,ntau,nchannel;
  vector<cntr::herm_matrix<double> > Delta,Delta_cc;
  vector<cntr::function<double> > J,Jcc;
  vector<double> kk;  // k(t)=kk[tstp+1], kk[0]=0
  double omega;
  
};




// -----------------------------------------------------------------------
// works only with HILB=n_channel_one_spin
template<class HILB>
ppsc::pp_ints_type get_pp_ints(hybridization_function_type & d,HILB & hil_) {

  // spin (u)p/(d)own and (c)reation/(a)nihilation operators
  int nchannel=hil_.nchannel_;
  int fermion=-1, fwd=+1, bwd=-1;
  ppsc::pp_ints_type pp_ints;
  for(int m=0;m<nchannel;m++){
      for(int s=0;s<2;s++){
          ppsc::operator_type ca = hil_.c_op_[hil_.flavor(m, s, 0)];
          ppsc::operator_type cc = hil_.c_op_[hil_.flavor(m, s, 1)];
          pp_ints.push_back(ppsc::pp_int_type(d.Delta[m],    cc, ca, fermion, fwd)); // spin up fwd
          pp_ints.push_back(ppsc::pp_int_type(d.Delta_cc[m], ca, cc, fermion, bwd)); // spin up bwd
     }
  }
  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::gf_verts_type get_gf_verts(HILB & hil_) {
  // compute all channel dependent "Greens functions",
  ppsc::gf_verts_type gf_verts;
  int nchannel=hil_.nchannel_;
  for(int m=0;m<nchannel;m++){
      for(int s=0;s<2;s++){
          ppsc::operator_type cua = hil_.c_op_[hil_.flavor(m, s, 0)];
          ppsc::operator_type cuc = hil_.c_op_[hil_.flavor(m, s, 1)];
          gf_verts.push_back(ppsc::gf_vert_type(0, 0, cua, cuc)); // spin up
       }
       // possibly use index ????
  }
  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {
  int itermax, iter_rtime, nt, ntau, kt, iter, tstp, order, nomp=1,nchannel,gf_out;

  int store_pp, read_eq_sym, read_rt_sym, read_state_from_file;
  
  double beta, h, errmax, dmfterr, dmfterr_equil,omega;
  bool matsubara_converged = false;
  
  try {
    // ---------------------------------------------------------------------
    // READ GENERAL INPUT (NOT YET MICROSCOPIC PARAMETERS)
    {
      if (argc < 2)
        throw("COMMAND LINE ARGUMENT MISSING");
      // scan the input file, double underscores to avoids mismatch
      find_param(argv[1], "__nt=", nt);
      find_param(argv[1], "__ntau=", ntau);
      find_param(argv[1], "__beta=", beta);
      find_param(argv[1], "__h=", h);
      find_param(argv[1], "__itermax=", itermax);
      find_param(argv[1], "__errmax=", errmax);
      find_param(argv[1], "__iter_rtime=", iter_rtime);
      find_param(argv[1], "__kt=", kt);
      // NOTE that there are no OCA diagrams for this problem (U=infty) !!
      //find_param(argv[1], "__order=", order);
      find_param(argv[1], "__store_pp=", store_pp);
      find_param(argv[1], "__read_eq_sym=", read_eq_sym);
      find_param(argv[1], "__read_rt_sym=", read_rt_sym);
      find_param(argv[1], "__gf_out=", gf_out);

      //if (order > 1) {
      //  cout << "using an OCA impurity solver" << endl;
      //  find_param(argv[1], "__nomp=", nomp);
      //} else {
      //nomp = 1;
      //}
    }
    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::n_channel_one_spin hilbert_space_type;
    typedef ppsc::hamiltonians::n_channel_siam<hilbert_space_type,NCHANNEL> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    //find_param(argv[1], "__nchannel=", nchannel);
    nchannel=NCHANNEL;
    //cout << "1..." << endl;
      
    hilbert_space_type hilbert_space;
    //cout << "2..." << endl;
    hilbert_space.init(nchannel);
    //cout << "3..." << endl;
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    //cout << "4..." << endl;
    { // -- Setup local Hamiltonian
      find_param_tvector(argv[1], "__epsf=", imp.hamiltonian.Epsf, nt);
      find_param_tvector(argv[1], "__bz=", imp.hamiltonian.Bz, nt);
      imp.update_hamiltonian();
    //cout << "5..." << endl;
    }
    // ---------------------------------------------------------------------
    // -- setup single particle greens functions and hybridizations
    // it will be Delta_{1/2} = J_{1/2}*Delta*J_{1/2}
    // since it is an impurity model, there is no self-consistency

    single_particle_greensfunction_type Gloc(nt,ntau,nchannel);
    Gloc.init_size();
    //cout << "6..." << endl;
    hybridization_function_type Delta(nt,ntau,nchannel);
    Delta.init_size();
    //cout << "7..." << endl;
    
    {
     // read in J:
      for(int m=0;m<nchannel;m++){
        char flag[20];
        vector<cdouble> jvector;
        sprintf(flag,"__J%d=",m);
        find_param_tvector(argv[1], flag, jvector, nt);
        for(int tstp=-1;tstp<=nt;tstp++){
            Delta.J[m][tstp]=jvector[tstp+1];
            //cout << flag << " at t= " << tstp << " = " << jvector[tstp+1] << endl;
            Delta.Jcc[m][tstp]=conj(jvector[tstp+1]);
        }
      }
      // read in kk
      find_param_tvector(argv[1],"__kk=",Delta.kk, nt);
        
       double cutoff;
       find_param(argv[1],"__D=",cutoff);
       find_param(argv[1],"__omega=",omega);
       Delta.init_box(beta,h,cutoff,omega);
       // Delta.init_bethe(beta,h);
    }

    Delta.Delta[0].print_to_file("d0_one.out");
    // ---------------------------------------------------------------------
    // MATSUBARA PART (EQUILIBRIUM INITIAL STATE)
    {
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      if(!read_state_from_file) {
	     imp.solve_atomic();
      } else {
	     gtmp.set_timestep(-1, Gloc.gup[0]); // this avoids one iteration
      }

      for (iter = 1; iter <= itermax; iter++) {
	// -- Construct interactions and verticies
	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	// -- Solve pseudo particle problem
	imp.update_diagrams(pp_ints, gf_verts);
	if(iter == 1 && read_eq_sym && !imp.has_symmetries()) {
	  std::string filename = "sym_eq.txt";
	  imp.read_symmetries(filename);
	}

	imp.pp_step(-1);

	// -- get spgf
	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(-1);
    Gloc.update(-1, gf_tstps);
	
    
	// -- Check error
	dmfterr_equil = cntr::distance_norm2(-1, gtmp, Gloc.gup[0]);
    gtmp.set_timestep(-1, Gloc.gup[0]);

	// -- Update Hybridization
	
	cout << "iter:  " << iter << " err: " << dmfterr_equil << endl;

        if (dmfterr_equil < errmax) {

	  if(!imp.has_symmetries()) {
	    imp.symmetry_reduction(-1);	    
	    std::string filename = "sym_eq.txt";
	    imp.write_symmetries(filename);
	  }
	  
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

    // ---------------------------------------------------------------------

    if(nt > 0) imp.clear_symmetries();

    if(read_rt_sym) {
      std::string filename = "sym_rt.txt";
      imp.read_symmetries(filename);
    }

    // ---------------------------------------------------------------------
    // 	START ... same iteration
    if (nt > 0 && matsubara_converged == true) {
      matsubara_converged = false;
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      imp.init_real_time();
      
      for (iter = 1; iter <= itermax; iter++) {

	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
	
	imp.update_diagrams(pp_ints, gf_verts);
	imp.pp_step(kt);

	for (int n = 0; n <= kt; n++) {
	  ppsc::gf_tstps_type gf_tstps = imp.get_spgf(n);
	  Gloc.update(n, gf_tstps);
	}
	
        dmfterr = cntr::distance_norm2(kt, gtmp, Gloc.gup[0]);
        gtmp.set_timestep(kt, Gloc.gup[0]);
	
        cout << "START: iter:  " << iter << " err: " << dmfterr << endl;
        if (dmfterr < errmax) {
          matsubara_converged = true;
          break;
        }
      }
    }


    // ---------------------------------------------------------------------
    // -- Reduce number of diagrams using symmetries

    if(matsubara_converged && nt > 0 && !imp.has_symmetries()) {
      imp.symmetry_reduction(kt);
      std::string filename = "sym_rt.txt";
      imp.write_symmetries(filename);
    }
    
    // ---------------------------------------------------------------------
    // 	REALTIME: ONLY FOR NT>0
    for (tstp = kt + 1; tstp <= nt; tstp++) {
      cout << "tstp:  " << tstp << endl;
      imp.extrapolate_timestep(tstp - 1);
      
      for (iter = 1; iter <= iter_rtime; iter++) {

	ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, imp.hamiltonian.hil_);
	ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);

	imp.update_diagrams(pp_ints, gf_verts);
	imp.pp_step(tstp);

	ppsc::gf_tstps_type gf_tstps = imp.get_spgf(tstp);
	Gloc.update(tstp, gf_tstps);
      }
    }

    // ---------------------------------------------------------------------
    // Kinetic energy

    //imp.calculate_pp_interaction_energy();
      
    //std::vector<double> Ekin_up = ppsc::get_kinetic_energy(Gloc_up, Delta_up, beta, h, kt);

    //std::vector<double> Ekin;
    //Ekin.resize(Ekin_up.size());
    
    //for( auto t : ppsc::range(0, Ekin.size()) ) Ekin[t] = 2*Ekin_up[t]; // two spin species

    // ---------------------------------------------------------------------
    // OBSERVABLES
    {
      
      //Gloc.gup[0].print_to_file("gup_C0.out");
        {
            FILE *out;
            out=fopen("obs.out","w");				
            for(int tstp=-1;tstp<=nt;tstp++){
                fprintf(out,"t: %d ",tstp);
                fprintf(out," n: %.10g",imp.hamiltonian.n_exp[tstp+1]);
                fprintf(out," m: %.10g",imp.hamiltonian.m_exp[tstp+1]);
                fprintf(out,"\n");
            }
            fclose(out);
        }
      if(gf_out){
      std:string filename = "data_ppsc.h5";
      hid_t file_id = open_hdf5_file(filename);
      hid_t group_id;

      imp.store(file_id, store_pp);
      Gloc.store(file_id);
      Delta.store(file_id);

      group_id = create_group(file_id, "bethe");
      //store_real_data_to_hid(group_id, "Ekin", Ekin.data(), Ekin.size());
      //store_real_data_to_hid(group_id, "Ekin_up", Ekin_up.data(), Ekin_up.size());
      store_double_attribute_to_hid(group_id, "dmfterr_equil", dmfterr_equil);
      close_group(group_id);
      close_hdf5_file(file_id);
      }
    }
    
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
