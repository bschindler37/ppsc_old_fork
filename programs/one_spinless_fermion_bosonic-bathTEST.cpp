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

#include "./ppsc/hamiltonians/single_band_hubbard.hpp"

// -----------------------------------------------------------------------

using namespace std;


#include <cassert>

// --------------------------------------------------------------------------------------
const char *string_hyborder(int order){
	std::ostringstream ostr;
	if(order==1) ostr << "nca";
	else {
		if(order==2) ostr << "oca";
		else { ostr << "order" << order; }
	}
	return (ostr.str()).c_str();
}

double n_bose(double e, double beta){
	return 1.0/(std::exp(beta*e)-1.0);
}

// calculate a(t)=h(t)*f(s) for t=0...tmax
void const_mult(cdouble *a, cdouble *h, cdouble *f, int s, int tmax){
	for(int n=0;n<=tmax;n++) { a[n]= h[n]*f[s]; }
}

// calculate a(t)=h(t)*f(t) for t=0...tmax
void time_mult(cdouble *a, cdouble *h, cdouble *f, int tmax){
	for(int n=0;n<=tmax;n++) { a[n]= h[n]*f[n]; }
}

// calculate a(t) = dt * int_0^t h(t-s) * f(s) ds for t=0...tmax (via trapezoidal integration)
void conv(cdouble *a, cdouble *h, cdouble *f, int tmax, double dt){
	cdouble trapz_int;
	a[0]=0.0; //set integral of measure zero equal to zero
	for(int n=1;n<=tmax;n++) {
		trapz_int=0.0;
		for(int m=1;m<n;m++) { trapz_int+=h[n-m]*f[m] ; }
		trapz_int+=(h[n]*f[0]+h[0]*f[n])/2.0;
		a[n]=dt*trapz_int;
	}
}

// calculate h * int_0^t h(t-s) * f(s) ds for given t<=tmax (via trapezoidal integration)
cdouble conv_timestep(int t, cdouble *h, cdouble *f, int tmax, double dt){
	cdouble trapz_int=0.0;
	for(int m=1;m<t;m++) { trapz_int+=h[t-m]*f[m] ; }
	trapz_int+=(h[t]*f[0]+h[0]*f[t])/2.0;
	return dt*trapz_int;
}


// determine value of oca diagram for time 'tstp' (assuming time translational invariance)
cdouble oca_diag(cdouble *gret, cdouble *f1, cdouble *f2, cdouble *f3, cdouble *f4, int tstp, double h) {
	cdouble cdiag, aux1[tstp+1], aux2[tstp+1];
	
	// constant multiplication: gret(t)*f1(0) for all t<=tstp
	const_mult(aux1,gret,f1,0,tstp);
	// time-wise multiplication: x(t)*f2(t) for all t<=tstp
	time_mult(aux2,aux1,f2,tstp);
	// convolution: h * int_0^t gret(t-s) * x(s) ds for all t<=tstp
	conv(aux1,gret,aux2,tstp,h);
	// time-wise multiplication: x(t)*f3(t) for all t<=tstp
	time_mult(aux2,aux1,f3,tstp);
	// convolution: h * int_0^t gret(t-s) * x(s) ds for _fixed_ t=tstp
	cdiag=conv_timestep(tstp,gret,aux2,tstp,h);
	
	return cdiag*f4[tstp]; // last, point-wise multiplication
}

cdouble diagrams_order_one(cdouble* gret, int tstp, double w0, double beta, double g, double h) {
	double nbose=n_bose(w0, beta);
	//return g*g* ( (1.0+nbose)*std::exp(-II*h*cdouble(tstp)*w0) + nbose*std::exp(II*h*cdouble(tstp)*w0) ) * gret[tstp] ;
	return g*g* ( (2.0*nbose+1.0)*std::cos(h*double(tstp)*w0) - II*std::sin(h*double(tstp)*w0) ) * gret[tstp];
	//factor 2=2^(order) due to indistinghuishability of hybridization line direction
}

cdouble diagrams_order_two(cdouble* gret, int tstp, double w0, double beta, double g, double h) {
	double nbose=n_bose(w0, beta), prefac=g*g/2.0;
	cdouble sigma_ret_tstp, fvert1[tstp+1], fvert2[tstp+1], fvert1_conj[tstp+1], fvert2_conj[tstp+1];
	// calculate vertex functions
	for(int n=0;n<=tstp;n++){
		fvert1[n]=std::sqrt(nbose)*std::exp(II*h*cdouble(n)*w0); fvert2[n]=std::sqrt(nbose+1)*std::exp(II*h*cdouble(n)*w0);
		fvert1_conj[n]=std::conj(fvert1[n]); fvert2_conj[n]=std::conj(fvert2[n]);
	}
	// sigma_ret_tstp is the sum of four diagram expressions; all of which have the same form but different vertex functions	
	sigma_ret_tstp=oca_diag(gret,fvert2,fvert2,fvert2_conj,fvert2_conj,tstp,h) + oca_diag(gret,fvert1_conj,fvert2,fvert1,fvert2_conj,tstp,h) + oca_diag(gret,fvert2,fvert1_conj,fvert2_conj,fvert1,tstp,h) + oca_diag(gret,fvert1_conj,fvert1_conj,fvert1,fvert1,tstp,h); 

	return 4.0*prefac*prefac*sigma_ret_tstp; //factor 4=2^(order) due to indistinghuishability of hybridization line direction
	
}

// update self energy at timestep 'tstp' for diagrams up to order 'order' for bosonic propagator with resonance freq. w0 & hybridization strength g
void update_sigma(GREEN &sigma, cntr::herm_pseudo<double> &G, int tstp, int order, double w0, double beta, double g, double h) {
	assert(sigma.nt()==G.nt());
	cdouble *g_ret, sigma_ret_tstp;
	g_ret = new cdouble[G.nt()+1];
	for(int j=0; j<= G.nt(); j++) { G.get_ret(j,0,g_ret[j]); }
	
	if(order>=1) sigma_ret_tstp = diagrams_order_one(g_ret,tstp,w0,beta,g,h);
	if(order>=2) sigma_ret_tstp += diagrams_order_two(g_ret,tstp,w0,beta,g,h);
	// ... and so on ...
	
	for(int n=tstp; n<=sigma.nt(); n++){ sigma.set_ret(n,n-tstp,sigma_ret_tstp); } 	// update (time-translational invariant) self energy for new value 'sigma_ret_tstp'
	
	delete[] g_ret;
}

// update self energy at timestep 'tstp' for diagrams up to order 'order' for bosonic propagator with resonance freq. w0 & hybridization strength g
void update_sigma(GREEN &sigma, GREEN &G, int tstp, int order, double w0, double beta, double g, double h) {
	assert(sigma.nt()==G.nt());
	cdouble *g_ret, sigma_ret_tstp;
	g_ret = new cdouble[G.nt()+1];
	for(int j=0; j<= G.nt(); j++) { G.get_ret(j,0,g_ret[j]); }
	
	if(order>=1) sigma_ret_tstp = diagrams_order_one(g_ret,tstp,w0,beta,g,h);
	if(order>=2) sigma_ret_tstp += diagrams_order_two(g_ret,tstp,w0,beta,g,h);
	// ... and so on ...
	
	for(int n=tstp; n<=sigma.nt(); n++){ sigma.set_ret(n,n-tstp,sigma_ret_tstp); } 	// update (time-translational invariant) self energy for new value 'sigma_ret_tstp'
	
	delete[] g_ret;
}
	
// update self energy at timestep 'tstp' for diagrams up to order 'order' for bosonic propagator with resonance freq. w0 & hybridization strength g
void update_sigma(cntr::herm_pseudo<double> &sigma, cntr::herm_pseudo<double> &G, int tstp, int order, double w0, double beta, double g, double h) {
	assert(sigma.nt()==G.nt());
	cdouble *g_ret, sigma_ret_tstp;
	g_ret = new cdouble[G.nt()+1];
	for(int j=0; j<= G.nt(); j++) { G.get_ret(j,0,g_ret[j]); }
	
	if(order>=1) sigma_ret_tstp = diagrams_order_one(g_ret,tstp,w0,beta,g,h);
	if(order>=2) sigma_ret_tstp += diagrams_order_two(g_ret,tstp,w0,beta,g,h);
	// ... and so on ...
	
	for(int n=tstp; n<=sigma.nt(); n++){ sigma.set_ret(n,n-tstp,sigma_ret_tstp); } 	// update (time-translational invariant) self energy for new value 'sigma_ret_tstp'
	
	delete[] g_ret;
}
	

// oscillation in time of ppGF due to pp_mu != 0
cdouble pp_oscillation(double beta, double eps, double mu_mats, int t1, int t2, double dt) {
    double x = 1.0 + std::exp(-beta*(eps-mu_mats));
    return std::pow(x,-II*dt*cdouble(t1-t2)/beta);
}
	
// oscillation in time of ppGF due to pp_mu != 0
cdouble pp_oscillation(double pp_mu, int t1, int t2, double dt) {
    return std::exp(+II*dt*cdouble(t1-t2)*pp_mu);
}
	
	
	
	
	








// -----------------------------------------------------------------------

#include "./ppsc/hilbert_spaces/hilbert_space_base.hpp"

namespace ppsc {
namespace hilbert_spaces {
  
class one_spinless_fermion : public hilbert_space_base {

public:
  
  typedef hilbert_space_base base_type;

  void init() {
    int norb = 1;
    int spin_degeneracy = 1;
    bool spin_conservation = true;
    base_type::init(norb, spin_degeneracy, spin_conservation);
  }
  
};

} // namespace hilbert_spaces
} // namespace ppsc

// -----------------------------------------------------------------------

namespace ppsc {
namespace hamiltonians {

// -----------------------------------------------------------------------
template<class HILB> class one_spinless_fermion {

public:

  // ---------------------------------------------------------------------
  one_spinless_fermion(int nt, HILB & hilbert_space) :
    nt(nt), mu(nt+2), Q_exp(nt+2),
    n_exp(nt+2), Eint_exp(nt+2), hil_(hilbert_space) {

    // -- Construct basic operators

    ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
    cc = hil_.c_op_[hil_.flavor(0, 0, 1)];

    n = cc * ca;

    Q = operator_type::Identity(hil_);
  }

  // ---------------------------------------------------------------------
  void get_hamiltonian(int tstp, operator_type & Htemp) {
    double mu_t = mu[tstp+1];
    Htemp = (eps-mu_t) * n;
    
  }

  // ---------------------------------------------------------------------
  void update_exp_vals(int tstp, operators_type & rho) {

    Q_exp[tstp + 1] = expectation_value(tstp, rho, Q).real();
    n_exp[tstp + 1] = expectation_value(tstp, rho, n).real();

    operator_type Ht;
    get_hamiltonian(tstp, Ht);
    Eint_exp[tstp + 1] = expectation_value(tstp, rho, Ht).real();
  }

  // ---------------------------------------------------------------------
  void local_obs(int tstp, double & n_out) {
    n_out = n_exp[tstp + 1];
  }

  // ---------------------------------------------------------------------
  void store(hid_t group_id) {

    store_real_data_to_hid(group_id, "mu", mu.data(), mu.size());
    store_double_attribute_to_hid(group_id, "eps", eps);

    store_real_data_to_hid(group_id, "Q_exp", Q_exp.data(), Q_exp.size());
    store_real_data_to_hid(group_id, "n_exp", n_exp.data(), n_exp.size());
    store_real_data_to_hid(group_id, "Eint_exp", Eint_exp.data(), Eint_exp.size());
  }

  int nt;
  HILB hil_;

  double eps; 	// need to be initialized manually (!)
  std::vector<double> mu; // need to be initialized manually (!)
  std::vector<double> Q_exp, n_exp, Eint_exp;

  operator_type ca, cc, n, Q;
};

} // end namespace hamiltonians
} // end namespace ppsc

// -----------------------------------------------------------------------
template<class HILB>
ppsc::pp_ints_type get_pp_ints(ppsc::gf_type & Delta,
			       ppsc::gf_type & Delta_cc,
			       HILB & hil_) {

  // (c)reation/(a)nihilation operators

  ppsc::operator_type ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cc = hil_.c_op_[hil_.flavor(0, 0, 1)];
  
  ppsc::operator_type n = cc * ca;
  
  int fermion=-1, boson=+1, fwd=+1, bwd=-1;

  ppsc::pp_ints_type pp_ints;
  
  pp_ints.push_back(ppsc::pp_int_type(Delta,    n, n, boson, fwd));
  pp_ints.push_back(ppsc::pp_int_type(Delta_cc, n, n, boson, bwd));
  
  return pp_ints;
}

// -----------------------------------------------------------------------
template<class HILB>
ppsc::gf_verts_type get_gf_verts(HILB & hil_) {

  ppsc::operator_type ca = hil_.c_op_[hil_.flavor(0, 0, 0)];
  ppsc::operator_type cc = hil_.c_op_[hil_.flavor(0, 0, 1)];

  ppsc::gf_verts_type gf_verts;
  gf_verts.push_back(ppsc::gf_vert_type(0, 0, ca, cc));
  return gf_verts;
}
 
// -----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  int nt, ntau, kt=5, order, nomp, itermax = 400, iter;
  double mu=0.0, beta, h, eps, gamma, gammatilde, Omega0, err=0.0, err2=0.0, errmax = 1e-10;
  CFUNC eps_fct;
  cdmatrix eps0(1,1);
  bool store_pp=true;	// enables storing of pp-gf and pp-selfenergy 


    // ---------------------------------------------------------------------
    // READ GENERAL INPUT (NOT YET MICROSCOPIC PARAMETERS)
    
    if (argc < 2) {
      std::cerr << "COMMAND LINE ARGUMENT (INPUT FILE) MISSING. ABORT PROGRAM..." << std::endl;
      return 1;
    }
      // scan the input file, double underscores to avoids mismatch
    find_param(argv[1], "__nt=", nt);
    find_param(argv[1], "__ntau=", ntau);
    find_param(argv[1], "__beta=", beta);
    find_param(argv[1], "__h=", h);
    find_param(argv[1], "__Omega0=", Omega0);
    find_param(argv[1], "__gamma=", gamma);
    find_param(argv[1], "__nomp=", nomp);
    find_param(argv[1], "__order=", order); // order = 1 for nca and = 2 for oca
          
    gammatilde = gamma*std::sqrt(2*Omega0);
    std::string order_string;
    if (order > 2){
      std::cerr << "ORDER MUST BE 1 (NCA) OR 2 (OCA). ABORT PROGRAM..." << std::endl;
      return 1;
    }
    if (order == 2)
    	order_string = "oca";
    else 
    	order_string = "nca";
    	
    // ---------------------------------------------------------------------
    // -- Setup pp calculator

    typedef ppsc::hilbert_spaces::one_spinless_fermion hilbert_space_type;
    typedef ppsc::hamiltonians::one_spinless_fermion<hilbert_space_type> hamiltonian_type;
    typedef ppsc::solver<hamiltonian_type> solver_type;
    
    hilbert_space_type hilbert_space;
    hilbert_space.init();
    
    solver_type imp(nt, ntau, beta, h, kt, nomp, hilbert_space, order);
    
    // -- Setup local Hamiltonian

    find_param(argv[1], "__mu_mats=", imp.hamiltonian.mu[0]);
    find_param(argv[1], "__eps=", imp.hamiltonian.eps);
    
    for ( int tstp = 0; tstp <= nt ; tstp++ )  { imp.hamiltonian.mu[tstp+1] = 0.0; }
    imp.update_hamiltonian();
    
    eps0(0,0) = imp.hamiltonian.eps;
    eps_fct = CFUNC(nt,1);
    for(int tstp=-1;tstp<=nt;tstp++) { eps_fct.set_value(tstp,eps0); }
    
    // ---------------------------------------------------------------------
    // -- setup single particle Greens functions
    
    cntr::herm_matrix<double> Sigma(nt,ntau,1,-1);
    cntr::herm_pseudo<double> G_pseudo(nt,ntau,1,-1);
    cntr::herm_matrix<double> Gloc(nt, ntau, 1, -1);    
    cntr::herm_matrix<double> Delta(nt, ntau, 1, +1);
    cntr::herm_matrix<double> Delta_cc(nt, ntau, 1, +1);
    
    // ---------------------------------------------------------------------
    // -- Setup Hybridization function, interaction vertices and green's function vertices

    for(int tstp = -1; tstp <= nt; tstp++) {	
      cntr::green_single_pole_XX_timestep(tstp, Delta, Omega0, beta, h);	// free bosonic propagator [-i < T_C X(t) X(t') >] for one timestep
      //Delta.smul(tstp, (gammatilde*gammatilde) / (2.0*Omega0) ); // line below is equivalent to this line
      Delta.smul(tstp, gamma*gamma );
      ppsc::set_bwd_from_fwd(tstp, Delta_cc, Delta);
    }
      
    ppsc::pp_ints_type pp_ints = get_pp_ints(Delta, Delta_cc, imp.hamiltonian.hil_);
    ppsc::gf_verts_type gf_verts = get_gf_verts(imp.hamiltonian.hil_);
   
    imp.update_diagrams(pp_ints, gf_verts);
       
    // ---------------------------------------------------------------------
    // -- Initialization (tstp = -1)
    
    ppsc::gf_tstps_type gf_tstps;
    
    imp.solve_atomic();
    
    imp.update_density_matrix(-1);
    imp.hamiltonian.update_exp_vals(-1, imp.rho);
    
    gf_tstps = imp.get_spgf(-1);	// extract local GF without hybridization (because routine 'pp_step' has not been called yet)
    Gloc.set_timestep(-1, gf_tstps[0]); 
    G_pseudo = imp.ppGfs[1];
    
    // ---------------------------------------------------------------------
    // -- extract Sigma_ret_bare & G_ret_bare from ppsc, and G_ret_bare from nessi routine
    
    for(int tstp=-1; tstp<=nt; tstp++ ){ imp.solve_dyson(tstp); }
    //imp.ppGfs[1].print_to_file("gf0_ppsc.txt");
    //for(int tstp=-1; tstp<=nt; tstp++){ imp.update_sigma(tstp); }
    //imp.ppSigmas[1].print_to_file("sf0_ppsc.txt");
    

    cntr::herm_matrix<double> G0(nt, ntau, 1, -1);   
    green_from_H(G0,mu,eps0,beta,h); // mu = eps0 = 0
    
    
    /* output
    cout.precision(10);
    cout << "# pp_mu = " << imp.pp_mu << endl;
    cout << "# tstp \t diff GF0.real \t diff GF0.imag" << endl;
    for(int tstp=0;tstp<=nt;tstp++){
    	cdouble ppg, g;
   	imp.ppGfs[1].get_ret(nt,tstp,ppg);
   	G0.get_ret(nt,tstp,g);
   	//ppg = ppg/pp_oscillation(beta,imp.hamiltonian.eps,imp.hamiltonian.mu[0],nt,tstp,h); // if mu=0 in 'green_from_H'
   	ppg = ppg/pp_oscillation(imp.pp_mu, nt, tstp, h);
   	cout << (int) (nt-tstp) << " " << (ppg-g).real() << " " << (ppg-g).imag() << endl;
    }
    */
    
    // ---------------------------------------------------------------------
    // -- MATSUBARA PART (EQUILIBRIUM INITIAL STATE)

    
    {
      
      cntr::herm_matrix<double> gtmp(-1, ntau, 1, -1);

      for (iter = 1; iter <= itermax; iter++) {

	// -- do pseudo-particle step (also updates expectation values)
	imp.pp_step(-1);

	// -- get spgf
	gf_tstps = imp.get_spgf(-1);
	Gloc.set_timestep(-1, gf_tstps[0]);
	
	// -- check error (difference to previous solution 'gtmp')
	err = cntr::distance_norm2(-1, gtmp, Gloc);
	
	// -- update 'gtmp'
        gtmp.set_timestep(-1, Gloc);

	// -- command line output
	//std::cout << "iter:  " << iter << " n_exp: " << imp.hamiltonian.n_exp[0] << " err: " << err << " pp_mu: " << imp.pp_mu << std::endl;

        if (err < errmax)
          break;

      }
      
      // -- error handling
      if (iter > itermax) {
        std::cerr << "WARNING: Matsubara not converged after " << itermax
             << "steps ... abort" << std::endl;
        std::cerr << "skip real-time calculation " << std::endl;
        return 1;
      }
      
    }

    
    // ---------------------------------------------------------------------
    // 	START: ONLY FOR NT>0
    
    integration::Integrator<double> integrator(kt);
    
    if (nt > 0) {
	
      cntr::herm_matrix<double> gtmp(kt, ntau, 1, -1);
      cntr::herm_matrix<double> Gtmp(kt, ntau, 1, -1);
      
      imp.init_real_time();
      
      for (iter = 0; iter <= itermax; iter++) {
	imp.pp_step(kt);

	for (int tstp = 0; tstp <= kt; tstp++) {
	   update_sigma(Sigma,G0,tstp,order,Omega0,beta,gamma,h); //update self energy
	   gf_tstps = imp.get_spgf(tstp);
	   Gloc.set_timestep(tstp, gf_tstps[0]);	
	}
        //cntr::dyson_start(G0,imp.pp_mu,eps_fct,Sigma,integrator,beta,h); // alternative to below
	cntr::dyson_start(G0, mu, eps_fct, Sigma, beta, h, kt);
        err2=distance_norm2(kt,G0,Gtmp); //calculate difference to previous iteration
        Gtmp.set_timestep(kt,G0); //update Gtemp_pseudo
	
        err = cntr::distance_norm2(kt, gtmp, Gloc);
        gtmp.set_timestep(kt, Gloc);	// copy only timeslice with tstp = kt to 'gtmp' as it is the latest time computed and the only one compared in the routine above
        //std::cout << "START: iter:  " << iter << " n_exp(kt): " << imp.hamiltonian.n_exp[kt+1] << " err: " << err << " pp_mu: " << imp.pp_mu << std::endl;
       if (err < errmax && err2 < errmax)	break;
    
      }
      std::cout << "iter = " << iter+1 << std::endl;
      
      // -- error handling
      if (iter > itermax) {
        std::cerr << "WARNING: Bootstrapping not converged after " << itermax
             << "steps ... abort" << std::endl;
        std::cerr << "skip real-time calculation for t>kt" << std::endl;
        return 1;
      }


    // ---------------------------------------------------------------------
    // 	REALTIME: ONLY FOR NT>0

      for (int tstp = kt + 1; tstp <= nt; tstp++) {
    
        cntr::herm_matrix<double> gtmp(tstp, ntau, 1, -1);
        cntr::herm_matrix<double> Gtmp(tstp, ntau, 1, -1);
        imp.extrapolate_timestep(tstp - 1);
      
        for (iter = 0; iter <= itermax; iter++) {
        	// ppsc
	  imp.pp_step(tstp);
	  gf_tstps = imp.get_spgf(tstp);
	  Gloc.set_timestep(tstp, gf_tstps[0]);	
          err = cntr::distance_norm2(tstp, gtmp, Gloc);
          gtmp.set_timestep(tstp, Gloc);	
	  //std::cout << "tstp:  " << tstp << " n_exp: " << imp.hamiltonian.n_exp[tstp+1] << " err: " << err << " pp_mu: " << imp.pp_mu << std::endl;
	  	// self-implemented Sigma (nessi)
	  update_sigma(Sigma,G0,tstp,order,Omega0,beta,gamma,h); //update self energy
	  //cntr::dyson_timestep(tstp,G0,imp.pp_mu,eps_fct,Sigma,integrator,beta,h); //solve Dyson equation; alternative to below
	  cntr::dyson_timestep(tstp, G0, mu, eps_fct, Sigma, beta, h, kt); //solve Dyson equation	
	
	  err2=distance_norm2(tstp,G0,Gtmp); //calculate difference to previous iteration
	  Gtmp.set_timestep(tstp,G0); //update Gtemp
	  
	  if (err < errmax && err2 < errmax)
            break;
          
        }
        std::cout << "iter = " << iter+1 << std::endl;
      
      // -- error handling
        if (iter > itermax) {
          std::cerr << "WARNING: Real-time evolution not converged after " << itermax
               << "steps ... abort" << std::endl;
          std::cerr << "skip remaining real-time calculation" << std::endl;
          return 1;
        }
      
      
      }

    }	// it nt > 0
   
   
    // ---------------------------------------------------------------------
    /* -- output I
    cntr::herm_matrix<double> Spseudo(nt, ntau, 1, FERMION);
    for(int tstp=0;tstp<=nt;tstp++) update_sigma(Spseudo,imp.ppGfs[1],tstp,order,Omega0,beta,gamma,h);
    cout.precision(10);
    for(int tstp=0;tstp<=nt;tstp++){
    	cdouble g1,g2,g3;
   	imp.ppGfs[1].get_ret(nt,tstp,g1);
    	Spseudo.get_ret(nt,tstp,g2);
    	imp.ppSigmas[1].get_ret(nt,tstp,g3);
    	cout << tstp;
    	cout << " " << g1.real() << " " << g1.imag();
        cout << " " << g2.real() << " " << g2.imag();
        cout << " " << g3.real() << " " << g3.imag();
        cout << endl;
    }
    */
    
    /* -- output II 
    cout.precision(10);
    for(int tstp=0;tstp<=nt;tstp++){
    	cdouble g1,g2;
   	G_pseudo.get_ret(nt,tstp,g1);	// self-implemented Sigma
    	imp.ppGfs[1].get_ret(nt,tstp,g2);	// from ppsc solver
    	cout << tstp;
    	cout << " " << g1.real() << " " << g1.imag();
        cout << " " << g2.real() << " " << g2.imag();
        //cout << " " << g3.real() << " " << g3.imag();
        cout << endl;
    }
    */
    
    
    // -- output III
    cout.precision(20);
    cout << "# pp_mu = " << imp.pp_mu << endl;
    cout << "# tstp \t GF.real \t GF.imag" << endl;
    for(int tstp=nt;tstp>=0;tstp--){
    	cdouble g1, g2;
   	imp.ppGfs[1].get_ret(nt,tstp,g1);
   	G0.get_ret(nt,tstp,g2);
   	g1 = g1/pp_oscillation(imp.pp_mu, nt, tstp, h);
   	cout << (int) (nt-tstp) << " " << g2.real() << " " << g2.imag() << endl;
    }
    
    //G0.print_to_file("g_oca.txt");
    
  return 0;
}



