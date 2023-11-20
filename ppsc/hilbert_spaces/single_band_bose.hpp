#ifndef _PPSC_HILBERT_SPACE_SINGLE_BAND_BOSE_HPP
#define _PPSC_HILBERT_SPACE_SINGLE_BAND_BOSE_HPP

// -----------------------------------------------------------------------
//
// Single band bosonic hilbert space
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

// May the power be with you...
// O(log_2(N)) templated power for (positive) integer types
// a' la Knuth
template<class T>
T TemplatePow(T x, T p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  // int tmp = TemplatePow<T>(x, p/2);
  int tmp = TemplatePow<T>(x, p >> 1);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}
  
// -----------------------------------------------------------------------

class single_band_bose : public hilbert_space_base {

public:
  
  typedef hilbert_space_base base_type;
  typedef double CoeffType;

  
  // ---------------------------------------------------------------------
  void init(int Nmax) {

    Nmax_ = Nmax;
    norb_ = 1;
    spin_degeneracy_ = 1;
    // spin_conservation_ = false;
    xi_ = 1;
    
    //base_type::init(norb, spin_degeneracy, spin_conservation);

    // hilbert space dimension Nmax^{nflavour}
    nflavor = norb_ * spin_degeneracy_;
    nh_ = std::pow(Nmax, nflavor);

    // Dense version all states are in one big merry sector
    ns_ = 1;

    // Initialize local class vectors
    hilbert_space_base::init_storage();
    init_states();
    init_c_cdag();
    init_vacuum(); 
  }
  
  // ---------------------------------------------------------------------

  void c(int orb,int spin,state &target,state psi, CoeffType & sign) const;
  void cdag(int orb,int spin,state &target,state psi, CoeffType & sign) const;
  // void c(int orb,int spin,state &target,state psi,int &sign) const;
  // void cdag(int orb,int spin,state &target,state psi, int & sign) const;

  // n(orb,spin)|psi> = sign |psi>
  void n(int orb,int spin,state psi, CoeffType &sign) const;
  // void n(int orb,int spin,state psi,int &sign) const;

  // -1 if state is nvalid (sign=0)
  int c_to_sector(int orb,int spin,int sector) const;
  int cdag_to_sector(int orb,int spin,int sector) const;

  // ---------------------------------------------------------------------
  // save c and cdag into a full Blockmatrix
  virtual void set_operator_matrix_c(int orb, int spin, cOp &op_c) {
    int s, dim1, j;
    CoeffType sign;
    state psi0, psi1;
    std::vector<int> tsector(ns_);
    for (s = 0; s < ns_; s++)
      tsector[s] = c_to_sector(orb, spin, s);
    op_c.InitZero(ssdim_, tsector); // sets correct dimensions
    for (s = 0; s < ns_; s++) {
      dim1 = ssdim_[s];
      if (op_c.to_sector_[s] != -1) {
        // std::cout << "C: s="<<s<<" -> ts="<<tsector[s]<<"
        // dim1="<<dim1<<std::endl;
        for (j = 0; j < dim1; j++) {
          psi0 = nstock_[s][j];
          c(orb, spin, psi1, psi0, sign);
          if (sign) {
            // std::cout << "psi0="<<psi0<<" idx="<< j << " -> psi1="<<psi1<<"
            // idx="<< ninv_[psi1] <<std::endl;
            op_c.M_[s](ninv_[psi1], j) += sign;
          }
        }
      }
    }
  }
  virtual void set_operator_matrix_cdag(int orb, int spin, cOp &op_cdag) {
    int s, dim1, j;
    CoeffType sign;
    state psi0, psi1;
    std::vector<int> tsector(ns_);
    for (s = 0; s < ns_; s++)
      tsector[s] = cdag_to_sector(orb, spin, s);
    op_cdag.InitZero(ssdim_, tsector); // sets correct dimensions
    for (s = 0; s < ns_; s++) {
      dim1 = ssdim_[s];
      if (op_cdag.to_sector_[s] != -1) {
        for (j = 0; j < dim1; j++) {
          psi0 = nstock_[s][j];
          cdag(orb, spin, psi1, psi0, sign);
          if (sign)
            op_cdag.M_[s](ninv_[psi1], j) += sign;
        }
      }
    }
  }

  // ---------------------------------------------------------------------
  // set up the c and cdag operators
  void init_c_cdag() {

    nf_ = norb_ * spin_degeneracy_ * 2;
    c_op_.resize(nf_);
    for (int orb = 0; orb < norb_; orb++) {
      for (int spin = 0; spin < spin_degeneracy_; spin++) {
        set_operator_matrix_c(   orb, spin, c_op_[flavor(orb, spin, 0)]);
        set_operator_matrix_cdag(orb, spin, c_op_[flavor(orb, spin, 1)]);
      }
    }
    
  }

  // ---------------------------------------------------------------------
  // set up vacuum
  void init_vacuum() {
    cdvector v1(1);
    int s0 = sector_[0];
    vacuum_.InitZero(ssdim_);
    v1(0) = 1.0;
    // just a check:
    // CNTR_ASSERT_EQ(NCA_SOLVER_ASSERT_0, 1, (int)vacuum_.M_[s0].rows(),
    //                __PRETTY_FUNCTION__)
    vacuum_.M_[sector_[0]] = v1;
  }

  // ---------------------------------------------------------------------
  void init_states() {
    int np, i, sector;

    // Loop over all states
    for(state psi = 0; psi < nh_; psi++){

      // count number of each spin flavor in psi, and find the sector
      // sector = (nspin(0)*norb + spin(1))*norb + nspin(2) ....  
      // (spin-conservation)
      // sector = nparticles (no spin-conservation)

      // Compute occupation numbers for each flavour
      size_t Ntot = 0;
      state psi1 = psi;
      std::vector<size_t> state_vec(nflavor);

      for(size_t flavor = 0; flavor < (size_t) nflavor; flavor++) {

	// <n> = numeric_cast<int>(psi / Nmax^flavor) % Nmax
	size_t n = (psi1 / TemplatePow(Nmax_, flavor)) % Nmax_;
	state_vec[flavor] = n;
	Ntot += n;

      }

      np = Ntot; // Number of particles
      sector = 0; // dense version has only one big sector

      i = ssdim_[sector];

      npart_[sector] = np;
      nstock_[sector][i] = psi;
      ninv_[psi] = i;
      sector_[psi] = sector;
      ssdim_[sector]++;
      sig_[sector] = 1;
    }

    // each element in nstock_ was original a vector of length nh_
    // resize to final (actually needed) size
    for(sector = 0; sector < ns_; sector++) {
      nstock_[sector].resize(ssdim_[sector]);
    }
  }

  // ---------------------------------------------------------------------
  
private:
  int nflavor;
  
};

// -----------------------------------------------------------------------
void single_band_bose::c(
  int orb, int spin, state &target, state psi, CoeffType &sign) const {

  size_t flavor = spin * norb_ + orb;
  //state psi1 = psi;
  
  assert((spin < spin_degeneracy_) && (orb < norb_));
  
  // Compute occupation in flavour
  size_t Nmax_pow_flavor = TemplatePow(Nmax_, flavor);
  size_t n = (psi / Nmax_pow_flavor) % Nmax_;
  
  target = psi - (n > 0) * Nmax_pow_flavor;
  sign = sqrt(n); // Bosonic coefficient

  /*
  if((1 << orb1) & psi) {
    
    target = (1 << orb1) xor psi; 

    int left_of = 0;
    // sign = 0;

    for(int l = 0 ; l < orb1; l++) {
      left_of += psi1 % 2;
      //sign += psi1 % 2;
      psi1 /= 2;
    } 

    sign = (left_of % 2 == 1 ? -1 : 1);

  } else {
    target = psi;
    sign = 0;
  }
  */
}

// -----------------------------------------------------------------------
void single_band_bose::cdag(
  int orb, int spin, state &target, state psi, CoeffType &sign) const {

  size_t flavor = spin * norb_ + orb;
  //state psi1 = psi;
  
  assert((spin < spin_degeneracy_) && (orb < norb_));
  
  // Compute occupation in flavour
  size_t Nmax_pow_flavor = TemplatePow(Nmax_, flavor);
  size_t n = (psi / Nmax_pow_flavor) % Nmax_;
  
  target = psi + (n < Nmax_-1) * Nmax_pow_flavor;
  //target = psi;

  sign = (n < Nmax_-1) * sqrt(n+1); // Bosonic coefficient

  //sign = 1.337;

  /*
  int orb1=spin*norb_+orb,l;
  state psi1=psi;

  assert(spin<spin_degeneracy_ && orb<norb_);

  if(  !((1<<orb1) & psi) ){
    target= (1 << orb1) | psi;     
    int left_of = 0;
    //sign=0;
    for(l=0;l<orb1;l++){
      left_of += psi1 % 2;
      //sign += psi1 % 2;
      psi1 /= 2;
    } 
    sign = (left_of % 2 == 1 ? -1 : 1);
    //sign = (sign % 2 == 1 ? -1 : 1);
  }else{
    target=psi;
    sign=0;
  }
  */
}

// -----------------------------------------------------------------------
void single_band_bose::n(int orb,int spin,state psi, CoeffType &sign) 
const{

  size_t flavor = spin * norb_ + orb;
  
  assert((spin < spin_degeneracy_) && (orb < norb_));
  
  // Compute occupation in flavour
  size_t Nmax_pow_flavor = TemplatePow(Nmax_, flavor);
  size_t n = (psi / Nmax_pow_flavor) % Nmax_;

  sign = n;

  /*
	  int orb1=spin*norb_+orb;
	  assert(spin<spin_degeneracy_ && orb<norb_);
	  if(  (1<<orb1) & psi){
	  	sign = 1;
	  }else{
	    sign=0;
	  }
  */
}

// -----------------------------------------------------------------------
int single_band_bose::c_to_sector(int orb,int spin,int sector)
const{
	  int tsector,dim,i;
	  state tstate,psi0;
	  single_band_bose::CoeffType sign;
	  tsector=-1;
	  sign=0;
	  dim=ssdim_[sector];
	  for(i=0;i<dim;i++){
	    psi0=nstock_[sector][i];
	    c(orb,spin,tstate,psi0,sign);
		if(sign){
		  tsector=sector_[tstate];
		  break;
		}
	  }
	  return tsector;
}

// -----------------------------------------------------------------------
int single_band_bose::cdag_to_sector(int orb,int spin,int sector)
const {
	  int tsector,dim,i;
	  state tstate,psi0;
	  single_band_bose::CoeffType sign;
	  tsector=-1;
	  sign=0;
	  dim=ssdim_[sector];
	  for(i=0;i<dim;i++){
	    psi0=nstock_[sector][i];
	    cdag(orb,spin,tstate,psi0,sign);
		if(sign){
		  tsector=sector_[tstate];
		  break;
		}
	  }
	  return tsector;
}

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_SINGLE_BAND_BOSE_HPP
// -----------------------------------------------------------------------
