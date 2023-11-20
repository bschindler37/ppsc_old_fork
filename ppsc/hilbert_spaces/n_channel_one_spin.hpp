
#ifndef _PPSC_HILBERT_SPACE_N_CHANNEL_ONE_SPIN_HPP
#define _PPSC_HILBERT_SPACE_N_CHANNEL_ONE_SPIN_HPP

// -----------------------------------------------------------------------
//
// Hilbert space with a doubly degenerate vacuum, in order to simulate a
// two-channel Kondo model
//
// Author: M. Eckstein, martin.eckstein@mpsd.cfel.de (2016)
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

class n_channel_one_spin {

public:

  // ------------------------------------------------------------------
  void init(int nchannel) {
    int m,sector;
    nchannel_ = nchannel;
    spin_degeneracy_ = 2;
    nh_ = nchannel+2;
    ns_ = nchannel+2;
    
    // set up sectors by hand:
    //    nchannel+2 states
    //    |m>    vacuum with flavor m=0...nchannel-1, Ndo=Nup=0  Sector 0
    //    |up>    vacuum with flavor 1, Ndo=0 Nup=1  Sector 2
    //    |do>    vacuum with flavor 1, Ndo=1 Nup=0  Sector 3

    ssdim_ = std::vector<int>(ns_, 1);
    npart_ = std::vector<int>(ns_);
    sig_ = std::vector<int>(ns_);
    for(m=0;m<nchannel_;m++){
        sector = m;
        npart_[sector] = 0;
        sig_[sector] = 1;
    }
    for(m=0;m<spin_degeneracy_;m++){
        sector = nchannel_+m;
        npart_[sector] = 1;
        sig_[sector] = -1;
    }
    // operator  c_{m,sigma} = f_{sigma} b_{m}^dagger is given by the matrix c_op[flavor(m,sigma,0)]
    // operator  c_{m,sigma}^dag = f_{sigma}^dag b_{m} is given by the matrix c_op[flavor(m,sigma,1)]
    init_c_cdag();
  }
  int cdag_to_sector(int m, int spin, int sector) {
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, m, nchannel_ - 1,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,__PRETTY_FUNCTION__)
    int tsector=-1;
    if(sector==m && spin==0) tsector=nchannel_;
    if(sector==m && spin==1) tsector=nchannel_+1;
    return tsector;
  }
  int c_to_sector(int m, int spin, int sector) {
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, m, nchannel_ - 1,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,__PRETTY_FUNCTION__)
    int tsector=-1;
    if(sector==nchannel_ && spin==0) tsector=m;
    if(sector==nchannel_+1 && spin==1) tsector=m;
    return tsector;
  }

  void set_operator_matrix_c(int m, int spin, cOp &op_c) {
    std::vector<int> tsector(ns_);
    for (int s = 0; s < ns_; s++) tsector[s] = c_to_sector(m, spin, s);
    op_c.InitZero(ssdim_, tsector); // sets correct dimensions
    // set sectors by hand
    if(spin==0) op_c.M_[nchannel_](0,0) = 1.0;
    if(spin==1) op_c.M_[nchannel_+1](0,0) = 1.0;
  }
  void set_operator_matrix_cdag(int m, int spin, cOp &op_cdag) {
    std::vector<int> tsector(ns_);
    for (int s = 0; s < ns_; s++) tsector[s] = cdag_to_sector(m, spin, s);
    op_cdag.InitZero(ssdim_, tsector); // sets correct dimensions
    // set sectors by hand
    if(spin==0) op_cdag.M_[m](0,0) = 1.0;
    if(spin==1) op_cdag.M_[m](0,0) = 1.0;
  }
  
  int flavor(int m, int spin, int dag) {
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, dag, 1, __PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, m, m_ - 1,
                        __PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,
                        __PRETTY_FUNCTION__)
    return dag * nchannel_ * spin_degeneracy_ + nchannel_ * spin + m;
  }
    void init_c_cdag() {
        nf_ = nchannel_ * spin_degeneracy_ * 2;
        c_op_.resize(nf_);
        for (int m = 0; m < nchannel_; m++) {
          for (int spin = 0; spin < spin_degeneracy_; spin++) {
            // std::cout << "set c[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
            set_operator_matrix_c(m, spin, c_op_[flavor(m, spin, 0)]);
            // std::cout << "set cdag[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
            set_operator_matrix_cdag(m, spin, c_op_[flavor(m, spin, 1)]);
          }
        }
      }
    
      int nchannel_,spin_degeneracy_,nh_,ns_;
      std::vector<int> ssdim_,sig_,npart_;
      int nf_;
      std::vector<cOp> c_op_;
};

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_N_CHANNEL_ONE_SPIN_HPP
// -----------------------------------------------------------------------
