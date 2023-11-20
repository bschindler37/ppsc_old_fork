
#ifndef _PPSC_HILBERT_SPACE_two_band_fermi_tJ_pos_POS_HPP
#define _PPSC_HILBERT_SPACE_two_band_fermi_tJ_pos_POS_HPP

// -----------------------------------------------------------------------
//
// Two band fermionic hilbert space restricted to states at half filling
// bosonic tJ with enabled spin and Hund interaction
//
// Author: Martin Eckstein (201?)
//         Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//         Denis Golez, denis.golez@gmail.com (2018)
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

  /*

    Almost identical to two_band_fermi_densdens,
    except that we integrate out all doublons, triplons

    The following type of terms in the Hamiltonian are allowed:
    *  density-density interaction terms
    *  hybridization diagonal in spin and orbital


    do  up
    B A B A    Ndo,Nup    Sector
    ---------------------------------------------
    0 0 0 0     0,0         0 Int out
    ---------------------------------------------
    0 0 0 1                 1 Int out
                0,1  ----------------------------
    0 0 1 0                 2 Int out
    ---------------------------------------------
    0 1 0 0                 3 Int out
                1,0  ----------------------------
    1 0 0 0                 4 Int out
    ---------------------------------------------
    0 0 1 1     0,2         5
    ---------------------------------------------
    1 1 0 0     2,0         6
    ---------------------------------------------
    0 1 0 1                 7 Int out 
                     ----------------------------
    1 0 1 0                 8 Int out
                1,1   ----------------------------
    0 1 1 0                 9
                     ----------------------------
    1 0 0 1                10
    ---------------------------------------------
    1 1 1 0                11 Int out 
               2,1   ----------------------------
    1 1 0 1                12 Int out
    ---------------------------------------------
    1 0 1 1                13 Int out
               1,2   ----------------------------
    0 1 1 1                14 Int out
    ---------------------------------------------
    1 1 1 1                15 Int out
    ---------------------------------------------

  */

// -----------------------------------------------------------------------
class two_band_fermi_tJ_pos : public hilbert_space_base {
public:
  std::vector<cOp> n_op_; // (***)
  std::vector<cOp> Sp_op_; // (***)
  std::vector<cOp> Sm_op_; // (***)
  // ------------------------------------------------------------------
  void init(void) {
    norb_ = 2;
    spin_degeneracy_ = 2;
    nh_ = 10;
    ns_ = 10;

    // ----------------------------------------------------------------
    init_storage();    
    ssdim_ = std::vector<int>(ns_, 1);
    
    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);
    // ----------------------------------------------------------------
    // determine all states by hand

    int sector;

    sector = 0;
    nstock_[sector][0] = 3;
    ninv_[3] = 0;
    sector_[3] = sector; // 0 0 1 1

    sector = 1;
    nstock_[sector][0] = 12;
    ninv_[12] = 0;
    sector_[12] = sector; // 1 1 0 0

    sector = 2;
    nstock_[sector][0] = 5;
    ninv_[5] = 0;
    sector_[5] = sector; // 0 1 0 1

    sector = 3;
    nstock_[sector][0] = 10;
    ninv_[10] = 0;
    sector_[10] = sector; // 1 0 1 0 

    sector = 4;
    nstock_[sector][0] = 6;
    ninv_[6] = 0;
    sector_[6] = sector; // 0 1 1 0

    sector = 5;
    nstock_[sector][0] = 9;
    ninv_[9] = 0;
    sector_[9] = sector; // 1 0 0 1

    sector = 6;
    nstock_[sector][0] = 14;
    ninv_[14] = 0;
    sector_[14] = sector; // 1 1 1 0 

    sector = 7;
    nstock_[sector][0] = 13;
    ninv_[13] = 0;
    sector_[13] = sector; // 1 1 0 1

    sector = 8;
    nstock_[sector][0] = 11;
    ninv_[11] = 0;
    sector_[11] = sector; // 1 0 1 1   

    sector = 9;
    nstock_[sector][0] = 7;
    ninv_[7] = 0;
    sector_[7] = sector; // 0 1 1 1

    sector_[0] = -1;  // 0000
    sector_[1] = -1;  // 0001
    sector_[2] = -1;  // 0010
    sector_[4] = -1;  // 0010
    sector_[8] = -1;  // 0010
    sector_[15] = -1; // 1111 
     

    // ----------------------------------------------------------------
    init_npart_sig();
    init_c_cdag();
    init_n();
    init_S();
    init_vacuum();
  }

  void init_storage() {
    ssdim_ = std::vector<int>(ns_, 0);
    npart_ = std::vector<int>(ns_);
    sig_ = std::vector<int>(ns_);
    ninv_ = std::vector<int>(16);
    sector_ = std::vector<int>(16);
    nstock_ = std::vector<std::vector<state>>(ns_, std::vector<state>(16));
  }

  // ------------------------------------------------------------------
  void init_npart_sig() {

    // There is an assumption that the number of particles
    // in a sector is constant! npart_.size() == ns_ !!
    // This is not the case for symmetry breaking,
    // i.e. |0> and |up,down> are in the same sector.

    for (int sector = 0; sector < ns_; sector++) {
      // count particles:
      int psi = nstock_[sector][0]; // Nb! using only one state from the sector
      int np;
      for (np = 0; psi; psi >>= 1)
        np += psi & 1; // Count the number of bits set (destructive on psi)

      npart_[sector] = np;
      sig_[sector] = 1 - 2 * (np % 2);
    }
  }

  // set n operators
  void init_n() {
    nf_ = norb_ * spin_degeneracy_ ;
    n_op_.resize(nf_);
    for (int orb = 0; orb < norb_; orb++) {
      for (int spin = 0; spin < spin_degeneracy_; spin++) {
        set_operator_matrix_n(orb, spin, n_op_[flavor(orb, spin, 0)]);
      }
    }
  }

  void init_S() {
    nf_ = norb_  ;
    Sp_op_.resize(nf_);
    for (int orb = 0; orb < norb_; orb++) {
        set_operator_matrix_Sp(orb,  Sp_op_[flavor(orb, 0, 0)]);
    }
    Sm_op_.resize(nf_);
    for (int orb = 0; orb < norb_; orb++) {
        set_operator_matrix_Sm(orb, Sm_op_[flavor(orb, 0, 0)]);
    }
  }

  void init_vacuum() {
    // sector 0 is the vacuum
    // cdvector v1(1);
    // vacuum_.InitZero(ssdim_);
    // v1(0) = 1.0;
    // vacuum_.M_[0] = v1;
  }

  void set_operator_matrix_n(int orb, int spin, cOp &op_n) {
    int s, dim1, j, sign;
    state psi0, psi1;
    std::vector<int> tsector(ns_);
    for (s = 0; s < ns_; s++)
      tsector[s] = n_to_sector(orb, spin,s);
    op_n.InitZero(ssdim_, tsector); // sets correct dimensions
    // std::cout << "ns " << ns_ << std::endl;
    for (s = 0; s < ns_; s++) {
      dim1 = ssdim_[s];
      if (op_n.to_sector_[s] != -1) {
        // std::cout << "n: s="<<s<<" -> ts="<<tsector[s]<<"dim1="<<dim1<<std::endl;
        for (j = 0; j < dim1; j++) {
          psi0 = nstock_[s][j];
          n(orb, spin, psi0, sign);
          if (sign) {
            // std::cout << "psi0="<<psi0<<" idx="<< j << " -> psi1="<<psi1<<"idx="<< ninv_[psi1] <<std::endl;
            op_n.M_[s](ninv_[psi0], j) += sign;
          }
        }
      }
    }
  }

  int n_to_sector(int orb, int spin,int sector) {
    int tsector, dim, i, sign;
    state tstate, psi0;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,
    //                     __PRETTY_FUNCTION__)
    tsector = -1;
    sign = 0;
    dim = ssdim_[sector];
    for (i = 0; i < dim; i++) {
      psi0 = nstock_[sector][i];
      n(orb, spin, psi0, sign);
      // std::cout << "n to sector " << sector << " " << orb << " " << spin << ": " << sign << std::endl;
      if (sign) {
    if( psi0 < sector_.size() ) {
          tsector = sector_[psi0];
    } else {
      tsector = -1;
    }
        break;
      }
    }
    return tsector;
  }

  void set_operator_matrix_Sm(int orb, cOp &op_c) {
    int s, dim1, j, sign;
    state psi0, psi1;
    std::vector<int> tsector(ns_);
    for (s = 0; s < ns_; s++)
      tsector[s] = Sm_to_sector(orb, s);
    op_c.InitZero(ssdim_, tsector); // sets correct dimensions
    for (s = 0; s < ns_; s++) {
      dim1 = ssdim_[s];
      if (op_c.to_sector_[s] != -1) {
        // std::cout << "C: s="<<s<<" -> ts="<<tsector[s]<<"
        // dim1="<<dim1<<std::endl;
        for (j = 0; j < dim1; j++) {
          psi0 = nstock_[s][j];
          Sm(orb, psi1, psi0, sign);
          if (sign) {
            // std::cout << "psi0="<<psi0<<" idx="<< j << " -> psi1="<<psi1<<"
            // idx="<< ninv_[psi1] <<std::endl;
            op_c.M_[s](ninv_[psi1], j) += sign;
          }
        }
      }
    }
  }

  int Sm_to_sector(int orb, int sector) {
    int tsector, dim, i, sign;
    state tstate, psi0;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    tsector = -1;
    sign = 0;
    dim = ssdim_[sector];
    for (i = 0; i < dim; i++) {
      psi0 = nstock_[sector][i];
      Sm(orb,  tstate, psi0, sign);
      if (sign) {
    if( tstate < sector_.size() ) {
          tsector = sector_[tstate];
    } else {
      tsector = -1;
    }
        break;
      }
    }
    return tsector;
  }
  
  void Sm(int orb, state &target, state psi, int &sign) {
    int orb1 = 0 * norb_ + orb, l;
    state psi1 = psi;
    state target1;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    if ((1 << orb1) & psi) {
      target1 = (1 << orb1) xor psi;
      sign = 0;
      for (l = 0; l < orb1; l++) {
        sign += psi1 % 2;
        psi1 /= 2;
      }
      sign = (sign % 2 == 1 ? -1 : 1);
    } else {
      target1 = psi;
      sign = 0;
    }
    int orb2 = 1 * norb_ + orb;
    state psi2 = target1;
    if (!((1 << orb2) & target1) && sign) { // The first and the second state should exist
      target = (1 << orb2) | target1;
      // sign = 0;
      for (l = 0; l < orb2; l++) {
        sign += psi2 % 2;
        psi2 /= 2;
      }
      sign = (sign % 2 == 1 ? -1 : 1);
    } else {
      target = target1;
      sign = 0;
    }
  }

  void set_operator_matrix_Sp(int orb, cOp &op_c) {
    int s, dim1, j, sign;
    state psi0, psi1;
    std::vector<int> tsector(ns_);
    for (s = 0; s < ns_; s++)
      tsector[s] = Sp_to_sector(orb, s);
    op_c.InitZero(ssdim_, tsector); // sets correct dimensions
    for (s = 0; s < ns_; s++) {
      dim1 = ssdim_[s];
      if (op_c.to_sector_[s] != -1) {
        // std::cout << "C: s="<<s<<" -> ts="<<tsector[s]<<"
        // dim1="<<dim1<<std::endl;
        for (j = 0; j < dim1; j++) {
          psi0 = nstock_[s][j];
          Sp(orb, psi1, psi0, sign);
          if (sign) {
            // std::cout << "psi0="<<psi0<<" idx="<< j << " -> psi1="<<psi1<<"
            // idx="<< ninv_[psi1] <<std::endl;
            op_c.M_[s](ninv_[psi1], j) += sign;
          }
        }
      }
    }
  }

  int Sp_to_sector(int orb, int sector) {
    int tsector, dim, i, sign;
    state tstate, psi0;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    tsector = -1;
    sign = 0;
    dim = ssdim_[sector];
    for (i = 0; i < dim; i++) {
      psi0 = nstock_[sector][i];
      Sp(orb, tstate, psi0, sign);
      if (sign) {
    if( tstate < sector_.size() ) {
          tsector = sector_[tstate];
    } else {
      tsector = -1;
    }
        break;
      }
    }
    return tsector;
  }
  
  void Sp(int orb, state &target, state psi, int &sign) {
    int orb1 = 1 * norb_ + orb, l;
    state psi1 = psi;
    state target1;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    if ((1 << orb1) & psi) {
      target1 = (1 << orb1) xor psi;
      sign = 0;
      for (l = 0; l < orb1; l++) {
        sign += psi1 % 2;
        psi1 /= 2;
      }
      sign = (sign % 2 == 1 ? -1 : 1);
    } else {
      target1 = psi;
      sign = 0;
    }

    int orb2 = 0 * norb_ + orb;
    state psi2 = target1;
    if (!((1 << orb2) & target1) && sign) { // The first and the second state should exist
      target = (1 << orb2) | target1;
      // sign = 0;  
      for (l = 0; l < orb2; l++) {
        sign += psi2 % 2;
        psi2 /= 2;
      }
      sign = (sign % 2 == 1 ? -1 : 1);
    } else {
      target = target1;
      sign = 0;
    }
  }

};
std::ostream &operator<<(
  std::ostream &os, two_band_fermi_tJ_pos & hilbert) {

  int nf = hilbert.spin_degeneracy_ * hilbert.norb_;
  hilbert_space_base::state psi;

  os << "--> ppsc::hilbert_spaces::hilbert_space_base:" << std::endl;
  os << "h.Nmax_ = " << hilbert.Nmax_ << std::endl;
  os << "h.nh_ = " << hilbert.nh_ << std::endl;
  os << "h.ns_ = " << hilbert.ns_ << std::endl;
  os << "h.norb_ = " << hilbert.norb_ << std::endl;
  os << "h.spin_degeneracy_ = " << hilbert.spin_degeneracy_ << std::endl;
  os << std::endl;

  // Loop over sectors
  for(int sector = 0; sector < hilbert.ns_; sector++) {

    os << "sector " << sector << ", nparticles: " << hilbert.npart_[sector]
       << ", dim: " << hilbert.ssdim_[sector]
       << ", sig: " << hilbert.sig_[sector]
       << std::endl;


    // Loop over state indices "i" in sector
    for(int i = 0; i < hilbert.ssdim_[sector]; i++) {

      psi = hilbert.nstock_[sector][i];
      os << "index " << i << "\t";

      psi = hilbert.nstock_[sector][i];
      os << "state " << psi << ", [0x" << std::bitset<8>(psi) << "] ";
      os << std::endl;
    }
  }

  for(int oidx = 0; oidx < hilbert.c_op_.size(); oidx++) {
    os << "operator idx " << oidx << std::endl;
    os << hilbert.c_op_[oidx];

    os << "dense repr" << std::endl;
    os << hilbert.c_op_[oidx].to_dense() << std::endl << std::endl;

  }

  cOp Ntot;
  Ntot.InitDiagZero(hilbert.ssdim_);
  
  std::vector<operator_type> n_op;
  n_op.resize(nf);

  for (int orb = 0; orb < hilbert.norb_; orb++) {
    for (int spin = 0; spin < hilbert.spin_degeneracy_; spin++) {
      hilbert.set_operator_matrix_n(orb,spin,n_op[hilbert.flavor(orb, spin, 0)]);
      os << "density operator for: orb, spin "
      << orb << ", " << spin << std::endl;
      std::cout << n_op[hilbert.flavor(orb, spin, 0)].to_dense() << std::endl;
    }
  }

  std::vector<operator_type> Sm;
  Sm.resize(nf);
  for (int orb = 0; orb < hilbert.norb_; orb++) {
      // std::cout << "set c[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
      hilbert.set_operator_matrix_Sm(orb, Sm[hilbert.flavor(orb, 0, 0)]);
      std::cout << "Spin- operator for: orb, spin "
      << orb << ", "  << std::endl;
      std::cout << Sm[hilbert.flavor(orb, 0, 0)].to_dense() << std::endl;
  }

  std::vector<operator_type> Sp;
  Sp.resize(nf);
  for (int orb = 0; orb < hilbert.norb_; orb++) {
      // std::cout << "set c[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
      hilbert.set_operator_matrix_Sp(orb, Sm[hilbert.flavor(orb, 0, 0)]);
      std::cout << "Spin+ operator for: orb, spin "
      << orb << ", " << std::endl;
      std::cout << Sm[hilbert.flavor(orb, 0, 0)].to_dense() << std::endl;
  }


  return os;
}

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_two_band_tJ_pos_HPP
// -----------------------------------------------------------------------
