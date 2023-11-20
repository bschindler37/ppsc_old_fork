
#ifndef _PPSC_HILBERT_SPACE_BASE_HPP
#define _PPSC_HILBERT_SPACE_BASE_HPP

// -----------------------------------------------------------------------
//
// Hilbert space constructions
//
// Author: Martin Eckstein (201?)
//
// Contributors:
//
// Hugo U. R. Strand, hugo.strand@gmail.com (2016)
//   Code reuse, L1 Nambu+Spin, L1 Bosons, ...
//
// -----------------------------------------------------------------------

#include <bitset>

// -----------------------------------------------------------------------

#include "ppsc/operator.hpp"

// -----------------------------------------------------------------------
namespace ppsc {
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
namespace hilbert_spaces {
// -----------------------------------------------------------------------

using namespace ppsc::operators; 	// use '::' for nested namespaces

// -----------------------------------------------------------------------
class hilbert_space_base {

public:

  typedef unsigned long state;

  int spin_degeneracy_;    // spin degeneracy=2S+1
  int norb_;               // number of orbitals (spin+orbital quantum numbers) [excluding spin !???]
  int ns_;                 // number of symmetry sectors  (***)
  unsigned long nh_;       // dimension of hilbert space
  unsigned long Nmax_;     // no. occupation number states (=1 fermions)
  int xi_;                 // statistics sign (-1 ferm, +1 bose)

  std::vector<int> ssdim_; // subspace dimensions in sectors   (***)
  std::vector<int> npart_; // number of partices in sectors
  std::vector<std::vector<state> > nstock_; // nstock_[s][i] = state i in sector s
  std::vector<int> ninv_;   // ninv_[psi] = index i of state psi in its sector
  std::vector<int> sector_; // sector_[psi] = sector number of state psi
  std::vector<int> sig_;    // even number of Ferimons: sign=+1, else -1 (***)
                            ///////
  // the annihilation and creation operators:
  int nf_; // nomber of flavors (including annihilation/creation character)
           // (***)
  std::vector<cOp> c_op_; 	// (***)  [remark Eva: class 'cOP' is defined in ppsc::operators]
  cVec vacuum_;		// cVec = Vector of standard complex double vectors (cdvector's), one cdvector for each symmetry sector

  // ---------------------------------------------------------------------
  void init(int norb, int spin_degeneracy = 2, bool spin_conservation = true) {
    int nspin, np, i, nflavor, sector, s, orb, spin;
    state psi, psi1;

    Nmax_ = 1;
    xi_ = -1;

    norb_ = norb;
    spin_degeneracy_ = spin_degeneracy;
    nflavor = norb * spin_degeneracy;
    nh_ = 1 << nflavor; // hilbert space dimension; each of the nflavor localized one-particle states can either exist or not (two options) in every many-particle state
    if (spin_conservation) {
      ns_ = 1;
      for (s = 0; s < spin_degeneracy; s++) {
        ns_ *= (norb + 1);	// there can be at most norb particles of the same spin in the model (plus the state where a particle of this spin does not exist at all) --> (norb + 1) many-particle states per spin
      }
    } else {
      ns_ = nflavor + 1;	// there are nflavor + 1 many-particle subspaces of different particle number (N=0, N=1, ... N=nflavor)
    }

    init_storage();
    init_states(spin_conservation = spin_conservation);
    init_c_cdag();
    init_vacuum();
  }

  void init_storage() {
    ssdim_ = std::vector<int>(ns_, 0);
    npart_ = std::vector<int>(ns_);
    sig_ = std::vector<int>(ns_);
    ninv_ = std::vector<int>(nh_);
    sector_ = std::vector<int>(nh_);
    nstock_ = std::vector<std::vector<state>>(ns_, std::vector<state>(nh_));
  }

  // set up all states
  void init_states(int spin_conservation = true) {
    // determine all states
    for (state psi = 0; psi < nh_; psi++) {
      // count number of each spin flavor in psi, and find the sector
      // sector = (nspin(0)*norb + spin(1))*norb + nspin(2) ....
      // (spin-conseration)
      // sector = nparticles (no spin-conservation)
      state psi1 = psi;
      int np = 0;
      int sector = 0;
      for (int s = 0; s < spin_degeneracy_; s++) {
        int nspin = 0;
        for (int orb = 0; orb < norb_; orb++) {
          nspin += psi1 % 2;	// read rightmost bit (0 or 1)
          psi1 /= 2;	// equivalent to psi1 = (psi1 >> 1)
        }
        np += nspin;
        sector = sector * (norb_ + 1) + nspin; 	// in base (norb_ + 1), sector has the digit-al representation [x(D) x(D-1) ... x(1) x(0)], where D=2S+1 is the spin degeneracy and x(i) is the number of particles with spin i of state psi (for each spin i, there can be 0, 1, ... norb_ particles of that kind)
      }
      if (!spin_conservation)
        sector = np;
      int i = ssdim_[sector];
      npart_[sector] = np;
      sig_[sector] = 1 - 2 * (np % 2);	// +1 (-1) if np is even (odd)
      nstock_[sector][i] = psi;
      ninv_[psi] = i;
      sector_[psi] = sector;
      ssdim_[sector]++;
    }
    for (int sector = 0; sector < ns_; sector++) {
      nstock_[sector].resize(ssdim_[sector]);
    }
  }

  // set up the c and cdag operators
  void init_c_cdag() {
    nf_ = norb_ * spin_degeneracy_ * 2;
    c_op_.resize(nf_);
    for (int orb = 0; orb < norb_; orb++) {
      for (int spin = 0; spin < spin_degeneracy_; spin++) {
        // std::cout << "set c[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
        set_operator_matrix_c(orb, spin, c_op_[flavor(orb, spin, 0)]);
        // std::cout << "set cdag[orb="<<orb<<",spin="<<spin<<"]" <<std::endl;
        set_operator_matrix_cdag(orb, spin, c_op_[flavor(orb, spin, 1)]);
      }
    }
  }

  // set up vacuum
  void init_vacuum() {
    cdvector v1(1);
    int s0 = sector_[0];
    vacuum_.InitZero(ssdim_);
    v1(0) = 1.0;
    // just a check:
    // CNTR_ASSERT_EQ(NCA_SOLVER_ASSERT_0, 1, (int)vacuum_.M_[s0].rows(),__PRETTY_FUNCTION__)
    vacuum_.M_[sector_[0]] = v1;
  }

  // map to determine the flavor index: (orb,spin)
  int flavor(int orb, int spin, int dag) {
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, dag, 1, __PRETTY_FUNCTION__)
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,__PRETTY_FUNCTION__)
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,__PRETTY_FUNCTION__)
    return dag * norb_ * spin_degeneracy_ + norb_ * spin + orb;
  }
  // the following operators are "simple" in occ-number basis:
  // Op(orb,spin)|psi> = sign |target>
  void c(int orb, int spin, state &target, state psi, int &sign) {
    int orb1 = spin * norb_ + orb, l;
    state psi1 = psi;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,__PRETTY_FUNCTION__)
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,__PRETTY_FUNCTION__)
    if ((1 << orb1) & psi) {
      target = (1 << orb1) xor psi;
      sign = 0;
      for (l = 0; l < orb1; l++) {
        sign += psi1 % 2;
        psi1 /= 2;
      }
      sign = (sign % 2 == 1 ? -1 : 1);
    } else {
      target = psi;
      sign = 0;
    }
  }
  void cdag(int orb, int spin, state &target, state psi, int &sign) {
    int orb1 = spin * norb_ + orb, l;
    state psi1 = psi;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,
    //                     __PRETTY_FUNCTION__)
    if (!((1 << orb1) & psi)) {
      target = (1 << orb1) | psi;
      sign = 0;
      for (l = 0; l < orb1; l++) {
        sign += psi1 % 2;
        psi1 /= 2;
      }
      sign = (sign % 2 == 1 ? -1 : 1);
    } else {
      target = psi;
      sign = 0;
    }
  }
  void n(int orb, int spin, state psi, int &sign) {
    int orb1 = spin * norb_ + orb;
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, orb, norb_ - 1,
    //                     __PRETTY_FUNCTION__)
    // CNTR_ASSERT_LESEQ_3(NCA_SOLVER_ASSERT_1, 0, spin, spin_degeneracy_ - 1,
    //                     __PRETTY_FUNCTION__)
    if ((1 << orb1) & psi) {
      sign = 1;
    } else {
      sign = 0;
    }
  }
  int c_to_sector(int orb, int spin, int sector) {
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
      c(orb, spin, tstate, psi0, sign);
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
  int cdag_to_sector(int orb, int spin, int sector) {
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
      cdag(orb, spin, tstate, psi0, sign);
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
  // save c and cdag into a full Blockmatrix
  virtual void set_operator_matrix_c(int orb, int spin, cOp &op_c) {
    int s, dim1, j, sign;
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
    int s, dim1, j, sign;
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
};

// -----------------------------------------------------------------------
std::ostream &operator<<(
  std::ostream &os, hilbert_space_base & hilbert) {

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

  for (int orb = 0; orb < hilbert.norb_; orb++) {
    for (int spin = 0; spin < hilbert.spin_degeneracy_; spin++) {
      auto ca = hilbert.c_op_[hilbert.flavor(orb, spin, 0)];
      auto cc = hilbert.c_op_[hilbert.flavor(orb, spin, 1)];
      auto n = cc * ca;
      os << "density operator for: orb, spin "
	 << orb << ", " << spin << std::endl;
      os << n.to_dense() << std::endl;
      Ntot = Ntot + n;
    }
  }

  os << std::endl;
  os << "total density operator: " << std::endl;
  os << Ntot << std::endl;
  os << Ntot.to_dense() << std::endl;


  return os;
}

// -----------------------------------------------------------------------
} // namespace hilbert_spaces
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
} // namespace ppsc
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
#endif // _PPSC_HILBERT_SPACE_BASE_HPP
// -----------------------------------------------------------------------
