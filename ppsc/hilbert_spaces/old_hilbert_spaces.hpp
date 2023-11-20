// ********************************************************************
class HilbertSpaceFermi_2orb_Nambu : public HilbertSpaceFermi {

public:
  void init() {

    norb_ = 2;
    spin_degeneracy_ = 2;
    nh_ = 16;
    ns_ = 9;

    // ----------------------------------------------------------------

    init_storage();

    ssdim_[0] = 4;
    ssdim_[1] = 2;
    ssdim_[2] = 2;
    ssdim_[3] = 2;
    ssdim_[4] = 1;
    ssdim_[5] = 1;
    ssdim_[6] = 2;
    ssdim_[7] = 1;
    ssdim_[8] = 1;

    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);

    // ----------------------------------------------------------------

    int sector;

    /*

Need to flip the original definition of orbital and spin order

    [ na_up, na_do, nb_up, nb_do ] --- to --- [ na_up, nb_up, na_do, nb_do ]

as assumed internally in the routines building the creation and annihilation
operators.

 0 = [0 0 0 0] -> [0 0 0 0] =  0
 1 = [0 0 0 1] -> [0 0 0 1] =  1
 2 = [0 0 1 0] -> [0 1 0 0] =  4
 3 = [0 0 1 1] -> [0 1 0 1] =  5
 4 = [0 1 0 0] -> [0 0 1 0] =  2
 5 = [0 1 0 1] -> [0 0 1 1] =  3
 6 = [0 1 1 0] -> [0 1 1 0] =  6
 7 = [0 1 1 1] -> [0 1 1 1] =  7
 8 = [1 0 0 0] -> [1 0 0 0] =  8
 9 = [1 0 0 1] -> [1 0 0 1] =  9
10 = [1 0 1 0] -> [1 1 0 0] = 12
11 = [1 0 1 1] -> [1 1 0 1] = 13
12 = [1 1 0 0] -> [1 0 1 0] = 10
13 = [1 1 0 1] -> [1 0 1 1] = 11
14 = [1 1 1 0] -> [1 1 1 0] = 14
15 = [1 1 1 1] -> [1 1 1 1] = 15

    */

    sector = 0;
    nstock_[sector][0] = 0;
    ninv_[0] = 0;
    sector_[0] = sector; //  0 = [0 0 0 0] -> [0 0 0 0] =  0
    nstock_[sector][1] = 5;
    ninv_[5] = 1;
    sector_[5] = sector; //  3 = [0 0 1 1] -> [0 1 0 1] =  5
    nstock_[sector][2] = 10;
    ninv_[10] = 2;
    sector_[10] = sector; // 12 = [1 1 0 0] -> [1 0 1 0] = 10
    nstock_[sector][3] = 15;
    ninv_[15] = 3;
    sector_[15] = sector; // 15 = [1 1 1 1] -> [1 1 1 1] = 15

    sector = 1;
    nstock_[sector][0] = 1;
    ninv_[1] = 0;
    sector_[1] = sector; //  1 = [0 0 0 1] -> [0 0 0 1] =  1
    nstock_[sector][1] = 11;
    ninv_[11] = 1;
    sector_[11] = sector; // 13 = [1 1 0 1] -> [1 0 1 1] = 11

    sector = 2;
    nstock_[sector][0] = 4;
    ninv_[4] = 0;
    sector_[4] = sector; //  2 = [0 0 1 0] -> [0 1 0 0] =  4
    nstock_[sector][1] = 14;
    ninv_[14] = 1;
    sector_[14] = sector; // 14 = [1 1 1 0] -> [1 1 1 0] = 14

    sector = 3;
    nstock_[sector][0] = 2;
    ninv_[2] = 0;
    sector_[2] = sector; //  4 = [0 1 0 0] -> [0 0 1 0] =  2
    nstock_[sector][1] = 7;
    ninv_[7] = 1;
    sector_[7] = sector; //  7 = [0 1 1 1] -> [0 1 1 1] =  7

    sector = 4;
    nstock_[sector][0] = 3;
    ninv_[3] = 0;
    sector_[3] = sector; //  5 = [0 1 0 1] -> [0 0 1 1] =  3

    sector = 5;
    nstock_[sector][0] = 6;
    ninv_[6] = 0;
    sector_[6] = sector; //  6 = [0 1 1 0] -> [0 1 1 0] =  6

    sector = 6;
    nstock_[sector][0] = 8;
    ninv_[8] = 0;
    sector_[8] = sector; //  8 = [1 0 0 0] -> [1 0 0 0] =  8
    nstock_[sector][1] = 13;
    ninv_[13] = 1;
    sector_[13] = sector; // 11 = [1 0 1 1] -> [1 1 0 1] = 13

    sector = 7;
    nstock_[sector][0] = 9;
    ninv_[9] = 0;
    sector_[9] = sector; // 9 = [1 0 0 1] -> [1 0 0 1] =  9

    sector = 8;
    nstock_[sector][0] = 12;
    ninv_[12] = 0;
    sector_[12] = sector; // 10 = [1 0 1 0] -> [1 1 0 0] = 12

    /*

    // Old states definition with orbital and spin index in the wrong order

    sector = 0;
    nstock_[sector][0] =  0; ninv_[ 0] = 0; sector_[ 0] = sector;   //  0 = 0000
    -> 0000 =  0
    nstock_[sector][1] =  3; ninv_[ 3] = 1; sector_[ 3] = sector;   //  3 = 1100
    -> 1010 =
    nstock_[sector][2] = 12; ninv_[12] = 2; sector_[12] = sector;   // 12 = 0011
    -> 0101 =
    nstock_[sector][3] = 15; ninv_[15] = 3; sector_[15] = sector;   // 15 = 1111
    -> 1111 = 15

    sector = 1;
    nstock_[sector][0] =  1; ninv_[ 1] = 0; sector_[ 1] = sector;   //  1 = 1000
    -> 1000 =  1
    nstock_[sector][1] = 13; ninv_[13] = 1; sector_[13] = sector;   // 13 = 1011
    -> 1101 =

    sector = 2;
    nstock_[sector][0] =  2; ninv_[ 2] = 0; sector_[ 2] = sector;   //  2 = 0100
    -> 0010 =
    nstock_[sector][1] = 14; ninv_[14] = 1; sector_[14] = sector;   // 0111

    sector = 3;
    nstock_[sector][0] =  4; ninv_[ 4] = 0; sector_[ 4] = sector;   // 0010
    nstock_[sector][1] =  7; ninv_[ 7] = 1; sector_[ 7] = sector;   // 1110

    sector = 4;
    nstock_[sector][0] =  5; ninv_[ 5] = 0; sector_[ 5] = sector;   // 1010

    sector = 5;
    nstock_[sector][0] =  6; ninv_[ 6] = 0; sector_[ 6] = sector;   // 0110

    sector = 6;
    nstock_[sector][0] =  8; ninv_[ 8] = 0; sector_[ 8] = sector;   // 0001
    nstock_[sector][1] = 11; ninv_[11] = 1; sector_[11] = sector;   // 1101

    sector = 7;
    nstock_[sector][0] =  9; ninv_[ 9] = 0; sector_[ 9] = sector;   // 1001

    sector = 8;
    nstock_[sector][0] = 10; ninv_[10] = 0; sector_[10] = sector;   // 0101
    */

    // ----------------------------------------------------------------

    init_npart_sig();
    init_c_cdag();
    init_vacuum();
  }

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

  void init_vacuum() {
    // sector 0 is the vacuum
    cdvector v1(1);
    vacuum_.InitZero(ssdim_);
    v1(0) = 1.0;
    vacuum_.M_[0] = v1;
  }
};

// ********************************************************************
class HilbertSpaceBose : public HilbertSpaceFermi {

public:
  void init(int Nmax = 4) {

    norb_ = 1;
    spin_degeneracy_ = 1;
    nh_ = Nmax;
    ns_ = 1;

    // ----------------------------------------------------------------

    init_storage();

    ssdim_[0] = nh_;

    for (int sector = 0; sector < ns_; sector++)
      nstock_[sector].resize(ssdim_[sector]);

    // ----------------------------------------------------------------

    int sector = 0;
    for (int psi = 0; psi < nh_; psi++) {
      nstock_[sector][psi] = psi;
      ninv_[psi] = psi;
      sector_[psi] = 0;
    }

    // ----------------------------------------------------------------

    init_npart_sig();
    init_c_cdag();
    init_vacuum();
  }

  void init_npart_sig() {
    int sector = 0;
    npart_[sector] = nh_;
    sig_[sector] = 1;
  }

  void init_vacuum() {
    // sector 0 is the vacuum
    cdvector v1(1);
    vacuum_.InitZero(ssdim_);
    v1(0) = 1.0;
    vacuum_.M_[0] = v1;
  }

  // save c and cdag into a full Blockmatrix
  void set_operator_matrix_c(int orb, int spin, cOp &op_c) {
    int sector = 0;
    std::vector<int> tsector(1, 0);
    op_c.InitZero(ssdim_, tsector);
    for (int n = 1; n < nh_; n++)
      op_c.M_[sector](n - 1, n) = std::sqrt(n);
  }

  void set_operator_matrix_cdag(int orb, int spin, cOp &op_cdag) {
    int sector = 0;
    std::vector<int> tsector(1, 0);
    op_cdag.InitZero(ssdim_, tsector);
    for (int n = 1; n < nh_; n++)
      op_cdag.M_[sector](n, n - 1) = std::sqrt(n);
  }
};
