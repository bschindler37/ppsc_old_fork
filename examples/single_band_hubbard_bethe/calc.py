""" Single band Hubbard model on the Bethe lattice.

Example simulation of U-quench.

Author: H. U.R. Strand (2020) """


import os


from ReadCNTR import write_input_file, Dummy


def make_calc(U0=3.0, U=5.0, beta=5.0):
    
    inp = Dummy()

    inp.nt = 100
    inp.ntau = 20
    inp.beta = beta
    inp.h = 0.02
    inp.itermax = 100
    inp.errmax = 1e-12
    inp.iter_rtime = 5
    inp.order = 1
    inp.kt = 5
    # --
    inp.U0 = U0
    inp.U = U
    inp.mu = 0.0
    inp.eps = 0.0
    # --
    inp.store_pp = 0
    inp.linear_mixing = 0.5
    inp.read_eq_sym = 0
    inp.read_rt_sym = 0
    inp.read_state_from_file = 0

    path = 'calc_U0%2.2E_U%2.2E_beta%2.2E' % (inp.U0, inp.U, inp.beta)
    input_filename = 'input.txt'
    os.mkdir(path)
    os.chdir(path)
    write_input_file(inp.__dict__, input_filename)

    cmd = 'single_band_hubbard_bethe.ex ' + input_filename
    print(cmd)
    os.system(cmd)
    os.chdir('../')


if __name__ == '__main__':

    os.system('rm -r calc_*')
    make_calc(U0=5.0, U=10.0, beta=1.0)    

