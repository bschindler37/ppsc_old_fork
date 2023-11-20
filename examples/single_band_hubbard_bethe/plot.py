""" Single band Hubbard model on the Bethe lattice.

Example plot of double occumpancy, density, and energy components
after a U-quench.

Author: H. U.R. Strand (2020) """


import glob
import numpy as np
import matplotlib.pyplot as plt


from ReadCNTR import read_input_file
from ReadCNTRhdf5 import read_imp_h5file


if __name__ == '__main__':

    paths = glob.glob('./calc*')
    print(paths)
    
    for path in paths:

        imp_filename = path + '/data_ppsc.h5'
        print('--> Loading:', imp_filename)
        d = read_imp_h5file(imp_filename)
        d.imp.parm =  read_input_file(path + '/input.txt')

        print(d.__dict__.keys())
        print(d.bethe.__dict__.keys())
        print(d.imp.__dict__.keys())

        d.imp.t = np.arange(-1, d.imp.nt + 1) * d.imp.h
        d.imp.n_exp = d.imp.nu_exp + d.imp.nd_exp
        
        d.Etot = d.imp.Eint_exp + d.bethe.Ekin

        plt.figure(figsize=(4, 10))
        subp = [3, 1, 1]
    
        plt.subplot(*subp); subp[-1] += 1
        plt.title(r'$\beta = %2.1f$' % d.imp.beta)
        plt.plot(d.imp.t[1:], d.imp.docc_exp[1:], '-',
                 label=r'$U_0=%2.0f, \, U=%2.0f$' \
                 % (d.imp.U[0], d.imp.U[-1]))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$\langle n^2 \rangle$')
        plt.legend(loc='best')
        
        plt.subplot(*subp); subp[-1] += 1
        plt.plot(d.imp.t[1:], d.Etot[1:], label=r'$E_{tot}$')
        plt.plot(d.imp.t[1:], d.bethe.Ekin[1:], label=r'$E_{kin}$')
        plt.plot(d.imp.t[1:], d.imp.Eint_exp[1:], label=r'$E_{int}$')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$E_{tot}$, $E_{kin}$, $E_{int}$')
        plt.legend(loc='best')

        plt.subplot(*subp); subp[-1] += 1
        plt.plot(d.imp.t[1:], np.abs(d.Etot[1:] - d.Etot[1]), label=r'$\Delta E_{tot}$')
        plt.plot(d.imp.t[1:], np.abs(d.imp.n_exp[1:] - d.imp.n_exp[0]), label=r'$\Delta n_{tot}$')
        plt.semilogy([], [])
        plt.xlabel(r'$t$')
        plt.legend(loc='best')
        
    plt.tight_layout()
    plt.savefig('figure_U_quench.pdf')
    plt.show()
            
