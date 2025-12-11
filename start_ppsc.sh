#!/bin/bash


#################################################################################################
#
#   some sample setup for a ppsc  calculation based on the ppsc strong coupling solver
#   
#   code:   git@bitbucket.org:computationalphysicsunifr/ppsc.git 
#   branch: ppcs_rixs
#
#
################################################################################################


################################################################################################
# program directory:
exe_dir=./build/programs

# working directory:
wdir=./data

################################################################################################
#   driven Hubbard model with a bethe density of states
#   scr code is ppsc/programs/impurity_one_non-int_bathsite.cpp

#   Parameters for simulation:
nt=1000  #number of timesteps
h=0.01   #timestep
ntau=500 #imaginary times (beta/ntau ~ imag timestep ... similar to h ) 
beta=5
# interaction defined as U(nup-1/2)*(ndo-1/2) + (epsd-mu) * (nup+ndo) 
U=8 # Hubbard U
mu=0 #half filling if epsd=0
epsd=0
eps_bath=0
J0=1.0 # hopping

# create a file containing the time-dependent vector potential (-> Peierls phase)
# the file has nt+2 entries where entry A(-1),A(0),...,A(nt), where A(-1)=0 is for imag branch
afile=${wdir}/A.txt
# a periodic drive with a ramp on over some periods
a_max=0
n_ramp=1
omega=3  # smaller than the Hubbard gap

awk -v h=$h -v nt=$nt -v a_max=$a_max -v n_ramp=$n_ramp -v w=$omega 'BEGIN{\
       pi=3.1415192;\
       tau=2*pi/w*n_ramp;\
       print 0;\
       for(j=0;j<=nt;j++){\
	  t=h*j;\
	  ramp=(t>tau ? 1 : sin(0.5*pi*t/tau)*sin(0.5*pi*t/tau));\
	  at=a_max*0.5*(1-cos(w*t))*ramp;\
	  printf("%.12f\n",at);\
       }}' > $afile

#   the following code produced an input file; uncomment and replace by your own
input_file=${wdir}/input.txt
cat <<EOL > "$input_file"
__nt= $nt
__ntau= $ntau
__beta= $beta
__h= $h
__itermax= 100
__errmax= 1e-6
__linear_mixing= 0.2
__iter_rtime= 5
__kt= 5
__order= 1
__store_pp= 1
__read_state_from_file= 0
__nomp= 1
__mu= $mu
__U= $U
__U0= $U
__eps= $epsd
__eps_bath= $eps_bath
__J0_real= $J0
__vector_potential= --$afile
__atomic_limit= 0
__text_output= 0
EOL

# further notes on params:
# * itermax,errmax control initial equilibrium simulation
# * iter_rtime controls the DMFT iteration at each timestep (fixed number of iterations)
# * kt is the integration order ... should be always maximum (5)
# * order is 1 for NCA and 2 for OCA. OCA is impractable for longer simulations. use NCA
# * store_pp must be "1" to generate a big hdf5-file with the full information on the impurity model
#   including pseudoparticle Greensfunction; this is needed to read in the information
#   for the calculation of RIXS later
# * U and eps could be time dependent (see program)
# * choosing U0 different from U would implement an interaction quench
# * text_output=1 will write some files in plain text format (in addition to hdf5)
#   these files are much larger, the option is only for those who are too lazy to
#   use scrpts to extract data from hdf5, should not be used in production runs on cluster

#   now you can run the code ...
out_prefix=./data/u8/
valence_simulation=${exe_dir}/impurity_one_non-int_bathsite.ex
echo $valence_simulation $input_file $out_prefix
$valence_simulation $input_file $out_prefix

#  the code should produce the following output
# * {out_prefix}data_ppsc.h5: contains all info on impurity model, including the hybridization function
# * {out_prefix}Gimp.out: Impurity greens function (if text_output= 1)
# * {out_prefix}obs.out: single-time observables (description in first line of code)
# * [possibly more to be added later ...]

# for example use gnuplot to check energy conservation:
# ($1*0.02):(8*$5+$4)  to plot the sum of (Ekin+Epot) ... the drift is due to energy absorption
# ($1*0.02):(8*$5+$4-$11) to plot (Ekin+Epot - int_0^t j*E ) which is constant
# note that  int_0^t j*E is computed wit trapez, which causes the remaining drift in (8*$5+$4-$11)
# if you reduce timestep to 0.01, (8*$5+$4) is almost unchnaged (h=0.02 is ok for the actual simaultion)
# while $11 changes ... anyway, intjE is colculated only as quick check

# you can also reduce the amplitude of the drive to a_max=0.5 ... the energy absorbtion should be strongly
# reduced as the fraquency is in the gap and the absoption is therefore highly nonlinear in a_max
