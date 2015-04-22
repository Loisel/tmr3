# Makefile for the calculation of the rates
# 
# make is used primarily because its simple parallelization mechanism
# make -j numberofcores
# 
# We use this to distribute the jobs for the calculation of the different rates
# e.g. in-/out-tunneling right/left lead and parallel/anti-parallel
# configuration of the leads
#
# This makes a total of 8 jobs per TMR trace

# The argv parameter contains the jobname. It will be overwritten
# by the call to make
argv 	= test

# We operate in the dat folder
fname 	= dat/$(argv)

# The script to calculate the rates
exec   	= rates.py

## FORTRAN part
# The runtime critical part of the code is written in Fortran
# We use f2py to interface python and Fortran.
# Compiler options
CC = gfortran
CFLAGS = -Wall -O3 -pg -g -ffree-line-length-none -lslatec -L/usr/local/lib/libslatec
OFLAGS = -c -fPIC -O3

# geval2 is the new version of the Fortran rate calculation routine.
# It supports multiple excited state shells.
geval2.so: geval2.f95 quadpack.o
	f2py -c -m geval2 --fcompiler=gfortran --f90flags="$(CFLAGS)" -lslatec -L/usr/local/lib/libslatec quadpack.o geval2.f95 only: gp_l gm_l init selfenergy :	

# the old version of the rate calculation code
geval.so: geval.f95 quadpack.o
	f2py -c -m geval --fcompiler=gfortran --f90flags="$(CFLAGS)" -lslatec -L/usr/local/lib/libslatec quadpack.o geval.f95 only: gp_l gm_l init selfenergy testcpsi :	

# we need to use the quadpack that contains the integration (quadrature) routine
quadpack.o: quadpack.f90
	$(CC) $(OFLAGS) quadpack.f90

## END FORTRAN part

# The targets to prepare the Fortran object files
geval2: geval2.so
geval: geval.so

## TMR calculations
# this is a tree of dependencies. In total, 8 rates are calculated.
# and tmr.py is called with the job name as parameter
tmr: rates
	./tmr.py $(argv)

# parallel and anti-parallel
rates: par apa

# ... in- and out-tunning in the parallel
par: par_p par_m

# ... in-tunneling to left and right lead
par_p: par_p0 par_p1

# ... and this is the actual call: 
# a: is the configuration
# i: is the in/out tunneling coefficient
# l: is the left/right lead
par_p0: geval2.so
	./$(exec) -a 0 -i 0 -l 0 $(argv) 
par_p1: geval2.so
	./$(exec) -a 0 -i 0 -l 1 $(argv) 

par_m: par_m0 par_m1

par_m0: geval2.so
	./$(exec) -a 0 -i 1 -l 0 $(argv) 

par_m1: geval2.so
	./$(exec) -a 0 -i 1 -l 1 $(argv) 

apa : apa_p apa_m

apa_p: apa_p0 apa_p1

apa_p0: geval2.so  
	./$(exec) -a 1 -i 0 -l 0 $(argv)

apa_p1: geval2.so
	./$(exec) -a 1 -i 0 -l 1 $(argv)

apa_m: apa_m0 apa_m1

apa_m0: geval2.so  
	./$(exec) -a 1 -i 1 -l 0 $(argv)

apa_m1: geval2.so
	./$(exec) -a 1 -i 1 -l 1 $(argv)
