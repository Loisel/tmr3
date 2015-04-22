
tmr3 Readme
~~~~~~~~~~~

tmr3 is the third version of a tool to calculate the conductance and the
TMR of CNT quantum dot contacted to ferromagnetic leads at (almost) zero bias.
The program was used to perform the numerical calculations for this paper:

[1] http://arxiv.org/abs/1502.02005  

Its structure arises, like most physic simulations, from needs of everyday
parameter scans and its "history", so many things might not be structured 
in a straightforward way. We added a section to this file to explain the 
building blocks and their interconnections.


Installation
------------

The runtime-critical part of the simulation is written in Fortran95.
From the source file geval2.f95 a python object file is compiled
via f2py. The parameters in the Makefile work for my system but
may have to be adjusted to compile correctly.
Compile the module via

  make geval2.so

and check for errors. 
The program depends on the digamma function routine

CPSI

provided by the SLATEC (http://www.netlib.org/slatec/).

Usage
-----

Let us first explain some simple use cases:

+ To reproduce the main result of the paper, i.e., Fig. 11 in [1], 
  execute the following line
  (given that you have 4 cores available, see numberofcores in lookup.py):

  python lookup.py

  All parameters should be set accordingly in the current version.

+ If you want to change some parameters, open lookup.py and search for the 
  definition of pdict. Here, you can specify parameter ranges for 
  different parameters.
  Executing lookup.py then generates all possible combinations of
  these parameters and starts to calculate it one by one using 4 cores 
  (number can be specified just below pdict).
  The current state of the calcaltion is saved in table.dat 
  (contains the list of parameter permutations) and 
  lnum.dat (contains the linenumber of the current parameterset in table.dat).

  IMPORTANT: If you start a new calculation, always remove lnum.dat
             and table.dat, otherwise the program tries to continue with 
             the previous calculations.

+ If you want to understand the program, read the paper and/or the
  thesis and then study the file cfg.py. 
  The important objects are generated there and if you understand
  them, the rest is quite simple.
  Then study geval.f95. It is written in Fortran and the actual rate 
  calculation, which is the most expensive part, is done there.


Structure:
----------

The structure of the simulation is as follows:

[lookup.py]
+ an interface for parameter scans
+ creates the job dir ./dat/jobname
   |
   |calls
   |
[make] ("Makefile")
+ triggers rate calculations 
+ when rates are ready, tmr is calculated
+ make is used to spawn multiple threads across the local cores
  basically it can also be replaced by a interface to a cluster
  qsub script. To this end, overwrite or extend the function
  submitJob in lookup.py
   |                                                |
   |calls                                           |
   |                                                |
[rates.py]                                          |calls
+ to calculate the rates.                           |
+ one run spawns 8 rates.py threads: 2 per lead,    |
  2 per config and 2 per in- or out-tunneling       |
   |                                 |              |
   |calls                            |provides      |
   |                                 |data          |
[geval2.f95]                         |for           | 
+ fortran module for rate integral   |              |
+ the "actual" computation           |              |
                                     |              |
                                     |              |
[tmr.py]
+ reads rates from ./dat/jobname/running_calc.h5
+ calculates density matrix and steady state solution for populations
+ current, conductance and tmr
+ generates plot ./dat/jobname.pdf
+ writes conductance and config to ./dat/simdata_new.h5
        |
        |uses
        |
[pop.py]
+ the module to calculate the populations 


Licence:
--------

the simulation is released under GNU GPLv3.
See LICENCE.txt.


Author:
-------

Alois Dirnaichner
alo.dir@gmail.com
