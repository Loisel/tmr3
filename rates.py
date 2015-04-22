#!/usr/bin/python
"""
This module is used to calculate the rates for a specific combination
of contact magnetization, in- or out-tunneling and left or right lead.

Usage:
  ./rates.py -a <config> -i <in/out> -l <right/left> <jobname>

Example:
  ./rates.py -a 0 -i 1 -l 0 lookup_G_scale=1.60e-04_E_0=1.40e-03_Pol=4.00e-01
"""

import scipy as sp
from scipy import integrate
import time
import cProfile
import sys
import getopt
import os

#import matplotlib.pyplot as plt

# Let us be a little picky about division or invalid operations that
# would cause NaNs to appear (does not look at the Fortran code, of course)
sp.seterr(divide="raise")
sp.seterr(invalid="raise")

import h5py

# the module to do the actual rate calculation
# this is a python .so file compiled from Fortran with f2py
import geval2 as hlpr

import cfg

def save_rates_hdf5(fname,dname,G):
    """
    This is a utility routine to save a rate to a h5 file.

    Args:
      fname: The filename of the h5 file.
      dname: The name of the dataset to save the rate in.
      G: The actual one-dimenstional numpy array containing
         values for all gate voltages.
    """
    # print "saving dataset {} to file {}".format(dname,fname)
    fileh = h5py.File(fname,"a")
    fileh.create_dataset(dname,data=G)
    fileh.close()


def calculate_rates(leadidx,gamidx,pm,jobname):
    """
    Plots the rates as a function of gate voltage.

    Args:
      leadidx: The index of the lead: 0 or 1 for source or drain
      gamidx: The index of the configuration: 0 for parallel, 1 for anti-p.
      pm: The index of the in- or out-tunneling character of the rate.
      jobname: The name of the job. The directory dat/jobname should exists.
    
    """

    # Lead mask with respect to spin and orbit: LEAD[0] is source...
    LEAD = sp.array([[1.,1.,0.,0.],[0.,0.,1.,1.]])

    # The temporary calculation file
    if not os.path.isdir(os.path.join(os.getcwd(),'dat',jobname)):
        raise Exception('Job directory {} does not exist.'.format(jobname))
    filename = './dat/{}/running_calc.h5'.format(jobname)

    # This is the number of rates for the given config.
    # It will determine the size of the density matrix.
    # It can be different for the different configurations
    NG = cfg.TLIST[gamidx].shape[0]

    # The array to contain the rates, cfg.conf['NV'] is the number of gate
    # voltage steps
    G = sp.zeros((cfg.conf["NV"],NG))

    t0 = 0.

    # For the status output
    PM = ['P','M']
    pvsap = ['parallel','anti-parallel']

    # we have to iterate over the bias voltages
    for Vb in cfg.conf["Vb_range"]:

        status = "Calculating {0}-tun. rates in lead {1} for {2} contacts at bias {3}.".format(PM[pm],leadidx,pvsap[gamidx],Vb)

        # this is the intialisation of the Fortran module
        hlpr.geval.init(leadidx,\
                        Vb,\
                        0,\
                        cfg.conf["W0_t"][gamidx],\
                        cfg.conf["W0_b"][gamidx],\
                        cfg.conf["WE_t"][gamidx],\
                        cfg.conf["WE_b"][gamidx],\
                        cfg.conf["BETA"],\
                        cfg.conf["ALPHA"],\
                        cfg.conf["MU"],\
                        cfg.TLIST[gamidx],\
                        cfg.upstates[gamidx],\
                        cfg.downstates[gamidx],\
                        cfg.statemap[gamidx],\
                        cfg.E[:,gamidx],\
                        cfg.conf["gamma"][gamidx],\
                        cfg.conf["Int_CO"],\
                        cfg.conf["Int_acc_abs"],\
                        cfg.conf["Int_acc_rel"],\
                        cfg.conf["Int_gkp_key"])
        
        # ... and calculate the rates
        for n,Vg in enumerate(cfg.conf["V_g"]):
            if n == 0:
                t0 = time.time()
            if pm == 1:        
                G[n] = hlpr.geval.gm_l(Vg)
            elif pm == 0:
                G[n] = hlpr.geval.gp_l(Vg)
            # G[n] = hlpr.geval.g_l_analytic(Vg,pm)
            if n == 0:
                t1 = time.time()
                print status
                print "Time estimate for one rate: {0}".format((t1-t0)*cfg.conf["NV"])
                sys.stdout.flush()
            
        print "Calculations finished, Vb {}. Runtime: {}".format(Vb,time.time()-t0)
        
        # then the rates are saved to the h5 file. The name of the dataset
        # is constructed as follows:
        # <config>_<in/out><leadindex>_V<bias>
        save_rates_hdf5(filename,'{}_{}{}_V{}'.format(['par','apa'][gamidx],PM[pm],leadidx,Vb),G)

    # this is just a debugging switch...
    # if switched on, the rates are plotted
    if False:
        import matplotlib.pyplot as plt
        Ax = plt.subplot(1,1,1)
        
        cm = plt.get_cmap('gist_rainbow')

        Ax.set_color_cycle([cm(1.*k/NG) for k in range(NG)])

        for i in sp.arange(4,NG-4):
            color = cm(1.*i/NG)
            LABEL = "{}->{}".format(cfg.int_to_state(cfg.TLIST[gamidx][i,1]),cfg.int_to_state(cfg.TLIST[gamidx][i,0]))
            plt.plot( cfg.conf["V_g"], G[:,i],label=LABEL)

        plt.title("{} tun. rates for {} contacts".format(["in","out"][pm],pvsap[gamidx]))

        lines =Ax.get_lines()
        labels = [l.get_label() for l in lines]
        leg = plt.figlegend(lines,labels,loc='upper right')
        plt.show()

def calculate_selfenergy(b,a,E_array,Vb,jobname):
    """
    Plots the real and imaginary part of the self energy
    as a function of gate voltage.
    This function is deprecated, since there is a standalone
    test module for the selfenergy tests (test_se.py).

    Args:
      b,a: the states that we calculate the SE for.
      E_array: the energy window for the running parameter of the SE.
      Vb: the bias voltage
      jobname: the name of the job in the ./dat directory.
    
    """
    import matplotlib.pyplot as plt
    import cfg

    # we have to parse the configuration file first
    cfg.parse_conf('cnt.conf')
    cfg.init()

    Vg = cfg.E[cfg.state_to_int(b)]-cfg.E[cfg.state_to_int(a)][0]

    LEAD = sp.array([[1.,1.,0.,0.],[0.,0.,1.,1.]])

    SE = sp.zeros((2,E_array.shape[0],2))

    t0 = 0.

    PM = ['P','M']
    pvsap = ['parallel','anti-parallel']

    directory = "dat/{0}".format(args[0])
    if not os.path.exists(directory):
        os.makedirs(directory)

    
    for gamidx in sp.arange(2):
        for lead in sp.arange(2):
            hlpr.geval.init(lead,\
                    Vb,\
                    Vg,\
                    cfg.conf["W0_t"][gamidx],\
                    cfg.conf["W0_b"][gamidx],\
                    cfg.conf["WE_t"][gamidx],\
                    cfg.conf["WE_b"][gamidx],\
                    cfg.conf["BETA"],\
                    cfg.conf["ALPHA"],\
                    cfg.conf["MU"],\
                    cfg.TLIST[gamidx],\
                    cfg.TARRAY_BUBBLES[gamidx],\
                    cfg.get_P(cfg.TARRAY_BUBBLES[gamidx]),\
                    cfg.CHARGE,\
                    cfg.E[:,gamidx],\
                    cfg.conf["gamma"][gamidx],\
                    cfg.conf["Int_CO"],\
                    cfg.conf["Int_acc_abs"],\
                    cfg.conf["Int_acc_rel"],\
                    cfg.conf["Int_gkp_key"])
        
            for nE,E in enumerate(E_array):
                SE[gamidx,nE] += hlpr.geval.selfenergy(E,cfg.state_to_int(b),cfg.state_to_int(a))
            
        print "Calculations finished. Runtime: {0}".format(time.time()-t0)

    
        f = open('dat/{0}/SE_{1}.dat'.format(args[0],['par','apa'][gamidx]), 'w')

        s = "# Self energy\n\
        # for {} contacts\n\
        # at bias {}\n".format(pvsap[gamidx],Vb)

        f.write(s)
    
        sp.savetxt(f,sp.column_stack((E_array,SE[gamidx,:,0],SE[gamidx,:,1])),fmt='%.18e', delimiter=' ')

        f.close()

    plt.plot(E_array,SE[0,:,0],'b',E_array,SE[0,:,1],'g')
    plt.plot(E_array,SE[1,:,0],'b--',E_array,SE[1,:,1],'g--')

    plt.show()
    


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    """
    Entry point for the rate calculation module.

    Usage:
      ./rates.py -a <config> -i <in/out> -l <right/left> <jobname>
    """
    leadidx = 0
    apar = 0

    pm = 0
    cfile = ""
    
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], 'hl:a:i:c:', ['help','lead=','apar=','inout=','config='])
            for o,a in opts:
                if o in ('-h','--help'):
                    usage()
                    exit()
                elif o in ('-l','--lead'):
                    leadidx = int(a)
                elif o in ('-a','--apar'):
                    apar = int(a)
                elif o in ('-i','--inout'):
                    pm = int(a)
                elif o in ('-c','--config'):
                    cfile = a
                else:
                    raise Usage('Invalid argument.')
        except getopt.error, msg:
             raise Usage(msg)
        # more code, unchanged
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

    # the configuration file is supposed to be in the job directory
    # with the name cnt.conf. It can also be specified via the -c|--config
    # parameter
    cfgpath = cfile if cfile else "dat/{0}/cnt.conf".format(args[0]) 

    # parsing and initializing the configuration file
    cfg.parse_conf(cfgpath)

    cfg.init()
    
    # ... and call the rate calculation routine
    calculate_rates(leadidx,apar,pm,args[0])

    
def usage():
	print "\n\
Calculate the rates for the Master equation of the CNT-QD.\n\
	\n\
	usage: rates.py [-hp] [-l|--lead lead] [-a|--apar apar] [-V bias] [-i|--inout inout] jname\n\
	\n\
	-h \t print help\n\
	-p \t enable plotting\n\
        --print \t prints the config\n\
	-l, --lead \t specify the lead (0 or 1)\n\
	-a, --apar \t specify the polarization of the contacts (parallel:0, antiparallel:1)\n\
	-i, --inout \t for a in (0,default) or out (1) tunneling rate\n\
	\n\
        -c, --config \t the path to the config file\n\
        \n\
       	-V \t the bias voltage (default is 0.)\n\
	jname is the basename of the resulting rate files.\n"

if __name__ == "__main__":
    sys.exit(main())
        
