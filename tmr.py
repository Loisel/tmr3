#!/usr/bin/python
"""
A module to calculate the current, the conductance and the TMR from
a set of rate arrays.
The rate arrays are supposed to be stored in a h5 file in the job directory.
The result is stored in a h5 file. The name of the dataset contains all
parameters. They are also stored as attributes in the dataset.
The conductance in the two lead configurations (parallel/anti-parallel)
are stored in arrays in the dataset.

Usage:
  ./tmr.py <jobname>
"""

import numpy as np
from numpy import linalg

import time
import sys
import getopt
import h5py
import os

# We are picky about possible floating point overflows
# to avoid calculating NaNs
np.seterr(divide="raise")
np.seterr(invalid="raise")

# A helper module to calculate the populations.
import pop

# The configuration module
import cfg

# path to the dat directory
datpath = "dat/"
# name of the temporary file where the rates are stored
ratefile = "running_calc.h5"
# name of the h5 file to store the conductance for the two configuration
# and the configuraion parameters.
hdffile = "simdata_new.h5"


def save_hdf5(fname,G_P,G_AP):
    """
    Store the conductance and the configuration to the h5 file.

    Args:
      fname: filename of the h5 file
      G_P: the conductance for leads with parallel magnetization
      G_AP: the conductance for leads with anti-parallel magnetization
    """
    print "Shape of GP {}".format(G_P.shape)
    
    fileh = h5py.File(fname,"a")

    # Note that the selection of parameters to construct the name of the
    # dataset should be chosen such that this string is unique!
    # That is, it should contain all running parameters.
    dset_name = "G={}_kbT={}_Ec={}_E0={}_Pol={}_PolOrb={}_SO={}_tau={}_DS={}_B_P={}_B_AP={}_B_ORB_P={}_B_ORB_AP={}_W_e={}_W_0={}".format(cfg.conf['G_scale'],cfg.conf['kBT'],cfg.conf['E_C'],cfg.conf['E_0'],cfg.conf['Pol'],cfg.conf['OrbPol'],cfg.conf['SO'],cfg.conf['tau_r'],cfg.conf['D_S_factor'],cfg.conf['B_P'],cfg.conf['B_AP'],cfg.conf['B_ORB_P'],cfg.conf['B_ORB_AP'],cfg.conf['W_E'],cfg.conf['W_0'])

    try:
        # we create the dataset
        dset = fileh.create_dataset(dset_name,data=np.vstack((G_P,G_AP)))

        # and store the config attributes
        dset.attrs['alpha'] = cfg.conf['ALPHA']
        dset.attrs['temperature'] = cfg.conf['kBT']
        dset.attrs['coupling'] = cfg.conf['G_scale']
        dset.attrs['electron_number'] = cfg.conf['N_0']
        dset.attrs['charging_energy'] = cfg.conf['E_C']
        dset.attrs['level_spacing'] = cfg.conf['E_0']
        dset.attrs['polarization_spin'] = cfg.conf['Pol']
        dset.attrs['polarization_orbit'] = cfg.conf['OrbPol']
        dset.attrs['spinorbit'] = cfg.conf['SO']
        dset.attrs['stonershift'] = cfg.conf['D_S_factor']
        dset.attrs['tau_r'] = cfg.conf['tau_r']
        dset.attrs['vg_min'] = cfg.conf['V_g_min']
        dset.attrs['vg_max'] = cfg.conf['V_g_max']
        dset.attrs['b_p'] = cfg.conf['B_P']
        dset.attrs['b_ap'] = cfg.conf['B_AP']
        dset.attrs['b_orb_p'] = cfg.conf['B_ORB_P']
        dset.attrs['b_orb_ap'] = cfg.conf['B_ORB_AP']
        dset.attrs['w_0'] = cfg.conf['W_0']    
        dset.attrs['w_e'] = cfg.conf['W_E']    
        dset.attrs['timestamp'] = time.time()
    except KeyError:
        # If the choice was not unique we complain but continue.
        print "Dataset exists."

    fileh.close()

def eval_DENKER(GM,GP,configuration):
    """
    Evaluate the density matrix kernel using the in- and out-tunneling rates.
    
    Args: 
      GM,GP: numpy arrays containing in- and out-tunneling rates
             in the order of cfg.TLIST.
      configuration: integer determining parallel (0) or anti-parallel(1)
                     configuration
    
    Returns:
      the density matrix as a square 2-d numpy array that is NP**2 in size,
      where NP is the number of states in the groundstatespace.
    """

    # we get a view on the transition list and, for simplicity, its transpose
    TLIST = cfg.TLIST[configuration]
    TLIST_T = np.transpose(TLIST)

    # from all transitions we extract all groundstates in the statespace
    # this is probably a complicated way to do it
    PLIST = list(set(TLIST_T[0]).union(TLIST_T[1]))
    # ... and sort it by index
    PLIST.sort()

    # the number of groundstates
    NP = len(PLIST)

    # let's create an empty density matrix
    ME = np.zeros((NP,NP))

    # we create a version of the transition list that does not contain
    # the indices in terms of the energy array (see cfg.py), but
    # in terms of the number in the state list (plist)
    # (the transition list can then be used to denote non-zero matrix elements)
    TMP = np.copy(TLIST)
    for idx,val in enumerate(PLIST):
        TMP[TLIST == val] = idx

    # We calculate diagonal elements of the density matrix:
    # TLIST_T[1] == num selects the correct in-tunneling rates for the
    # state with label num
    # have a look at numpy.where to understand this line
    
    for idx,num in enumerate(PLIST):
        ME[idx,idx] = -np.sum(np.where(TLIST_T[1] == num,GP,0.)) - np.sum(np.where(TLIST_T[0] == num,GM,0.))

    # for the off diagonal elements we can directly use the generated TMP
    # transition list
    
    for k,tup in enumerate(TMP):
        ME[tup[0],tup[1]] = GP[k]
        ME[tup[1],tup[0]] = GM[k]
#        print "tup: {} and matrix element {}".format(tup,ME[tuple(tup)])

    return ME

def eval_CURKER(GM,GP,configuration):
    """
    Evaluate the current kernel using the in- and out-tunneling rates.

    Args:
      GM,GP: numpy arrays containing in- and out-tunneling rates
             in the order of cfg.TLIST.
      configuration: integer determining parallel (0) or anti-parallel(1)
                     configuration
    
    Returns:
      the current kernel as a 1-d numpy array.
    """

    # We get a view on the transition list and its transpose
    TLIST = cfg.TLIST[configuration]
    TLIST_T = np.transpose(TLIST)

    # ... and extract the list of groundstates (see also eval_DENKER)
    PLIST = list(set(TLIST_T[0]).union(TLIST_T[1]))
    PLIST.sort()

    # this determines the size of the statespace
    NP = len(PLIST)
    CUR = np.zeros(NP)
    
    # Note that the current kernel can be calculated by summing the diagonal elements
    # of the density matrix with opposite sign
    # compare eval_DENKER
    
    for idx,num in enumerate(PLIST):

        CUR[idx] = np.sum(np.where(TLIST_T[1] == num,GP,0.)) - np.sum(np.where(TLIST_T[0] == num,GM,0.))
    return CUR


def current(GP,GM,POP,configuration):
    """
    Calculate the current using the rates and populations.

    Args:
      GP, GM: np-arrays containing in- and out-tunneling rates.
      POP: np-array for the populations
      configuration: integer determining parallel (0) or anti-parallel(1)
                     configuration

    Returns:
      current as a float. 
    """

    # We calculate the current kernel
    CURKER = eval_CURKER(GM,GP,configuration)

    # and vector-multiply it with the population vector
    I = -np.sum(cfg.conf["ELE"]*np.dot( CURKER, POP))

    return I
    



def eval_tmr(fname,plotname,pop):
    """
    Calculates the TMR by evaluating conductance through
    parallel and anti-parallel polarized contacts.

    Args:
      fname: the h5 file to load the rates from.
      plotname: A name for the pdf output to produce.
      pop: If True, we plot the populations, too.
    
    """

    # We prepare the current and conductance vectors for different
    # values of gate and bias voltage
    
    C_p = np.zeros((cfg.conf['NV'],cfg.conf['NVb']))
    C_ap = np.zeros((cfg.conf['NV'],cfg.conf['NVb']))
    G_p =  np.zeros((cfg.conf['NV'],cfg.conf['NVb']-1))
    G_ap = np.zeros((cfg.conf['NV'],cfg.conf['NVb']-1))
    dVb = cfg.conf['Vb_range'][1]- cfg.conf['Vb_range'][0]

    # the population vectors, for all values of gate and bias
    POP_p =  np.zeros((cfg.conf['NVb'],cfg.conf['NV'],cfg.N_GS[0]))
    POP_ap = np.zeros((cfg.conf['NVb'],cfg.conf['NV'],cfg.N_GS[1]))

    # We iterate over two bias values first
    for nV,Vb in enumerate(cfg.conf["Vb_range"]):
        # now the rates are loaded from the h5 file
        # note that the label of the specific rate arrays are fixed
        with h5py.File(fname) as file:
            GP0_p = np.array(file['par_P0_V{}'.format(Vb)])
            GP0_ap = np.array(file['apa_P0_V{}'.format(Vb)])
            GP1_p = np.array(file['par_P1_V{}'.format(Vb)])
            GP1_ap = np.array(file['apa_P1_V{}'.format(Vb)])
            GM0_p = np.array(file['par_M0_V{}'.format(Vb)])
            GM0_ap = np.array(file['apa_M0_V{}'.format(Vb)])
            GM1_p = np.array(file['par_M1_V{}'.format(Vb)])
            GM1_ap = np.array(file['apa_M1_V{}'.format(Vb)])

        # for the density kernel, we sum all rates over both leads
        DENKER_p = np.array([eval_DENKER(GM0_p[n]+GM1_p[n],GP0_p[n]+GP1_p[n],0)for n in range(cfg.conf["NV"])])

        DENKER_ap = np.array([eval_DENKER(GM0_ap[n]+GM1_ap[n],GP0_ap[n]+GP1_ap[n],1)for n in range(cfg.conf["NV"])])

        # the populations are calculated from the density kernel by an asymptotic
        # approximation scheme
        POP_ap[nV] = np.array([pop.asymptotic_ssp(DENKER_ap[n]) for n in range(cfg.conf["NV"])])
        POP_p[nV] = np.array([pop.asymptotic_ssp(DENKER_p[n]) for n in range(cfg.conf["NV"])])
        
        # note that the current is calculated from the rates in one of the leads only
        C_p[:,nV] = np.array([ current(GP0_p[n],GM0_p[n],POP_p[nV,n],0) for n in np.arange(cfg.conf["NV"]) ])
        C_ap[:,nV] = np.array([ current(GP0_ap[n],GM0_ap[n],POP_ap[nV,n],1) for n in np.arange(cfg.conf["NV"]) ])
        
        
    # the numerical derivative gives the conductance
    G_p = np.diff(C_p).flatten()/dVb
    G_ap = np.diff(C_ap).flatten()/dVb

    # we save the conductance traces to a h5 file specified as a global variable
    # hdffile in the path datpath
    # It is possible that the dataset already exists. In this case, we issue a warning.
    try:
        save_hdf5("{}{}".format(datpath,hdffile),G_p,G_ap)
    except RuntimeError:
        print "Unable to save to {}, maybe there is already a dataset with similar parameters...".format(hdffile)

    # the tmr and conductance graphs are plotted to a pdf file for review.
    plot_tmr_pdf(G_p,G_ap,plotname)

    # if the pop flag is set, we also plot the population for one bias value
    if pop:
        plot_population([POP_p[0],POP_ap[0]],os.path.splitext(plotname)[0]+"_POP.pdf")

        
def plot_tmr_pdf(C_p,C_ap,fname):
    """
    A helper routine to plot the conductance and TMR to a pdf file in the datpath.

    Args:
      C_p, C_ap: the parallel and anti-parallel conductance.
      fname: the filename to plot to
    """
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    # we plot the conductance graph on top, p and ap with different colors
    Axes1 = plt.subplot(2,1,1)

    Axes1.set_xticklabels([])
    plt.ylabel("Conductance (e^2/h)")

    plt.title("Conductance at zero bias")

    # parallel is plotted in red, and anti-parallel as blue dashed line
    plt.plot( cfg.conf["V_g"],C_p,'r',cfg.conf["V_g"], C_ap, 'b--')

    # on the second panel, the TMR is plotted
    Axes2 = plt.subplot(2,1,2)

    plt.xlabel("gate voltage (V)")
    plt.ylabel("TMR")
    plt.title("TMR")

    plt.ylim((-0.3,1.5))
    
    TMR = np.zeros(cfg.conf["NV"])

    for i in range(cfg.conf["NV"]):
        try:
            TMR[i] = C_p[i]/C_ap[i]-1.
        except ZeroDivisionError:
            print "Zero Division, returning null."
            TMR[i] = 0.
    plt.plot( cfg.conf["V_g"], TMR)
    plt.savefig(fname, bbox_inches='tight')



def plot_population(POP, fname):
    """
    Calculates and plots selected populations of the quantum dot
    with gate voltage. The edge states N=-1 and 5 are neglected.

    Args:
      POP: a list with the two population vectors
           for parallel and anti-parallel configurations
      fname: the filename to plot to

    """
    import matplotlib.pyplot as plt

    NV = cfg.conf["NV"]

    
    print "Calculating populations..."

    # We plot the populations for both configurations
    # the parallel populations on top
    # the anti-parallel on bottom
    Ax = [plt.subplot(2,1,1),plt.subplot(2,1,2)]

    cm = plt.get_cmap('gist_rainbow')

    PopPlots = [1,4,8,12,17,18]
    NP = len(PopPlots)

    for gamidx in range(2):
        TLIST = cfg.TLIST[gamidx]
        TLIST_T = np.transpose(TLIST)
    
        PLIST = list(set(TLIST_T[0]).union(TLIST_T[1]))
        PLIST.sort()

        # we cycle through the linecolors to distinguish the different
        # groundstates
        Ax[gamidx].set_color_cycle([cm(1.*k/NP) for k in range(NP)])

        for i in PopPlots:
            color = cm(1.*i/NP)
            LABEL = "P_{}".format(cfg.int_to_state(PLIST[i]))
            Ax[gamidx].plot( cfg.conf["V_g"], POP[gamidx][:,i],label=LABEL)


        lines =Ax[gamidx].get_lines()
        labels = [l.get_label() for l in lines]
        leg = plt.figlegend(lines,labels,loc='upper right')

    plt.savefig(fname)
    plt.show()

        

    

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    """
    Interface routine to call the tmr module.
    Example:
      ./tmr.py <jobname>

    In principle, there were routines to plot rates, populations,
    conductances etc. but apart from the population plotting,
    none of the use cases was needed anymore. 
    
    """
    POP = False

    # The default config file is called cnt.conf
    cfile = "cnt.conf"
    rlist = [0.,]
    
    
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hc:P", ["help","config=","pop"])
        except getopt.error, msg:
             raise Usage(msg)
        for o,a in opts:
            if o in ('-h','--help'):
                usage()
                exit()
            elif o in ('-c','--config'):
                cfile = a
            elif o in ('-P','--pop'):
                POP = True
            else:
                raise Usage('Invalid argument.')

        # we parse the config and initialize it
        cfg.parse_conf("dat/{0}/{1}".format(args[0],cfile))
        cfg.init()


        h5file = "{}{}/{}".format(datpath,args[0],ratefile)
        pdffile = "{}{}.pdf".format(datpath,args[0])
        print "Try to open {}".format(h5file)
        eval_tmr(h5file,pdffile,POP)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2


def usage():
    print "This is a tool to process rate files.\n\
    \n\
    usage: tmr.py [-hP] [--pop] jobname\n\
    \n\
    --pop or -P: Plot the populations.\n\
    \n\
    jobname: The name of the directory for the rate files.\n\
    \n\
    The script searches for files dat/jobname/running_calc.h5\n\
    and dat/jobname/cnt.conf"
        
if __name__ == "__main__":
    sys.exit(main())
        
