"""
Configuration module for the TMR calculation.
To understand the guts of the program, it is essential to understand
the objects which are defined here.
Important objects are

- the array of ground state transitions which
determines the rates that are calculated
- the list of groundstates
- the list of excited states
- the energy functional and and the list of ground state energies

It is usually used in a two step way, parsing a config file and initializing
the configuration.

Example:
  cfg.parse_conf('cnt.conf')
  cfg.init()

"""


import scipy as sp
import numpy as np

# Array masks for spin up/down for the states of the system
# and for the polarization of the leads    

SPINU = np.array([1.,0.,1.,0.])
SPIND = np.array([0.,1.,0.,1.])

SO = np.array([1.,-1.,-1.,1.])
S = np.array([1.,-1.,1.,-1.])
KK = np.array([-1.,-1.,1.,1.])

# global constants

conf = {}

# the list of transitions for the two configurations

TLIST = [[],[]]

# The array of ground state energies in the vicinity of one shell
# it is chosen such that we can cover the whole shell and -1 and +5 charging states.
# the rule is the following: be n a state vector specifying
# a dot configuration, then the number a = n[0]*4**3 + n[1]*4**2 + n[2]*4 + n[3]
# is a unique number between 0 and 255 labelling one of these states
# we call this number the STATENUMBER
E = np.zeros((256,2))

# a "mask" (not in the numpy sense) indicating transitions by non-zero elements
# for groundstate transitions
TARRAY_TRANSPORT = np.zeros((2,256,256))

# the charge for the possible groundstates
CHARGE = np.zeros(256)

# the list of groundstates labelled by the STATENUMBER
GROUNDSTATES = [[],[]]
N_GS = [0,0]

def parse_conf( infile ):
    """
    Reads a configfile and feeds the parameters into a config
    dictionary.

    Args:
      infile: the path to the config file
    """
    global conf
    with open( infile, 'rt') as file:
            for line in file:
                    line = line.strip()
                    if not line == '' and not line.startswith("#"):
                            name, nums = line.split('=',1)
                            conf[name.strip()] = float(nums.strip())

    


def print_conf():
    """
    Print the configuration.
    """
    STRING = "CNT configuration:\n"
    for k in conf.keys():
        if isinstance(conf[k],float) or isinstance(conf[k],int):
            STRING += "{}\t{}\n".format(k,conf[k])

    return STRING

def state_to_int(n):
    """
    Convert a state vector (a numpy array of size 4) to a unique
    STATENUMBER between 0 and 255.
    
    Example:
      In: state_to_int([1,1,0,1])
      Out: 166

    Args:
      n: np.array of shape (1,4)
    Returns:
      an integer between 0 and 255
    """
    # the state vector is always given with respect to the
    # main shell (there can be entries n[x] = -1)
    n += np.ones(4)
    return n[0]*4**3 + n[1]*4**2 + n[2]*4 + n[3]

def int_to_state(num):
    """
    The reverse operation of state_to_int. Given a STATENUMBER
    it returns a state vector.

    Args:
      num: integer
    Returns:
      np.array of shape (1,4)
    """
    num = int(num)
    n = np.zeros(4)
    div = 3
    for div in np.linspace(3,0,4):
        n[div] = num%4
        num = num/4
    # substract one on each position to give the electron number
    # relative to the main shell
    return n-np.ones(4)


def init():
    """
    Initializes the system parameters and states based on a previously
    parsed configuration.
    """
    global TARRAY_TRANSPORT
    global E
    global CHARGE
    global GROUNDSTATES
    global conf

    # the upstates and downstates arrays contain lists of energies
    # of states that can be reached from a specific groundstate
    # by in- or out-tunneling
    # the statemap is a map between these entries and the STATENUMBER
    # of the groundstate
    global upstates
    global downstates
    global statemap

    upstates = list()
    downstates = list()
    statemap = list()
    
    TARRAY_TRANSPORT = np.zeros((2,256,256))
    E = np.zeros((256,2))
    CHARGE = np.zeros(256)

    # the definition of the Stoner shift relative to the size of the bandwidth
    # the control parameter of this shift is D_S_factor, see Koller (2012)
    conf["D_S"] = conf["D_S_factor"]*conf["W_0"]/2.
    # we replace the bandwidth by a spin-split version
    conf["W0"] = np.array([conf["W_0"]-conf["D_S"] , conf["W_0"]+conf["D_S"]])
    # and then distinguish between the top and the bottom of the bandwidth
    # for the two configurations and the left and right lead
    # so W0_t[parallel/anti-parallel][left/right] gives a 4-vector of the
    # upper bandwidth boundary (depending on valley and spin)
    conf["W0_t"] = np.array([ [[conf["W0"][0],conf["W0"][1],conf["W0"][0],conf["W0"][1]],\
                              [conf["W0"][0],conf["W0"][1],conf["W0"][0],conf["W0"][1]]],\
                             [[conf["W0"][1],conf["W0"][0],conf["W0"][1],conf["W0"][0]],\
                              [conf["W0"][0],conf["W0"][1],conf["W0"][0],conf["W0"][1]]] ])


    conf["W0_b"] = np.array([ [[conf["W0"][1],conf["W0"][0],conf["W0"][1],conf["W0"][0]],\
                              [conf["W0"][1],conf["W0"][0],conf["W0"][1],conf["W0"][0]]],\
                             [[conf["W0"][0],conf["W0"][1],conf["W0"][0],conf["W0"][1]],\
                              [conf["W0"][1],conf["W0"][0],conf["W0"][1],conf["W0"][0]]] ])
    # We generate also a seperate bandwidth for excited states
    # this is just a additional possibility of fine-tuning.
    # It should coincide with the number of excited state transitions that is allowed
    # e.g. if NX is the number of excited states this should be around NX*E_0
    # where E_0 is the shell spacing
    conf["WE"] = np.array([conf["W_E"]-conf["D_S"] , conf["W_E"]+conf["D_S"]])

    conf["WE_t"] = np.array([ [[conf["WE"][0],conf["WE"][1],conf["WE"][0],conf["WE"][1]],\
                              [conf["WE"][0],conf["WE"][1],conf["WE"][0],conf["WE"][1]]],\
                             [[conf["WE"][1],conf["WE"][0],conf["WE"][1],conf["WE"][0]],\
                              [conf["WE"][0],conf["WE"][1],conf["WE"][0],conf["WE"][1]]] ])


    conf["WE_b"] = np.array([ [[conf["WE"][1],conf["WE"][0],conf["WE"][1],conf["WE"][0]],\
                              [conf["WE"][1],conf["WE"][0],conf["WE"][1],conf["WE"][0]]],\
                             [[conf["WE"][0],conf["WE"][1],conf["WE"][0],conf["WE"][1]],\
                              [conf["WE"][1],conf["WE"][0],conf["WE"][1],conf["WE"][0]]] ])


    conf["BETA"] = 1./conf["kBT"]

    conf["Vb_range"] = np.linspace(-conf["Vb"],conf["Vb"],conf["NVb"])

    conf["NV"] = int(conf["NV"])

    # In the following lines we specify the couplings.
    # They are defined according to the paper and everything is
    # multiplied by a coupling scale G_scale (corresponding to \gamma_0 in the paper)
    
    # tau_r is the parameter "a" of the paper
    conf["tau_s"] = ((1.+conf["tau_r"])/2.)*np.ones(4)*conf["G_scale"]
    conf["tau_d"] = ((1.-conf["tau_r"])/2.)*np.ones(4)*conf["G_scale"]

    # Pol is P
    s_up = (1. + conf["Pol"])/2.
    s_down = (1. - conf["Pol"])/2.

    # OrbPol is a possible valley dependent coupling or polarization
    # it is not used in the paper
    so_up = (1. + conf["OrbPol"])/2.
    so_down = (1. - conf["OrbPol"])/2.
    
    
    conf["DOS_P"] = np.array([s_up*so_up , s_down*so_up , s_up*so_down , s_down*so_down])

    conf["DOS_AP"] = np.array([s_down*so_up , s_up*so_up , s_down*so_down , s_up*so_down])

    # the whole coupling array is now structured as follows
    # gamma[parallel/anti-parallel][source/drain] is a (1,4) vector with the couplings
    conf["gamma"] = np.array([[conf["DOS_P"]*conf["tau_s"],conf["DOS_P"]*conf["tau_d"]],\
                              [conf["DOS_P"]*conf["tau_s"],conf["DOS_AP"]*conf["tau_d"]]])
    # this is the effective Zeeman fields
    conf["B"] = np.array([conf["B_P"],conf["B_AP"]])
    # ... that we can in principle also apply for orbital splitting
    conf["B_ORB"] = np.array([conf["B_ORB_P"],conf["B_ORB_AP"]])

    # generate the gate range
    # the boundary values are chose such that we "see" the main shell
    conf["V_g_min"] = (get_state_energy([0,0,0,0],conf['N_0'],1)-get_state_energy([0,0,-1,0],conf['N_0'],1))/conf["ALPHA"]
    conf["V_g_max"] = (get_state_energy([1,2,1,1],conf['N_0'],1)-get_state_energy([1,1,1,1],conf['N_0'],1))/conf["ALPHA"]
    conf["V_g"] = np.linspace(conf["V_g_min"],conf["V_g_max"],conf["NV"])

    # we pre-calculate all possible groundstate energies
    E = get_E()

    for config in np.arange(2):
        # here the ground states are constructed depending on the parameters in the
        # different configurations
        GROUNDSTATES[config] = get_groundstates(config)

        # and the number of groundstates is saved
        N_GS[config] = len([item for sublist in GROUNDSTATES[config] for item in sublist])

        # for the groundstate array we generate excited state arrays
        u,d,smap = generate_excitedstate_array(GROUNDSTATES[config],config)
        # the structure of upstates and downstates is
        # upstates[parallel/anti-parallel][number-of-shell][spin/valley (0-3)][num]
        # where num links to STATENUMBER by the statemap
        # the number-of-shell is given by twice the number
        # of excited states plus two, 2*NX + 2.
        # Note that this arrays contain all the necessary information for all
        # possible tunneling events for all states in the ground state space,
        # given the coupling array gamma and the bandwidth W_E which
        # do only depend on spin and valley numbers.
        upstates.append(u)
        downstates.append(d)
        statemap.append(smap)
        
    TARRAY_TRANSPORT = get_TARRAY_TRANSPORT()

    CHARGE = get_CHARGE()

    return 0

def get_selected_groundstates(configuration):
    """
    a convenient method to generate a "tailored" ground state
    space that can depend on the configuration for experimenting with it.

    Args:
      configuration: integer (0 for parallel / 1 for anti-parallel)
    Returns:
      A list of lists of state vectors, one list for each charging state
    """
    if configuration == 0:
        groundstates = [
            [[0,-1,0,0],[0,0,0,-1],[-1,0,0,0],[0,0,-1,0]],
            [[0,0,0,0]],
            [[0,1,0,0],[0,0,0,1]],
            [[0,1,0,1]],
            [[0,2,0,1],[0,1,0,2]],
            [[1,2,0,1],[0,2,1,1],[1,1,0,2],[0,1,1,2]],
            [[1,2,1,1],[1,1,1,2],[2,1,1,1],[1,1,2,1]]]
    else:
        groundstates = [
            [[0,-1,0,0],[0,0,0,-1],[-1,0,0,0],[0,0,-1,0]],
            [[0,0,0,0]],
            [[0,1,0,0],[0,0,0,1]],
            [[0,1,0,1]],
            [[0,2,0,1],[0,1,0,2]],
            [[1,2,0,1],[0,2,1,1],[1,1,0,2],[0,1,1,2]],
            [[1,2,1,1],[1,1,1,2],[2,1,1,1],[1,1,2,1]]]

    return groundstates

def get_groundstates(configuration):
    """
    The method to generate the basic set of groundstates.
    This replaces a more "natural" way of generating all possible
    groundstates automatically, which, by the size of the state space
    does not really make a difference.

    Args:
      configuration: the configuration of the leads (0/1)
    Returns:
      A list of lists of state vectors, one list for each charging state
    """
    groundstates = [
        [[0,-1,0,0],[0,0,0,-1],[-1,0,0,0],[0,0,-1,0]],
        [[0,0,0,0]],
        [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],
        [[1,1,0,0],[1,0,1,0],[1,0,0,1],[0,1,1,0],[0,1,0,1],[0,0,1,1]],
        [[1,1,1,0],[1,1,0,1],[1,0,1,1],[0,1,1,1]],
        [[1,1,1,1]],
        [[2,1,1,1],[1,1,2,1],[1,2,1,1],[1,1,1,2]]]

    # for charge,state_array in enumerate(groundstates):
    #     energies = np.array([get_state_energy(np.array(state),conf['N_0'],configuration) for state in state_array])
    #     min_energy = np.amin(energies)
    #     gs_idx = np.nonzero(energies <= min_energy+conf['GS_cutoff'])[0]
    #     groundstates[charge] = [state_array[idx] for idx in gs_idx]

         
    return groundstates

def generate_excitedstate_array(groundstates,configuration):
    """
    Given a groundstate list of lists, this method generates the up-
    and downstates arrays.

    Args:
      groundstates: A list of lists of state vectors specifying all possible
                    groundstates
      configuration: two possible configurations (0/1)

    Returns:
      A list of two 3-dimensional arrays (one per configuration) with energies.
      The energies are calculated with respect to the given groundstate.
      The arrays are structured as
      my_upstates[parallel/anti-parallel][number-of-shell][spin/valley (0-3)][num]
      where the number-of-shell is given by twice the number
      of excited states plus two, 2*NX + 2.
    """
    my_upstates = list()
    my_downstates = list()
    statemap = np.zeros((256))
    ns = 0

    # for each charging state
    for line in groundstates:
        #print line
        # and for each state in the charging state
        for state in line:
            # we compute an array of excited states, append it to the
            # return arrays
            u,d = get_excited_states(np.array(state),configuration)
            my_upstates.append(u)
            my_downstates.append(d)
            # and keep the index in the statemap to connect STATENUMBER and
            # index
            statemap[state_to_int(state)] = ns
            ns += 1

    return np.dstack((my_upstates)),np.dstack((my_downstates)),statemap
    
def get_TARRAY_TRANSPORT():
    """
    Generates all the possible transitions between the groundstates and saves them
    to TLIST. It returns the TARRAY mask that has non-zero entries (i,j) where a
    transition between a state with STATENUMBER i and STATENUMBER j is allowed.

    Returns:
      numpy.array of size (2,256,256), i.e., one array per configuration.
    """
    global TLIST
    myTARRAY = np.zeros((2,256,256))

    for config in np.arange(2):
        groundstates = GROUNDSTATES[config]

        #if config == 0:
        #    groundstates = get_selected_groundstates()
        #else:
        #    groundstates = get_groundstates(config)

        # we first treat the "extremal" cases, the empty and the full shell
        for state in groundstates[0]:
            for trans_state in groundstates[1]:
                myTARRAY[config,state_to_int(state),state_to_int(trans_state)] = 1.
        for state in groundstates[6]:
            for trans_state in groundstates[5]:
                myTARRAY[config,state_to_int(state),state_to_int(trans_state)] = 1.

        # and all others in between
        for charge in np.arange(1,6):
            for state in groundstates[charge]:
                for trans_state in groundstates[charge-1]:
                    if np.sum(np.absolute(np.array(state)-np.array(trans_state))) == 1:
                        myTARRAY[config,state_to_int(state),state_to_int(trans_state)] = 1.
                for trans_state in groundstates[charge+1]:
                    if np.sum(np.absolute(np.array(state)-np.array(trans_state))) == 1:
                        myTARRAY[config,state_to_int(state),state_to_int(trans_state)] = 1.

        T = np.transpose(np.nonzero(myTARRAY[config]))

        # we also write the information in the TLIST list of transitions
        # (that is ordered). It will determine the list of rates.
        for transition in T:
            # note that the transition list is structured such that the first
            # entry in each line is always the entry with the higher charging state (b)
            # and the second the one which is connected by out-tunneling (a)
            # i.e. TLIST[parallel/anti-parallel][num of transition][b/a]
            if get_CHARGE()[transition[0]] > get_CHARGE()[transition[1]]:
                TLIST[config].append(transition)

    TLIST = np.array(TLIST)

    return myTARRAY



def get_E():
    """
    Prepares an array for all possible groundstate energies ordered by STATENUMBER.

    Returns:
      numpy.array of shape (256,2), one column per configuration
    """
    global E
    myE = np.zeros((256,2))

    for num in range(256):
        n = np.array([int_to_state(num),])
        myE[num,0] = get_state_energy(n,conf['N_0'],0)

        n = np.array([int_to_state(num),])
        myE[num,1] = get_state_energy(n,conf['N_0'],1)

    return myE



def get_CHARGE():
    """
    Returns charging states or electron numbers for all possible groundstates ordered
    by STATENUMBER.
    """
    return np.array([np.sum(int_to_state(num)) for num in range(256) ])



def get_P(T):
    """
    Returns an array with P[b,a] = +1 if the transition a to b is an in-tunneling event
    P[b,a] = -1 if a to b is an out-tunneling event. Thereby is the ordering
    according to the TARRAY, so a non-zero entry (i,j) corresponds to a transition
    from state with STATENUMBER j to STATENUMBER i.

    Args:
      T: A (256,256) array, usually a TARRAY

    """
    myP = np.zeros((256,256))
    select = np.where(T)

    for num in np.arange(select[0].size):
        a = select[0][num]
        b = select[1][num]
        myP[b,a] = CHARGE[b]-CHARGE[a]
    return myP


def get_excited_states(state,config):
    """
    Given a state, generates an array of accessible excited state energies,
    both for in- and out-tunneling of a single electron. 

    The result contains the energies of reachable
    states from the given state.
    Thereby represents the line with linenumber NX always the shell where the
    groundstate resides.
    This means, an entry (2,1) for a state [0,0,1,0] contains the energy difference
    between the state [0,0,1,0] plus an extra electron in the higher shell with
    spin down in the first valley and the state [0,0,1,0] (NX=1).
    
    Args:
      state: a state vector (1,4)
      config: the configuration of the leads (p/ap = 0/1)

    Returns:
      two arrays of shape (2*NX+2,4). 
    """
    # If we recieve a state of the next shell,
    # or a state from the underlying shell n-1,
    # we shift the shell number
    N = conf['N_0']
    if np.sum(state) >= 4:
        state = state-np.ones(4)
        N = N+1
    elif np.sum(state) < 0:
        state = state+np.ones(4)
        N = N-1
                
    # shells that can be accessed. This is put by hand in the config.
    # It could also be calculated from W_E, the excited state bandwidth,
    # but this way more fine tuning is allowed.
    n_s = conf['NX']
    
    # We create the groundstate in the representation with the nearby
    # accessible shells.
    # This structure has the shape (2*NX+2).
    # If you try to really understand what is happening here, please
    # make sure you know what gs contains and why.
    gs = np.vstack((np.ones((n_s+1,4)),state,np.zeros((n_s,4)),))
    in_e = -np.ones((gs.shape[0],gs.shape[1]))
    out_e = -np.ones((gs.shape[0],gs.shape[1]))
    in_qn = -np.ones((gs.shape[0],gs.shape[1]))
    out_qn = -np.ones((gs.shape[0],gs.shape[1]))
    #print gs
    # now we construct excited states
    # lists will contain all states, energies and the tunneling quantum number
    out_transitions = list()
    in_transitions = list()

    # we have to decide whether the empty shell state corresponds to shell
    # n or shell n-1. Depending on this, we start from 0 or from 1 in the
    # gs structure.
    # Without excited states, NX=0 and the empty state has the representation
    # [[1,1,1,1],[0,0,0,0]]. To collect all accessible states we collect all states
    # with one electron less on the lower (out-tunneling) and one additional electron
    # on the upper (in-tunneling) shell.
    if np.sum(state) == 0:
        shell0 = 0
    else:
        shell0 = 1

    # Line by line, we replace 0 by 1 or the other way round and calculate the
    # (relative) energy of this state
    # The energy is saved in the returned array.
    for exc_shell in range(shell0,gs.shape[0]):
        #print "shell {}".format(exc_shell)
        for n in range(4):
            if gs[exc_shell,n] == 1:
                tmp = np.copy(gs)
                tmp[exc_shell,n] = 0
                out_transitions.append(tmp)
                out_e[exc_shell,n] = get_state_energy(tmp,N,config)
            else:
                tmp = np.copy(gs)
                tmp[exc_shell,n] = 1
                in_transitions.append(tmp)
                in_e[exc_shell,n] = get_state_energy(tmp,N,config)

    #print in_transitions
    #print out_transitions
    #print in_trans_de
    #print out_trans_de

    return in_e,out_e
                
                
def get_state_energy(state,N,config):
    """
    Is the Hamiltonian of the system. For a state vector relative to a shell N
    for a configuration config it calculates the energy of this state.
    It can handle simple (1,4) state vectors (as used for the groundstates)
    as well as more complicated
    representations as the ones used by get_excited_states.

    Args:
      state: (1,4) state vector or (2*NX+2,4) state representation
      N: the shell number
      config: configuration of the leads (p/ap = 0/1)

    Returns:
      energy, float
    """
    # we can call this the old fashioned way with a list of four elements
    # from this we generate a (2,4) array in the spirit of get_excited_states
    if isinstance(state,list) and len(state)==4:
        state = np.array([state,])
    if np.any(state == -1):
        N = N-1
        state += np.ones(4)
    elif np.any(state == 2):
        N = N+1
        state -= np.ones(4)

    # the number of filled shells
    n_0 = N - state.shape[0]/2
    # the charging energy is easy
    charging_energy = 0.5*conf['E_C']*(4*n_0 + np.sum(state))**2
    #print "Charging energy: {}".format(charging_energy)
    # we create a shell energy ladder
    shell_array = (np.vstack([np.ones((4))*(n_0+n+1)*conf['E_0'] for n in range(state.shape[0])]))
    shell_energy = 2.*n_0*(n_0+1)*conf['E_0'] + np.sum(state*shell_array)
    #print "Shell energy: {}".format(shell_energy)

    spinorbit = np.dot(n_0+np.sum(state,axis=0),SO)*conf['SO']
    magneticfield = np.dot(n_0+np.sum(state,axis=0),S)*conf['MU_B']*conf['B'][config]
    #if config == 0:
    #    print "state {} b {}".format(n_0+np.sum(state,axis=0),magneticfield)
    orbitalsplit = np.dot(n_0+np.sum(state,axis=0),KK)*conf['B_ORB'][config]

    return charging_energy+shell_energy+spinorbit+magneticfield+orbitalsplit

    
    
