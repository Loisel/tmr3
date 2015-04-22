#!/usr/bin/python

"""
The interface module for the tmr3 code.

In this file you can specify parameter ranges. The parameter names are
taken from the cnt.conf config file.
A table with all possible combinations of the parameters is stored,
with a configuration per line, in a file that is either a default file
(table.dat) or specified by the user.
The module's main routine parses the table file line by line, starting
from 0 or from the line number in the lnum.dat file.

The module uses a system call to span multiple (4) make threads that calculate
the rates.
When all rates are calculated, the tmr.py module is called and collects
the rates to evaluate the conductance and the TMR.

It is also very instructive to have a look at the Makefile that is
called by this module.

Example:
  python lookup.py -c -t mytable.dat

"""

import numpy as np
from collections import OrderedDict
import shutil
import os
import sys
import subprocess
import getopt

# This is a module to calculate the cartesian product of the input
# parameter space
from cartesian import cartesian

# The ordered dictionary with the parameters
# A cartesian product is done to span the possible parameterspace
# All combinations are considered.
# Note that for x parameters 'a' and for y parameters 'b'
# the program does x*y calculations

pdict = OrderedDict([('G_scale',np.array([0.8e-4])),
                     ('E_0',np.array([1.4e-3])),
                     ('Pol',np.array([0.4])),
                     ('OrbPol',np.array([0])),
                     ('E_C',np.array([6.0e-3])),
                     ('SO',np.array([0.])),
                     ('tau_r',np.array([-0.7])),
                     ('B_P',np.array([-12e-5])),
                     ('B_AP',np.array([-16e-5])),
                     ('W_E',np.array([4.2e-3])),
                     ('W_0',np.array([1.e-1])),
                     ('B_ORB_P',np.array([0])),
                     ('B_ORB_AP',np.array([0])),
                     ('kBT',np.array([5e-5]))])

# path to the dat file where the results are saved and the config file
tmrhome = "./"
datdir = tmrhome+"dat/"
defaultconfig = tmrhome+"cnt.conf"

# the number of (local) cores where the calculations are performed
numberofcores = 4



def generateTable(fname):
    """
    Generates the table file with all possible combinations of the input
    parameters in pdict.

    Saves the file to fname.

    Args:
      fname: the filename to save the table to.
    """

    table = cartesian((pdict[k] for k in pdict))
    
    np.savetxt(fname,table)


def submitJob(myDict):
    """
    Submits the jobs with the parameters in a dictionary.

    Args:
      myDict: a dictionary in the form {'parameter':value}
    """

    print "Submitting job with config: {}".format(myDict)

    # we create a directory to save the rates in
    jobname = "lookup_"+'_'.join("{}={:.2e}".format(key,value) for (key,value) in myDict.items())
    jobhome = datdir+jobname+"/"

    if not os.path.exists(jobhome): os.makedirs(jobhome)

    print "Saving the rates to {}".format(jobhome)

    # copy the config to the folder
    shutil.copy(defaultconfig,jobhome)

    # add the configuration in the dictionary to the config
    with open(jobhome+defaultconfig,'a') as configfile:
        for (key,value) in myDict.items():
            configfile.write('{} = {}\n'.format(key,value))

    # call make with -j for multiple threads and argv=jobname
    subprocess.call(["make", "-j", "{}".format(numberofcores), "argv={}".format(jobname), "tmr"])

    # cleanup the jobhome
    shutil.rmtree(jobhome)
    

def main(argv=None):
    """
    The main routine generates a table file if there is none and parses
    it to successively submit jobs with all possible parameter combinations.
    Thereby the calculation starts with a specified line number
    in the tablefile to resume calculations that were interrupted.
    
    Example:
      python lookup.py -t table.dat -l 5

    Args:
      tablefilename (-t --table): the filename of the tablefile.
                                  if not given, table.dat is used.
      linenumber (-l --line): number of the line in the table file to start
                              the calculation. First line is 0.
    """
    

    # the file with the line number in the table file to start the calculations
    linenumberfile = "lnum.dat"
    tablefile = None
    cont = False

    # parse the command line options
    if argv is None:
        argv = sys.argv
    try:
        opts, args = getopt.getopt(argv[1:], "hct:", ["help","continue","table"])
    except getopt.error, msg:
        raise Usage(msg)
    for o,a in opts:
        if o in ('-h','--help'):
            usage()
            exit()
        elif o in ('-t','--table'):
            tablefile = a
        elif o in ('-c', '--continue'):
            cont = True
        else:
            raise Usage('Invalid argument.')

    # We generate the tablefile if there is none specified
    if not tablefile:
        tablefile = "table.dat"
        print "Generating Tablefile"
        generateTable(tablefile)

    # if the continue switch is used, we start by line number
    # that can be found in the lnum dat file
    # Otherwise we start from zero.
    if cont:
        try:
            with open(linenumberfile, 'r') as lnf:
                line = int(lnf.readline())
                print "Starting from linenumber {}".format(line)
        except:
            print "No file with linenumber found and no linenumber given. Can not continue."
            exit
    else:
        line = 0


    # open the tablefile and iterate through the lines
    with open(tablefile) as tfile:
        table = tfile.readlines()

        while(True):
            myDict = OrderedDict()

            # we assume that the line in the table file are ordered
            # in the same way as pdict is and copy the content
            # to a temporary dictionary
            for num,key in enumerate(pdict):
                #print table[line]
                myDict[key] = float(table[line].split()[num])
            # the job is submitted with the current line from the table
            submitJob(myDict)
            line += 1

            # when we are done, we save the linenumber to a lnum.dat
            with open(linenumberfile, 'w') as lnf:
                lnf.write(str(line))


        
    
def usage():
    print "\n\
This is a tool to generate and populate a lookup table\n\
and perform the calculation of the TMR accordingly.\n\
\n\
usage: lookup.py [-t|--table tablefile] [-c|--continue]\n\
\n\
-t|--table tablefile specifies a file containing a lookup table.\n\
-c|--continue continue calculations in the default or specified table\n\
              file with the linenumber found in lnum.dat."
        
if __name__ == "__main__":
    sys.exit(main())
