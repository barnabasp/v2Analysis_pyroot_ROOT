import sys
import os
import math
import ROOT
import numpy as np

import time


ROOT.gInterpreter.ProcessLine('#include "example_py_so/particle_tree.h"')
ROOT.gSystem.Load('./example_py_so/lib/libcommon.dylib')


def main(argv):
    if len(argv) < 3:
        print("Usage: {} <input file name> <output file name> <max events=-1>".format(argv[0]))
        return -1

    infilename = argv[1]
    outfilename = argv[2]
#    infilename = "data.root"
#    outfilename = "test.root"
    print("Writing to", outfilename)

    Nmaxevt = -1
    if len(argv) >= 4:
        Nmaxevt = int(argv[3])
    if Nmaxevt < 1:
        Nmaxevt = -1


    ROOT.gStyle.SetOptFit(1111)

    # Create a histogram to be filled by looping through all events
    pTdist = ROOT.TH1F("pTdist", "pT distribution", 100, 0, 2)

    p = ROOT.particle_tree(infilename)
    if(p.fChain):
        print("Tree initialized")
    else:
        print("Tree not found")
        sys.exit(1)
    
    # Determine how many events to run on
    Nevents = p.fChain.GetEntries()
    if 0 < Nmaxevt < Nevents:
        Nevents = Nmaxevt
    print("Will run on {} events (out of {}).".format(Nevents, p.fChain.GetEntries()))
    # Loop through events in the given dataset
    for ievent in range(Nevents):
    #  tree.GetEntry(ievent)
      p.GetEntry(ievent)
      # Monitor progress through output
      if ievent > 0 and ievent % 1000 == 0:
          print(".", end="", flush=True)
      if ievent > 0 and ievent % 10000 == 0:
          print("Analyzing event #", ievent)

      # Loop through all particles of the given event
      for ipart in range(p.Ntracks):
          pT = math.sqrt(p.px[ipart] * p.px[ipart] + p.py[ipart] * p.py[ipart])
          pTdist.Fill(pT)
    # Write all histograms to the output ROOT file
    f = ROOT.TFile(outfilename, "RECREATE")
    if not f.IsWritable():
        print("File", outfilename, "was not opened!")
    else:
        print("Analysis done, writing histos to", outfilename)
    f.cd()
    pTdist.Write()

    f.Write()
    f.Close()

    return 0

sys.exit(main(sys.argv))
