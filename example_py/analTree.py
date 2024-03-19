import sys
import os
import math
import ROOT
import numpy as np

import time

def main(argv):
#    if len(argv) < 3:
#        print("Usage: {} <input file name> <output file name> <max events=-1>".format(argv[0]))
#        return -1

#    infilename = argv[1]
#    outfilename = argv[2]
    infilename = "data.root"
    outfilename = "test.root"
    print("Writing to", outfilename)

    Nmaxevt = -1
    if len(argv) >= 4:
        Nmaxevt = int(argv[3])
    if Nmaxevt < 1:
        Nmaxevt = -1


    ROOT.gStyle.SetOptFit(1111)

    # Create a histogram to be filled by looping through all events
    pTdist = ROOT.TH1F("pTdist", "pT distribution", 100, 0, 2)

    histFile = ROOT.TFile.Open(infilename, "READ")
    tree = histFile.Get("particle_tree")
    # Determine how many events to run on
    Nevents = tree.GetEntries()
    if 0 < Nmaxevt < Nevents:
        Nevents = Nmaxevt
    print("Will run on {} events (out of {}).".format(Nevents, tree.GetEntries()))

    # Loop through events in the given dataset
    for ievent in range(Nevents):
      tree.GetEntry(ievent)
      # Monitor progress through output
      if ievent > 0 and ievent % 1000 == 0:
          print(".", end="", flush=True)
      if ievent > 0 and ievent % 10000 == 0:
          print("Analyzing event #", ievent)

      # Loop through all particles of the given event
      for ipart in range(tree.Ntracks):
          pT = math.sqrt(tree.px[ipart] * tree.px[ipart] + tree.py[ipart] * tree.py[ipart])
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
