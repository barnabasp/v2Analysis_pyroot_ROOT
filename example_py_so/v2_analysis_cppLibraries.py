import sys
import os
import math
import ROOT
import numpy as np

import time


ROOT.gInterpreter.ProcessLine('#include "particle_tree.h"')
ROOT.gSystem.Load('./lib/libcommon.dylib')



NCENT = 6
NPT = 18
centLims = np.array([0,10,20,30,40,60,100])
pTLims = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0])

phidist = [[ROOT.TH1D() for _ in range(NCENT)] for _ in range(NPT)]
calcEventPlaneResDist = [[ROOT.TH1D() for _ in range(NCENT)] for _ in range(NPT)]
v2values = [ROOT.TGraph() for _ in range(NCENT)]
v2errs = [ROOT.TGraphErrors() for _ in range(NCENT)]

def whichCent(centr):
    ibin = np.searchsorted(centLims, centr, side='right') - 1
    return ibin

def whichpT(pTval):
    ibin = 0
    if pTval > pTLims[NPT]:
        return NPT - 1
    if pTval < pTLims[0]:
        return 0
    while pTval > pTLims[ibin + 1]:
        ibin += 1
    return ibin


def main(argv):
    t0 = time.time()
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
    centdist = ROOT.TH1F("centdist", "Centrality distribution", 100, 0, 100)
    zvertexdist = ROOT.TH1F("zvertexdist", "z Vertex distribution", 80, -40, 40)
    reacpldist = ROOT.TH1F("reacplandist", "Reaction plane distribution", 100, -2.0, 2.0)
    v2val = ROOT.TH1F("v2val", "True v2 values", 1000, 0.0, 0.1)
    RxNP = ROOT.TGraph("RxNP.csv","%lg %lg"," \\t,;")

    for icent in range(NCENT):
        v2values[icent] = ROOT.TGraph(NPT)
        v2values[icent].SetName(f"CentGraph{icent}")
        v2errs[icent] = ROOT.TGraphErrors(NPT)
        v2errs[icent].SetName(f"CentErrsGraph{icent}")
        for ipT in range(NPT):
            phidist[ipT][icent] = ROOT.TH1D(f"phidist_cent{icent}_pT{ipT}", "PhiDistribution;phi;Entries", 200,
                                             -np.pi / 2.0, np.pi / 2.0)
            calcEventPlaneResDist[ipT][icent] = ROOT.TH1D(f"phiCalcRP_RPDist_cent{icent}_pT{ipT}",
                                                   "cos(2*(Measured #varphi - RP)) Distribution; cos(2*(Measured #varphi - RP)); Entries",
                                                   200, -np.pi / 2.0, np.pi / 2.0)
    p = ROOT.particle_tree(infilename)
    if(p.fChain):
        print("Tree initialized")
    else:
        print("Tree not found")
        sys.exit(1)
    
    histFile = ROOT.TFile.Open(infilename, "READ")
    tree = histFile.Get("particle_tree")
    # Determine how many events to run on
    Nevents = p.fChain.GetEntries()
    if 0 < Nmaxevt < Nevents:
        Nevents = Nmaxevt
    print("Will run on {} events (out of {}).".format(Nevents, p.fChain.GetEntries()))
    correction = []
    corrCentBin = 0
    y = []
    for i in range(RxNP.GetN()+1):
        x = RxNP.GetPointX(i)
        if(corrCentBin == whichCent(x)):
            y.append(RxNP.GetPointY(i))
        else:
            correction.append(np.average(y))
            y = []
        corrCentBin = whichCent(x)
    print(correction)
    # Loop through events in the given dataset
    for ievent in range(Nevents):
    #  tree.GetEntry(ievent)
      p.GetEntry(ievent)
      # Monitor progress through output
      if ievent > 0 and ievent % 1000 == 0:
          print(".", end="", flush=True)
      if ievent > 0 and ievent % 10000 == 0:
          print("Analyzing event #", ievent)

      zvertexdist.Fill(p.Zvertex)

      cent = p.Centrality
      centdist.Fill(cent)
      centBin = whichCent(cent)
 


      reactionPlane = p.ReactionPlane
      reacpldist.Fill(reactionPlane)

      sumofSin = np.array(NPT*[0])
      sumofCos = np.array(NPT*[0])

      # Loop through all particles of the given event
      for ipart in range(p.Ntracks):
          pT = math.sqrt(p.px[ipart] * p.px[ipart] + p.py[ipart] * p.py[ipart])
          pTdist.Fill(pT)
          pTBin = whichpT(pT)

          #calculate azimuthal angle for given particle (in lab)
          phi0 = math.atan2(p.py[ipart], p.px[ipart])
          #calculate azimuthal angle for given particle (in reaction plane)
          phi = phi0 - reactionPlane
          # fix to [-pi / 2, pi / 2]
          while phi > math.pi / 2:
              phi -= math.pi
          while phi < -math.pi / 2:
              phi += math.pi
          #fill to histogram
          phidist[pTBin][centBin].Fill(phi)
          #Manual way to calculate the Reaction plane resolution
          sumofSin[pTBin] += math.sin(2 * phi0)
          sumofCos[pTBin] += math.cos(2 * phi0)
      for ipT in range(NPT):
          #calculate eventplane value for all - manual way
          phiCalcRP = math.atan2(sumofCos[ipT] * 2, sumofSin[ipT])
          phiCalcRP_RP = 2 * (phiCalcRP - reactionPlane)
          calcEventPlaneRes = math.cos(phiCalcRP_RP)
          #fill in histogram
          calcEventPlaneResDist[ipT][centBin].Fill(calcEventPlaneRes)
    t1 = time.time()
    total = t1-t0
    print()
    #define fit function for n = 2
    fitFunc = ROOT.TF1("fitFunc","[0]+[1]*cos(2*x)", -math.pi/2,  math.pi/2)
    for ipT in range(NPT):
        for icent in range(NCENT): 
            print("cent:", icent, "\tpt =", ipT)
            #fit phi to get parameters for v2
            #phidist[ipT][icent].Fit(fitFunc,"n")
            fitPtr = ROOT.TFitResultPtr(phidist[ipT][icent].Fit(fitFunc, "S"))
            #get parameters, constant and coefficient
            a = fitPtr.Parameter(0)
            b = fitPtr.Parameter(1)
            #get errors
            aerr = fitPtr.ParError(0)
            berr = fitPtr.ParError(1)
            #get covariance between parameters
            abcov = fitPtr.CovMatrix(1, 0)

            v2_raw = b/(2.0*a)
            print("Raw v2 is:", v2_raw)
            #R
            evPlaneRes = np.sqrt(abs(calcEventPlaneResDist[ipT][icent].GetMean()))
            evPlaneRes = correction[icent]
            print("event plane resolution:", evPlaneRes)

            #error propagation
            v2error_raw = abs(b/(2.0*a)) * np.sqrt((aerr**2) / (a**2) + (berr**2) / (b**2) - 2 * abcov / (a*b))

            v2_true = v2_raw/evPlaneRes
            v2error_true = v2error_raw / evPlaneRes

            print("true v2:", v2_true, "+-", v2error_true)
            v2val.Fill(v2_true)
            v2values[icent].SetPoint(icent,pTLims[ipT],v2_true)
            v2errs[icent].SetPoint(ipT,pTLims[ipT]+0.05,v2_true)
            v2errs[icent].SetPointError(ipT,0,v2error_true)
    t2 = time.time()
    total_2 = t2-t1
    print(total)
    print(total_2)
    # Write all histograms to the output ROOT file
    f = ROOT.TFile(outfilename, "RECREATE")
    if not f.IsWritable():
        print("File", outfilename, "was not opened!")
    else:
        print("Analysis done, writing histos to", outfilename)
    f.cd()
    pTdist.Write()
    centdist.Write()
    zvertexdist.Write()
    reacpldist.Write()
    v2val.Write()
    for icent in range(NCENT):
        v2values[icent].Write()
        v2errs[icent].Write()
        for ipT in range(NPT):
            phidist[ipT][icent].Write()
            calcEventPlaneResDist[ipT][icent].Write()

    f.Write()
    f.Close()

    return 0

sys.exit(main(sys.argv))
