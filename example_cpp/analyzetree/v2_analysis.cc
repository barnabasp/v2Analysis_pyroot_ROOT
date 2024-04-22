#define particle_tree_cxx
#include "particle_tree.h"
#include <iostream>
#include <cmath>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TWebFile.h>
#include <vector>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <numeric>


#define NCENT 6
#define NPT 18

const double centLims[NCENT+1] = {0,10,20,30,40,60,100};
const double pTLims[NPT+1] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};

TH1D *phidist[NPT][NCENT];
TH1D *phirawphidist[NPT][NCENT];
TGraph *v2values[NCENT];
TGraphErrors *v2errs[NCENT];
using namespace std;

int whichCent(int centr)
{
//if the interval is not 0 - 100
//    if( centr > centLims[NCENT] ) return NCENT-1;
//    if( centr < centLims[0] ) return 0;

    int ibin=0;
    while( centr > centLims[ibin+1] ) ibin++;
    return ibin;
}
int whichpT(double pTval)
{

    int ibin=0;
    if( pTval > pTLims[NPT] ) return NPT-1;
    if( pTval < pTLims[0] ) return 0;
    while( pTval > pTLims[ibin+1] ) ibin++;
    return ibin;
}

// -----------------------------------------------------------------------


int main(int argc, const char** argv)
{
    if(argc<2)
    {
        cerr << "Usage: " << argv[0] << " <input file name> <output file name> <max events=-1>" << endl;
        return -1;
    }
    string infilename(argv[1]);
    string outfilename(argv[2]);
    cerr << "Writing to " << outfilename.c_str() << endl;
    int Nmaxevt = -1;
    if(argc>=4) Nmaxevt = atoi(argv[3]);
    if(Nmaxevt<1) Nmaxevt = -1;

    //CREATE A HISTOGRAM, AND FILL IT BY LOOPING THROUGH ALL EVENTS
    TH1 *pTdist = new TH1F("pTdist","pT distribution",100,0,2);
    TH1 *centdist = new TH1F("centdist", "Centrality distribution", 100, 0, 100);
    TH1 *zvertexdist = new TH1F("zvertexdist", "z Vertex distribution", 80, -40, 40);
    TH1 *reacpldist = new TH1F("reacplandist", "Reaction plane distribution", 100, -2.0, 2.0);
    TH1 *v2val = new TH1F("v2val", "True v2 values", 1000, 0.0, 0.1);
    TGraph *RxNP = new TGraph("RxNP.csv","%lg %lg"," \\t,;");

    for(int icent = 0; icent < NCENT; icent++)
    {
        v2values[icent] = new TGraph(NPT);
        v2values[icent]->SetName(Form("CentGraph%i",icent));
        v2errs[icent] = new TGraphErrors(NPT);
        v2errs[icent]->SetName(Form("CentErrsGraph%i",icent));
        for(int ipT = 0; ipT < NPT; ipT++)
        {
            phidist[ipT][icent] = new TH1D(Form("phidist_cent%i_pT%i",icent, ipT), "PhiDistribution;phi;Entries", 200, -M_PI/2.0,M_PI/2.0);
            phirawphidist[ipT][icent] = new TH1D(Form("phiCalcRP_RPDist_cent%i_pT%i",icent, ipT), "cos(2*(Measured #varphi - RP)) Distribution; cos(2*(Measured #varphi - RP)); Entries", 200, -M_PI/2.0,M_PI/2.0);
        }
    }
    std::vector<double> correction;
    int corrCentBin = 0;
    std::vector<double> y;

    for(int i = 0; i < RxNP->GetN(); i++)
    {
      int x = RxNP->GetPointX(i);
      if(corrCentBin == whichCent(x))
      {
        y.push_back(RxNP->GetPointY(i));
      }
      else
      {
        double sum = std::accumulate(y.begin(), y.end(), 0.0); 

        // average of the vector elements 
        double avg = sum / y.size(); 
        correction.push_back(avg);
        y.clear();
      }
      corrCentBin = whichCent(x);
    }
    

    particle_tree p(infilename.c_str());
    if(p.fChain) cerr << "Tree initialized" << endl;
    else { cerr << "No tree found." << endl; return 1; }
    long unsigned int Nevents = p.fChain->GetEntries();
    if(Nmaxevt>0 && Nmaxevt<(int)Nevents) Nevents=Nmaxevt;
    cerr << "Will run on " << Nevents << " events (out of " << p.fChain->GetEntries()  << ")." << endl;

    for(long unsigned int ievent=0;ievent<Nevents;ievent++)
    {
        if(ievent>0&&ievent%1000==0) cerr << ".";
        if(ievent>0&&ievent%10000==0) cerr << "Analyzing event #" << ievent << endl;
        p.GetEntry(ievent);

        zvertexdist->Fill(p.Zvertex);

        //centrality determination
        int cent = p.Centrality;
        centdist->Fill(cent);
        int centBin = whichCent(cent);



        double reactionPlane = p.ReactionPlane;
        reacpldist->Fill(reactionPlane);
        //variables for reactionplanePhi
        double sumofSin[NPT] = {0};
        double sumofCos[NPT] = {0};

        //LOOP THROUGH ALL PARTICLES OF THE GIVEN EVENT
        for(int ipart=0;ipart<p.Ntracks;ipart++)
        {
            double pT = sqrt(p.px[ipart]*p.px[ipart]+p.py[ipart]*p.py[ipart]);
            int pTBin = whichpT(pT);
            pTdist->Fill(pT);

            //calculating raw phi from py/px 
            double rawPhi = std::atan2(p.py[ipart],p.px[ipart]);
            double phi = rawPhi-reactionPlane;
            while(phi > M_PI/2) phi -= M_PI;
            while(phi < -M_PI/2) phi += M_PI;
//            cerr << "phi = " <<phi << "for pt, cent" << pTBin << " " << centBin <<endl;
            phidist[pTBin][centBin]->Fill(phi);


            sumofSin[pTBin] += sin(2*rawPhi);
            sumofCos[pTBin] += cos(2*rawPhi);
        }
        for(int ipT = 0; ipT < NPT; ipT++)
        { 
            double realPhi = atan2(sumofCos[ipT]*2,sumofSin[ipT]);
            //while(realPhi > M_PI/2) realPhi -= M_PI;
            //while(realPhi < -M_PI/2) realPhi += M_PI;
            double phirpDiff = 2*(realPhi-reactionPlane);
//            while(phirpDiff > M_PI/2) phirpDiff -= M_PI;
//            while(phirpDiff < -M_PI/2) phirpDiff += M_PI;
            double phiRawPhi = cos(phirpDiff);

            phirawphidist[ipT][centBin]->Fill(phiRawPhi);
        }
    }
    cerr << endl;

    TF1* fitFunc = new TF1("fitFunc","[0]+[1]*cos(2*x)", -TMath::Pi()/2,  TMath::Pi()/2);
    for(int ipT = 0; ipT < NPT; ipT++)
    {
        for(int icent = 0; icent < NCENT; icent++)
        { 
            cerr << "cent: " <<icent << "\tpt = "<<ipT<< endl;
            TFitResultPtr r = phidist[ipT][icent]->Fit(fitFunc, "S");

            // constant term
            double A = r->Parameter(0);
            // coefficient for 2 * cos(2 * x)
            double B = r->Parameter(1);
            // get errors
            double AErr = r->ParError(0);
            double BErr = r->ParError(1);
            // get covariance
            double ABCov = r->CovMatrix(1, 0);

            // elliptic flow (v2) and correction via reaction plane resolution

            double v2_raw = B/(2.0*A);
            cout << "Raw v2 is: " << v2_raw << endl;
            //sigm
            double evPlaneRes = sqrt(abs(phirawphidist[ipT][icent]->GetMean()));
            evPlaneRes = correction[icent];
            cout << "event plane resolution: "<< evPlaneRes << endl;

            double v2_true = v2_raw/evPlaneRes;



            // error estimation through error propagation
            double v2ErrNonCorr = std::abs(B / A) * std::sqrt((AErr * AErr) / (A * A) + (BErr * BErr) / (B * B) - 2 * ABCov / (A * B));
            double v2Err = std::abs(B / A) * std::sqrt((AErr * AErr) / (A * A) + (BErr * BErr) / (B * B) - 2 * ABCov / (A * B)) / evPlaneRes;
            double v2_truerr = v2_true*v2Err;
            cout << "true v2: " << v2_true << "+-" << v2_truerr<< endl;

            /*
            phidist[ipT][icent]->Fit(fitFunc,"n");
            double a = fitFunc->GetParameter(0);
            double b = fitFunc->GetParameter(1);
            double aerr = fitFunc->GetParError(0);
            double berr = fitFunc->GetParError(1);
            double error = (aerr/a + berr/b);
            */


            v2val->Fill(v2_true);
            v2values[icent]->SetPoint(icent,pTLims[ipT],v2_true);
            v2errs[icent]->SetPoint(ipT,pTLims[ipT]+0.05,v2_true);
            v2errs[icent]->SetPointError(ipT,0,v2_truerr);
            
        }
    }

    

    //WRITE ALL HISTOGRAMS TO THE OUTPUT ROOT FILE
    TFile *f = new TFile(outfilename.c_str(),"RECREATE");
    if(!f->IsWritable()) cerr << "File " << outfilename.c_str() << " was not opened!" << endl;
    else cerr << "Analysis done, writing histos to " << outfilename.c_str() << endl;
    f->cd();


    pTdist->Write();
    centdist->Write();
    zvertexdist->Write();
    reacpldist->Write();
    v2val->Write();
    for(int icent = 0; icent < NCENT; icent++)
    {
        v2values[icent]->Write();
        v2errs[icent]->Write();
        for(int ipT = 0; ipT < NPT; ipT++)
        {
            phidist[ipT][icent]->Write();
            phirawphidist[ipT][icent]->Write();
        }
    }
    f->Write();
    f->Close();

  //BELOW ARE SOME EXAMPLES TO DRAW STUFF WITHOUT THE NEED TO LOOP THROUGH EVENTS
  //p.fChain->Draw("ReactionPlane","ReactionPlane>-9000");
  //c1->Print("figs/ReactionPlane.png");
  //p.fChain->Draw("Centrality","Centrality>-9000");
  //c1->Print("figs/Centrality.png");
  //p.fChain->Draw("Zvertex","Zvertex>-9000");
  //c1->Print("figs/Zvertex.png");
  //p.fChain->Draw("Ntracks");
  //c1->Print("figs/Ntracks.png");
  return 0; 
}
// -----------------------------------------------------------------------
