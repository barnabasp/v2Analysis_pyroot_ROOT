using namespace std;
#define NCENT 6
#define NPT 18

void Plot_analyzetree(const char* filename="/Users/barnabasp/Documents/Code/Nehezion_spec/analyzetree/test_5M.root", const char* figdir="figs")
{
    TFile *f = new TFile(filename);

    TCanvas *c1 = new TCanvas();

    TH1 *ptdist = (TH1F*)f->Get("pTdist");
    TH1 *zvertexdist = (TH1F*)f->Get("zvertexdist");
    TH1 *centdist = (TH1F*)f->Get("centdist");
    TH1 *reacpldist = (TH1F*)f->Get("reacplandist");
    TH1 *v2val = (TH1D*)f->Get("v2val");
    //TGraph *v2 = (TGraph*)f->Get("Graph");


    ptdist->SetTitle("N(p_{T})");
    ptdist->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ptdist->GetYaxis()->SetTitle("N(p_{T}) [c/GeV]");
    ptdist->Draw("e");
    c1->Print(Form("%s/ptdist.pdf",figdir));
    c1->Clear();

    centdist->SetTitle("Centrality distribution");
    centdist->GetXaxis()->SetTitle("Centrality");
    centdist->GetYaxis()->SetTitle("");
    centdist->Draw("e");
    c1->Print(Form("%s/centdist.pdf",figdir));
    c1->Clear();

    zvertexdist->SetTitle("z Vertex distribution");
    zvertexdist->GetXaxis()->SetTitle("ZVertex");
    zvertexdist->GetYaxis()->SetTitle("");
    zvertexdist->Draw("e");
    c1->Print(Form("%s/zvertexdist.pdf",figdir));
    c1->Clear();

    reacpldist->SetTitle("Reaction Plane distribution");
    reacpldist->GetXaxis()->SetTitle("Reaction Plane");
    reacpldist->GetYaxis()->SetTitle("");
    reacpldist->Draw("e");
    c1->Print(Form("%s/reacpldist.pdf",figdir));
    c1->Clear();
    
    v2val->Draw("e");
    c1->Print(Form("%s/v2val.pdf",figdir));
    c1->Clear();

    TF1* fitFunc = new TF1("fitFunc","[0]+[1]*cos(2*x)", -TMath::Pi()/2,  TMath::Pi()/2);
    for(int ipT = 0; ipT < NPT; ipT++)
    {
        for(int icent = 0; icent < NCENT; icent++)
        {
            TH1 *phidist = (TH1D*)f->Get(Form("phidist_cent%i_pT%i",icent,ipT));
            TH1 *phiRawphidist = (TH1D*)f->Get(Form("phiCalcRP_RPDist_cent%i_pT%i",icent,ipT));
            cerr << icent << ipT <<endl;
            double yMax = phidist->GetMaximum();
            phidist->GetYaxis()->SetRangeUser(yMax-12e3,yMax+5e3);
            phidist->Fit(fitFunc,"n");
            double p0 = fitFunc->GetParameter(0);
            double p1 = fitFunc->GetParameter(1);
            phidist->Draw("e");
            fitFunc->Draw("same");
//            c1->Print(Form("%s/phidist_cent%i_pT%i.pdf",figdir,icent,ipT));
            c1->Clear();
            phiRawphidist->Draw("e");
//            c1->Print(Form("%s/phirawphidist_cent%i_pT%i.pdf",figdir,icent,ipT));
            c1->Clear();
        }
    }
    
    TLegend* lptcent = new TLegend(0.15,0.55,0.30,0.75);

    const double centLims[NCENT+1] = {0,10,20,30,40,60,100};
    int color[NCENT] = {2,3,4,5,6,7};
    const char* labels[NCENT] = {"0 - 10 %","10 - 20 %","20 - 30 %","30 - 40%","40 - 60 %","60 - 100 %"};
    
    for(int icent = 0; icent < NCENT; icent++)
    {
        TGraphErrors *v2err = (TGraphErrors*)f->Get(Form("CentErrsGraph%i",icent));
        v2err->SetMarkerStyle(4);
        v2err->SetMarkerColor(color[icent]);
        v2err->SetMarkerSize(1.1);
        if(icent == 0)
        {
            v2err->GetXaxis()->SetTitle("pT (GeV/c)");
            v2err->GetXaxis()->SetRangeUser(0.2,4.7);
            v2err->GetYaxis()->SetTitle("v_{2}");
            //v2err->GetYaxis()->SetRangeUser(0.0,0.25);
            v2err->SetTitle("v_{2}(p_{T}) values for different centrality bins ");
            v2err->Draw("AP");
        }
        else v2err->Draw("P");
        lptcent->AddEntry(v2err,labels[icent],"P");
        if(icent == 5) 
        {
            lptcent->SetBorderSize(0);
            lptcent->Draw("same");
            c1->Print(Form("%s/v2err.pdf",figdir));
            c1->Clear();
        }
    }
/*
    v2->Draw("a2p");
    v2->SetMarkerStyle(20);
    v2->SetMarkerSize(1);
    c1->Print(Form("%s/v2.pdf",figdir));
    c1->Clear();
*/
    

}
