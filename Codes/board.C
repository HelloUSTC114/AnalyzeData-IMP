#define board_cxx
#include "board.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void board::Loop()
{
    //   In a ROOT session, you can do:
    //      root> .L board.C
    //      root> board t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    // by  b_branchname->GetEntry(ientry); //read only this branch

    auto c = new TCanvas("c", "c", 1);
    auto h = new TH1D("h", "h", 2e2, -1e2, 1e2);
    auto h2 = new TH1D("h2", "h2", 2e2, -4e2, 4e2);

    auto tg = new TGraph();
    int tgPointCounter = 0;

    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;

    double time1 = -1, time2 = -1;
    int oddCounter = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // std::cout << jentry << std::endl;
        // if (Cut(ientry) < 0) continue;
        h->Fill(TDCTime[FiredCh[0]] - TDCTime[FiredCh[1]]);

        if (TDCTime[32] != 0)
        {
            if (time1 < 0)
                time1 = TDCTime[32];
            else
            {
                double dev = ((time2 - time1) - 1e9 + 347527);
                h2->Fill(dev);
                time1 = time2;
                time2 = TDCTime[32];
                // std::cout << ((time2 - time1) - 1e9 + 347527) << std::endl;
                if (dev > -4e2 && dev < 4e2)
                {
                    tg->SetPoint(tgPointCounter, tgPointCounter, dev);
                    tgPointCounter++;
                }
                if (TMath::Abs(time2 - time1) > 1.1e9)
                {
                    std::cout << oddCounter++ << '\t' << time1 << '\t' << time2 << '\t' << time2 - time1 << '\t' << std::endl;
                }
            }
        }
    }
    h->Draw();
    c->SaveAs("test.pdf");
    h2->Draw();
    c->SaveAs("test2.pdf");
    tg->SetTitle(";Time/s;T_{meas}-T_{0}/ns");
    tg->Draw("AZL");
    c->SaveAs("test3.pdf");
}
