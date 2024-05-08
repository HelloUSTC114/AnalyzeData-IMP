#include "cor.C"

#include <iostream>
#include <iomanip>
#include <Math/MinimizerOptions.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TGraph.h>
#include <TProfile.h>

std::string gDataFolder = "./";
auto c2 = new TCanvas("c", "c", 1);

// Fit function with 2 linear functions
double FitFun(double *x, double *par)
{
    // par[0]:t, par[1]:b, par[2]:k1, par[3]:k2, k1,k2>0
    if (x[0] <= par[0])
        return par[1] + par[2] * (x[0] - par[0]);
    if (x[0] > par[0])
        return par[1] - par[2] * (x[0] - par[0]);
    // return par[1] - par[3] * (x[0] - par[0]);
    return 0;
}
// auto fitfun = new TF1("fitfun", FitFun, -2000, 2000, 4);
auto fitfun = new TF1("fitfun", FitFun, -2000, 2000, 3);

void Align()
{
    cor cor1("1stCorrection.root");

    // Set Fit Function
    TF1 *fitfunArray[64]{0};
    for (int strip = 0; strip < 64; strip++)
    {
        auto fitfun = new TF1("fitfun", FitFun, -2000, 2000, 3);
        fitfunArray[strip] = fitfun;
        fitfun->SetParNames("t", "b", "k1", "k2");
        fitfun->SetParLimits(2, 0, 100);
        fitfun->SetParLimits(3, 0, 100);
    }

    // Open File
    auto tree = cor1.fChain;
    gSystem->Exec("mkdir Alignment");
    auto fileW = new TFile(Form("%s/Alignment.root", gDataFolder.c_str()), "recreate");
    tree->Draw("calcPos3[0]-derivedPos[0]>>hDiff(200,-20,20)", "decodeFlag[0]>0", "COLZ");
    auto hDiff = (TH1F *)gDirectory->Get("hDiff");
    auto fGaus = new TF1("fGaus", "gaus", -20, 20);
    hDiff->Fit(fGaus, "Q");
    auto hDiffMean = fGaus->GetParameter(1);
    std::cout << "Mean: " << hDiffMean << std::endl;

    // Do Alignment
    TH2F *hArray[64]{0};
    TGraph *tgArray[64]{0};
    TH2F *hAll = new TH2F("hAll", "All", 100, 0, 640, 100, 0, 10);
    gStyle->SetOptFit(11);
    for (int strip = 0; strip < 64; strip++)
    {
        hArray[strip] = new TH2F(Form("h%d", strip), Form("Strip %d;x/mm;Q_{rel}", strip), 100, 11 * (strip - 1) - 10, 11 * (strip) + 30, 100, 0, 5);
        tgArray[strip] = new TGraph();
        tgArray[strip]->SetName(Form("tg%d", strip));
        tgArray[strip]->SetTitle(Form("Graph %d;x/mm;Q_{rel}", strip));

        // char cDraw[1024];
        // char cCut[1024];
        // sprintf(cDraw, "(firedStrip1[0]==%d)?(weight1a[0]+weight1b[0]):((firedStrip1[0]==%d&&decodeFlag[0]>1)?(weight2a[0]+weight2b[0]):0):derivedPos[0]>>h%d", strip, strip - 1, strip);
        // sprintf(cCut, "decodeFlag[0]>0&&(firedStrip1[0]==%d||(firedStrip1[0]==%d&&decodeFlag[0]>1))", strip, strip - 1);
        // auto entries = tree->Draw(cDraw, cCut, "COLZ");
        // std::cout << "Strip: " << strip << " Entries: " << entries << std::endl;
    }

    for (int entry = 0; entry < cor1.fChain->GetEntries(); entry++)
    {
        cor1.GetEntry(entry);

        int strip1 = cor1.firedStrip1[0];
        int strip2 = strip1 + 1;

        double pos = cor1.derivedPos[0];
        if (isnan(pos) || isnan(cor1.weight1a[0]) || isnan(cor1.weight1b[0]) || isnan(cor1.weight2a[0]) || isnan(cor1.weight2b[0]))
            continue;
        // if (pos < strip1 * 11 - 10 || pos > strip1 * 11 + 20)
        if (pos < -10 || pos > 800)
            continue;
        if (TMath::Abs(cor1.calcPos3[0] - hDiffMean - pos) > 5)
            continue;
        hAll->Fill(pos, cor1.weight1a[0] + cor1.weight1b[0] + cor1.weight2a[0] + cor1.weight2b[0]);

        if (strip1 < 0 || strip1 > 63)
            continue;
        if (cor1.decodeFlag[0] < 1)
            continue;
        tgArray[strip1]->SetPoint(tgArray[strip1]->GetN(), pos, cor1.weight1a[0] + cor1.weight1b[0]);
        hArray[strip1]->Fill(pos, cor1.weight1a[0] + cor1.weight1b[0]);

        if (strip2 < 0 || strip2 > 63)
            continue;
        if (cor1.decodeFlag[0] < 2)
            continue;
        tgArray[strip2]->SetPoint(tgArray[strip2]->GetN(), pos, cor1.weight2a[0] + cor1.weight2b[0]);
        hArray[strip2]->Fill(pos, cor1.weight2a[0] + cor1.weight2b[0]);
    }

    double tPre = 0;
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer()
    for (int strip = 0; strip < 64; strip++)
    {
        // fitfun->SetParameters(strip * 11, 2, 2.0 / 10, 2.0 / 10);
        fitfunArray[strip]->SetParameters(strip * 11, 2, 2.0 / 10);

        // tgArray[strip]->GetXaxis()->SetRangeUser(11 * (strip)-15 - hDiffMean, 11 * (strip) + 15 - hDiffMean);
        // tgArray[strip]->GetYaxis()->SetRangeUser(0, 3);
        // tgArray[strip]->Fit(fitfun, "QL", "", strip * 11 - 11 - hDiffMean, strip * 11 + 11 - hDiffMean);
        // tgArray[strip]->Draw("AZP");
        // c2->SaveAs(Form("Alignment/Strip-%d.jpg", strip));

        hArray[strip]->Draw("COLZ");
        hArray[strip]->GetXaxis()->SetRangeUser(11 * (strip)-15 - hDiffMean, 11 * (strip) + 15 - hDiffMean);
        hArray[strip]->GetYaxis()->SetRangeUser(0, 3);
        auto hpf = hArray[strip]->ProfileX();
        hpf->GetXaxis()->SetRangeUser(11 * (strip)-15 - hDiffMean, 11 * (strip) + 15 - hDiffMean);
        hpf->GetYaxis()->SetRangeUser(0, 3);
        hpf->Fit(fitfunArray[strip], "QL", "", strip * 11 - 11 - hDiffMean, strip * 11 + 11 - hDiffMean);
        c2->SaveAs(Form("Alignment/pStrip-%d.jpg", strip));

        // hArray[strip]->FitSlicesY();
        // auto h1 = (TH1D *)gDirectory->Get(Form("h%d_1", strip));
        // h1->GetXaxis()->SetRangeUser(11 * (strip)-15 - hDiffMean, 11 * (strip) + 15 - hDiffMean);
        // h1->GetYaxis()->SetRangeUser(0, 3);
        // h1->Fit(fitfunArray[strip], "QL", "", strip * 11 - 11 - hDiffMean, strip * 11 + 11 - hDiffMean);
        // c2->SaveAs(Form("Alignment/h1Strip-%d.jpg", strip));

        tgArray[strip]->GetXaxis()->SetRangeUser(11 * (strip)-15 - hDiffMean, 11 * (strip) + 15 - hDiffMean);
        tgArray[strip]->GetYaxis()->SetRangeUser(0, 3);
        tgArray[strip]->Draw("AZP");
        fitfunArray[strip]->Draw("same");
        c2->SaveAs(Form("Alignment/Strip-%d.jpg", strip));

        double tNow = fitfunArray[strip]->GetParameter(0);
        std::cout << std::setprecision(10) << "Strip: " << strip << " t: " << tNow << " t-tPre: " << tNow - tPre << std::endl;
        tPre = tNow;
    }

    // First iteration
    TGraph *tgArray2[64]{0};
    for (int strip = 0; strip < 64; strip++)
    {
        tgArray2[strip] = new TGraph();
        tgArray2[strip]->SetName(Form("tg2%d", strip));
        tgArray2[strip]->SetTitle(Form("Graph 2 %d;x/mm;Q_{rel}", strip));

        for (int i = 0; i < tgArray[strip]->GetN(); i++)
        {
            double x, y;
            tgArray[strip]->GetPoint(i, x, y);
            double y2 = y - fitfunArray[strip]->Eval(x);
            tgArray2[strip]->SetPoint(tgArray2[strip]->GetN(), x, y2);
        }

        tgArray2[strip]->GetXaxis()->SetRangeUser(11 * (strip)-15 - hDiffMean, 11 * (strip) + 15 - hDiffMean);
        tgArray2[strip]->GetYaxis()->SetRangeUser(-3, 3);
        tgArray2[strip]->Draw("AZP");
        c2->SaveAs(Form("Alignment/tg2Strip-%d.jpg", strip));
    }

    for (int strip = 0; strip < 64; strip++)
    {
        hArray[strip]->Write();
        tgArray[strip]->Write();
        if (hArray[strip] == NULL)
            continue;

        // hArray[strip]->Draw("COLZ same");
    }
    hAll->Draw("COLZ");
    hAll->Write();
    c2->SaveAs("Alignment/All.jpg");
    delete fileW;
}