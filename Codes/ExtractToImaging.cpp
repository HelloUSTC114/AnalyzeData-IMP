#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <string>
#include <iostream>
#include <TProfile.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TMath.h>
#include <iomanip>
#include <fstream>
#include <TH3.h>
#include <TProfile2D.h>

std::string gDataFolder = "./";
int boardX[4] = {6, 2, 4, 0};
int boardY[4] = {7, 3, 5, 1};

void ExtractToImaging1()
{
    auto fileIn = new TFile(Form("%s/GetPos.root", gDataFolder.c_str()));
    auto treeIn = (TTree *)(fileIn->Get("data"));

    double calcPos2[8];
    int decodeFlag[8];
    double T[8];
    treeIn->SetBranchAddress("calcPos2", calcPos2);
    treeIn->SetBranchAddress("decodeFlag", decodeFlag);

    auto fileOut = new TFile(Form("%s/Imaging1.root", gDataFolder.c_str()), "RECREATE");
    auto treeOut = new TTree("CRY", "CRY data");
    treeOut->Branch("X0", &calcPos2[boardX[0]], "X0/D");
    treeOut->Branch("Y0", &calcPos2[boardY[0]], "Y0/D");
    treeOut->Branch("X1", &calcPos2[boardX[1]], "X1/D");
    treeOut->Branch("Y1", &calcPos2[boardY[1]], "Y1/D");
    treeOut->Branch("X2", &calcPos2[boardX[2]], "X2/D");
    treeOut->Branch("Y2", &calcPos2[boardY[2]], "Y2/D");
    treeOut->Branch("X3", &calcPos2[boardX[3]], "X3/D");
    treeOut->Branch("Y3", &calcPos2[boardY[3]], "Y3/D");
    for (int i = 0; i < 4; i++)
        treeOut->Branch(Form("T%d", i), &T[i], Form("T%d/D", i));
    double ke;
    treeOut->Branch("KE", &ke, "KE/D");

    for (int entry = 0; entry < treeIn->GetEntries(); entry++)
    {
        treeIn->GetEntry(entry);
        bool flag = 1;
        for (int i = 0; i < 8; i++)
        {
            // std::cout << decodeFlag[i] << " ";
            if (decodeFlag[i] < 0)
                flag = 0;
        }
        // std::cout << std::endl;
        if (!flag)
            continue;

        for (int i = 0; i < 4; i++)
            T[i] = 1;
        ke = 1; // GeV

        treeOut->Fill();
    }
    fileOut->cd();
    treeOut->Write();
    fileOut->Close();

    fileIn->Close();
}

void ExtractToImaging2()
{
    auto fileIn = new TFile(Form("%s/1stCorrection.root", gDataFolder.c_str()));
    auto treeIn = (TTree *)(fileIn->Get("cor"));
    auto fGaus = new TF1("fGaus", "gaus", -100, 100);

    double calcPos3[8];
    int decodeFlag[8];
    double T[8];
    treeIn->SetBranchAddress("calcPos3", calcPos3);
    treeIn->SetBranchAddress("decodeFlag", decodeFlag);

    // Read the mean of the deviation
    double dev1XMean, dev1YMean, dev2XMean, dev2YMean;
    auto hDevLayer1X = (TH1F *)(gFile->Get("hDevLayer1X"));
    auto hDevLayer1Y = (TH1F *)(gFile->Get("hDevLayer1Y"));
    auto hDevLayer2X = (TH1F *)(gFile->Get("hDevLayer2X"));
    auto hDevLayer2Y = (TH1F *)(gFile->Get("hDevLayer2Y"));
    hDevLayer1X->Fit(fGaus, "Q");
    hDevLayer1X->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1.82 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1.82 * fGaus->GetParameter(2));
    dev1XMean = fGaus->GetParameter(1);
    hDevLayer1Y->Fit(fGaus, "Q");
    hDevLayer1Y->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1.82 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1.82 * fGaus->GetParameter(2));
    dev1YMean = fGaus->GetParameter(1);
    hDevLayer2X->Fit(fGaus, "Q");
    hDevLayer2X->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1.82 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1.82 * fGaus->GetParameter(2));
    dev2XMean = fGaus->GetParameter(1);
    hDevLayer2Y->Fit(fGaus, "Q");
    hDevLayer2Y->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1.82 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1.82 * fGaus->GetParameter(2));
    dev2YMean = fGaus->GetParameter(1);

    auto fileOut = new TFile(Form("%s/Imaging2.root", gDataFolder.c_str()), "RECREATE");
    auto treeOut = new TTree("CRY", "CRY data");
    treeOut->Branch("X0", &calcPos3[boardX[0]], "X0/D");
    treeOut->Branch("Y0", &calcPos3[boardY[0]], "Y0/D");
    treeOut->Branch("X1", &calcPos3[boardX[1]], "X1/D");
    treeOut->Branch("Y1", &calcPos3[boardY[1]], "Y1/D");
    treeOut->Branch("X2", &calcPos3[boardX[2]], "X2/D");
    treeOut->Branch("Y2", &calcPos3[boardY[2]], "Y2/D");
    treeOut->Branch("X3", &calcPos3[boardX[3]], "X3/D");
    treeOut->Branch("Y3", &calcPos3[boardY[3]], "Y3/D");
    for (int i = 0; i < 4; i++)
        treeOut->Branch(Form("T%d", i), &T[i], Form("T%d/D", i));
    double ke;
    treeOut->Branch("KE", &ke, "KE/D");

    for (int entry = 0; entry < treeIn->GetEntries(); entry++)
    {
        treeIn->GetEntry(entry);
        bool flag = 1;
        for (int i = 0; i < 8; i++)
        {
            // std::cout << decodeFlag[i] << " ";
            if (decodeFlag[i] < 0)
                flag = 0;
        }
        // calcPos3[2] -= dev1XMean;
        // calcPos3[3] -= dev1YMean;
        // calcPos3[6] -= dev2XMean;
        // calcPos3[7] -= dev2YMean;

        // std::cout << std::endl;
        if (!flag)
            continue;

        for (int i = 0; i < 4; i++)
            T[i] = 1;
        ke = 1; // GeV

        treeOut->Fill();
    }
    fileOut->cd();
    treeOut->Write();
    fileOut->Close();

    fileIn->Close();
}

void ExtractToImaging(std::string sDataFolder = "../Sum/")
{
    gDataFolder = sDataFolder;
    // ExtractToImaging1();
    ExtractToImaging2();
}