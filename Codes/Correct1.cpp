
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

std::string gDataFolder = "../Sum/";
auto c = new TCanvas("c", "c", 1);
auto fGaus = new TF1("gaus", "gaus", -100, 100);
TFile *file;
TTree *cor;
double gPitch = 11; // unit: mm

const int gBoardCount = 8;

//  X direction is defined as even boards, which means top layers, while Y direction is defined as odd boards, which means bottom layers
double Z_POS[4] = {0., 918.0, 1814., 2737.};
int boardX[4] = {6, 2, 4, 0};
int boardY[4] = {7, 3, 5, 1};
const double layerDistance = 25.0; // unit: mm
double zArrayX[4], zArrayY[4];

// Graphs used in loops
TGraph *gDevX = new TGraph();
TGraph *gDevY = new TGraph();
auto fpolX = new TF1("fpolX", "pol1", -2000, 2000);
auto fpolY = new TF1("fpolY", "pol1", -2000, 2000);

TGraph *tgXTemp = new TGraph();
TGraph *tgYTemp = new TGraph();

void Correction()
{
    auto file = new TFile(Form("%s/GetPos.root", gDataFolder.c_str()));
    auto tree = (TTree *)(file->Get("data"));

    auto fileW = new TFile(Form("%s/1stCorrection.root", gDataFolder.c_str()), "recreate");
    auto treeW = new TTree("cor", "corrrection pos");

    double calcPos2[gBoardCount];
    int firedStrip1[gBoardCount];
    double calcPos3[gBoardCount];
    double eta[gBoardCount];
    int decodeFlag[gBoardCount];
    int decodeType[gBoardCount];
    int evenOdd[gBoardCount]; // judge even odd for strip number, even is 1, odd is 0
    double tanX, tanY;        // X is 0, 2, 4, Y is 1, 3, 5, calculate first order tan theta
    double corDev[gBoardCount], corEta[gBoardCount];
    double weight1[gBoardCount], weight1a[gBoardCount], weight1b[gBoardCount];
    double weight2[gBoardCount], weight2a[gBoardCount], weight2b[gBoardCount];
    double boardHG[gBoardCount][32];
    double boardLG[gBoardCount][32];
    double boardTDC[gBoardCount][33];

    double tdcTime[gBoardCount];
    double lastSeg[gBoardCount];
    double nextSeg[gBoardCount];
    double interval[gBoardCount];

    tree->SetBranchAddress("calcPos2", calcPos2);
    tree->SetBranchAddress("firedStrip1", firedStrip1);
    tree->SetBranchAddress("decodeFlag", decodeFlag);
    tree->SetBranchAddress("decodeType", decodeType);
    tree->SetBranchAddress("weight1", weight1);
    tree->SetBranchAddress("weight1a", weight1a);
    tree->SetBranchAddress("weight1b", weight1b);
    tree->SetBranchAddress("weight2", weight2);
    tree->SetBranchAddress("weight2a", weight2a);
    tree->SetBranchAddress("weight2b", weight2b);
    tree->SetBranchAddress("boardHG", boardHG);
    tree->SetBranchAddress("boardLG", boardLG);
    tree->SetBranchAddress("boardTDC", boardTDC);
    tree->SetBranchAddress("tdcTime", tdcTime);
    tree->SetBranchAddress("lastSeg", lastSeg);
    tree->SetBranchAddress("nextSeg", nextSeg);
    tree->SetBranchAddress("interval", interval);

    treeW->Branch("calcPos2", calcPos2, Form("calcPos2[%d]/D", gBoardCount));
    treeW->Branch("firedStrip1", firedStrip1, Form("firedStrip1[%d]/I", gBoardCount));
    treeW->Branch("weight1", weight1, Form("weight1[%d]/D", gBoardCount));
    treeW->Branch("weight1a", weight1a, Form("weight1a[%d]/D", gBoardCount));
    treeW->Branch("weight1b", weight1b, Form("weight1b[%d]/D", gBoardCount));
    treeW->Branch("weight2", weight2, Form("weight2[%d]/D", gBoardCount));
    treeW->Branch("weight2a", weight2a, Form("weight2a[%d]/D", gBoardCount));
    treeW->Branch("weight2b", weight2b, Form("weight2b[%d]/D", gBoardCount));

    treeW->Branch("calcPos3", calcPos3, Form("calcPos3[%d]/D", gBoardCount));
    treeW->Branch("eta", eta, Form("eta[%d]/D", gBoardCount));
    treeW->Branch("decodeFlag", decodeFlag, Form("decodeFlag[%d]/I", gBoardCount));
    treeW->Branch("decodeType", decodeType, Form("decodeType[%d]/I", gBoardCount));
    treeW->Branch("corDev", corDev, Form("corDev[%d]/D", gBoardCount));
    treeW->Branch("corEta", corEta, Form("corEta[%d]/D", gBoardCount));
    treeW->Branch("evenOdd", evenOdd, Form("evenOdd[%d]/I", gBoardCount));
    treeW->Branch("tanX", &tanX, "tanX/D");
    treeW->Branch("tanY", &tanY, "tanY/D");

    treeW->Branch("boardHG", boardHG, Form("boardHG[%d][32]/D", gBoardCount));
    treeW->Branch("boardLG", boardLG, Form("boardLG[%d][32]/D", gBoardCount));
    treeW->Branch("boardTDC", boardTDC, Form("boardTDC[%d][33]/D", gBoardCount));

    treeW->Branch("tdcTime", tdcTime, Form("tdcTime[%d]/D", gBoardCount));
    treeW->Branch("lastSeg", lastSeg, Form("lastSeg[%d]/D", gBoardCount));
    treeW->Branch("nextSeg", nextSeg, Form("nextSeg[%d]/D", gBoardCount));
    treeW->Branch("interval", interval, Form("interval[%d]/D", gBoardCount));

    // Position derived from other 3 detectors
    double derivedPos[gBoardCount];
    double fitSigma[gBoardCount];
    treeW->Branch("derivedPos", derivedPos, Form("derivedPos[%d]/D", gBoardCount));
    treeW->Branch("fitSigma", fitSigma, Form("fitSigma[%d]/D", gBoardCount));

    // TVector3 for each Point
    TVector2 vx[4], vy[4];
    TVector3 v[4];
    TVector3 vin, vout;
    TVector3 vScat;
    treeW->Branch("v0", "TVector3", &v[0]);
    treeW->Branch("v1", "TVector3", &v[1]);
    treeW->Branch("v2", "TVector3", &v[2]);
    treeW->Branch("v3", "TVector3", &v[3]);
    treeW->Branch("vin", "TVector3", &vin);
    treeW->Branch("vout", "TVector3", &vout);
    treeW->Branch("vScat", "TVector3", &vScat);

    // Scattering angle
    double thetaX, thetaY, theta;
    treeW->Branch("thetaX", &thetaX, "thetaX/D");
    treeW->Branch("thetaY", &thetaY, "thetaY/D");
    treeW->Branch("theta", &theta, "theta/D");

    //  X direction is defined as even boards, which means top layers, while Y direction is defined as odd boards, which means bottom layers
    for (int layer = 0; layer < 4; layer++)
    {
        zArrayX[layer] = Z_POS[layer] + layerDistance / 2.0;
        zArrayY[layer] = Z_POS[layer] - layerDistance / 2.0;
    }

    // Calculate deviation for layer 1, layer 2
    auto hDevLayer1X = new TH1F("hDevLayer1X", "Deviation for Layer 1X", 1000, -20, 20);
    auto hDevLayer1Y = new TH1F("hDevLayer1Y", "Deviation for Layer 1Y", 1000, -20, 20);
    auto hDevLayer2X = new TH1F("hDevLayer2X", "Deviation for Layer 2X", 1000, -20, 20);
    auto hDevLayer2Y = new TH1F("hDevLayer2Y", "Deviation for Layer 2Y", 1000, -20, 20);

    auto hDevLayer1XCor = new TH1F("hDevLayer1XCor", "Deviation for Layer 1X", 1000, -20, 20);
    auto hDevLayer1YCor = new TH1F("hDevLayer1YCor", "Deviation for Layer 1Y", 1000, -20, 20);
    auto hDevLayer2XCor = new TH1F("hDevLayer2XCor", "Deviation for Layer 2X", 1000, -20, 20);
    auto hDevLayer2YCor = new TH1F("hDevLayer2YCor", "Deviation for Layer 2Y", 1000, -20, 20);

    for (int entry = 0; entry < tree->GetEntries(); entry++)
    {
        tree->GetEntry(entry);
        if (entry % 10000 == 0)
            std::cout << std::setprecision(10) << "Deviation " << entry << std::endl;
        // Calculate deviation for layer 1, layer 2
        double x0 = calcPos2[boardX[0]], x1 = calcPos2[boardX[1]], x2 = calcPos2[boardX[2]], x3 = calcPos2[boardX[3]];
        double y0 = calcPos2[boardY[0]], y1 = calcPos2[boardY[1]], y2 = calcPos2[boardY[2]], y3 = calcPos2[boardY[3]];
        double z0 = Z_POS[0], z1 = Z_POS[1], z2 = Z_POS[2], z3 = Z_POS[3];

        double dev1X = 0, dev1Y = 0;
        double calcPosLayerX1, calcPosLayerY1;
        calcPosLayerX1 = x0 + (x3 - x0) * (z1 - z0) / (z3 - z0);
        calcPosLayerY1 = y0 + (y3 - y0) * (z1 - z0) / (z3 - z0);
        dev1X = x1 - calcPosLayerX1;
        dev1Y = y1 - calcPosLayerY1;
        hDevLayer1X->Fill(dev1X);
        hDevLayer1Y->Fill(dev1Y);

        double dev2X = 0, dev2Y = 0;
        double calcPosLayerX2, calcPosLayerY2;
        calcPosLayerX2 = x0 + (x3 - x0) * (z2 - z0) / (z3 - z0);
        calcPosLayerY2 = y0 + (y3 - y0) * (z2 - z0) / (z3 - z0);
        dev2X = x2 - calcPosLayerX2;
        dev2Y = y2 - calcPosLayerY2;
        hDevLayer2X->Fill(dev2X);
        hDevLayer2Y->Fill(dev2Y);
    }

    auto fGaus = new TF1("gaus", "gaus", -20, 20);
    hDevLayer1X->Fit(fGaus, "Q", "", hDevLayer1X->GetMean() - 1 * hDevLayer1X->GetRMS(), hDevLayer1X->GetMean() + 1 * hDevLayer1X->GetRMS());
    double dev1XMean = fGaus->GetParameter(1);
    hDevLayer1Y->Fit(fGaus, "Q", "", hDevLayer1Y->GetMean() - 1 * hDevLayer1Y->GetRMS(), hDevLayer1Y->GetMean() + 1 * hDevLayer1Y->GetRMS());
    double dev1YMean = fGaus->GetParameter(1);
    hDevLayer2X->Fit(fGaus, "Q", "", hDevLayer2X->GetMean() - 1 * hDevLayer2X->GetRMS(), hDevLayer2X->GetMean() + 1 * hDevLayer2X->GetRMS());
    double dev2XMean = fGaus->GetParameter(1);
    hDevLayer2Y->Fit(fGaus, "Q", "", hDevLayer2Y->GetMean() - 1 * hDevLayer2Y->GetRMS(), hDevLayer2Y->GetMean() + 1 * hDevLayer2Y->GetRMS());
    double dev2YMean = fGaus->GetParameter(1);

    for (int entry = 0; entry < tree->GetEntries(); entry++)
    {
        tree->GetEntry(entry);
        if (entry % 10000 == 0)
            std::cout << std::setprecision(10) << "Correction " << entry << std::endl;

        // Calculate Angle for X, Y
        gDevX->Set(0);
        gDevY->Set(0);
        for (int i : {0, 3})
        {
            gDevX->SetPoint(i, Z_POS[i], calcPos2[boardX[i]]);
            gDevY->SetPoint(i, Z_POS[i], calcPos2[boardY[i]]);
        }
        gDevX->SetPoint(1, Z_POS[1], calcPos2[boardX[1]] - dev1XMean);
        gDevY->SetPoint(1, Z_POS[1], calcPos2[boardY[1]] - dev1YMean);
        gDevX->SetPoint(2, Z_POS[2], calcPos2[boardX[2]] - dev2XMean);
        gDevY->SetPoint(2, Z_POS[2], calcPos2[boardY[2]] - dev2YMean);

        gDevX->Fit(fpolX, "Q");
        gDevY->Fit(fpolY, "Q");

        tanX = fpolX->GetParameter(1);
        tanY = fpolY->GetParameter(1);

        for (int layer = 0; layer < gBoardCount; layer++)
        {
            eta[layer] = (weight2a[layer] + weight2b[layer]) / (weight1a[layer] + weight1b[layer] + weight2a[layer] + weight2b[layer]);
            if ((int)(calcPos2[layer] / 10) % 2 == 0)
            {
                evenOdd[layer] = 0;
                corEta[layer] = 10.0 * (eta[layer] - 0.5);
            }
            else
            {
                evenOdd[layer] = 1;
                corEta[layer] = 10.0 * (0.5 - eta[layer]);
            }
            if (layer % 2 == 0)
                corDev[layer] = corEta[layer] * tanX;
            else
                corDev[layer] = corEta[layer] * tanY;

            calcPos3[layer] = calcPos2[layer] * 1.1 + corDev[layer];
        }

        calcPos3[boardX[1]] -= dev1XMean;
        calcPos3[boardY[1]] -= dev1YMean;
        calcPos3[boardX[2]] -= dev2XMean;
        calcPos3[boardY[2]] -= dev2YMean;
        // Calculate deviation for layer 1, layer 2
        double x0 = calcPos3[boardX[0]], x1 = calcPos3[boardX[1]], x2 = calcPos3[boardX[2]], x3 = calcPos3[boardX[3]];
        double y0 = calcPos3[boardY[0]], y1 = calcPos3[boardY[1]], y2 = calcPos3[boardY[2]], y3 = calcPos3[boardY[3]];
        double z0 = Z_POS[0], z1 = Z_POS[1], z2 = Z_POS[2], z3 = Z_POS[3];
        double z0X = zArrayX[0], z1X = zArrayX[1], z2X = zArrayX[2], z3X = zArrayX[3];
        double z0Y = zArrayY[0], z1Y = zArrayY[1], z2Y = zArrayY[2], z3Y = zArrayY[3];

        double dev1X = 0, dev1Y = 0;
        double calcPosLayerX1, calcPosLayerY1;
        calcPosLayerX1 = x0 + (x3 - x0) * (z1 - z0) / (z3 - z0);
        calcPosLayerY1 = y0 + (y3 - y0) * (z1 - z0) / (z3 - z0);
        dev1X = x1 - calcPosLayerX1;
        dev1Y = y1 - calcPosLayerY1;
        hDevLayer1XCor->Fill(dev1X);
        hDevLayer1YCor->Fill(dev1Y);

        double dev2X = 0, dev2Y = 0;
        double calcPosLayerX2, calcPosLayerY2;
        calcPosLayerX2 = x0 + (x3 - x0) * (z2 - z0) / (z3 - z0);
        calcPosLayerY2 = y0 + (y3 - y0) * (z2 - z0) / (z3 - z0);
        dev2X = x2 - calcPosLayerX2;
        dev2Y = y2 - calcPosLayerY2;
        hDevLayer2XCor->Fill(dev2X);
        hDevLayer2YCor->Fill(dev2Y);

        // Calculate derived position
        double xArray[4] = {x0, x1, x2, x3};
        double yArray[4] = {y0, y1, y2, y3};
        double zArray[4] = {z0, z1, z2, z3};

        for (int layer = 0; layer < 4; layer++)
        {
            vx[layer].Set(zArrayX[layer], xArray[layer]);
            vy[layer].Set(zArrayY[layer], yArray[layer]);
        }
        for (int module = 0; module < 2; module++)
        {
            int layer1 = module * 2;
            int layer2 = module * 2 + 1;

            auto vecX = vx[layer2] - vx[layer1];
            auto vdX1 = vx[layer1] + (zArray[layer1] - zArrayX[layer1]) / (zArrayX[layer2] - zArrayX[layer1]) * vecX;
            auto vdX2 = vx[layer2] + (zArray[layer2] - zArrayX[layer2]) / (zArrayX[layer2] - zArrayX[layer1]) * vecX;

            auto vecY = vy[layer2] - vy[layer1];
            auto vdY1 = vy[layer1] + (zArray[layer1] - zArrayY[layer1]) / (zArrayY[layer2] - zArrayY[layer1]) * vecY;
            auto vdY2 = vy[layer2] + (zArray[layer2] - zArrayY[layer2]) / (zArrayY[layer2] - zArrayY[layer1]) * vecY;

            v[layer1].SetXYZ(vdX1.Y(), vdY1.Y(), zArray[layer1]);
            v[layer2].SetXYZ(vdX2.Y(), vdY2.Y(), zArray[layer2]);
        }
        vin = v[1] - v[0];
        vout = v[3] - v[2];
        vScat = vout.Unit() - vin.Unit();
        theta = vin.Angle(vout);
        thetaX = TMath::ATan(vin.X() / vin.Z()) - TMath::ATan(vout.X() / vout.Z());
        thetaY = TMath::ATan(vin.Y() / vin.Z()) - TMath::ATan(vout.Y() / vout.Z());

        for (int i = 0; i < 4; i++)
        {
            tgXTemp->SetPoint(tgXTemp->GetN(), zArray[i], xArray[i]);
            tgYTemp->SetPoint(tgYTemp->GetN(), zArray[i], yArray[i]);
        }
        tgXTemp->Fit(fpolX, "Q");
        tgYTemp->Fit(fpolY, "Q");

        for (int layer = 0; layer < 4; layer++)
        {
            tgXTemp->Set(0);
            tgYTemp->Set(0);
            // for (int i = 0; i < 4; i++)
            // {
            //     if (i == layer)
            //         continue;
            //     tgXTemp->SetPoint(tgXTemp->GetN(), zArray[i], xArray[i]);
            //     tgYTemp->SetPoint(tgYTemp->GetN(), zArray[i], yArray[i]);
            // }

            double xDerived = 0, yDerived = 0;
            // tgXTemp->Fit(fpolX, "Q");
            xDerived = fpolX->Eval(zArray[layer]);
            derivedPos[boardX[layer]] = xDerived;
            fitSigma[boardX[layer]] = fpolX->GetParError(0) + fpolX->GetParError(1) * zArray[layer];

            // tgYTemp->Fit(fpolY, "Q");
            yDerived = fpolY->Eval(zArray[layer]);
            derivedPos[boardY[layer]] = yDerived;
            fitSigma[boardY[layer]] = fpolY->GetParError(0) + fpolY->GetParError(1) * zArray[layer];
        }

        treeW->Fill();
    }

    fileW->cd();
    treeW->Write();
    hDevLayer1X->Write();
    hDevLayer1Y->Write();
    hDevLayer2X->Write();
    hDevLayer2Y->Write();

    hDevLayer1XCor->Write();
    hDevLayer1YCor->Write();
    hDevLayer2XCor->Write();
    hDevLayer2YCor->Write();

    fileW->Close();
}

double dev1XMean, dev1YMean, dev2XMean, dev2YMean;

void DebugCor()
{
    auto c2 = new TCanvas("c2", "c2", 720, 720);
    c2->cd();
    gStyle->SetOptStat(0);

    bool flagScattering = true;
    if (flagScattering)
    {
        // Scattering angle for X direction
        cor->Draw(Form("TMath::Abs( (calcPos3[0]-(calcPos3[4]))/923.0 - ((calcPos3[2])-calcPos3[6])/918.0 ):640-(calcPos3[3]/1.1):(calcPos3[2]/1.1+247*tanX)>>hX(256,0,640,256, 0, 640,100,0.0,0.1)"), "decodeFlag[0]>0&&decodeFlag[2]>0&&decodeFlag[4]>0&&decodeFlag[6]>0", "");
        auto hX = (TH3F *)(gFile->Get("hX"));
        hX->SetTitle("Scattering angle for X direction;X/mm;Y/mm;Scattering angle");
        auto hX_pxy = hX->Project3DProfile("xy");
        hX_pxy->SetTitle("Scattering angle for X direction;X/mm;Y/mm;Scattering angle");
        // hX_pxy->GetZaxis()->SetRangeUser(0, 0.04);
        // hX_pxy->GetXaxis()->SetRangeUser(0, 480);
        hX_pxy->GetXaxis()->SetRangeUser(0, 480);
        hX_pxy->GetYaxis()->SetRangeUser(0, 480);
        hX_pxy->Draw("COLZ");
        c2->SaveAs(Form("%s/hX.jpg", gDataFolder.c_str()));

        // Scattering angle for Y direction
        cor->Draw(Form("TMath::Abs((calcPos3[1]-(calcPos3[5]))/923.0-((calcPos3[3])-calcPos3[7])/918.0):640-(calcPos3[3]/1.1):(calcPos3[2]/1.1+247*tanY)>>hY(256,0,640,256, 0, 640,100,0.0,0.1)"), "decodeFlag[0]>0&&decodeFlag[2]>0&&decodeFlag[4]>0&&decodeFlag[6]>0", "");
        auto hY = (TH3F *)(gFile->Get("hY"));
        hY->SetTitle("Scattering angle for Y direction;X/mm;Y/mm;Scattering angle");
        auto hY_pxy = hY->Project3DProfile("xy");
        hY_pxy->SetTitle("Scattering angle for Y direction;X/mm;Y/mm;Scattering angle");
        // hY_pxy->GetZaxis()->SetRangeUser(0, 0.04);
        // hY_pxy->GetXaxis()->SetRangeUser(0, 480);
        hY_pxy->GetXaxis()->SetRangeUser(0, 480);
        hY_pxy->GetYaxis()->SetRangeUser(0, 480);
        hY_pxy->Draw("COLZ");
        c2->SaveAs(Form("%s/hY.jpg", gDataFolder.c_str()));

        hX_pxy->Rebin2D(2, 2);
        hY_pxy->Rebin2D(2, 2);
        double threshold = 0.005;
        // for (int binX = 0; binX < hY_pxy->GetNbinsX(); binX++)
        //     for (int binY = 0; binY < hY_pxy->GetNbinsY(); binY++)
        //     {
        //         if (hY_pxy->GetBinContent(binX, binY) < 0.01)
        //             hY_pxy->SetBinContent(binX, binY, 0);
        //         if (hX_pxy->GetBinContent(binX, binY) < 0.01)
        //             hX_pxy->SetBinContent(binX, binY, 0);
        //     }
        hX_pxy->GetZaxis()->SetRangeUser(threshold, 0.01);
        hY_pxy->GetZaxis()->SetRangeUser(threshold, 0.01);

        // hX_pxy->GetXaxis()->SetRangeUser(0, 480);
        hX_pxy->GetXaxis()->SetRangeUser(0, 480);
        hX_pxy->GetYaxis()->SetRangeUser(0, 480);
        // hY_pxy->GetXaxis()->SetRangeUser(0, 480);
        hY_pxy->GetXaxis()->SetRangeUser(0, 480);
        hY_pxy->GetYaxis()->SetRangeUser(0, 480);

        hX_pxy->Draw("COLZ");
        c2->SaveAs(Form("%s/hX2.jpg", gDataFolder.c_str()));
        hY_pxy->Draw("COLZ");
        c2->SaveAs(Form("%s/hY2.jpg", gDataFolder.c_str()));
    }

    bool flagDeviation = true;
    if (flagDeviation)
    {

        // Get position resolution
        // // cor->Draw("TMath::Abs(calcPos3[0]-derivedPos[0]):calcPos3[0]>>hDiffX(200,0,640,200,0,20)", "decodeFlag[0]>0", "COLZ");
        // cor->Draw("TMath::Abs(calcPos3[0]-derivedPos[0]):calcPos3[1]:calcPos3[0]>>hDiffX(128,0,640,128,0,640,100,-30,30)", "decodeFlag[0]>0", "");
        // auto hDiffX = (TH3F *)(gFile->Get("hDiffX"));
        // hDiffX->SetTitle("Position resolution for X direction;X/mm;Position resolution/mm");
        // auto hDiffX_pxy = hDiffX->Project3DProfile("xy");
        // std::cout << hDiffX << '\t' << hDiffX_pxy << std::endl;
        // hDiffX_pxy->SetTitle("Position resolution for X direction;X/mm;Position resolution/mm");
        // hDiffX_pxy->Draw("COLZ");
        // c2->SaveAs(Form("%s/hDiffX.jpg", gDataFolder.c_str()));

        gStyle->SetOptStat(1);
        gStyle->SetOptFit(11);

        for (int layer = 0; layer < 4; layer++)
        {
            char cDraw[1024];
            char cCut[1024];
            c->cd();

            int xboardNo = 2 * layer;
            sprintf(cDraw, "(calcPos3[%d]-derivedPos[%d])>>hDiffX(1000,-10,10)", xboardNo, xboardNo);
            sprintf(cCut, "decodeFlag[%d]>1&&fitSigma[%d]<10", xboardNo, xboardNo);
            cor->Draw(cDraw, cCut);
            auto hDiffX = (TH1F *)(gFile->Get("hDiffX"));
            // hDiffX->SetTitle(Form("Board %d;#Delta_{x}/mm;Counts", xboardNo));
            hDiffX->SetTitle(Form("Layer %d X,Board %d;#Delta_{x}/mm;Counts", layer, xboardNo));
            hDiffX->Fit(fGaus, "Q");
            hDiffX->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1. * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1. * fGaus->GetParameter(2));
            hDiffX->Draw();
            fGaus->SetRange(fGaus->GetParameter(1) - 5 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 5 * fGaus->GetParameter(2));
            fGaus->Draw("same");
            c->SaveAs(Form("%s/Layer%d-hDiffX.jpg", gDataFolder.c_str(), layer));

            int yboardNo = 2 * layer + 1;
            sprintf(cDraw, "(calcPos3[%d]-derivedPos[%d])>>hDiffY(1000,-10,10)", yboardNo, yboardNo);
            sprintf(cCut, "decodeFlag[%d]>1&&fitSigma[%d]<10", yboardNo, yboardNo);
            cor->Draw(cDraw, cCut);
            auto hDiffY = (TH1F *)(gFile->Get("hDiffY"));
            // hDiffY->SetTitle("Position resolution for Y direction;Position resolution/mm;Counts");
            // hDiffY->SetTitle(Form("Board %d;#Delta_{x}/mm;Counts", yboardNo));
            hDiffY->SetTitle(Form("Layer %d Y,Board %d;#Delta_{x}/mm;Counts", layer, yboardNo));
            hDiffY->Fit(fGaus, "Q");
            hDiffY->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1. * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1. * fGaus->GetParameter(2));
            hDiffY->Draw();
            fGaus->SetRange(fGaus->GetParameter(1) - 5 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 5 * fGaus->GetParameter(2));
            fGaus->Draw("same");
            c->SaveAs(Form("%s/Layer%d-hDiffY.jpg", gDataFolder.c_str(), layer));
        }
    }
}

void Correct1(std::string sDataFolder = "../Sum/")
{
    gDataFolder = sDataFolder;
    // TestCaclA0();
    // return;
    std::ifstream fin(Form("%s/1stCorrection.root", gDataFolder.c_str()));
    if (!fin.is_open())
        Correction();

    file = new TFile(Form("%s/1stCorrection.root", gDataFolder.c_str()));
    cor = (TTree *)(gFile->Get("cor"));

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

    DebugCor();
    return;
}