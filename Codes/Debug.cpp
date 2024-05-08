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

#include "board.C"

std::string gDataFolder = "../";
auto c = new TCanvas("c", "c", 1);

// Each Board information
int gBoardNo[6];
const int gBoardCount = 6;
board *gBoard[6];

void GenerateBoardMap()
{
    for (int i = 0; i < 6; i++)
    {
        gBoardNo[i] = i;
        std::cout << Form("%s/Board%d-Aligned.root", gDataFolder.c_str(), i) << std::endl;
        gBoard[i] = new board(Form("%s/Board%d-Aligned.root", gDataFolder.c_str(), i));
    }
}

void DebugMatchEntries()
{
    auto file = new TFile(Form("%s/MatchEntries.root", gDataFolder.c_str()));
    auto tree = (TTree *)(file->Get("match"));

    for (int i = 0; i < 6; i++)
    {
        std::cout << gBoard[i]->fChain->GetEntries() << std::endl;
    }

    int board1 = 0, board2 = 2;
    char cbuf1[512], cbuf2[512];
    sprintf(cbuf1, "(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9>>h(100000,-1000,1000)", board1, board1, board1, board2, board2, board2);
    // sprintf(cbuf1, "(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9:matchTime[%d]/1e9>>h(100000,0, 50e3,100,-0.4e6,0.4e6)", board1, board1, board1, board2, board2, board2, board1);
    // sprintf(cbuf1, "(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9:Entry$>>h(100000,0, 50e3,100,-0.4e6,0.4e6)", board1, board1, board1, board2, board2, board2);
    // sprintf(cbuf2, "matchFlag[%d]&&matchFlag[%d]&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)>100", board1, board2, board1, board1, board1, board2, board2, board2);
    sprintf(cbuf2, "matchFlag[%d]&&matchFlag[%d]", board1, board2);
    auto rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
    std::cout << rtn << std::endl;

    sprintf(cbuf1, "(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9>>h(100000,-1000,1000)", board1, board1, board1, board2, board2, board2);
    sprintf(cbuf2, "matchFlag[%d]&&matchFlag[%d]&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000", board1, board2, board1, board1, board1, board2, board2, board2);
    rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
    std::cout << rtn << std::endl;

    c->SaveAs("Debug/DebugMatchEntries.jpg");
}

TFile *gTSFile;
TTree *gTSTree;
bool gTSInit = 0;
double ts[6];
uint8_t tsFlag[6];
auto fGaus = new TF1("gaus", "gaus", -100, 100);
void DebugTSInit()
{
    gTSFile = new TFile(Form("%s/TS.root", gDataFolder.c_str()));
    gTSTree = (TTree *)(gTSFile->Get("tsTree"));
    gTSTree->SetBranchAddress("ts", ts);
    gTSTree->SetBranchAddress("tsFlag", tsFlag);
    // gTSTree->Draw("(ts[4]-ts[1])/1e9:ts[1]/1e9", "tsFlag[1]&&tsFlag[4]");
    // (TGraph *)(gTSFile->Get("Graph;1"));

    auto tgTemp = new TGraph();
    int counter = 0;
    for (int idx = 0; idx < gTSTree->GetEntries(); idx++)
    {
        gTSTree->GetEntry(idx);
        if (!tsFlag[1] || !tsFlag[4])
            continue;
        tgTemp->SetPoint(counter++, ts[1] / 1e9, (ts[4] - ts[1]) / 1e9);
    }
    tgTemp->SetTitle("board 4 - board 1 dev;real time/s;#Delta_{t}/s");
    tgTemp->Draw("AZP");
    c->SaveAs(Form("Debug/board%d-%dDev.jpg", 1, 4));
    gTSInit = 1;
}

void DebugTS(int board = 0)
{
    if (!gTSInit)
        DebugTSInit();
    auto tg = new TGraph();
    int pointCount = 0;

    gTSTree->GetEntry(1000);
    double ts1 = ts[board];
    gTSTree->GetEntry(2000);
    double ts2 = ts[board];
    double interval = (ts2 - ts1) / (2000 - 1000);

    for (int i = 0; i < gTSTree->GetEntries(); i++)
    {
        double tsPre, tsNow;
        gTSTree->GetEntry(i - 1);
        tsPre = ts[board];
        gTSTree->GetEntry(i);
        tsNow = ts[board];
        tg->SetPoint(pointCount++, tsNow / interval, (tsNow - tsPre) / interval - 1);
    }
    tg->SetMarkerStyle(2);
    tg->SetTitle(Form("T0 TS Interval for board %d;t/s;tsNext-tsPre/s;", board));
    tg->GetXaxis()->SetRangeUser(2, 10000);
    tg->GetYaxis()->SetRangeUser(-1000e-9, 1000e-9);
    tg->Draw("AZPL");
    c->SaveAs(Form("Debug/board%d.jpg", board));
}

void DebugPos()
{
    auto file = new TFile(Form("%s/GetPos.root", gDataFolder.c_str()));
    auto tree = (TTree *)(file->Get("data"));

    char cbuf1[1024], cbuf2[1024];
    int board1 = 0, board2 = 2, board3 = 4;
    // int board1 = 1, board2 = 3, board3 = 5;
    Long64_t rtn;
    c->SetLogz(1);

    // X
    if (1)
    {
        sprintf(cbuf1, "calcPos2[%d]-calcPos2[%d]:calcPos2[%d]-calcPos2[%d]>>h2D(600,-150,150, 100,-150,150)", board1, board2, board2, board3);
        sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3);
        rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
        std::cout << rtn << std::endl;
        auto h2D = (TH2F *)(file->Get("h2D"));
        h2D->SetTitle(Form("Position measurement;pos%d-pos%d/mm;pos%d-pos%d/mm", board2, board3, board1, board2));
        c->SaveAs(Form("Debug/Pos-%d%d%d-2D1.jpg", board1, board2, board3));
    }

    // Spatial resolution, with Calibration
    gStyle->SetOptFit(111);
    sprintf(cbuf1, "calcPos2[%d]/2.0-calcPos2[%d]+calcPos2[%d]/2.0>>hCali(1000,-15,5)", board1, board2, board3);
    sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&TMath::Abs(calcPos2[%d]-calcPos2[%d])<90", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3, board1, board3);
    rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
    std::cout << rtn << std::endl;
    auto hCali = (TH1F *)(file->Get("hCali"));
    hCali->SetTitle(Form("Position measurement;#Delta_{x}/mm;Counts"));
    hCali->Fit(fGaus, "", "", -10, -5);
    fGaus->SetRange(-20, 0);
    fGaus->Draw("same");
    c->SaveAs(Form("Debug/Pos-%d%d%d-Cali1.jpg", board1, board2, board3));

    // Compare with No calibration
    sprintf(cbuf1, "calcPos[%d]/2.0-calcPos[%d]+calcPos[%d]/2.0>>hNoCali(1000,-15,5)", board1, board2, board3);
    sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&TMath::Abs(calcPos[%d]-calcPos[%d])<90", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3, board1, board3);
    rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
    std::cout << rtn << std::endl;
    auto hNoCali = (TH1F *)(file->Get("hNoCali"));
    hNoCali->SetTitle(Form("Position measurement;#Delta_{x}/mm;Counts"));
    hNoCali->Fit(fGaus, "", "", -10, -5);
    fGaus->SetRange(-20, 0);
    fGaus->Draw("same");
    c->SaveAs(Form("Debug/Pos-%d%d%d-NoCali1.jpg", board1, board2, board3));

    // Y
    board1 = 1, board2 = 3, board3 = 5;
    if (1)
    {
        sprintf(cbuf1, "calcPos2[%d]-calcPos2[%d]:calcPos2[%d]-calcPos2[%d]>>h2D(600,-150,150, 100,-150,150)", board1, board2, board2, board3);
        sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3);
        rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
        std::cout << rtn << std::endl;
        auto h2D = (TH2F *)(file->Get("h2D"));
        h2D->SetTitle(Form("Position measurement;pos%d-pos%d/mm;pos%d-pos%d/mm", board2, board3, board1, board2));
        c->SaveAs(Form("Debug/Pos-%d%d%d-2D1.jpg", board1, board2, board3));
    }

    fGaus = new TF1("gaus", "gaus", -100, 100);
    // Spatial resolution, with Calibration
    gStyle->SetOptFit(111);
    sprintf(cbuf1, "calcPos2[%d]/2.0-calcPos2[%d]+calcPos2[%d]/2.0>>hCali(1000,-16,16)", board1, board2, board3);
    sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&TMath::Abs(calcPos2[%d]-calcPos2[%d])<90", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3, board1, board3);
    rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
    std::cout << rtn << std::endl;
    hCali = (TH1F *)(file->Get("hCali"));
    hCali->SetTitle(Form("Position measurement;#Delta_{x}/mm;Counts"));
    hCali->Fit(fGaus, "", "", -2, 1);
    fGaus->SetRange(-8, 8);
    fGaus->Draw("same");
    c->SaveAs(Form("Debug/Pos-%d%d%d-Cali1.jpg", board1, board2, board3));

    // Compare with No calibration
    sprintf(cbuf1, "calcPos[%d]/2.0-calcPos[%d]+calcPos[%d]/2.0>>hNoCali(1000,-16,16)", board1, board2, board3);
    sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&TMath::Abs(calcPos[%d]-calcPos[%d])<90", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3, board1, board3);
    rtn = tree->Draw(cbuf1, cbuf2, "COLZ");
    std::cout << rtn << std::endl;
    hNoCali = (TH1F *)(file->Get("hNoCali"));
    hNoCali->SetTitle(Form("Position measurement;#Delta_{x}/mm;Counts"));
    hNoCali->Fit(fGaus, "", "", -2, 1);
    fGaus->SetRange(-8, 8);
    fGaus->Draw("same");
    c->SaveAs(Form("Debug/Pos-%d%d%d-NoCali1.jpg", board1, board2, board3));

    sprintf(cbuf1, "calcPos2[%d]/2.0-calcPos2[%d]+calcPos2[%d]/2.0:calcPos2[%d]-calcPos2[%d]>>h2(600,-150,150, 100,-10,10)", board1, board2, board3, board3, board1);
    sprintf(cbuf2, "TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&TMath::Abs((matchTime[%d]-lastSeg[%d])/interval[%d]*1e9-(matchTime[%d]-lastSeg[%d])/interval[%d]*1e9)<1000&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&decodeFlag[%d]>1&&TMath::Abs(calcPos2[%d]-calcPos2[%d])<150", board1, board1, board1, board2, board2, board2, board2, board2, board2, board3, board3, board3, board1, board2, board3, board1, board3);
}

void DebugPos2()
{
    auto file = new TFile(Form("%s/GetPos.root", gDataFolder.c_str()));
    auto tree = (TTree *)(file->Get("data"));

    char cbuf1[1024], cbuf2[1024];
    int board1 = 0, board2 = 2, board3 = 4;
    // int board1 = 1, board2 = 3, board3 = 5;
    Long64_t rtn;
    // c->SetLogz(1);
    auto fGaus = new TF1("gaus", "gaus", -100, 100);
    // Spatial resolution, with Calibration
    // gStyle->SetOptFit(111);
    gStyle->SetOptStat(0);

    auto fline = new TF1("pol1", "pol1", 0, 1000);

    // X
    // Check deviation vs x2, calcPos
    // tree->Draw("calcPos[2]-(calcPos[0]+calcPos[4])/2-4.632:calcPos[2]>>(12800,0,640,100,-10,10)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos[2]-(calcPos[0]+calcPos[4])/2.0-4.632)<10", "COLZ");
    // c->cd();
    // tree->Draw("(calcPos[0]+calcPos[4])/2+4.632-floor(calcPos[2]/10)*10:calcPos[2]>>h1(12800,0,640,100,-20,20)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos[2]-(calcPos[0]+calcPos[4])/2.0-4.632)<10", "COLZ");
    // auto h1 = (TH2F *)(gFile->Get("h1"));
    // h1->GetXaxis()->SetRangeUser(210, 220);
    // fline->SetParameters(-210, 1);
    // fline->Draw("same");
    // for (int strip = 0; strip < 64; strip++)
    // {
    //     sprintf(cbuf1, "calcPos[2]-(calcPos[0]+calcPos[4])/2-4.632:calcPos[2]>>(200,%d,%d,100,-10,10)", strip * 10, strip * 10 + 10);
    //     sprintf(cbuf2, "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>0&&decodeFlag[4]>0&&(calcPos[2]>%d&calcPos[2]<%d)&&TMath::Abs(calcPos[2]-(calcPos[0]+calcPos[4])/2.0-4.632)<10", strip * 10, strip * 10 + 10);
    //     tree->Draw(cbuf1, cbuf2, "COLZ");
    // }

    // Check deviation vs x2, calcPos2
    // tree->Draw("calcPos2[2]-(calcPos2[0]+calcPos2[4])/2-4.632:calcPos2[2]>>(12800,0,640,100,-10,10)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");
    // tree->Draw("(calcPos2[0]+calcPos2[4])/2+4.632-floor(calcPos2[2]/10)*10:calcPos2[2]>>h12(12800,0,640,100,-10,20)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");
    tree->Draw("(calcPos2[0]+calcPos2[4])/2+4.632-calcPos2[2]:calcPos2[2]>>h12(12800,0,480,100,-10,10)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");
    auto h12 = (TH2F *)(gFile->Get("h12"));
    h12->SetTitle(";x_{2}/mm;#DeltaX/mm");
    c->SaveAs("Debug/dx-x.pdf");
    h12->GetXaxis()->SetRangeUser(310, 320);
    c->SaveAs("Debug/dx-x2.pdf");

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    tree->Draw("(calcPos2[0]+calcPos2[4])/2+4.632-calcPos2[2]>>hd(1000,-10,10)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");
    auto hd = (TH1 *)(gFile->Get("hd"));
    hd->Fit(fGaus, "Q", "", -2, 2);
    fGaus->SetRange(-10, 10);
    fGaus->Draw("same");
    c->SaveAs("Debug/DeltaX.pdf");

    // h12->GetXaxis()->SetRangeUser(210, 220);
    // fline->SetParameters(-210, 1);
    // fline->Draw("same");
    // for (int strip = 0; strip < 64; strip++)
    // {
    //     sprintf(cbuf1, "calcPos2[2]-(calcPos2[0]+calcPos2[4])/2-4.632:calcPos2[2]>>(200,%d,%d,100,-10,10)", strip * 10, strip * 10 + 10);
    //     sprintf(cbuf2, "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>0&&decodeFlag[4]>0&&(calcPos2[2]>%d&calcPos2[2]<%d)&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", strip * 10, strip * 10 + 10);
    //     tree->Draw(cbuf1, cbuf2, "COLZ");
    // }

    tree->Draw("(calcPos2[1]+calcPos2[5])/2-calcPos2[3]:calcPos2[3]>>h12Y(12800,0,480,100,-10,10)", "!(decodeFlag[1]==1&&decodeFlag[5]==1)&&decodeFlag[1]>0&&decodeFlag[3]>1&&decodeFlag[5]>0&&TMath::Abs(calcPos2[3]-(calcPos2[1]+calcPos2[5])/2.0-4.632)<10", "COLZ");
    auto h12Y = (TH2F *)(gFile->Get("h12Y"));
    h12Y->SetTitle(";x_{2}/mm;#DeltaX/mm");
    c->SaveAs("Debug/dx-xY.pdf");
    h12Y->GetXaxis()->SetRangeUser(310, 320);
    c->SaveAs("Debug/dx-x2Y.pdf");

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    tree->Draw("(calcPos2[1]+calcPos2[5])/2-calcPos2[3]>>hdY(1000,-10,10)", "!(decodeFlag[1]==1&&decodeFlag[5]==1)&&decodeFlag[1]>0&&decodeFlag[3]>1&&decodeFlag[5]>0&&TMath::Abs(calcPos2[3]-(calcPos2[1]+calcPos2[5])/2.0-4.632)<10", "COLZ");
    auto hdY = (TH1 *)(gFile->Get("hdY"));
    hdY->Fit(fGaus, "Q", "", -2, 2);
    fGaus->SetRange(-10, 10);
    fGaus->Draw("same");
    c->SaveAs("Debug/DeltaXY.pdf");
}
// cor->Draw("(calcPos2[0]+corDev[0]+calcPos2[4]+corDev[4])/2+4.632-((calcPos2[2]+corDev[2])/10)*10:(calcPos2[2]+corDev[2])>>h12(12800,0,640,100,-20,20)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");
// cor->Draw("(calcPos2[0]+corDev[0]+calcPos2[4]+corDev[4])/2+4.632-((calcPos2[2]+corDev[2])):(calcPos2[2]+corDev[2])>>h12(12800,0,640,100,-20,20)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");
// cor->Draw("(calcPos2[0]+corDev[0]+calcPos2[4]+corDev[4])/2+4.632-((calcPos2[2]+corDev[2])):(calcPos2[2]+corDev[2])>>h12(100,310,320,100,-10,10)", "!(decodeFlag[0]==1&&decodeFlag[4]==1)&&decodeFlag[0]>0&&decodeFlag[2]>1&&decodeFlag[4]>0&&TMath::Abs(calcPos2[2]-(calcPos2[0]+calcPos2[4])/2.0-4.632)<10", "COLZ");

void Debug()
{
    // GenerateBoardMap();
    // DebugMatchEntries();
    // for (int board = 0; board < 6; board++)
    //     DebugTS(board);

    // DebugPos();
    // DebugPos2();
}