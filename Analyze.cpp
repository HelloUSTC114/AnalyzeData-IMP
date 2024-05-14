#include <TCanvas.h>
#include "Codes/cor.C"
#include "Codes/Decode.cpp"

auto c = new TCanvas("c", "c", 1);
int gXNormalStrip = 24;
int gYNormalStrip = 40;
const int gBoardCount = 8;
const int gStripCount = 64;
double gDeltaTa[gBoardCount][gStripCount];
double gDeltaTb[gBoardCount][gStripCount];

const double layerDistance = 25.0; // unit: mm
double Z_POS[4] = {0., 918.0, 1814., 2737.};
int boardX[4] = {6, 2, 4, 0};
int boardY[4] = {7, 3, 5, 1};

TFile *f = new TFile("Sum/1stCorrection.root");
TTree *cor = (TTree *)f->Get("cor");
// typedef john::cor classcor;
auto corfile_ptr = new john::cor(cor);
auto &corfile = *corfile_ptr;
int gTotalEntries = corfile.fChain->GetEntries();
// int gTotalEntries = 300000;

double FitStripCenter(double *x, double *par)
{
    double center = par[0];
    double amp = par[1];

    const double gStripHalfLength = 11;
    double xmin = center - gStripHalfLength;
    double xmax = center + gStripHalfLength;

    double distance = TMath::Abs(x[0] - center);

    return (1.0 - distance / gStripHalfLength) * amp;
}

// 5. Draw time difference for each board
void TDCCorrectionForEachBoard()
{
    // Loop for all entries in cor
    auto fileW = new TFile("Sum/tdcTimeDiff.root", "recreate");
    TH3F *htaXY[gBoardCount / 2];
    TH3F *htbXY[gBoardCount / 2];
    TH3F *htabXY[gBoardCount];

    TProfile2D *htaXY_pyx[gBoardCount / 2];
    TProfile2D *htbXY_pyx[gBoardCount / 2];
    TProfile2D *htabXY_pyx[gBoardCount];

    for (int i = 0; i < gBoardCount / 2; i++)
    {
        htaXY[i] = new TH3F(Form("htaXY%d", i), Form("TDC Time Difference for Group A Board %d;x/mm;y/mm;#Delta_t/ns", i), 48, 0, 48, 48, 16, 64, 5000, -50, 50);
        htbXY[i] = new TH3F(Form("htbXY%d", i), Form("TDC Time Difference for Group B Board %d;x/mm;y/mm;#Delta_t/ns", i), 48, 0, 48, 48, 16, 64, 5000, -50, 50);
    }
    for (int i = 0; i < gBoardCount; i++)
    {
        htabXY[i] = new TH3F(Form("htabXY%d", i), Form("TDC Time Difference for Group A Board %d;x/mm;y/mm;#Delta_t/ns", i), 48, 0, 48, 48, 16, 64, 5000, -50, 50);
    }

    // for (int entry = 0; entry < corfile.fChain->GetEntries(); entry++)
    // int gTotalEntries = 30000;
    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile.GetEntry(entry);
        // Do something

        if (corfile.decodeFlag[0] == -3 || corfile.decodeFlag[1] == -3)
            continue;

        int strip[gBoardCount];
        double ta[gBoardCount], tb[gBoardCount];

        // Find valid strip, ta, tb for Board 0
        for (int board = 0; board < gBoardCount; board++)
        {
            int strip0 = corfile.firedStrip1[board];
            int firedSiPM = corfile.processedCh[board][0];
            int strip1 = strip0 + 1;
            int a0 = UserDefine::GetFiberCode(strip0, UserDefine::a);
            int b0 = UserDefine::GetFiberCode(strip0, UserDefine::b);
            int a1 = UserDefine::GetFiberCode(strip1, UserDefine::a);
            int b1 = UserDefine::GetFiberCode(strip1, UserDefine::b);

            int validStrip;
            int validCha, validChb;
            bool flag;
            if (a0 == firedSiPM || b0 == firedSiPM)
            {
                validStrip = strip0;
                validCha = a0;
                validChb = b0;
            }
            else if (a1 == firedSiPM || b1 == firedSiPM)
            {
                validStrip = strip1;
                validCha = a1;
                validChb = b1;
            }
            else
            {
                flag = false;
                continue;
            }
            strip[board] = validStrip;
            ta[board] = corfile.processedTDC[board][validCha];
            tb[board] = corfile.processedTDC[board][validChb];
        }

        for (int group = 0; group < gBoardCount / 2; group++)
        {
            int x0strip = strip[group * 2];
            int x1strip = strip[group * 2 + 1];
            htaXY[group]->Fill(x0strip, x1strip, ta[group * 2] - ta[group * 2 + 1]);
            htbXY[group]->Fill(x0strip, x1strip, tb[group * 2] - tb[group * 2 + 1]);

            htabXY[group * 2]->Fill(x0strip, x1strip, ta[group * 2] - tb[group * 2]);
            htabXY[group * 2 + 1]->Fill(x0strip, x1strip, ta[group * 2 + 1] - tb[group * 2 + 1]);
        }
    }

    gStyle->SetOptStat(0);
    fileW->cd();
    for (int group = 0; group < gBoardCount / 2; group++)
    {
        htaXY[group]->Draw("COLZ");
        auto h01a = htaXY[group]->Project3DProfile("yx");
        htaXY_pyx[group] = h01a;
        h01a->SetTitle(Form("TDC Time Difference for A fibers in Group %d;x strip;y strip;#Delta_t/ns", group));
        h01a->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-x-y-a-%d.jpg", group));
        htaXY[group]->Write();
        h01a->Write();

        htbXY[group]->Draw("COLZ");
        auto h01b = htbXY[group]->Project3DProfile("yx");
        htbXY_pyx[group] = h01b;
        h01b->SetTitle(Form("TDC Time Difference for B fibers in Group %d;x strip;y strip;#Delta_t/ns", group));
        h01b->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-x-y-b-%d.jpg", group));
        htbXY[group]->Write();
        h01b->Write();
    }
    for (int board = 0; board < gBoardCount; board++)
    {
        htabXY[board]->Draw("COLZ");
        auto htab = htabXY[board]->Project3DProfile("yx");
        htabXY_pyx[board] = htab;
        htab->SetTitle(Form("TDC Time Difference for Group ab in Board %d;x strip;y strip;#Delta_t/ns", board));
        htab->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-ab-%d.jpg", board));
        htabXY[board]->Write();
        htab->Write();
    }

    // Normalize delta t
    for (int group = 0; group < gBoardCount / 2; group++)
    {

        for (int strip = 0; strip < 48; strip++)
        {
            gDeltaTa[group * 2][strip] = htaXY_pyx[group]->GetBinContent(htaXY_pyx[group]->FindBin(strip, gYNormalStrip)) - htaXY_pyx[group]->GetBinContent(htaXY_pyx[group]->FindBin(gXNormalStrip, gYNormalStrip));
            gDeltaTb[group * 2][strip] = htbXY_pyx[group]->GetBinContent(htbXY_pyx[group]->FindBin(strip, gYNormalStrip)) - htbXY_pyx[group]->GetBinContent(htbXY_pyx[group]->FindBin(gXNormalStrip, gYNormalStrip));
        }

        for (int strip = 16; strip < 64; strip++)
        {
            gDeltaTa[group * 2 + 1][strip] = -htaXY_pyx[group]->GetBinContent(htaXY_pyx[group]->FindBin(gXNormalStrip, strip)) + htaXY_pyx[group]->GetBinContent(htaXY_pyx[group]->FindBin(gXNormalStrip, gYNormalStrip));
            gDeltaTb[group * 2 + 1][strip] = -htbXY_pyx[group]->GetBinContent(htbXY_pyx[group]->FindBin(gXNormalStrip, strip)) + htbXY_pyx[group]->GetBinContent(htbXY_pyx[group]->FindBin(gXNormalStrip, gYNormalStrip));
        }
    }

    // Delta t correction
    TH3F *htaXYcor[gBoardCount / 2];
    TH3F *htbXYcor[gBoardCount / 2];
    TH3F *htabXYcor[gBoardCount];

    TProfile2D *htaXYcor_pyx[gBoardCount / 2];
    TProfile2D *htbXYcor_pyx[gBoardCount / 2];
    TProfile2D *htabXYcor_pyx[gBoardCount];

    for (int i = 0; i < gBoardCount / 2; i++)
    {
        htaXYcor[i] = new TH3F(Form("htaXYcor%d", i), Form("TDC Time Difference for Group A Board %d;x/mm;y/mm;#Delta_t/ns", i), 48, 0, 48, 48, 16, 64, 5000, -50, 50);
        htbXYcor[i] = new TH3F(Form("htbXYcor%d", i), Form("TDC Time Difference for Group B Board %d;x/mm;y/mm;#Delta_t/ns", i), 48, 0, 48, 48, 16, 64, 5000, -50, 50);
    }
    for (int i = 0; i < gBoardCount; i++)
    {
        htabXYcor[i] = new TH3F(Form("htabXYcor%d", i), Form("TDC Time Difference for Group A Board %d;x/mm;y/mm;#Delta_t/ns", i), 48, 0, 48, 48, 16, 64, 5000, -50, 50);
    }
    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile.GetEntry(entry);
        // Do something

        if (corfile.decodeFlag[0] == -3 || corfile.decodeFlag[1] == -3)
            continue;

        int strip[gBoardCount];
        double ta[gBoardCount], tb[gBoardCount];

        // Find valid strip, ta, tb for Board 0
        for (int board = 0; board < gBoardCount; board++)
        {
            int strip0 = corfile.firedStrip1[board];
            int firedSiPM = corfile.processedCh[board][0];
            int strip1 = strip0 + 1;
            int a0 = UserDefine::GetFiberCode(strip0, UserDefine::a);
            int b0 = UserDefine::GetFiberCode(strip0, UserDefine::b);
            int a1 = UserDefine::GetFiberCode(strip1, UserDefine::a);
            int b1 = UserDefine::GetFiberCode(strip1, UserDefine::b);

            int validStrip;
            int validCha, validChb;
            bool flag;
            if (a0 == firedSiPM || b0 == firedSiPM)
            {
                validStrip = strip0;
                validCha = a0;
                validChb = b0;
            }
            else if (a1 == firedSiPM || b1 == firedSiPM)
            {
                validStrip = strip1;
                validCha = a1;
                validChb = b1;
            }
            else
            {
                flag = false;
                continue;
            }
            strip[board] = validStrip;
            ta[board] = corfile.processedTDC[board][validCha];
            tb[board] = corfile.processedTDC[board][validChb];
        }

        for (int group = 0; group < gBoardCount / 2; group++)
        {
            int x0strip = strip[group * 2];
            int x1strip = strip[group * 2 + 1];
            htaXYcor[group]->Fill(x0strip, x1strip, ta[group * 2] - ta[group * 2 + 1] - (gDeltaTa[group * 2][x0strip] - gDeltaTa[group * 2 + 1][x1strip]));
            htbXYcor[group]->Fill(x0strip, x1strip, tb[group * 2] - tb[group * 2 + 1] - (gDeltaTb[group * 2][x0strip] - gDeltaTb[group * 2 + 1][x1strip]));

            htabXYcor[group * 2]->Fill(x0strip, x1strip, ta[group * 2] - tb[group * 2] - (gDeltaTa[group * 2][x0strip] - gDeltaTb[group * 2][x0strip]));
            htabXYcor[group * 2 + 1]->Fill(x0strip, x1strip, ta[group * 2 + 1] - tb[group * 2 + 1] - (gDeltaTa[group * 2 + 1][x1strip] - gDeltaTb[group * 2 + 1][x1strip]));
        }
    }

    for (int group = 0; group < gBoardCount / 2; group++)
    {
        htaXYcor[group]->Draw("COLZ");
        auto h01a = htaXYcor[group]->Project3DProfile("yx");
        htaXYcor_pyx[group] = h01a;
        h01a->SetTitle(Form("TDC Time Difference for A fibers in Group %d;x strip;y strip;#Delta_t/ns", group));
        h01a->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-x-y-a-cor-%d.jpg", group));
        htaXYcor[group]->Write();
        h01a->Write();

        htbXYcor[group]->Draw("COLZ");
        auto h01b = htbXYcor[group]->Project3DProfile("yx");
        htbXYcor_pyx[group] = h01b;
        h01b->SetTitle(Form("TDC Time Difference for B fibers in Group %d;x strip;y strip;#Delta_t/ns", group));
        h01b->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-x-y-b-cor-%d.jpg", group));
        htbXYcor[group]->Write();
        h01b->Write();
    }
    for (int board = 0; board < gBoardCount; board++)
    {
        htabXYcor[board]->Draw("COLZ");
        auto htab = htabXYcor[board]->Project3DProfile("yx");
        htabXYcor_pyx[board] = htab;
        htab->SetTitle(Form("TDC Time Difference for Group ab in Board %d;x strip;y strip;#Delta_t/ns", board));
        htab->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-ab-cor-%d.jpg", board));
        htabXYcor[board]->Write();
        htab->Write();
    }

    fileW->Close();
    delete fileW;
    std::ofstream fout("Sum/tdcTimeDiff.txt");
    for (int i = 0; i < gBoardCount; i++)
    {
        for (int j = 0; j < gStripCount; j++)
        {
            fout << gDeltaTa[i][j] << '\t' << gDeltaTb[i][j] << std::endl;
        }
    }
    fout.close();
}

// 6. Draw time difference for system
void TDCCorrectionForSystem()
{
    std::ifstream fin("Sum/tdcTimeDiff.txt");
    if (!fin.is_open())
    {
        std::cout << "Error: Cannot open Sum/tdcTimeDiff.txt" << std::endl;

        TDCCorrectionForEachBoard();
    }
    else
    {
        for (int i = 0; i < gBoardCount; i++)
        {
            for (int j = 0; j < gStripCount; j++)
            {
                fin >> gDeltaTa[i][j] >> gDeltaTb[i][j];
                std::cout << gDeltaTa[i][j] << '\t' << gDeltaTb[i][j] << std::endl;
            }
        }
        fin.close();
    }

    TH1F *hxt = new TH1F("hxt", "TDC Time Difference ;#Delta t/ns", 1000, 60, 100);
    TH2F *hxtx = new TH2F("hxtx", "TDC Time Difference ;#DeltaX/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);
    TH2F *hxty = new TH2F("hxty", "TDC Time Difference ;#DeltaY/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);
    TH2F *hxttheta = new TH2F("hxttheta", "TDC Time Difference ;#theta;#Delta t/ns", 1000, 0, 15, 1000, 60, 100);
    TH3F *hxtxy = new TH3F("hxtxy", "TDC Time Difference ;x/mm;y/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);
    TH2F *hxtxyd = new TH2F("hxtxyd", "TDC Time Difference ;#Deltax-#Deltay/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);

    TH1F *hyt = new TH1F("hyt", "TDC Time Difference ;#Delta t/ns", 1000, 50, 90);
    TH2F *hytx = new TH2F("hytx", "TDC Time Difference ;#DeltaX/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 50, 90);
    TH2F *hyty = new TH2F("hyty", "TDC Time Difference ;#DeltaY/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 50, 90);
    TH2F *hyttheta = new TH2F("hyttheta", "TDC Time Difference ;#theta;#Delta t/ns", 1000, 0, 15, 1000, 50, 90);
    TH3F *hytxy = new TH3F("hytxy", "TDC Time Difference ;x/mm;y/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 480, -480 * 1.1, 480 * 1.1, 1000, 50, 90);
    TH2F *hytxyd = new TH2F("hytxyd", "TDC Time Difference ;#Deltax-#Deltay/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);

    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile.GetEntry(entry);
        // Do something

        if (corfile.decodeFlag[0] == -3 || corfile.decodeFlag[1] == -3)
            continue;

        int strip[gBoardCount];
        double ta[gBoardCount], tb[gBoardCount];

        // Find valid strip, ta, tb for Board 0
        for (int board = 0; board < gBoardCount; board++)
        {
            int strip0 = corfile.firedStrip1[board];
            int firedSiPM = corfile.processedCh[board][0];
            int strip1 = strip0 + 1;
            int a0 = UserDefine::GetFiberCode(strip0, UserDefine::a);
            int b0 = UserDefine::GetFiberCode(strip0, UserDefine::b);
            int a1 = UserDefine::GetFiberCode(strip1, UserDefine::a);
            int b1 = UserDefine::GetFiberCode(strip1, UserDefine::b);

            int validStrip;
            int validCha, validChb;
            bool flag;
            if (a0 == firedSiPM || b0 == firedSiPM)
            {
                validStrip = strip0;
                validCha = a0;
                validChb = b0;
            }
            else if (a1 == firedSiPM || b1 == firedSiPM)
            {
                validStrip = strip1;
                validCha = a1;
                validChb = b1;
            }
            else
            {
                flag = false;
                continue;
            }
            strip[board] = validStrip;
            ta[board] = corfile.processedTDC[board][validCha] - gDeltaTa[board][validStrip];
            tb[board] = corfile.processedTDC[board][validChb] - gDeltaTb[board][validStrip];
        }

        double theta = TMath::ATan(TMath::Sqrt(TMath::Power(corfile.tanX, 2) + TMath::Power(corfile.tanY, 2))) / TMath::Pi() * 180;
        double phi = TMath::ATan(TMath::Abs(corfile.tanX / corfile.tanY)) / TMath::Pi() * 180;
        double xdev = corfile.calcPos3[0] - corfile.calcPos3[6];
        double ydev = corfile.calcPos3[1] - corfile.calcPos3[7];
        double zdev = Z_POS[0] - Z_POS[3];
        double r = TMath::Sqrt(xdev * xdev + ydev * ydev + zdev * zdev);
        double tdev = r / 299.792458;
        // tdev = 10;
        hxt->Fill(ta[0] - tb[6] + tdev);
        hxtx->Fill(xdev, ta[0] - tb[6] + tdev);
        hxty->Fill(ydev, ta[0] - tb[6] + tdev);
        hxttheta->Fill(theta, ta[0] - tb[6] + tdev);
        hxtxy->Fill(xdev, ydev, ta[0] - tb[6] + tdev);
        hxtxyd->Fill(xdev - ydev, ta[0] - tb[6] + tdev);

        hyt->Fill(ta[1] - tb[7] + tdev);
        hytx->Fill(xdev, ta[1] - tb[7] + tdev);
        hyty->Fill(ydev, ta[1] - tb[7] + tdev);
        hyttheta->Fill(theta, ta[1] - tb[7] + tdev);
        hytxy->Fill(xdev, ydev, ta[1] - tb[7] + tdev);
        hytxyd->Fill(xdev - ydev, ta[1] - tb[7] + tdev);
    }

    hxttheta->Draw("COLZ");
    auto hxtxyp = hxtxy->Project3DProfile("yx");
    hxtxyp->SetTitle("TDC Time Difference ;#Delta x/mm;#Delta y/mm;t/ns");
    hxtxyp->GetZaxis()->SetRangeUser(70, 90);
    hxtxyp->Draw("COLZ");

    auto fpol1 = new TF1("fpol1", "pol1", -1000, 1000);
    auto hxtxp = hxtx->ProfileX();
    hxtxp->GetYaxis()->SetRangeUser(60, 90);
    hxtxp->Fit(fpol1, "", "", -400, 400);
    double p1xx = fpol1->GetParameter(1);
    auto hxtyp = hxty->ProfileX();
    hxtyp->Fit(fpol1, "", "", -400, 400);
    hxtyp->GetYaxis()->SetRangeUser(60, 90);
    double p1xy = fpol1->GetParameter(1);
    double p1x = (TMath::Abs(p1xx) + TMath::Abs(p1xy)) / 2.0;

    auto hytxp = hytx->ProfileX();
    hytxp->GetYaxis()->SetRangeUser(60, 90);
    hytxp->Fit(fpol1, "", "", -400, 400);
    double p1yx = fpol1->GetParameter(1);
    auto hytyp = hyty->ProfileX();
    hytyp->Fit(fpol1, "", "", -400, 400);
    hytyp->GetYaxis()->SetRangeUser(60, 90);
    double p1yy = fpol1->GetParameter(1);
    double p1y = (TMath::Abs(p1yx) + TMath::Abs(p1yy)) / 2.0;

    double p1 = (p1x + p1y) / 2.0;

    gStyle->SetOptFit(111);
    hxtxyd->Fit(fpol1, "", "", -400, 400);
    hxtxyd->Draw("COLZ");
    c->SaveAs("Sum/tdcTimeDiff-tx-xyd.jpg");
    double p1dx = TMath::Abs(fpol1->GetParameter(1));

    hytxyd->Fit(fpol1, "", "", -400, 400);
    hytxyd->Draw("COLZ");
    c->SaveAs("Sum/tdcTimeDiff-ty-yxd.jpg");
    double p1dy = TMath::Abs(fpol1->GetParameter(1));
    gStyle->SetOptFit(0000);

    std::cout << p1xx << '\t' << p1xy << '\t' << p1x << '\t' << p1dx << std::endl;
    std::cout << p1yx << '\t' << p1yy << '\t' << p1y << '\t' << p1dy << std::endl;

    p1 = (p1dx + p1dy) / 2.0;

    // return;

    TH1F *hxt2 = new TH1F("hxt2", "TDC Time Difference ;#Delta t/ns", 1000, 60, 100);
    TH2F *hxtx2 = new TH2F("hxtx2", "TDC Time Difference ;#DeltaX/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);
    TH2F *hxty2 = new TH2F("hxty2", "TDC Time Difference ;#DeltaY/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);
    TH2F *hxttheta2 = new TH2F("hxttheta2", "TDC Time Difference ;#theta;#Delta t/ns", 1000, 0, 15, 1000, 60, 100);
    TH3F *hxtxy2 = new TH3F("hxtxy2", "TDC Time Difference ;x/mm;y/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);
    TH2F *hxtxyd2 = new TH2F("hxtxyd2", "TDC Time Difference ;#Deltax-#Deltay/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);

    TH1F *hyt2 = new TH1F("hyt2", "TDC Time Difference ;#Delta t/ns", 1000, 50, 90);
    TH2F *hytx2 = new TH2F("hytx2", "TDC Time Difference ;#DeltaX/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 50, 90);
    TH2F *hyty2 = new TH2F("hyty2", "TDC Time Difference ;#DeltaY/mm;#Delta t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 50, 90);
    TH2F *hyttheta2 = new TH2F("hyttheta2", "TDC Time Difference ;#theta;#Delta t/ns", 1000, 0, 15, 1000, 50, 90);
    TH3F *hytxy2 = new TH3F("hytxy2", "TDC Time Difference ;x/mm;y/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 480, -480 * 1.1, 480 * 1.1, 1000, 50, 90);
    TH2F *hytxyd2 = new TH2F("hytxyd2", "TDC Time Difference ;#Deltax-#Deltay/mm;t/ns", 480, -480 * 1.1, 480 * 1.1, 1000, 60, 100);

    TH2F *hlxt = new TH2F("hlxt", "TDC Time Difference ;l/mm;#Delta t/ns", 100, 2737, 2837, 100, 60, 100);
    TH2F *hlyt = new TH2F("hlyt", "TDC Time Difference ;l/mm;#Delta t/ns", 100, 2737, 2837, 100, 50, 90);

    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile.GetEntry(entry);
        // Do something

        if (corfile.decodeFlag[0] == -3 || corfile.decodeFlag[1] == -3)
            continue;

        int strip[gBoardCount];
        double ta[gBoardCount], tb[gBoardCount];

        // Find valid strip, ta, tb for Board 0
        for (int board = 0; board < gBoardCount; board++)
        {
            int strip0 = corfile.firedStrip1[board];
            int firedSiPM = corfile.processedCh[board][0];
            int strip1 = strip0 + 1;
            int a0 = UserDefine::GetFiberCode(strip0, UserDefine::a);
            int b0 = UserDefine::GetFiberCode(strip0, UserDefine::b);
            int a1 = UserDefine::GetFiberCode(strip1, UserDefine::a);
            int b1 = UserDefine::GetFiberCode(strip1, UserDefine::b);

            int validStrip;
            int validCha, validChb;
            bool flag;
            if (a0 == firedSiPM || b0 == firedSiPM)
            {
                validStrip = strip0;
                validCha = a0;
                validChb = b0;
            }
            else if (a1 == firedSiPM || b1 == firedSiPM)
            {
                validStrip = strip1;
                validCha = a1;
                validChb = b1;
            }
            else
            {
                flag = false;
                continue;
            }
            strip[board] = validStrip;
            ta[board] = corfile.processedTDC[board][validCha] - gDeltaTa[board][validStrip];
            tb[board] = corfile.processedTDC[board][validChb] - gDeltaTb[board][validStrip];
        }

        double theta = TMath::ATan(TMath::Sqrt(TMath::Power(corfile.tanX, 2) + TMath::Power(corfile.tanY, 2))) / TMath::Pi() * 180;
        double phi = TMath::ATan(TMath::Abs(corfile.tanX / corfile.tanY)) / TMath::Pi() * 180;
        double xdev = corfile.calcPos3[0] - corfile.calcPos3[6];
        double ydev = corfile.calcPos3[1] - corfile.calcPos3[7];
        double zdev = Z_POS[0] - Z_POS[3];
        double r = TMath::Sqrt(xdev * xdev + ydev * ydev + zdev * zdev);
        double tdev = r / 299.792458;
        // double tdev = p1 * (xdev - ydev);
        tdev += p1 * (xdev - ydev);

        hxt2->Fill(ta[0] - tb[6] + tdev);
        hxtx2->Fill(xdev, ta[0] - tb[6] + tdev);
        hxty2->Fill(ydev, ta[0] - tb[6] + tdev);
        hxttheta2->Fill(theta, ta[0] - tb[6] + tdev);
        hxtxy2->Fill(xdev, ydev, ta[0] - tb[6] + tdev);
        hxtxyd2->Fill(xdev - ydev, ta[0] - tb[6] + tdev);

        hyt2->Fill(ta[1] - tb[7] + tdev);
        hytx2->Fill(xdev, ta[1] - tb[7] + tdev);
        hyty2->Fill(ydev, ta[1] - tb[7] + tdev);
        hyttheta2->Fill(theta, ta[1] - tb[7] + tdev);
        hytxy2->Fill(xdev, ydev, ta[1] - tb[7] + tdev);
        hytxyd2->Fill(xdev - ydev, ta[1] - tb[7] + tdev);

        hlxt->Fill(r, ta[0] - tb[6] + tdev);
        hlyt->Fill(r, ta[1] - tb[7] + tdev);
    }

    hlxt->Draw("COLZ");
    auto hlxtp = hlxt->ProfileX();
    hlxtp->Draw();
    hlxtp->GetYaxis()->SetRangeUser(80, 80.5);
    hlxtp->Fit(fpol1, "", "", 2737, 2787);
    c->SaveAs("Sum/tdcTimeDiff-l-x.jpg");

    hlyt->Draw("COLZ");
    auto hlytp = hlyt->ProfileX();
    hlytp->Draw();
    hlytp->GetYaxis()->SetRangeUser(67, 67.5);
    hlytp->Fit(fpol1, "", "", 2737, 2787);
    c->SaveAs("Sum/tdcTimeDiff-l-y.jpg");

    auto fileW = new TFile("Sum/tdcCor.root", "recreate");
    fileW->cd();
    hxt->Write();
    hxtx->Write();
    hxty->Write();
    hxttheta->Write();
    hxtxy->Write();
    hxtxyd->Write();

    hyt->Write();
    hytx->Write();
    hyty->Write();
    hyttheta->Write();
    hytxy->Write();
    hytxyd->Write();

    hxt2->Write();
    hxtx2->Write();
    hxty2->Write();
    hxttheta2->Write();
    hxtxy2->Write();
    hxtxyd2->Write();

    hyt2->Write();
    hytx2->Write();
    hyty2->Write();
    hyttheta2->Write();
    hytxy2->Write();
    hytxyd2->Write();

    hlxt->Write();
    hlyt->Write();

    fileW->Close();
    delete fileW;
}

// 7. Draw time difference for each board, only use strip
void TDCStripCorrelation()
{
    TH2F *hsatal[64];
    TH2F *hsatar[64];
    TH2F *hsbtal[64];
    TH2F *hsbtar[64];
    for (int i = 0; i < 64; i++)
    {
        hsatal[i] = new TH2F(Form("hsatal%d", i), Form("TDC Time Difference for A fibers in Strip %d;Strip %d Amplitude;tl-t0/ns", i, i - 1), 100, 0, 3, 1000, -30, 30);
        hsatar[i] = new TH2F(Form("hsatar%d", i), Form("TDC Time Difference for A fibers in Strip %d;Strip %d Amplitude;tr-t0/ns", i, i + 1), 100, 0, 3, 1000, -30, 30);
        hsbtal[i] = new TH2F(Form("hsbtal%d", i), Form("TDC Time Difference for B fibers in Strip %d;Strip %d Amplitude;tl-t0/ns", i, i - 1), 100, 0, 3, 1000, -30, 30);
        hsbtar[i] = new TH2F(Form("hsbtar%d", i), Form("TDC Time Difference for B fibers in Strip %d;Strip %d Amplitude;tr-t0/ns", i, i + 1), 100, 0, 3, 1000, -30, 30);
    }

    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile.GetEntry(entry);
        // Do something

        // if (corfile.decodeFlag[0] == -3)
        //     continue;

        int strip[gBoardCount];
        double ta[gBoardCount], tb[gBoardCount];

        // Find valid strip, ta, tb for Board 0
        // for (int board = 0; board < gBoardCount; board++)
        for (int board = 0; board < 1; board++)
        {
            if (corfile.weight2[board] == 0)
                continue;

            int strip0 = corfile.firedStrip1[board];
            int firedSiPM = corfile.processedCh[board][0];
            int strip1 = strip0 + 1;
            int a0 = UserDefine::GetFiberCode(strip0, UserDefine::a);
            int b0 = UserDefine::GetFiberCode(strip0, UserDefine::b);
            int a1 = UserDefine::GetFiberCode(strip1, UserDefine::a);
            int b1 = UserDefine::GetFiberCode(strip1, UserDefine::b);

            int validStrip, neighborStrip;
            int validCha, validChb, neighborCha, neighborChb;
            double validAmpa, neighborAmpa, validAmpb, neighborAmpb;
            double validTima, neighborTima, validTimb, neighborTimb;
            bool flag;
            // Valid strip is strip0, or weight 1
            if (a0 == firedSiPM || b0 == firedSiPM)
            {
                validStrip = strip0;
                neighborStrip = strip1;
                validCha = a0;
                validChb = b0;
                neighborCha = a1;
                neighborChb = b1;

                validAmpa = corfile.weight1a[board];
                neighborAmpa = corfile.weight2a[board];
                validAmpb = corfile.weight1b[board];
                neighborAmpb = corfile.weight2b[board];
                flag = 1;
            }
            else if (a1 == firedSiPM || b1 == firedSiPM)
            // Valid strip is strip1, or weight 2
            {
                validStrip = strip1;
                neighborStrip = strip0;
                validCha = a1;
                validChb = b1;
                neighborCha = a0;
                neighborChb = b0;

                validAmpa = corfile.weight2a[board];
                neighborAmpa = corfile.weight1a[board];
                validAmpb = corfile.weight2b[board];
                neighborAmpb = corfile.weight1b[board];

                flag = 0;
            }
            else
                continue;

            validTima = corfile.processedTDC[board][validCha];
            validTimb = corfile.processedTDC[board][validChb];
            neighborTima = corfile.processedTDC[board][neighborCha];
            neighborTimb = corfile.processedTDC[board][neighborChb];

            if (flag)
            {
                hsatar[validStrip]->Fill(neighborAmpa, neighborTima - validTima);
                hsbtar[validStrip]->Fill(neighborAmpb, neighborTimb - validTima);

                hsatal[neighborStrip]->Fill(validAmpa, validTima - neighborTima);
                hsbtal[neighborStrip]->Fill(validAmpb, validTimb - neighborTimb);
            }
            else
            {
                hsatal[validStrip]->Fill(neighborAmpa, neighborTima - validTima);
                hsbtal[validStrip]->Fill(neighborAmpb, neighborTimb - validTima);

                hsatar[neighborStrip]->Fill(validAmpa, validTima - neighborTima);
                hsbtar[neighborStrip]->Fill(validAmpb, validTimb - neighborTimb);
            }
        }
    }

    // Write to file
    auto fileW = new TFile("Sum/tdcTimeDiff-strip.root", "recreate");
    fileW->cd();
    for (int i = 0; i < 64; i++)
    {
        hsatal[i]->Write();
        hsatar[i]->Write();
        hsbtal[i]->Write();
        hsbtar[i]->Write();

        // Draw histos
        hsatal[i]->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-strip-a-l-%d.jpg", i));
        hsatar[i]->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-strip-a-r-%d.jpg", i));
        hsbtal[i]->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-strip-b-l-%d.jpg", i));
        hsbtar[i]->Draw("COLZ");
        c->SaveAs(Form("Sum/tdcTimeDiff-strip-b-r-%d.jpg", i));
    }
    fileW->Close();
    delete fileW;
}

// 8. Draw strip alignment, Write to Sum/stripAlign.root
void StripAlignStep1()
{
    auto fileW = new TFile("Sum/stripAlign.root", "recreate");

    TGraph *tgXTemp = new TGraph();
    TGraph *tgYTemp = new TGraph();
    TF1 *fpolX = new TF1("fpolX", "pol1", -1000, 1000);
    TF1 *fpolY = new TF1("fpolY", "pol1", -1000, 1000);

    auto tree = new TTree("align", "Strip Alignment");
    double calcPos3[gBoardCount];
    double derivedPos[gBoardCount];
    int strip1[gBoardCount];
    int strip2[gBoardCount];
    double weight1[gBoardCount], weight2[gBoardCount];
    double weight1a[gBoardCount], weight1b[gBoardCount], weight2a[gBoardCount], weight2b[gBoardCount];
    tree->Branch("calcPos3", calcPos3, Form("calcPos3[%d]/D", gBoardCount));
    tree->Branch("derivedPos", derivedPos, Form("derivedPos[%d]/D", gBoardCount));
    tree->Branch("strip1", strip1, Form("strip1[%d]/I", gBoardCount));
    tree->Branch("strip2", strip2, Form("strip2[%d]/I", gBoardCount));
    tree->Branch("weight1", weight1, Form("weight1[%d]/D", gBoardCount));
    tree->Branch("weight2", weight2, Form("weight2[%d]/D", gBoardCount));
    tree->Branch("weight1a", weight1a, Form("weight1a[%d]/D", gBoardCount));
    tree->Branch("weight1b", weight1b, Form("weight1b[%d]/D", gBoardCount));
    tree->Branch("weight2a", weight2a, Form("weight2a[%d]/D", gBoardCount));
    tree->Branch("weight2b", weight2b, Form("weight2b[%d]/D", gBoardCount));

    double zArrayX[4], zArrayY[4];
    for (int layer = 0; layer < 4; layer++)
    {
        zArrayX[layer] = Z_POS[layer] + layerDistance / 2.0;
        zArrayY[layer] = Z_POS[layer] - layerDistance / 2.0;
    }

    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile.GetEntry(entry);
        // Do something

        // Calculate deviation for layer 1, layer 2
        double x0 = corfile.calcPos3[boardX[0]], x1 = corfile.calcPos3[boardX[1]], x2 = corfile.calcPos3[boardX[2]], x3 = corfile.calcPos3[boardX[3]];
        double y0 = corfile.calcPos3[boardY[0]], y1 = corfile.calcPos3[boardY[1]], y2 = corfile.calcPos3[boardY[2]], y3 = corfile.calcPos3[boardY[3]];
        double z0 = Z_POS[0], z1 = Z_POS[1], z2 = Z_POS[2], z3 = Z_POS[3];
        double z0X = zArrayX[0], z1X = zArrayX[1], z2X = zArrayX[2], z3X = zArrayX[3];
        double z0Y = zArrayY[0], z1Y = zArrayY[1], z2Y = zArrayY[2], z3Y = zArrayY[3];

        // Calculate derived position
        double xArray[4] = {x0, x1, x2, x3};
        double yArray[4] = {y0, y1, y2, y3};
        double zArray[4] = {z0, z1, z2, z3};

        for (int layer = 0; layer < 4; layer++)
        {
            tgXTemp->Set(0);
            tgYTemp->Set(0);
            for (int i = 0; i < 4; i++)
            {
                if (i == layer)
                    continue;
                tgXTemp->SetPoint(tgXTemp->GetN(), zArray[i], xArray[i]);
                tgYTemp->SetPoint(tgYTemp->GetN(), zArray[i], yArray[i]);
            }

            double xDerived = 0, yDerived = 0;
            tgXTemp->Fit(fpolX, "Q");
            xDerived = fpolX->Eval(zArray[layer]);
            derivedPos[boardX[layer]] = xDerived;

            tgYTemp->Fit(fpolY, "Q");
            yDerived = fpolY->Eval(zArray[layer]);
            derivedPos[boardY[layer]] = yDerived;
        }

        for (int board = 0; board < gBoardCount; board++)
        {
            calcPos3[board] = corfile.calcPos3[board];
            strip1[board] = corfile.firedStrip1[board];
            strip2[board] = strip1[board] + 1;
            weight1[board] = corfile.weight1[board];
            weight2[board] = corfile.weight2[board];
            weight1a[board] = corfile.weight1a[board];
            weight1b[board] = corfile.weight1b[board];
            weight2a[board] = corfile.weight2a[board];
            weight2b[board] = corfile.weight2b[board];

            if (corfile.decodeFlag[board] == 1)
            {
                strip2[board] = -1;
                continue;
            }
            if (strip2[board] >= 64)
            {
                strip2[board] = -1;
                std::cout << board << '\t' << strip2[board] << '\t' << weight2[board] << std::endl;
                continue;
            }
        }
        tree->Fill();
    }
    // Write to file
    fileW->cd();
    tree->Write();
    fileW->Close();
    delete fileW;
}

// 9. Read Sum/striptAlign.root file, generate histos & graphs, write into Sum/Strip-Align/saPlots.root
void StripAlignStep2()
{
    auto fAlign = new TFile("Sum/stripAlign.root");
    if (fAlign->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/stripAlign.root file, try to create." << std::endl;
        StripAlignStep1();
        fAlign = new TFile("Sum/stripAlign.root");
    }
    auto treeAlign = (TTree *)fAlign->Get("align");
    if (!treeAlign)
    {
        std::cout << "Error: Cannot find align tree" << std::endl;
        return;
    }

    double derivedPos[gBoardCount];
    int strip1[gBoardCount];
    int strip2[gBoardCount];
    double weight1[gBoardCount], weight2[gBoardCount];
    double weight1a[gBoardCount], weight1b[gBoardCount], weight2a[gBoardCount], weight2b[gBoardCount];
    treeAlign->SetBranchAddress("derivedPos", derivedPos);
    treeAlign->SetBranchAddress("strip1", strip1);
    treeAlign->SetBranchAddress("strip2", strip2);
    treeAlign->SetBranchAddress("weight1", weight1);
    treeAlign->SetBranchAddress("weight2", weight2);
    treeAlign->SetBranchAddress("weight1a", weight1a);
    treeAlign->SetBranchAddress("weight1b", weight1b);
    treeAlign->SetBranchAddress("weight2a", weight2a);
    treeAlign->SetBranchAddress("weight2b", weight2b);

    auto fileW = new TFile("Sum/Strip-Align/saPlots.root", "recreate");

    TH2F *hsna[gBoardCount][64];
    TH2F *hsnb[gBoardCount][64];
    TH2F *hsnab[gBoardCount][64];

    TGraph *tgsna[gBoardCount][64];
    TGraph *tgsnb[gBoardCount][64];
    TGraph *tgsnab[gBoardCount][64];

    for (int board = 0; board < gBoardCount; board++)
    {
        for (int strip = 0; strip < 64; strip++)
        {
            hsna[board][strip] = new TH2F(Form("hsna%d_%d", board, strip), Form("Strip a Alignment for Board %d Strip %d;x_{fit}/mm;Q_{a}/ns", board, strip), 100, strip * 11 - 20, strip * 11 + 20, 100, 0, 1.5);
            hsnb[board][strip] = new TH2F(Form("hsnb%d_%d", board, strip), Form("Strip a Alignment for Board %d Strip %d;x_{fit}/mm;Q_{a}/ns", board, strip), 100, strip * 11 - 20, strip * 11 + 20, 100, 0, 1.5);
            hsnab[board][strip] = new TH2F(Form("hsnab%d_%d", board, strip), Form("Strip a Alignment for Board %d Strip %d;x_{fit}/mm;Q_{a}/ns", board, strip), 100, strip * 11 - 20, strip * 11 + 20, 100, 0, 3);

            tgsna[board][strip] = new TGraph();
            tgsnb[board][strip] = new TGraph();
            tgsnab[board][strip] = new TGraph();
            tgsna[board][strip]->SetName(Form("tgsna%d_%d", board, strip));
            tgsnb[board][strip]->SetName(Form("tgsnb%d_%d", board, strip));
            tgsnab[board][strip]->SetName(Form("tgsnab%d_%d", board, strip));
        }
    }

    for (int entry = 0; entry < treeAlign->GetEntries(); entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        treeAlign->GetEntry(entry);
        for (int board = 0; board < gBoardCount; board++)
        {

            hsna[board][strip1[board]]->Fill(derivedPos[board], weight1a[board]);
            hsnb[board][strip1[board]]->Fill(derivedPos[board], weight1b[board]);
            hsnab[board][strip1[board]]->Fill(derivedPos[board], weight1a[board] + weight1b[board]);

            if (TMath::Abs(derivedPos[board] - strip1[board] * 11) < 15)
            {
                tgsna[board][strip1[board]]->SetPoint(tgsna[board][strip1[board]]->GetN(), derivedPos[board], weight1a[board]);
                tgsnb[board][strip1[board]]->SetPoint(tgsnb[board][strip1[board]]->GetN(), derivedPos[board], weight1b[board]);
                tgsnab[board][strip1[board]]->SetPoint(tgsnab[board][strip1[board]]->GetN(), derivedPos[board], weight1a[board] + weight1b[board]);
            }

            if (strip2[board] < 0 || strip2[board] > 63)
                continue;

            hsna[board][strip2[board]]->Fill(derivedPos[board], weight2a[board]);
            hsnb[board][strip2[board]]->Fill(derivedPos[board], weight2b[board]);
            hsnab[board][strip2[board]]->Fill(derivedPos[board], weight2a[board] + weight2b[board]);

            if (TMath::Abs(derivedPos[board] - strip2[board] * 11) < 15)
            {
                tgsna[board][strip2[board]]->SetPoint(tgsna[board][strip2[board]]->GetN(), derivedPos[board], weight2a[board]);
                tgsnb[board][strip2[board]]->SetPoint(tgsnb[board][strip2[board]]->GetN(), derivedPos[board], weight2b[board]);
                tgsnab[board][strip2[board]]->SetPoint(tgsnab[board][strip2[board]]->GetN(), derivedPos[board], weight2a[board] + weight2b[board]);
            }
        }
    }
    // Write to file
    fileW->cd();
    for (int board = 0; board < gBoardCount; board++)
    {
        for (int strip = 0; strip < 64; strip++)
        {
            hsna[board][strip]->Write();
            hsnb[board][strip]->Write();
            hsnab[board][strip]->Write();
            tgsna[board][strip]->Write(Form("tgsna%d_%d", board, strip));
            tgsnb[board][strip]->Write(Form("tgsnb%d_%d", board, strip));
            tgsnab[board][strip]->Write(Form("tgsnab%d_%d", board, strip));
        }
    }
    fileW->Close();
    delete fileW;
    delete fAlign;
}

// 10. Read Sum/Strip-Align/saPlots.root file
void StripAlignStep3()
{
    auto fAlign = new TFile("Sum/Strip-Align/saPlots.root");
    if (fAlign->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/Strip-Align/saPlots.root file, try to create." << std::endl;
        // StripAlignStep1();
        fAlign = new TFile("Sum/Strip-Align/saPlots.root");
    }

    TH2F *hsna[gBoardCount][64];
    TH2F *hsnb[gBoardCount][64];
    TH2F *hsnab[gBoardCount][64];
    TGraph *tgsna[gBoardCount][64];
    TGraph *tgsnb[gBoardCount][64];
    TGraph *tgsnab[gBoardCount][64];

    auto fstrip = new TF1("fstrip", FitStripCenter, -1000, 1000, 2);
    fstrip->SetParName(0, "center");
    fstrip->SetParName(1, "amp");

    double stripAlignA[gBoardCount][64]; // Group A
    double stripAlignB[gBoardCount][64]; // Group B
    double stripAlignC[gBoardCount][64]; // Center

    double stripAmpA[gBoardCount][64]; // Group A
    double stripAmpB[gBoardCount][64]; // Group B
    double stripAmpC[gBoardCount][64]; // Center

    gSystem->Exec("mkdir Sum/Strip-Align/");
    for (int board = 0; board < gBoardCount; board++)
    // for (int board = 0; board < 1; board++)
    {
        gSystem->Exec(Form("mkdir Sum/Strip-Align/Board%d/", board));
        for (int strip = 0; strip < 64; strip++)
        // for (int strip = 0; strip < 1; strip++)
        {
            auto htsna = (TH2F *)fAlign->Get(Form("hsna%d_%d", board, strip));
            auto htsnb = (TH2F *)fAlign->Get(Form("hsnb%d_%d", board, strip));
            auto htsnab = (TH2F *)fAlign->Get(Form("hsnab%d_%d", board, strip));

            auto tgsna = (TGraph *)fAlign->Get(Form("tgsna%d_%d", board, strip));
            auto tgsnb = (TGraph *)fAlign->Get(Form("tgsnb%d_%d", board, strip));
            auto tgsnab = (TGraph *)fAlign->Get(Form("tgsnab%d_%d", board, strip));

            hsna[board][strip] = htsna;
            hsnb[board][strip] = htsnb;
            hsnab[board][strip] = htsnab;

            if (!htsna || !htsnb || !htsnab)
            {
                std::cout << "Error: Cannot find histos" << std::endl;
                return;
            }
            if (!tgsna || !tgsnb || !tgsnab)
            {
                std::cout << "Error: Cannot find graphs" << std::endl;
                return;
            }

            gStyle->SetOptFit(111);

            fstrip->SetParLimits(0, strip * 11 - 30, strip * 11 + 30);
            fstrip->SetParLimits(1, 0.5, 3);

            fstrip->SetParameter(0, strip * 11);
            fstrip->SetParameter(1, 1);
            tgsna->GetXaxis()->SetRangeUser(strip * 11 - 30, strip * 11 + 30);
            tgsna->Fit(fstrip, "QN", "", strip * 11 - 9, strip * 11 + 9);
            tgsna->Fit(fstrip, "Q", "", fstrip->GetParameter(0) - 9, fstrip->GetParameter(0) + 9);
            // htsna->Draw("COLZ");
            // fstrip->Draw("same");
            tgsna->Draw("AZP");
            c->SaveAs(Form("Sum/Strip-Align/Board%d/stripAlign-a-%d-%d.jpg", board, board, strip));
            stripAlignA[board][strip] = fstrip->GetParameter(0);
            stripAmpA[board][strip] = fstrip->GetParameter(1);

            fstrip->SetParameter(0, strip * 11);
            fstrip->SetParameter(1, 1);
            tgsnb->GetXaxis()->SetRangeUser(strip * 11 - 30, strip * 11 + 30);
            tgsnb->Fit(fstrip, "QN", "", strip * 11 - 9, strip * 11 + 9);
            tgsnb->Fit(fstrip, "Q", "", fstrip->GetParameter(0) - 9, fstrip->GetParameter(0) + 9);
            // htsnb->Draw("COLZ");
            // fstrip->Draw("same");
            tgsnb->Draw("AZP");
            c->SaveAs(Form("Sum/Strip-Align/Board%d/stripAlign-b-%d-%d.jpg", board, board, strip));
            stripAlignB[board][strip] = fstrip->GetParameter(0);
            stripAmpB[board][strip] = fstrip->GetParameter(1);

            fstrip->SetParameter(0, strip * 11);
            fstrip->SetParameter(1, 2);
            tgsnab->GetXaxis()->SetRangeUser(strip * 11 - 30, strip * 11 + 30);
            tgsnab->Fit(fstrip, "QN", "", strip * 11 - 9, strip * 11 + 9);
            tgsnab->Fit(fstrip, "Q", "", fstrip->GetParameter(0) - 9, fstrip->GetParameter(0) + 9);
            // htsnab->Draw("COLZ");
            // fstrip->Draw("same");
            tgsnab->Draw("AZP");
            c->SaveAs(Form("Sum/Strip-Align/Board%d/stripAlign-ab-%d-%d.jpg", board, board, strip));
            stripAlignC[board][strip] = fstrip->GetParameter(0);
            stripAmpC[board][strip] = fstrip->GetParameter(1);
        }
    }

    // Output fit results
    std::ofstream fout("Sum/Strip-Align/stripAlign.txt");
    for (int board = 0; board < gBoardCount; board++)
    {
        for (int strip = 0; strip < 64; strip++)
        {
            fout << board << '\t' << strip << '\t' << stripAlignA[board][strip] << '\t' << stripAmpA[board][strip] << '\t' << stripAlignB[board][strip] << '\t' << stripAmpB[board][strip] << '\t' << stripAlignC[board][strip] << '\t' << stripAmpC[board][strip] << std::endl;
        }
    }
    fout.close();
}

void Analyze()
{
    // Open Sum/1stCorrection.root file
    if (f->IsZombie())
    {
        cout << "Error: Cannot open Sum/1stCorrection.root file" << endl;
        return;
    }
    // Get cor tree
    if (!cor)
    {
        cout << "Error: Cannot find cor tree" << endl;
        return;
    }

    // Analyze
    // 1. Draw time difference between t0 and t3
    if (0)
    {
        gStyle->SetOptStat(1);
        gStyle->SetOptFit(10);
        cor->Draw("(tdcTime[6]-lastSeg[6])/interval[6]*1e9-(tdcTime[0]-lastSeg[0])/interval[0]*1e9>>ht06(1000,0,-45)", "", "");
        TH1F *ht06 = (TH1F *)gDirectory->Get("ht06");
        ht06->SetTitle("TDC Time Difference;t0 - t3 / ns;Counts");
        ht06->Draw();
        ht06->Fit("gaus", "", "", -80, -50);
        c->SaveAs("t06.pdf");
    }
    // 2. Draw position difference between t0 and t3

    // 3. Draw theta, phi distribution
    auto fCos = new TF1("fCos", "[0]*(TMath::Sin(x/180*TMath::Pi())*TMath::Power(TMath::Cos(x/180*TMath::Pi()),2))", 0, 90);
    if (0)
    {
        gStyle->SetOptStat(1);
        gStyle->SetOptFit(111);
        cor->Draw("TMath::ATan(TMath::Sqrt(TMath::Power(tanX,2)+TMath::Power(tanY,2)))/TMath::Pi()*180:TMath::ATan(TMath::Abs(tanX/tanY))/TMath::Pi()*180>>hr50(180,0,90,150,0,15)", "TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))<50", "COLZ", 100e5);
        TH2F *hr50 = (TH2F *)gDirectory->Get("hr50");
        hr50->SetTitle("#theta vs #phi at r<50 mm;#phi/degree;#theta/degree");
        hr50->Draw("COLZ");
        c->SaveAs("Sum/theta-phi-r50.jpg");
        for (int i = 0; i < 18; i++)
        {
            auto proj = hr50->ProjectionY(Form("proj%d", i), i * 10, i * 10 + 9);
            proj->SetTitle(Form("Projection of #theta at #phi = %d degree;#theta/degree;Counts", i * 5));
            proj->Draw();
            if (TMath::Abs(i - 9) > 2)
                proj->Fit(fCos, "", "", 0, 8);
            else
                proj->Fit(fCos, "", "", 0, 10);
            c->SaveAs(Form("Sum/theta-phi-r50-%d.jpg", i));
        }

        // cor->Draw("TMath::ATan(TMath::Sqrt(TMath::Power(tanX,2)+TMath::Power(tanY,2)))/TMath::Pi()*180:TMath::ATan(TMath::Abs(tanX/tanY))/TMath::Pi()*180>>hr100(180,0,90,300,0,15)", "TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))<100&&TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))>50", "COLZ", 100e5);
        // TH2F *hr100 = (TH2F *)gDirectory->Get("hr100");
        // hr100->SetTitle("#theta vs #phi at 100 mm > r >50 mm ;#phi/degree;#theta/degree");
        // hr100->Draw("COLZ");
        // c->SaveAs("Sum/theta-phi-r100.jpg");

        // cor->Draw("TMath::ATan(TMath::Sqrt(TMath::Power(tanX,2)+TMath::Power(tanY,2)))/TMath::Pi()*180:TMath::ATan(TMath::Abs(tanX/tanY))/TMath::Pi()*180>>hr150(180,0,90,300,0,15)", "TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))<150&&TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))>100", "COLZ", 100e5);
        // TH2F *hr150 = (TH2F *)gDirectory->Get("hr150");
        // hr150->SetTitle("#theta vs #phi at 150 mm > r >100 mm ;#phi/degree;#theta/degree");
        // hr150->Draw("COLZ");
        // c->SaveAs("Sum/theta-phi-r150.jpg");

        // cor->Draw("TMath::ATan(TMath::Sqrt(TMath::Power(tanX,2)+TMath::Power(tanY,2)))/TMath::Pi()*180:TMath::ATan(TMath::Abs(tanX/tanY))/TMath::Pi()*180>>hr1000(180,0,90,300,0,15)", "TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))<1000&&TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))>150", "COLZ", 100e5);
        // TH2F *hr1000 = (TH2F *)gDirectory->Get("hr1000");
        // hr1000->SetTitle("#theta vs #phi at 1000 mm > r >150 mm ;#phi/degree;#theta/degree");
        // hr1000->Draw("COLZ");
        // c->SaveAs("Sum/theta-phi-r1000.jpg");
    }

    // 4. Draw time difference
    if (0)
    {
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        // cor->Draw("((tdcTime[0]-lastSeg[0])/interval[0]-(tdcTime[1]-lastSeg[1])/interval[1])*1e9:calcPos2[1]*1.1:calcPos2[0]*1.1>>ht01(480,0,528,480,176,704,100, 0,30)");
        // cor->Draw("((tdcTime[0]-lastSeg[0])/interval[0]-(tdcTime[1]-lastSeg[1])/interval[1])*1e9:calcPos2[1]*1.1:calcPos2[0]*1.1>>ht01(480,0,528,480,176,704,100, 0,30)");
        cor->Draw("processedTDC[0][processedCh[0][0]]-processedTDC[1][processedCh[1][0]]:calcPos2[1]*1.1:calcPos2[0]*1.1>>ht01(480,0,528,480,176,704,100, 0,30)");
        TH3F *ht01 = (TH3F *)gDirectory->Get("ht01");
        ht01->SetTitle("TDC Time Difference;x0/mm;x1/mm;#Delta_t/ns");
        auto ht01pyx = ht01->Project3DProfile("yx");
        ht01pyx->SetTitle("TDC Time Difference;x0/mm;x1/mm;#Delta_t/ns");
        ht01pyx->Draw("COLZ");
        c->SaveAs("Sum/tdcTimeDiff-x0-x1.jpg");
    }

    // TDCCorrectionForEachBoard();
    // TDCCorrectionForSystem();
    // StripAlignStep1();
    StripAlignStep2();
    StripAlignStep3();
}