#include <functional>
#include <TCanvas.h>
#include "cor.C"
#include "Decode.cpp"
std::string gRepoFolder = "../Sum";
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

TFile *f;
TTree *cor;
// typedef john::cor classcor;
john::cor *corfile;
int gTotalEntries;
// int gTotalEntries = 300000;

class DoStripAlign
{
public:
private:
    TTree *treeIn;
    TTree *treeOut;
};

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

struct EntryData
{
    double calcPos[gBoardCount];
    double weight1[gBoardCount], weight2[gBoardCount];
    double weight1a[gBoardCount], weight1b[gBoardCount], weight2a[gBoardCount], weight2b[gBoardCount];
    int strip1[gBoardCount];
    int strip2[gBoardCount];
    int decodeFlag[gBoardCount];
    int decodeType[gBoardCount];
} gEntryData;

// 0. Give calcPos3/calcPos4, weight, strip, for entry
void GetEntryData(int entry, double *calcPos, EntryData *data = &gEntryData)
{
    // Only calcPos determined by input
    memcpy(data->calcPos, calcPos, sizeof(double) * gBoardCount);

    corfile->GetEntry(entry);
    memcpy(data->weight1, corfile->weight1, sizeof(double) * gBoardCount);
    memcpy(data->weight2, corfile->weight2, sizeof(double) * gBoardCount);
    memcpy(data->weight1a, corfile->weight1a, sizeof(double) * gBoardCount);
    memcpy(data->weight1b, corfile->weight1b, sizeof(double) * gBoardCount);
    memcpy(data->weight2a, corfile->weight2a, sizeof(double) * gBoardCount);
    memcpy(data->weight2b, corfile->weight2b, sizeof(double) * gBoardCount);
    memcpy(data->strip1, corfile->firedStrip1, sizeof(int) * gBoardCount);
    for (int board = 0; board < gBoardCount; board++)
        if (data->decodeFlag[board] <= 1 || data->strip1[board] >= 64)
            data->strip2[board] = -1;
        else
            data->strip2[board] = data->strip1[board] + 1;
    memcpy(data->decodeFlag, corfile->decodeFlag, sizeof(int) * gBoardCount);
    memcpy(data->decodeType, corfile->decodeType, sizeof(int) * gBoardCount);
}

// 1. Calculate derived position for each board using left layers' position, Write to Sum/Strip-Align/stripAlign.root
void CalcDerivedPos()
{
    auto fileW = new TFile(Form("%s/Strip-Align/stripAlign.root", gRepoFolder.c_str()), "recreate");

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
        corfile->GetEntry(entry);
        // Do something

        // Calculate deviation for layer 1, layer 2
        double x0 = corfile->calcPos3[boardX[0]], x1 = corfile->calcPos3[boardX[1]], x2 = corfile->calcPos3[boardX[2]], x3 = corfile->calcPos3[boardX[3]];
        double y0 = corfile->calcPos3[boardY[0]], y1 = corfile->calcPos3[boardY[1]], y2 = corfile->calcPos3[boardY[2]], y3 = corfile->calcPos3[boardY[3]];
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
            calcPos3[board] = corfile->calcPos3[board];
            strip1[board] = corfile->firedStrip1[board];
            strip2[board] = strip1[board] + 1;
            weight1[board] = corfile->weight1[board];
            weight2[board] = corfile->weight2[board];
            weight1a[board] = corfile->weight1a[board];
            weight1b[board] = corfile->weight1b[board];
            weight2a[board] = corfile->weight2a[board];
            weight2b[board] = corfile->weight2b[board];

            if (corfile->decodeFlag[board] == 1)
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

void CalcDerivedPosIter(std::function<double *(int)> func = [](int entry)
                        { corfile->GetEntry(entry);return corfile->calcPos3; })
{
    auto fileW = new TFile(Form("%s/Strip-Align/stripAlign.root", gRepoFolder.c_str()), "recreate");

    TGraph *tgXTemp = new TGraph();
    TGraph *tgYTemp = new TGraph();
    TF1 *fpolX = new TF1("fpolX", "pol1", -1000, 1000);
    TF1 *fpolY = new TF1("fpolY", "pol1", -1000, 1000);

    auto tree = new TTree("align", "Strip Alignment");

    double calcPos[8];
    double derivedPos[gBoardCount];
    tree->Branch("derivedPos", derivedPos, Form("derivedPos[%d]/D", gBoardCount));

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
        // Do something
        memcpy(calcPos, func(entry), sizeof(double) * 8);

        // Calculate deviation for layer 1, layer 2
        double x0 = calcPos[boardX[0]], x1 = calcPos[boardX[1]], x2 = calcPos[boardX[2]], x3 = calcPos[boardX[3]];
        double y0 = calcPos[boardY[0]], y1 = calcPos[boardY[1]], y2 = calcPos[boardY[2]], y3 = calcPos[boardY[3]];
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

        tree->Fill();
    }
    // Write to file
    fileW->cd();
    tree->Write();
    fileW->Close();
    delete fileW;
}

// 2. Read Sum/striptAlign.root file, generate Charge-derived position histos & graphs, write into Sum/Strip-Align/saPlots.root
void PlotQX()
{
    auto fAlign = new TFile(Form("%s/Strip-Align/stripAlign.root", gRepoFolder.c_str()));
    if (fAlign->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/Strip-Align/stripAlign.root file, try to create." << std::endl;
        CalcDerivedPos();
        fAlign = new TFile(Form("%s/Strip-Align/stripAlign.root", gRepoFolder.c_str()));
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

    auto fileW = new TFile(Form("%s/Strip-Align/saPlots.root", gRepoFolder.c_str()), "recreate");

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

void PlotQXIter(std::function<double *(int)> func = [](int entry)
                { corfile->GetEntry(entry);return corfile->calcPos3; })
{
    auto fAlign = new TFile(Form("%s/Strip-Align/stripAlign.root", gRepoFolder.c_str()));
    if (fAlign->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/Strip-Align/stripAlign.root file, try to create." << std::endl;
        CalcDerivedPos();
        fAlign = new TFile(Form("%s/Strip-Align/stripAlign.root", gRepoFolder.c_str()));
    }
    auto treeAlign = (TTree *)fAlign->Get("align");
    if (!treeAlign)
    {
        std::cout << "Error: Cannot find align tree" << std::endl;
        return;
    }

    double derivedPos[gBoardCount];
    treeAlign->SetBranchAddress("derivedPos", derivedPos);

    auto fileW = new TFile(Form("%s/Strip-Align/saPlots.root", gRepoFolder.c_str()), "recreate");

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
        corfile->GetEntry(entry);
        double calcPos[gBoardCount];
        memcpy(calcPos, func(entry), sizeof(double) * gBoardCount);

        GetEntryData(entry, calcPos);

        for (int board = 0; board < gBoardCount; board++)
        {

            hsna[board][gEntryData.strip1[board]]->Fill(derivedPos[board], gEntryData.weight1a[board]);
            hsnb[board][gEntryData.strip1[board]]->Fill(derivedPos[board], gEntryData.weight1b[board]);
            hsnab[board][gEntryData.strip1[board]]->Fill(derivedPos[board], gEntryData.weight1a[board] + gEntryData.weight1b[board]);

            if (TMath::Abs(derivedPos[board] - gEntryData.strip1[board] * 11) < 15)
            {
                tgsna[board][gEntryData.strip1[board]]->SetPoint(tgsna[board][gEntryData.strip1[board]]->GetN(), derivedPos[board], gEntryData.weight1a[board]);
                tgsnb[board][gEntryData.strip1[board]]->SetPoint(tgsnb[board][gEntryData.strip1[board]]->GetN(), derivedPos[board], gEntryData.weight1b[board]);
                tgsnab[board][gEntryData.strip1[board]]->SetPoint(tgsnab[board][gEntryData.strip1[board]]->GetN(), derivedPos[board], gEntryData.weight1a[board] + gEntryData.weight1b[board]);
            }

            if (gEntryData.strip2[board] < 0 || gEntryData.strip2[board] > 63)
                continue;

            hsna[board][gEntryData.strip2[board]]->Fill(derivedPos[board], gEntryData.weight2a[board]);
            hsnb[board][gEntryData.strip2[board]]->Fill(derivedPos[board], gEntryData.weight2b[board]);
            hsnab[board][gEntryData.strip2[board]]->Fill(derivedPos[board], gEntryData.weight2a[board] + gEntryData.weight2b[board]);

            if (TMath::Abs(derivedPos[board] - gEntryData.strip2[board] * 11) < 15)
            {
                tgsna[board][gEntryData.strip2[board]]->SetPoint(tgsna[board][gEntryData.strip2[board]]->GetN(), derivedPos[board], gEntryData.weight2a[board]);
                tgsnb[board][gEntryData.strip2[board]]->SetPoint(tgsnb[board][gEntryData.strip2[board]]->GetN(), derivedPos[board], gEntryData.weight2b[board]);
                tgsnab[board][gEntryData.strip2[board]]->SetPoint(tgsnab[board][gEntryData.strip2[board]]->GetN(), derivedPos[board], gEntryData.weight2a[board] + gEntryData.weight2b[board]);
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

// 3. Fit Q-X relation to get strip alignment, write to Sum/Strip-Align/stripAlign.txt
void CalcAlignment()
{
    auto fAlign = new TFile(Form("%s/Strip-Align/saPlots.root", gRepoFolder.c_str()));
    if (fAlign->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/Strip-Align/saPlots.root file, try to create." << std::endl;
        // CalcDerivedPos();
        fAlign = new TFile(Form("%s/Strip-Align/saPlots.root", gRepoFolder.c_str()));
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
            c->SaveAs(Form("%s/Strip-Align/Board%d/stripAlign-a-%d-%d.jpg", gRepoFolder.c_str(), board, board, strip));
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
            c->SaveAs(Form("%s/Strip-Align/Board%d/stripAlign-b-%d-%d.jpg", gRepoFolder.c_str(), board, board, strip));
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
            c->SaveAs(Form("%s/Strip-Align/Board%d/stripAlign-ab-%d-%d.jpg", gRepoFolder.c_str(), board, board, strip));
            stripAlignC[board][strip] = fstrip->GetParameter(0);
            stripAmpC[board][strip] = fstrip->GetParameter(1);
        }
    }

    // Output fit results
    std::ofstream fout(Form("%s/Strip-Align/stripAlign.txt", gRepoFolder.c_str()));
    for (int board = 0; board < gBoardCount; board++)
    {
        for (int strip = 0; strip < 64; strip++)
        {
            fout << board << '\t' << strip << '\t' << stripAlignA[board][strip] << '\t' << stripAmpA[board][strip] << '\t' << stripAlignB[board][strip] << '\t' << stripAmpB[board][strip] << '\t' << stripAlignC[board][strip] << '\t' << stripAmpC[board][strip] << std::endl;
        }
    }
    fout.close();
}

// 4. Read Sum/stripAlign.txt file, calculate new position for each strip, write to Sum/Strip-Align/stripPosition.txt
void CalcStripPosition()
{
    std::ifstream fin(Form("%s/Strip-Align/stripAlign.txt", gRepoFolder.c_str()));
    if (!fin.is_open())
    {
        std::cout << "Error: Cannot open Sum/Strip-Align/stripAlign.txt file" << std::endl;
        return;
    }

    double stripAlignA[gBoardCount][64]; // Group A
    double stripAlignB[gBoardCount][64]; // Group B
    double stripAlignC[gBoardCount][64]; // Center
    double stripAmpA[gBoardCount][64];   // Group A
    double stripAmpB[gBoardCount][64];   // Group B
    double stripAmpC[gBoardCount][64];   // Center

    while (!fin.eof())
    {
        int board, strip;
        double alignA, ampA, alignB, ampB, alignC, ampC;
        fin >> board >> strip >> alignA >> ampA >> alignB >> ampB >> alignC >> ampC;
        stripAlignA[board][strip] = alignA;
        stripAmpA[board][strip] = ampA;
        stripAlignB[board][strip] = alignB;
        stripAmpB[board][strip] = ampB;
        stripAlignC[board][strip] = alignC;
        stripAmpC[board][strip] = ampC;
    }
    fin.close();

    TGraph *tg[gBoardCount];
    for (int board = 0; board < gBoardCount; board++)
    {
        tg[board] = new TGraph();
        tg[board]->SetName(Form("tg%d", board));
    }

    // Calculate strip offset
    double stripOffset[gBoardCount];
    for (int layer = 0; layer < gBoardCount / 2; layer++)
    {
        double offsetX = 0, offsetY = 0;
        int countX = 0, countY = 0;

        for (int stripx = 2; stripx < 46; stripx++)
        {
            double pitch = stripAlignC[layer * 2][stripx + 1] - stripAlignC[layer * 2][stripx];
            if (TMath::Abs(pitch - 11) > 1)
                continue;
            offsetX += stripAlignC[layer * 2][stripx] - stripx * 11;
            countX++;
        }
        offsetX /= countX;

        for (int stripy = 18; stripy < 62; stripy++)
        {
            double pitch = stripAlignC[layer * 2 + 1][stripy + 1] - stripAlignC[layer * 2 + 1][stripy];
            if (TMath::Abs(pitch - 11) > 1)
                continue;
            offsetY += stripAlignC[layer * 2 + 1][stripy] - stripy * 11;
            countY++;
        }
        offsetY /= countY;

        stripOffset[layer * 2] = offsetX;
        stripOffset[layer * 2 + 1] = offsetY;
    }

    // Calculate new position
    double stripPosition[gBoardCount][64];
    double distance[gBoardCount][64];
    for (int layer = 0; layer < gBoardCount / 2; layer++)
    {
        for (int stripx = 0; stripx < 48; stripx++)
            stripPosition[layer * 2][stripx] = stripAlignC[layer * 2][stripx] - stripOffset[layer * 2];

        for (int stripy = 16; stripy < 64; stripy++)
            stripPosition[layer * 2 + 1][stripy] = stripAlignC[layer * 2 + 1][stripy] - stripOffset[layer * 2 + 1];
    }

    std::ofstream fout(Form("%s/Strip-Align/stripPosition.txt", gRepoFolder.c_str()));
    for (int board = 0; board < gBoardCount; board++)
    {
        std::cout << "Board: " << board << '\t' << stripOffset[board] << std::endl;
        bool flagX = (board % 2 == 0);
        for (int strip = 0; strip < 64; strip++)
        {
            if (flagX && strip > 47)
                continue;
            if (!flagX && strip < 16)
                continue;

            // If offset is too large, set position as strip * 11
            double offset = TMath::Abs(stripPosition[board][strip] - strip * 11);
            if (offset > 2)
            {
                std::cout << "Board: " << board << '\t' << strip << '\t' << stripPosition[board][strip] << '\t' << strip * 11 << '\t' << offset << std::endl;
                stripPosition[board][strip] = strip * 11;
            }
            tg[board]->SetPoint(tg[board]->GetN(), strip, stripPosition[board][strip] - strip * 11);
            fout << board << '\t' << strip << '\t' << stripPosition[board][strip] << std::endl;
        }
        std::cout << std::endl;
    }
    fout.close();
}

// 5. Give Calculated hit position using new strip position
void CalcHitPosition()
{
    std::ifstream fin;

    fin.open(Form("%s/Strip-Align/stripPosition.txt", gRepoFolder.c_str()));
    if (!fin.is_open())
    {
        std::cout << "Error: Cannot open Sum/Strip-Align/stripPosition.txt file" << std::endl;
        return;
    }

    double stripPosition[gBoardCount][64];
    while (!fin.eof())
    {
        int board, strip;
        double pos;
        fin >> board >> strip >> pos;
        stripPosition[board][strip] = pos;
    }
    fin.close();

    // Create new file
    auto fileW = new TFile(Form("%s/Strip-Align/stripAlign-Verify.root", gRepoFolder.c_str()), "recreate");

    auto tree = new TTree("pos4", "calculate position 4th time");
    double calcPos3[gBoardCount];
    double calcPos4[gBoardCount];
    double derivedPos[gBoardCount];
    double fitSigma[gBoardCount];
    int decodeFlag[gBoardCount];
    int decodeType[gBoardCount];

    tree->Branch("calcPos3", calcPos3, Form("calcPos3[%d]/D", gBoardCount));
    tree->Branch("calcPos4", calcPos4, Form("calcPos4[%d]/D", gBoardCount));
    tree->Branch("derivedPos", derivedPos, Form("derivedPos[%d]/D", gBoardCount));
    tree->Branch("fitSigma", fitSigma, Form("fitSigma[%d]/D", gBoardCount));
    tree->Branch("decodeFlag", decodeFlag, Form("decodeFlag[%d]/I", gBoardCount));
    tree->Branch("decodeType", decodeType, Form("decodeType[%d]/I", gBoardCount));

    double zArrayX[4], zArrayY[4];
    for (int layer = 0; layer < 4; layer++)
    {
        zArrayX[layer] = Z_POS[layer] + layerDistance / 2.0;
        zArrayY[layer] = Z_POS[layer] - layerDistance / 2.0;
    }
    TGraph *tgXTemp = new TGraph();
    TGraph *tgYTemp = new TGraph();
    TF1 *fpolX = new TF1("fpolX", "pol1", -1000, 1000);
    TF1 *fpolY = new TF1("fpolY", "pol1", -1000, 1000);

    for (int entry = 0; entry < gTotalEntries; entry++)
    {
        if (entry % 1000 == 1000 - 1)
            std::cout << "Entry: " << entry << std::endl;
        corfile->GetEntry(entry);

        // Calculate calcPos4
        for (int board = 0; board < gBoardCount; board++)
        {
            if (corfile->decodeFlag[board] < 2 || corfile->fitSigma[board] > 10)
                continue;
            double devTemp = corfile->calcPos2[board] - corfile->firedStrip1[0] * 10;
            int strip1 = corfile->firedStrip1[board];
            if (TMath::Abs(devTemp) < 1e-3)
                calcPos4[board] = stripPosition[board][strip1];
            else
            {
                double pos1 = stripPosition[board][strip1];
                double pos2 = stripPosition[board][strip1 + 1];
                double weight1a = corfile->weight1a[board];
                double weight1b = corfile->weight1b[board];
                double weight2a = corfile->weight2a[board];
                double weight2b = corfile->weight2b[board];
                calcPos4[board] = (pos1 * (weight1a + weight1b) + pos2 * (weight2a + weight2b)) / (weight1a + weight1b + weight2a + weight2b);
            }
            calcPos4[board] += corfile->calcPos3[board] - corfile->calcPos2[board] * 1.1;
        }

        // Decode flag and type
        for (int board = 0; board < gBoardCount; board++)
        {
            decodeFlag[board] = corfile->decodeFlag[board];
            decodeType[board] = corfile->decodeType[board];
        }

        // Calculate deviation for layer 1, layer 2
        double x0 = calcPos4[boardX[0]], x1 = calcPos4[boardX[1]], x2 = calcPos4[boardX[2]], x3 = calcPos4[boardX[3]];
        double y0 = calcPos4[boardY[0]], y1 = calcPos4[boardY[1]], y2 = calcPos4[boardY[2]], y3 = calcPos4[boardY[3]];
        double z0 = Z_POS[0], z1 = Z_POS[1], z2 = Z_POS[2], z3 = Z_POS[3];
        double z0X = zArrayX[0], z1X = zArrayX[1], z2X = zArrayX[2], z3X = zArrayX[3];
        double z0Y = zArrayY[0], z1Y = zArrayY[1], z2Y = zArrayY[2], z3Y = zArrayY[3];

        // Calculate derived position
        double xArray[4] = {x0, x1, x2, x3};
        double yArray[4] = {y0, y1, y2, y3};
        double zArray[4] = {z0, z1, z2, z3};

        tgXTemp->Set(0);
        tgYTemp->Set(0);
        for (int i = 0; i < 4; i++)
        {
            tgXTemp->SetPoint(tgXTemp->GetN(), zArray[i], xArray[i]);
            tgYTemp->SetPoint(tgYTemp->GetN(), zArray[i], yArray[i]);
        }
        tgXTemp->Fit(fpolX, "Q");
        tgYTemp->Fit(fpolY, "Q");

        for (int layer = 0; layer < 4; layer++)
        {
            double xDerived = 0, yDerived = 0;
            xDerived = fpolX->Eval(zArray[layer]);
            derivedPos[boardX[layer]] = xDerived;
            fitSigma[boardX[layer]] = fpolX->GetParError(0) + fpolX->GetParError(1) * zArray[layer];

            yDerived = fpolY->Eval(zArray[layer]);
            derivedPos[boardY[layer]] = yDerived;
            fitSigma[boardY[layer]] = fpolY->GetParError(0) + fpolY->GetParError(1) * zArray[layer];
        }
        tree->Fill();
    }

    tree->Write();
    fileW->Close();
    delete fileW;
}

// 6. Verify the strip alignment, check position deviation
void VerifyAlign()
{
    auto fileR = new TFile(Form("%s/Strip-Align/stripAlign-Verify.root", gRepoFolder.c_str()));
    if (fileR->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/stripAlign-Verify.root file" << std::endl;
        return;
    }

    auto fGaus = new TF1("fGaus", "gaus", -20, 20);
    auto tree = (TTree *)fileR->Get("pos4");

    gStyle->SetOptStat(1);
    gStyle->SetOptFit(11);
    for (int layer = 0; layer < 4; layer++)
    {
        char cDraw[1024];
        char cCut[1024];
        c->cd();

        int xboardNo = 2 * layer;
        sprintf(cDraw, "(calcPos4[%d]-derivedPos[%d])>>hDiffX(1000,-10,10)", xboardNo, xboardNo);
        sprintf(cCut, "decodeFlag[%d]>1&&fitSigma[%d]<10", xboardNo, xboardNo);
        tree->Draw(cDraw, cCut);
        auto hDiffX = (TH1F *)(gFile->Get("hDiffX"));
        // hDiffX->SetTitle(Form("Board %d;#Delta_{x}/mm;Counts", xboardNo));
        hDiffX->SetTitle(Form("Layer %d X,Board %d;#Delta_{x}/mm;Counts", layer, xboardNo));
        hDiffX->Fit(fGaus, "Q");
        hDiffX->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1. * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1. * fGaus->GetParameter(2));
        hDiffX->Draw();
        fGaus->SetRange(fGaus->GetParameter(1) - 5 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 5 * fGaus->GetParameter(2));
        fGaus->Draw("same");
        c->SaveAs(Form("%s/Strip-Align/Layer%d-hDiffX.jpg", gRepoFolder.c_str(), layer));

        int yboardNo = 2 * layer + 1;
        sprintf(cDraw, "(calcPos4[%d]-derivedPos[%d])>>hDiffY(1000,-10,10)", yboardNo, yboardNo);
        sprintf(cCut, "decodeFlag[%d]>1&&fitSigma[%d]<10", yboardNo, yboardNo);
        tree->Draw(cDraw, cCut);
        auto hDiffY = (TH1F *)(gFile->Get("hDiffY"));
        // hDiffY->SetTitle("Position resolution for Y direction;Position resolution/mm;Counts");
        // hDiffY->SetTitle(Form("Board %d;#Delta_{x}/mm;Counts", yboardNo));
        hDiffY->SetTitle(Form("Layer %d Y,Board %d;#Delta_{x}/mm;Counts", layer, yboardNo));
        hDiffY->Fit(fGaus, "Q");
        hDiffY->Fit(fGaus, "Q", "", fGaus->GetParameter(1) - 1. * fGaus->GetParameter(2), fGaus->GetParameter(1) + 1. * fGaus->GetParameter(2));
        hDiffY->Draw();
        fGaus->SetRange(fGaus->GetParameter(1) - 5 * fGaus->GetParameter(2), fGaus->GetParameter(1) + 5 * fGaus->GetParameter(2));
        fGaus->Draw("same");
        c->SaveAs(Form("%s/Strip-Align/Layer%d-hDiffY.jpg", gRepoFolder.c_str(), layer));
    }

    fileR->Close();
    delete fileR;
    return;
}

// 7. Iteration
void StripAlignIter(std::string sRepoFolder = "../Sum")
{
    gRepoFolder = sRepoFolder;

    auto fileR = new TFile(Form("%s/Strip-Align/stripAlign-Verify.root", gRepoFolder.c_str()));
    if (fileR->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/stripAlign-Verify.root file" << std::endl;
        return;
    }

    auto fGaus = new TF1("fGaus", "gaus", -20, 20);
    auto tree = (TTree *)fileR->Get("pos4");
    double *calcPos4 = new double[gBoardCount];

    tree->SetBranchAddress("calcPos4", calcPos4);

    auto lambda0 = ([](int entry)
                    { corfile->GetEntry(entry);return corfile->calcPos3; });
    auto lambda1 = ([&tree, &calcPos4](int entry)
                    { tree->GetEntry(entry);return calcPos4; });

    // CalcDerivedPosIter(lambda1);
    // PlotQXIter(lambda1);
    // CalcAlignment();
    CalcStripPosition();
    // CalcHitPosition();
    // VerifyAlign();

    delete[] calcPos4;
}

void StripAlign(std::string sRepoFolder = "../Sum")
{
    gRepoFolder = sRepoFolder;

    f = new TFile(Form("%s/1stCorrection.root", gRepoFolder.c_str()));
    cor = (TTree *)f->Get("cor");
    // typedef john::cor classcor;
    corfile = new john::cor(cor);
    gTotalEntries = corfile->fChain->GetEntries();
    // Open Sum/1stCorrection.root file
    if (f->IsZombie())
    {
        std::cout << "Error: Cannot open Sum/1stCorrection.root file" << std::endl;
        return;
    }
    // Get cor tree
    if (!cor)
    {
        std::cout << "Error: Cannot find cor tree" << std::endl;
        return;
    }

    // CalcDerivedPos();
    // PlotQX();
    // CalcAlignment();
    // CalcStripPosition();
    // CalcHitPosition();
    // VerifyAlign();

    StripAlignIter(sRepoFolder);
}