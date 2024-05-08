#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>
#include <string>
#include <memory>

// #include "board.C"
#include "Calibration.cpp"

const int gBoardCount = 8;
std::string gDataFolder = "../";

TFile *gMatchFile = NULL;
TTree *gMatchTree = NULL;

double gTSInterval[gBoardCount];
double gLastSeg[gBoardCount]; // segment time into pieces, length is about 1 s
double gNextSeg[gBoardCount]; // Same use with gLastTS

int gMatchCounter;                  // restore match counter
int gMatchedBoard[256];             // restore matched channel
int gMatchFlag[gBoardCount];        // restore whether this board has matched index
ULong64_t gMatchEntry[gBoardCount]; // restore entries
double gMatchTime[gBoardCount];     // restore matched time

uint64_t gCurrentEntries[gBoardCount]{0}; // Current reading entries
double gT0Delay[gBoardCount];             // Restore T0Delay for each board, set board 0 at 0

std::vector<double> gEventTime[gBoardCount];
std::vector<int> gEventEntry[gBoardCount];

double gTDCCut[gBoardCount] = {};

bool OpenMatchFile(std::string sFileName = gDataFolder + "/" + "MatchEntries.root")
{
    gMatchFile = new TFile(sFileName.c_str());
    gMatchTree = (TTree *)(gMatchFile->Get("match"));
    if (!gMatchTree)
        return false;
    gMatchTree->SetBranchAddress("counter", &gMatchCounter);
    gMatchTree->SetBranchAddress("matchedBoard", gMatchedBoard);
    gMatchTree->SetBranchAddress("matchFlag", gMatchFlag);
    gMatchTree->SetBranchAddress("matchEntry", gMatchEntry);
    gMatchTree->SetBranchAddress("matchTime", gMatchTime);
    gMatchTree->SetBranchAddress("lastSeg", gLastSeg);
    gMatchTree->SetBranchAddress("nextSeg", gNextSeg);
    gMatchTree->SetBranchAddress("interval", gTSInterval);

    return true;
}

// Each Board information
int gBoardNo[gBoardCount];
board *gBoard[gBoardCount];
BoardCalibration *gBoardCali[gBoardCount];

void GenerateBoardMap(std::string sCaliFolder = "../", std::string sDataFolder = "../")
{
    for (int i = 0; i < gBoardCount; i++)
    {
        gBoardNo[i] = i;

        gBoardCali[i] = new BoardCalibration;
        gBoardCali[i]->DoCalibration(gBoardNo[i], sCaliFolder, sDataFolder);
        gBoardCali[i]->GetStripCaliFactor();
        gBoardCali[i]->ReadDataFile(sDataFolder);

        // gBoard[i] = new board(Form("%s/Board%d-Aligned.root", gDataFolder.c_str(), i));
        gBoard[i] = gBoardCali[i]->GetBoard();
    }
}

int gDecodeFlag[gBoardCount];
int gDecodeType[gBoardCount];
int gFiredStrip1[gBoardCount];
int gFiredStrip2[gBoardCount];
int gFiredCount[gBoardCount];
double gWeight1[gBoardCount];
double gWeight2[gBoardCount];
double gWeight1A[gBoardCount];
double gWeight2A[gBoardCount];
double gWeight1B[gBoardCount];
double gWeight2B[gBoardCount];
double gCalcPos[gBoardCount];
double gCalcPos2[gBoardCount];

double gTDCTime[gBoardCount];

double gBoardHG[gBoardCount][32];
double gBoardLG[gBoardCount][32];
double gBoardTDC[gBoardCount][33];

#include <TSystem.h>
void GetPos(std::string sCaliFolder = "../", std::string sDataFolder = "../")
{
    gDataFolder = sDataFolder;
    OpenMatchFile(Form("%s/MatchEntries.root", gDataFolder.c_str()));
    GenerateBoardMap(sCaliFolder, sDataFolder);

    auto file = new TFile(Form("%s/GetPos.root", gDataFolder.c_str()), "recreate");
    auto tree = new TTree("data", "Calculated position");
    tree->Branch("decodeFlag", gDecodeFlag, Form("decodeFlag[%d]/I", gBoardCount));
    tree->Branch("decodeType", gDecodeType, Form("decodeType[%d]/I", gBoardCount));
    tree->Branch("firedStrip1", gFiredStrip1, Form("firedStrip1[%d]/I", gBoardCount));
    tree->Branch("firedStrip2", gFiredStrip2, Form("firedStrip2[%d]/I", gBoardCount));
    tree->Branch("firedCount", gFiredCount, Form("firedCount[%d]/I", gBoardCount));
    tree->Branch("weight1", gWeight1, Form("weight1[%d]/D", gBoardCount));
    tree->Branch("weight2", gWeight2, Form("weight2[%d]/D", gBoardCount));
    tree->Branch("weight1a", gWeight1A, Form("weight1a[%d]/D", gBoardCount));
    tree->Branch("weight2a", gWeight2A, Form("weight2a[%d]/D", gBoardCount));
    tree->Branch("weight1b", gWeight1B, Form("weight1b[%d]/D", gBoardCount));
    tree->Branch("weight2b", gWeight2B, Form("weight2b[%d]/D", gBoardCount));
    tree->Branch("calcPos", gCalcPos, Form("calcPos[%d]/D", gBoardCount));
    tree->Branch("calcPos2", gCalcPos2, Form("calcPos2[%d]/D", gBoardCount));

    tree->Branch("tdcTime", gTDCTime, Form("tdcTime[%d]/D", gBoardCount));
    tree->Branch("matchEntry", gMatchEntry, Form("matchEntry[%d]/I", gBoardCount));
    tree->Branch("matchTime", gMatchTime, Form("matchTime[%d]/D", gBoardCount));
    tree->Branch("lastSeg", gLastSeg, Form("lastSeg[%d]/D", gBoardCount));
    tree->Branch("nextSeg", gNextSeg, Form("nextSeg[%d]/D", gBoardCount));
    tree->Branch("interval", gTSInterval, Form("interval[%d]/D", gBoardCount));
    tree->Branch("matchFlag", gMatchFlag, Form("matchFlag[%d]/I", gBoardCount));

    tree->Branch("boardHG", gBoardHG, Form("boardHG[%d][32]/D", gBoardCount));
    tree->Branch("boardLG", gBoardLG, Form("boardLG[%d][32]/D", gBoardCount));
    tree->Branch("boardTDC", gBoardTDC, Form("boardTDC[%d][33]/D", gBoardCount));

#ifdef SAVERAW
    char cbuf1[256], cbuf2[256];
    for (int board = 0; board < gBoardCount; board++)
    {
        sprintf(cbuf1, "board%dhg", board);
        sprintf(cbuf2, "board%dhg[32]/D", board);
        tree->Branch(cbuf1, gBoard[board]->HGamp, cbuf2);
        sprintf(cbuf1, "board%dlg", board);
        sprintf(cbuf2, "board%dlg[32]/D", board);
        tree->Branch(cbuf1, gBoard[board]->LGamp, cbuf2);
        sprintf(cbuf1, "board%dtdc", board);
        sprintf(cbuf2, "board%dtdc[33]/D", board);
        tree->Branch(cbuf1, gBoard[board]->TDCTime, cbuf2);
        sprintf(cbuf1, "board%dfiredCount", board);
        sprintf(cbuf2, "board%dfiredCount/b", board);
        tree->Branch(cbuf1, &gBoard[board]->FiredCount, cbuf2);
        sprintf(cbuf1, "board%dfiredCh", board);
        sprintf(cbuf2, "board%dfiredCh[board%dfiredCount]/b", board, board);
        // std::cout << cbuf1 << '\t' << cbuf2 << std::endl;
        tree->Branch(cbuf1, gBoard[board]->FiredCh, cbuf2);
    }
#endif

    auto debugTree = new TTree("debug", "debug for failure decoding");
    int dBoardNo;
    int dDecodeFlag;
    int dnChannel = 0;
    int dFiredChannel[8];
    double dFiredAmp[8];
    debugTree->Branch("boardNo", &dBoardNo, "boardNo/I");
    debugTree->Branch("decodeFlag", &dDecodeFlag, "decodeFlag/I");
    debugTree->Branch("nChannel", &dnChannel, "nChannel/I");
    debugTree->Branch("firedCh", dFiredChannel, "firedCh[nChannel]/I");
    debugTree->Branch("firedAmp", dFiredAmp, "firedAmp/D");

    for (int entry = 0; entry < gMatchTree->GetEntries(); entry++)
    // for (int entry = 0; entry < 100; entry++)
    {
        // if (entry > 300)
        //     break;
        if (entry % 1000 == 0)
            std::cout << entry << std::endl;
        gMatchTree->GetEntry(entry);
        for (int board = 0; board < gBoardCount; board++)
        {
            gBoard[board]->GetEntry(gMatchEntry[board]);
            // copy data from board gBoard to the tree
            memcpy(gBoardHG[board], gBoard[board]->HGamp, sizeof(double) * 32);
            memcpy(gBoardLG[board], gBoard[board]->LGamp, sizeof(double) * 32);
            memcpy(gBoardTDC[board], gBoard[board]->TDCTime, sizeof(double) * 33);

            auto firedCount = gFiredCount[board] = gBoard[board]->FiredCount;
            gTDCTime[board] = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            // if (firedCount > gBoardCount)
            //     continue;

            std::vector<std::pair<int, double>> firedMap0;
            std::vector<std::pair<int, double>> firedMap;
            // auto &firedMap = firedMap0;
            for (int i = 0; i < firedCount; i++)
            {
                int sipm = gBoard[board]->FiredCh[i];
                double hgEffValue = gBoardCali[board]->CalcSiPMEffValue(sipm);
                if (hgEffValue > 0.05)
                    firedMap0.push_back(std::pair<int, double>(sipm, hgEffValue));
            }
            std::sort(firedMap0.begin(), firedMap0.end(), [](std::pair<int, double> a, std::pair<int, double> b)
                      { return a.second > b.second; });

            if (firedMap0.size() > 0)
                gTDCTime[board] = gBoard[board]->TDCTime[firedMap0[0].first];

            for (int i = 0; i < firedMap0.size(); i++)
            {
                auto sipm = firedMap0[i].first;
                auto hgEffValue = firedMap0[i].second;
                if (TMath::Abs(gBoard[board]->TDCTime[sipm] - gTDCTime[board]) < 30)
                    firedMap.push_back(std::pair<int, double>(sipm, hgEffValue));
            }

            // std::cout << "FiredCount: " << firedCount << "\t filtered count:" << firedMap.size() << std::endl;

            // Decode for the first try
            if (firedMap.size() < 4 && firedMap.size() >= 2)
            {
                gDecodeType[board] = 1;
                double weight;
                gDecodeFlag[board] = UserDefine::Decode2SiPM(firedMap[0].first, firedMap[1].first, gFiredStrip1[board], weight, gBoard[board]->HGamp);
                if (gDecodeFlag[board] < 0)
                {
                    std::cout << "Error in LOOP2: " << entry << '\t' << board << '\t' << firedCount << '\t' << (int)gDecodeFlag[board] << '\t' << firedMap[0].first << '\t' << firedMap[0].second << '\t' << firedMap[1].first << '\t' << gFiredStrip1[board] / 16 << '\t' << gWeight1[board] << std::endl
                              << std::endl;
                    dBoardNo = board;
                    dDecodeFlag = 1;
                    dnChannel = 2;
                    for (int idx = 0; idx < dnChannel; idx++)
                    {
                        dFiredChannel[idx] = firedMap[idx].first;
                        dFiredAmp[idx] = firedMap[idx].second;
                    }
                    debugTree->Fill();
                }

                gBoardCali[board]->CalcStripEffValue(gFiredStrip1[board], gWeight1A[board], gWeight1B[board]);
                gFiredStrip2[board] = -1;
                gWeight2[board] = gWeight2A[board] = gWeight2B[board] = 0;
                gCalcPos[board] = gFiredStrip1[board] * 10;
                gCalcPos2[board] = gFiredStrip1[board] * 10;
            }
            else if (firedMap.size() >= 4)
            {
                int strip1, strip2;
                gDecodeFlag[board] = UserDefine::Decode4SiPM(firedMap[0].first, firedMap[1].first, firedMap[2].first, firedMap[3].first, gBoard[board]->HGamp, gWeight1[board], gWeight2[board], strip1, strip2);
                if (gDecodeFlag[board] > 0)
                {
                    gFiredStrip1[board] = strip1;
                    gFiredStrip2[board] = strip2;
                    gBoardCali[board]->CalcStripEffValue(gFiredStrip1[board], gWeight1A[board], gWeight1B[board]);
                    gBoardCali[board]->CalcStripEffValue(gFiredStrip2[board], gWeight2A[board], gWeight2B[board]);
                    gDecodeType[board] = 2;

                    gCalcPos[board] = (gFiredStrip1[board] * gWeight1[board] + gFiredStrip2[board] * gWeight2[board]) / (gWeight1[board] + gWeight2[board]);
                    gCalcPos[board] *= 10;
                    gCalcPos2[board] = (gFiredStrip1[board] * (gWeight1A[board] + gWeight1B[board]) + gFiredStrip2[board] * (gWeight2A[board] + gWeight2B[board])) / (gWeight1A[board] + gWeight1B[board] + gWeight2A[board] + gWeight2B[board]);
                    gCalcPos2[board] *= 10;
                }

                if (gDecodeFlag[board] < 0)
                {
                    std::cout << std::setprecision(12) << "Error in LOOP4: " << entry << '\t' << board << '\t' << firedCount << '\t' << (int)gDecodeFlag[board] << '\t' << std::endl;
                    for (int i = 0; i < firedMap.size(); i++)
                        std::cout << std::setprecision(12) << firedMap[i].first << '\t' << firedMap[i].second << '\t' << gBoard[board]->TDCTime[firedMap[i].first] - gTDCTime[board] << std::endl;
                    std::cout << std::endl
                              << std::setprecision(12)
                              << gFiredStrip1[board] / 16 << '\t' << gFiredStrip2[board] / 16 << '\t' << gWeight1[board] << '\t' << gWeight2[board] << std::endl
                              << std::endl;
                    dBoardNo = board;
                    dDecodeFlag = 2;
                    dnChannel = 4;
                    for (int idx = 0; idx < dnChannel; idx++)
                    {
                        dFiredChannel[idx] = firedMap[idx].first;
                        dFiredAmp[idx] = firedMap[idx].second;
                    }
                    debugTree->Fill();

                    // Decode using 2 SiPM
                    double weight;
                    int decodeFlag2;
                    decodeFlag2 = UserDefine::Decode2SiPM(firedMap[0].first, firedMap[1].first, gFiredStrip1[board], weight, gBoard[board]->HGamp);
                    if (decodeFlag2 < 0)
                    {
                        std::cout << "Try to Decode in LOOP2 error: " << entry << '\t' << board << '\t' << firedCount << '\t' << (int)decodeFlag2 << '\t' << firedMap[0].first << '\t' << firedMap[0].second << '\t' << firedMap[1].first << '\t' << gFiredStrip1[board] / 16 << '\t' << gWeight1[board] << std::endl
                                  << std::endl;
                        gDecodeFlag[board] = -3;
                        gDecodeType[board] = 3;
                    }
                    else
                    {
                        gBoardCali[board]->CalcStripEffValue(gFiredStrip1[board], gWeight1A[board], gWeight1B[board]);
                        gFiredStrip2[board] = -1;
                        gWeight2[board] = gWeight2A[board] = gWeight2B[board] = 0;
                        gCalcPos[board] = gFiredStrip1[board] * 10;
                        gCalcPos2[board] = gFiredStrip1[board] * 10;
                        gDecodeFlag[board] = -2;
                    }
                }

                if ((gFiredStrip2[board] > 0) && (gFiredStrip1[board] / 16 != gFiredStrip2[board] / 16) && (gFiredStrip1[board] != gFiredStrip2[board] - 1))
                    std::cout << "Group Error: " << gFiredStrip1[board] / 16 << '\t' << gFiredStrip2[board] / 16 << '\t' << gWeight1[board] << '\t' << gWeight2[board] << std::endl;
            }

            // if only find one fired strip
            if (gDecodeFlag[board] == 1 || gDecodeFlag[board] == -2)
            {
                int firedStrip = gFiredStrip1[board];
                // sipms for the left strip, right strip
                int sipmla, sipmlb, sipmra, sipmrb;
                double ampla, amplb, ampra, amprb;
                if (firedStrip > 0)
                {
                    sipmla = UserDefine::GetFiberCode(firedStrip - 1, UserDefine::a);
                    ampla = gBoardCali[board]->CalcSiPMEffValue(sipmla);
                    sipmlb = UserDefine::GetFiberCode(firedStrip - 1, UserDefine::b);
                    amplb = gBoardCali[board]->CalcSiPMEffValue(sipmlb);
                }
                else
                {
                    sipmla = sipmlb = -1;
                    ampla = amplb = 0;
                }

                if (firedStrip < 63)
                {
                    sipmra = UserDefine::GetFiberCode(firedStrip + 1, UserDefine::a);
                    ampra = gBoardCali[board]->CalcSiPMEffValue(sipmra);
                    sipmrb = UserDefine::GetFiberCode(firedStrip + 1, UserDefine::b);
                    amprb = gBoardCali[board]->CalcSiPMEffValue(sipmrb);
                }
                else
                {
                    sipmra = sipmrb = -1;
                    ampra = amprb = 0;
                }

                bool leftValidFlag = (ampla + amplb) > 0.05 && (ampla + amplb) > (ampra + amprb);
                bool rightValidFlag = (ampra + amprb) > 0.05 && (ampra + amprb) > (ampla + amplb);

                if (leftValidFlag)
                {
                    gDecodeType[board] = 4;
                    gFiredStrip2[board] = firedStrip - 1;
                    gWeight2[board] = gBoard[board]->HGamp[sipmla] + gBoard[board]->HGamp[sipmlb];
                    gWeight2A[board] = ampla;
                    gWeight2B[board] = amplb;

                    gCalcPos[board] = (firedStrip * gWeight1[board] + (firedStrip - 1) * gWeight2[board]) / (gWeight1[board] + gWeight2[board]);
                    gCalcPos[board] *= 10;
                    gCalcPos2[board] = (firedStrip * (gWeight1A[board] + gWeight1B[board]) + (firedStrip - 1) * (gWeight2A[board] + gWeight2B[board])) / (gWeight1A[board] + gWeight1B[board] + gWeight2A[board] + gWeight2B[board]);
                    gCalcPos2[board] *= 10;
                }
                else if (rightValidFlag)
                {
                    gDecodeType[board] = 5;
                    gFiredStrip2[board] = firedStrip + 1;
                    gWeight2[board] = gBoard[board]->HGamp[sipmra] + gBoard[board]->HGamp[sipmrb];
                    gWeight2A[board] = ampra;
                    gWeight2B[board] = amprb;

                    gCalcPos[board] = (firedStrip * gWeight1[board] + (firedStrip + 1) * gWeight2[board]) / (gWeight1[board] + gWeight2[board]);
                    gCalcPos[board] *= 10;
                    gCalcPos2[board] = (firedStrip * (gWeight1A[board] + gWeight1B[board]) + (firedStrip + 1) * (gWeight2A[board] + gWeight2B[board])) / (gWeight1A[board] + gWeight1B[board] + gWeight2A[board] + gWeight2B[board]);
                    gCalcPos2[board] *= 10;
                }
                else
                {
                    gDecodeType[board] = 6;
                }
            }
        }

        tree->Fill();
    }
    tree->Write();
    debugTree->Write();
    file->Close();
}