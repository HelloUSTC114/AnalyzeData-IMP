#include <vector>
#include <map>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

#include "board.C"
#include "tsTree.C"
std::string gDataFolder = "../";
uint64_t gTSReadingEntry = 0;
double gCurrentTS[6];
double gTSInterval[6];

// Matched information
double gLastSeg[6]; // segment time into pieces, length is about 1 s
double gNextSeg[6]; // Same use with gLastTS

tsTree *gTS = NULL;

bool JudgeBoardFlag(bool *flags, int boardNumber)
{
    for (int i = 0; i < boardNumber; i++)
    {
        if (!flags[i])
            return false;
    }
    return true;
}

bool FindStamp(int startEntry, double *stampArray, int *devFromStart)
{
    gTS->fChain->GetEntry(startEntry);

    bool findFlag[6]{0};

    for (int board = 0; board < 6; board++)
    {
        findFlag[board] = 0;
        devFromStart[board] = 0;
    }
    for (int entry = startEntry; !JudgeBoardFlag(findFlag, 6); entry++)
    {
        if (entry >= gTS->fChain->GetEntries())
            return false;
        gTS->GetEntry(entry);
        for (int board = 0; board < 6; board++)
        {
            if (findFlag[board])
                continue;

            if (gTS->tsFlag[board])
            {
                findFlag[board] = true;
                stampArray[board] = gTS->ts[board];
            }
            else
            {
                devFromStart[board]++;
            }
        }
    }
    gTS->fChain->GetEntry(startEntry);
    return true;
}

bool UpdateInterval(int startEntry)
{
    double tsPre[6];
    int devCount[6];
    bool tsIntervalFlag[6]{0};

    for (int entry = startEntry; !JudgeBoardFlag(tsIntervalFlag, 6); entry++)
    {
        if (entry >= gTS->fChain->GetEntries())
            return false;
        gTS->GetEntry(entry);
#ifdef DEBUG_INTERVAL
        std::cout << entry << '\t';
        std::cout << "Calculating first time stamp interval" << std::endl;
#endif
        for (int board = 0; board < 6; board++)
        {
            if (gTS->tsFlag[board])
            {
                double tsDev = gTS->ts[board] - tsPre[board];
                tsDev /= devCount[board] + 1;
                if (TMath::Abs(tsDev - 1e9) < 0.01e9)
                {
                    gTSInterval[board] = tsDev;
                    tsIntervalFlag[board] = 1;
                }
#ifdef DEBUG_INTERVAL
                std::cout << gTS->ts[board] / 1e9 << '\t' << tsDev / 1e9 - 1 << '\t';
#endif
                devCount[board] = 0;
                tsPre[board] = gTS->ts[board];
            }
            else
            {
                devCount[board]++;
#ifdef DEBUG_INTERVAL
                std::cout << "Void \t"
                          << "Void \t";
#endif
            }
        }
#ifdef DEBUG_INTERVAL
        std::cout
            << std::endl;
#endif
    }
    return true;
}

void InitiateTS()
{
    gTS = new tsTree(Form("%s/TS.root", gDataFolder.c_str()));
    UpdateInterval(0);

    gTS->GetEntry(0);
    std::cout << "calculated initiated time interval: " << std::endl;
    gTSReadingEntry = 0;
    for (int board = 0; board < 6; board++)
    {
        gLastSeg[board] = gTS->ts[board];
        gCurrentTS[board] = gTS->ts[board];
        std::cout << "Board index: " << board << '\t' << "First Time stamp(s): " << '\t' << gCurrentTS[board] / 1e9 << '\t' << " Time Interval(s): " << gTSInterval[board] / 1e9 << '\t' << std::endl;

        int startIdx = 30;
        int endIdx = gTS->fChain->GetEntries() - 100;
        // int endIdx = 100;
        // int endIdx = 1020;
        double ts1, ts2;
        for (; startIdx < gTS->fChain->GetEntries(); startIdx++)
        {
            gTS->fChain->GetEntry(startIdx);
            if (gTS->tsFlag[board])
                break;
        }
        ts1 = gTS->ts[board];

        for (; endIdx > 0; endIdx--)
        {
            gTS->fChain->GetEntry(endIdx);
            if (gTS->tsFlag[board])
                break;
        }
        ts2 = gTS->ts[board];
        int clockCount = 0;
        double tsPre, tsNow;
        gTS->GetEntry(startIdx);
        tsPre = gTS->ts[board];
        for (int idx = startIdx + 1; idx <= endIdx; idx++)
        {
            gTS->GetEntry(idx);
            tsNow = gTS->ts[board];
            double dev = tsNow - tsPre;
            clockCount += round(dev / 1e9);
            // std::cout << dev << '\t' << round(dev / 1e9) << std::endl;
            tsPre = tsNow;
        }
    }
}

#include <TMath.h>
bool UpdateSegPoint(int startEntry)
{
    double stamp[6];
    int dev[6];
    bool flag;
    flag = FindStamp(startEntry, stamp, dev);
    if (!flag)
        return false;
    for (int board = 0; board < 6; board++)
    {
        gLastSeg[board] = stamp[board] - dev[board] * gTSInterval[board];
    }
    flag = FindStamp(startEntry + 1, stamp, dev);
    if (flag)
        for (int board = 0; board < 6; board++)
            gNextSeg[board] = stamp[board] - dev[board] * gTSInterval[board];
    else
        for (int board = 0; board < 6; board++)
            gNextSeg[board] = gLastSeg[board] + gTSInterval[board];

    // std::cout << "Segment: " << startEntry << '\t';
    // for (int board = 0; board < 6; board++)
    //     std::cout << gLastSeg[board] / 1e9 << '\t' << gNextSeg[board] / 1e9 << '\t';
    // std::cout << std::endl;
    return true;
}

bool UpdateSeg(int startEntry)
{
    return (UpdateInterval(startEntry) && UpdateSegPoint(startEntry));
}

// Each Board information
int gBoardNo[6];
const int gBoardCount = 6;
board *gBoard[6];

void GenerateBoardMap()
{
    for (int i = 0; i < 6; i++)
    {
        gBoardNo[i] = i;
        gBoard[i] = new board(Form("%s/Board%d-Aligned.root", gDataFolder.c_str(), i));
    }
}

TFile *gMatchFile = NULL;
TTree *gMatchTree = NULL;

int gMatchCounter;        // restore match counter
int gMatchedBoard[256];   // restore matched channel
int gMatchFlag[6];        // restore whether this board has matched index
ULong64_t gMatchEntry[6]; // restore entries
double gMatchTime[6];     // restore matched time

uint64_t gCurrentEntries[6]{0}; // Current reading entries
double gT0Delay[6];             // Restore T0Delay for each board, set board 0 at 0

std::vector<double> gEventTime[6];
std::vector<int> gEventEntry[6];

void InitMatchFile()
{
    gMatchFile = new TFile(Form("%s/MatchEntries.root", gDataFolder.c_str()), "recreate");
    gMatchTree = new TTree("match", "matched entries");
    gMatchTree->Branch("counter", &gMatchCounter, "counter/I");
    gMatchTree->Branch("matchedBoard", gMatchedBoard, "matchedBoard[counter]/I");
    gMatchTree->Branch("matchFlag", gMatchFlag, "matchFlag[6]/I");
    gMatchTree->Branch("matchEntry", gMatchEntry, "matchEntry[6]/l");
    gMatchTree->Branch("matchTime", gMatchTime, "matchTime[6]/D");
    gMatchTree->Branch("lastSeg", gLastSeg, "lastSeg[6]/D");
    gMatchTree->Branch("nextSeg", gNextSeg, "nextSeg[6]/D");
    gMatchTree->Branch("interval", gTSInterval, "interval[6]/D");
    gMatchTree->AutoSave();
}

int JudgeEventTime(double eventTime, int boardNo)
{
    if (eventTime < gLastSeg[boardNo])
        return -1;
    else if (eventTime < gNextSeg[boardNo])
        return 0;
    else
        return 1;
};
bool FindAllEventsInSeg(int segEntry, int *startEntries, int *endEntries)
{
    if (!UpdateSeg(segEntry))
        return false;
    for (int board = 0; board < 6; board++)
    {
        int idx = gCurrentEntries[board];
        gBoard[board]->GetEntry(idx);
        auto eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
        if (JudgeEventTime(eventTime, board) < 0)
        {
            for (; JudgeEventTime(eventTime, board) < 0; idx++)
            {
                if (idx < 0 || idx >= gBoard[board]->fChain->GetEntries())
                    break;
                gBoard[board]->GetEntry(idx);
                eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            }
            startEntries[board] = idx;
            for (; JudgeEventTime(eventTime, board) == 0; idx++)
            {
                if (idx < 0 || idx >= gBoard[board]->fChain->GetEntries())
                    break;
                gBoard[board]->GetEntry(idx);
                eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            }
            endEntries[board] = idx - 1;
        }
        else if (JudgeEventTime(eventTime, board) > 0)
        {
            for (; JudgeEventTime(eventTime, board) <= 0; idx--)
            {
                if (idx < 0 || idx >= gBoard[board]->fChain->GetEntries())
                    break;
                gBoard[board]->GetEntry(idx);
                eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            }
            endEntries[board] = idx;
            for (; JudgeEventTime(eventTime, board) == 0; idx--)
            {
                if (idx < 0 || idx >= gBoard[board]->fChain->GetEntries())
                    break;
                gBoard[board]->GetEntry(idx);
                eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            }
            startEntries[board] = idx + 1;
        }
        else if (JudgeEventTime(eventTime, board) == 0)
        {
            for (idx = gCurrentEntries[board]; JudgeEventTime(eventTime, board) == 0; idx--)
            {
                if (idx < 0 || idx >= gBoard[board]->fChain->GetEntries())
                    break;
                gBoard[board]->GetEntry(idx);
                eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            }
            startEntries[board] = idx + 1;
            for (idx = gCurrentEntries[board]; JudgeEventTime(eventTime, board) == 0; idx++)
            {
                if (idx < 0 || idx >= gBoard[board]->fChain->GetEntries())
                    break;
                gBoard[board]->GetEntry(idx);
                eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
            }
            startEntries[board] = idx - 1;
        }

        gCurrentEntries[board] = startEntries[board];
    }
    // Match: find how many
    // for (int)
    return true;
}

bool MatchEventsInSeg(int *startEntries, int *endEntries)
{
    int eventsCount[6];
    int judgeCounter[6]; // how many judge left
    int idxCount[6];     // current reading time
    int maxEventsCount = 0;
    int maxBoard = 0;
    for (int board = 0; board < 6; board++)
    {
        eventsCount[board] = endEntries[board] - startEntries[board] + 1;
        if (eventsCount[board] > maxEventsCount)
        {
            maxEventsCount = eventsCount[board];
            maxBoard = board;
        }

        gEventEntry[board].clear();
        gEventTime[board].clear();
        for (int entry = startEntries[board]; entry <= endEntries[board]; entry++)
        {
            gBoard[board]->GetEntry(entry);
            gEventEntry[board].push_back(entry);
            gEventTime[board].push_back(gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]]);
        }
    }
    for (int board = 0; board < 6; board++)
    {
        judgeCounter[board] = maxEventsCount - eventsCount[board];
        if (judgeCounter[board] < 0)
            judgeCounter[board] = 0;
        idxCount[board] = 0;
    }

    for (int idx0 = 0; idx0 < maxEventsCount; idx0++)
    {
        gMatchCounter = 0;
        for (int board = 0; board < 6; board++)
        {
            int &idx = idxCount[board];
            double maxBoardTime = gEventTime[maxBoard][idx];
            bool skipFlag = 0;
            double timeNow;
            double ctimeNow;
            double cmaxBoardTime;
            if (idx >= gEventEntry[board].size())
            {
                skipFlag = 1;
            }
            else
            {
                timeNow = gEventTime[board][idx];
                ctimeNow = (timeNow - gLastSeg[board]) / gTSInterval[board];
                cmaxBoardTime = (maxBoardTime - gLastSeg[maxBoard]) / gTSInterval[maxBoard];
                skipFlag = !(TMath::Abs(ctimeNow - cmaxBoardTime) < 1000e3);
            }
            // Judge range: 2000 ns
            // skipFlag = !(TMath::Abs((ctimeNow - cmaxBoardTime) * 1e9) < 2000);
            if (judgeCounter[board] > 0)
            {
                if (skipFlag)
                {
                    gMatchFlag[board] = 0;
                    judgeCounter[board]--;
                }
                else
                {
                    gMatchedBoard[gMatchCounter++] = board;
                    gMatchEntry[board] = gEventEntry[board][idx];
                    gMatchFlag[board] = 1;
                    gMatchTime[board] = gEventTime[board][idx];
                    idx++;
                }
            }
            else
            {
                gMatchedBoard[gMatchCounter++] = board;
                gMatchEntry[board] = gEventEntry[board][idx];
                gMatchFlag[board] = !skipFlag;
                gMatchTime[board] = gEventTime[board][idx];
                idx++;
            }
        }
        if (gMatchTree)
        {
            gMatchTree->Fill();
        }
    }
    return true;
}

void GenerateMatchKeyTree()
{
    InitiateTS();
    GenerateBoardMap();
    InitMatchFile();
    // UpdateSeg(1000);
    int startEntries[6], endEntries[6];
    double preTime[6];
    for (int entry = 0; entry < 1000000; entry++)
    // for (int entry = 10; entry < 11; entry++)
    {
        auto flag = FindAllEventsInSeg(entry, startEntries, endEntries);
        if (!flag)
            break;
        if (entry % 100 == 0)
            std::cout << entry << std::endl;
        MatchEventsInSeg(startEntries, endEntries);
        // for (int board = 0; board < 6; board++)
        // {
        //     // std::cout << (endEntries[board] - startEntries[board] + 1) << '\t';
        //     // std::cout << startEntries[board] << '-' << endEntries[board] << '\t' << (endEntries[board] - startEntries[board] + 1) << '\t';
        //     continue;
        //     std::cout << board << std::endl;
        //     std::cout << gLastSeg[board] / 1e9 << '\t' << gNextSeg[board] / 1e9 << '\t' << gTSInterval[board] / 1e9 << '\t' << (gNextSeg[board] - gLastSeg[board]) / (gTSInterval[board]) - 1 << std::endl;
        //     std::cout << startEntries[board] << '\t' << endEntries[board] << std::endl;
        //     // for (int entries = startEntries[board]; entries <= endEntries[board]; entries++)
        //     // {
        //     //     gBoard[board]->GetEntry(entries);
        //     //     auto eventTime = gBoard[board]->TDCTime[gBoard[board]->FiredCh[0]];
        //     //     // std::cout << eventTime / 1e9 << '\t';
        //     //     // std::cout << (eventTime - gLastSeg[board]) / gTSInterval[board] << '\t';
        //     //     double correctedTime = (eventTime - gLastSeg[board]) / gTSInterval[board] * 1e9;
        //     //     std::cout << Form("%.1f", correctedTime - preTime[board]) << '\t';

        //     //     preTime[board] = correctedTime;
        //     // }
        //     // std::cout << std::endl;
        //     // std::cout << std::endl;
        // }
        // std::cout << std::endl;
    }
    gMatchTree->Write();
    gMatchFile->Close();
}

void MatchBoards()
{
    GenerateMatchKeyTree();
    // InitiateTS();
}