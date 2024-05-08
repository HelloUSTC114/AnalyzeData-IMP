#include "Decode.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>

int UserDefine::GetFiberCode(int stripNo, SiPMAB ab)
{
    if (stripNo < 0 || stripNo > 63)
        return -1;
    int stripGroup = stripNo / 16;
    int sipmGroup = 3 - stripGroup;
    int sipmNo = sipmGroup * 8;
    int stripNoInside = stripNo % 16;
    if (stripNoInside == 0)
        if (ab == a)
            sipmNo += 6;
        else
            sipmNo += 2;
    else if (stripNoInside == 1)
        if (ab == a)
            sipmNo += 7;
        else
            sipmNo += 3;
    else if (stripNoInside == 2)
        if (ab == a)
            sipmNo += 4;
        else
            sipmNo += 0;
    else if (stripNoInside == 3)
        if (ab == a)
            sipmNo += 5;
        else
            sipmNo += 1;
    else if (stripNoInside == 4)
        if (ab == a)
            sipmNo += 6;
        else
            sipmNo += 3;
    else if (stripNoInside == 5)
        if (ab == a)
            sipmNo += 7;
        else
            sipmNo += 0;
    else if (stripNoInside == 6)
        if (ab == a)
            sipmNo += 5;
        else
            sipmNo += 2;
    else if (stripNoInside == 7)
        if (ab == a)
            sipmNo += 4;
        else
            sipmNo += 3;
    else if (stripNoInside == 8)
        if (ab == a)
            sipmNo += 6;
        else
            sipmNo += 0;
    else if (stripNoInside == 9)
        if (ab == a)
            sipmNo += 7;
        else
            sipmNo += 1;
    else if (stripNoInside == 10)
        if (ab == a)
            sipmNo += 4;
        else
            sipmNo += 2;
    else if (stripNoInside == 11)
        if (ab == a)
            sipmNo += 5;
        else
            sipmNo += 0;
    else if (stripNoInside == 12)
        if (ab == a)
            sipmNo += 6;
        else
            sipmNo += 1;
    else if (stripNoInside == 13)
        if (ab == a)
            sipmNo += 7;
        else
            sipmNo += 2;
    else if (stripNoInside == 14)
        if (ab == a)
            sipmNo += 5;
        else
            sipmNo += 3;
    else if (stripNoInside == 15)
        if (ab == a)
            sipmNo += 4;
        else
            sipmNo += 1;
    else
        sipmNo = -1;
    return sipmNo;
}

int UserDefine::Decode2SiPM(int ch1NoSort, int ch2NoSort, int &firedStrip1, double &weight1, double *amp)
{
    int sipmGroup = ch1NoSort / 8;
    int ch1 = ((ch1NoSort % 8) < (ch2NoSort % 8) ? (ch1NoSort % 8) : (ch2NoSort % 8)) % 8;
    int ch2 = ((ch1NoSort % 8) > (ch2NoSort % 8) ? (ch1NoSort % 8) : (ch2NoSort % 8)) % 8;

    if (ch1 == 0 && ch2 == 4)
    {
        firedStrip1 = 2;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 4];
    }
    else if (ch1 == 0 && ch2 == 5)
    {
        firedStrip1 = 11;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 5];
    }
    else if (ch1 == 0 && ch2 == 6)
    {
        firedStrip1 = 8;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 6];
    }
    else if (ch1 == 0 && ch2 == 7)
    {
        firedStrip1 = 5;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 7];
    }
    else if (ch1 == 1 && ch2 == 4)
    {
        firedStrip1 = 15;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 4];
    }
    else if (ch1 == 1 && ch2 == 5)
    {
        firedStrip1 = 3;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 5];
    }
    else if (ch1 == 1 && ch2 == 6)
    {
        firedStrip1 = 12;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 6];
    }
    else if (ch1 == 1 && ch2 == 7)
    {
        firedStrip1 = 9;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 7];
    }
    else if (ch1 == 2 && ch2 == 4)
    {
        firedStrip1 = 10;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 4];
    }
    else if (ch1 == 2 && ch2 == 5)
    {
        firedStrip1 = 6;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 5];
    }
    else if (ch1 == 2 && ch2 == 6)
    {
        firedStrip1 = 0;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 6];
    }
    else if (ch1 == 2 && ch2 == 7)
    {
        firedStrip1 = 13;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 7];
    }
    else if (ch1 == 3 && ch2 == 4)
    {
        firedStrip1 = 7;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 4];
    }
    else if (ch1 == 3 && ch2 == 5)
    {
        firedStrip1 = 14;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 5];
    }
    else if (ch1 == 3 && ch2 == 6)
    {
        firedStrip1 = 4;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 6];
    }
    else if (ch1 == 3 && ch2 == 7)
    {
        firedStrip1 = 1;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 7];
    }
    else
    {
        int ch1Sort = (ch1NoSort > ch2NoSort ? ch2NoSort : ch1NoSort);
        int ch2Sort = (ch1NoSort > ch2NoSort ? ch1NoSort : ch2NoSort);
        std::cout << "Decoding Error in Algorithm2: " << ch1Sort << '\t' << amp[ch1Sort] << '\t' << ch2Sort << '\t' << amp[ch2Sort] << std::endl;
        return -1;
    }
    int stripGroup = 3 - sipmGroup;
    firedStrip1 += 16 * stripGroup;
    return 1;
}

int UserDefine::Decode4SiPM(int ch1NoSort, int ch2NoSort, int ch3NoSort, int ch4NoSort, double *amp, double &weight1, double &weight2, int &firedStrip1, int &firedStrip2)
{
    int sipmGroup = ch1NoSort / 8;
    int chIdx[4];
    chIdx[0] = ch1NoSort;
    chIdx[1] = ch2NoSort;
    chIdx[2] = ch3NoSort;
    chIdx[3] = ch4NoSort;

    std::sort(chIdx, chIdx + 4, [](int a, int b)
              { return a < b; });
    int ch1 = chIdx[0];
    int ch2 = chIdx[1];
    int ch3 = chIdx[2];
    int ch4 = chIdx[3];
    // std::cout << "Input: " << ch1 << '\t' << ch2 << '\t' << ch3 << '\t' << ch4 << std::endl;

    ch1 = ch1 % 8;
    ch2 = ch2 % 8;
    ch3 = ch3 % 8;
    ch4 = ch4 % 8;

    // Abnormal case: boarder of strips
    if (chIdx[0] == 10 && chIdx[1] == 14 && chIdx[2] == 17 && chIdx[3] == 20)
    {
        Decode2SiPM(10, 14, firedStrip2, weight2, amp);
        Decode2SiPM(17, 20, firedStrip1, weight1, amp);
        return 2;
    }
    if (chIdx[0] == 2 && chIdx[1] == 6 && chIdx[2] == 9 && chIdx[3] == 12)
    {
        Decode2SiPM(2, 6, firedStrip2, weight2, amp);
        Decode2SiPM(9, 12, firedStrip1, weight1, amp);
        return 2;
    }
    if (chIdx[0] == 18 && chIdx[1] == 22 && chIdx[2] == 25 && chIdx[3] == 28)
    {
        Decode2SiPM(18, 22, firedStrip2, weight2, amp);
        Decode2SiPM(25, 28, firedStrip1, weight1, amp);
        return 2;
    }

    // Normal cases:
    if (ch1 == 2 && ch2 == 3 && ch3 == 6 && ch4 == 7)
    {
        firedStrip1 = 0;
        firedStrip2 = 1;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 6];
        weight2 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 7];

        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 3 && ch3 == 4 && ch4 == 7)
    {
        firedStrip1 = 1;
        firedStrip2 = 2;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 7];
        weight2 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 4];

        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 1 && ch3 == 4 && ch4 == 5)
    {
        firedStrip1 = 2;
        firedStrip2 = 3;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 4];
        weight2 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 5];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 1 && ch2 == 3 && ch3 == 5 && ch4 == 6)
    {
        firedStrip1 = 3;
        firedStrip2 = 4;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 5];
        weight2 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 6];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 3 && ch3 == 6 && ch4 == 7)
    {
        firedStrip1 = 4;
        firedStrip2 = 5;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 6];
        weight2 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 7];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 2 && ch3 == 5 && ch4 == 7)
    {
        firedStrip1 = 5;
        firedStrip2 = 6;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 7];
        weight2 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 5];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 2 && ch2 == 3 && ch3 == 4 && ch4 == 5)
    {
        firedStrip1 = 6;
        firedStrip2 = 7;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 5];
        weight2 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 4];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 3 && ch3 == 4 && ch4 == 6)
    {
        firedStrip1 = 7;
        firedStrip2 = 8;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 4];
        weight2 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 6];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 1 && ch3 == 6 && ch4 == 7)
    {
        firedStrip1 = 8;
        firedStrip2 = 9;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 6];
        weight2 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 7];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 1 && ch2 == 2 && ch3 == 4 && ch4 == 7)
    {
        firedStrip1 = 9;
        firedStrip2 = 10;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 7];
        weight2 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 4];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 2 && ch3 == 4 && ch4 == 5)
    {
        firedStrip1 = 10;
        firedStrip2 = 11;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 4];
        weight2 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 5];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 0 && ch2 == 1 && ch3 == 5 && ch4 == 6)
    {
        firedStrip1 = 11;
        firedStrip2 = 12;
        weight1 = amp[8 * sipmGroup + 0] + amp[8 * sipmGroup + 5];
        weight2 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 6];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 1 && ch2 == 2 && ch3 == 6 && ch4 == 7)
    {
        firedStrip1 = 12;
        firedStrip2 = 13;
        weight1 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 6];
        weight2 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 7];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 2 && ch2 == 3 && ch3 == 5 && ch4 == 7)
    {
        firedStrip1 = 13;
        firedStrip2 = 14;
        weight1 = amp[8 * sipmGroup + 2] + amp[8 * sipmGroup + 7];
        weight2 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 5];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else if (ch1 == 1 && ch2 == 3 && ch3 == 4 && ch4 == 5)
    {
        firedStrip1 = 14;
        firedStrip2 = 15;
        weight1 = amp[8 * sipmGroup + 3] + amp[8 * sipmGroup + 5];
        weight2 = amp[8 * sipmGroup + 1] + amp[8 * sipmGroup + 4];
        // pos = 1.1 * (weight1 * firedStrip1 + weight2 * firedStrip2) / (weight1 + weight2);
    }
    else
    {
        // std::cout << "Decoding Error in Algorithm: " << ch1 << '\t' << ch1NoSort << '\t' << ch2 << '\t' << ch2NoSort << '\t' << ch3 << '\t' << ch3NoSort << '\t' << ch4 << '\t' << ch4NoSort << std::endl;
        // std::cout << "Decoding Error in Algorithm: " << chIdx[0] << '\t' << chIdx[1] << '\t' << chIdx[2] << '\t' << chIdx[3] << std::endl;
        // std::cout << ch1NoSort << '\t' << ch2NoSort << '\t' << ch3NoSort << '\t' << ch4NoSort << std::endl;
        return -1;
    }
    int stripGroup = 3 - sipmGroup;
    firedStrip1 += 16 * stripGroup;
    firedStrip2 += 16 * stripGroup;
    // std::cout << (int)firedStrip1 << '\t' << (int)firedStrip2 << std::endl;
    return 3;
}

struct UserDefine::PeakFitValue
{
    double first = 0;  // mean value
    double second = 0; // constant value
    double sigma = 0;  // sigma value
    PeakFitValue(double x, double y, double z) : first(x), second(y), sigma(z){};
    PeakFitValue() : first(0), second(0), sigma(0){};
};

int UserDefine::FindPeaks(TH1D *h, std::vector<PeakFitValue> &outPeaks, int nRebins, double maxPeakSigma, double threshold)
{
    // TH1D hCopy(*h);
    // TH1D *hInput = &hCopy;
    TH1D *hInput = h;

    static TSpectrum *sp = new TSpectrum;
    if (nRebins > 0)
        hInput->Rebin(1 << nRebins);
    hInput->GetXaxis()->SetRangeUser(3000, 12000);

    sp->Search(hInput, maxPeakSigma, "", threshold);
    hInput->Draw();
    sp->Draw("same");

    std::vector<PeakFitValue> tempPeaks;
    for (int i = 0; i < sp->GetNPeaks(); i++)
        tempPeaks.push_back(PeakFitValue(sp->GetPositionX()[i], sp->GetPositionY()[i], 0));
    std::sort(tempPeaks.begin(), tempPeaks.end(), [](const PeakFitValue &temp1, const PeakFitValue &temp2)
              { return temp1.first < temp2.first; });

    if (tempPeaks.size() <= 1)
    {
        std::cout << "Error: only found 1 peak. return false." << std::endl;
        return 0;
    }

    for (int peakNo = 0; peakNo < tempPeaks.size(); peakNo++)
    {
        double dev2Next = 0;
        if (peakNo == 0)
        {
            dev2Next = tempPeaks[1].first - tempPeaks[0].first;
        }
        else
        {
            dev2Next = tempPeaks[peakNo].first - tempPeaks[peakNo - 1].first;
        }

        double currentPeak = tempPeaks[peakNo].first;
        static int count = 0;
        // auto f0 = new TF1(Form("gausn%d", count++), "gausn", 0, 65536);
        // TF1 &f = *f0;
        TF1 f(Form("gausn%d", count++), "gausn", 0, 65536);

        if (peakNo == 0)
            hInput->Fit(&f, "Q", "", currentPeak - dev2Next * 0.3, currentPeak + dev2Next * 0.3);
        else
            hInput->Fit(&f, "+Q", "", currentPeak - dev2Next * 0.3, currentPeak + dev2Next * 0.3);
        double fitCurrentPeak = f.GetParameter(1);
        double fitPeakAmplitude = f.GetParameter(0);
        if (fitCurrentPeak <= currentPeak - dev2Next * 0.3 || fitCurrentPeak >= currentPeak + dev2Next * 0.3)
        {
            tempPeaks[peakNo].first = -1;
            tempPeaks[peakNo].second = -1;
            std::cout << "Error: Peak " << peakNo << " is out of range: [" << currentPeak - dev2Next * 0.3 << ", " << currentPeak + dev2Next * 0.3 << "]. Return false." << std::endl;
            tempPeaks.resize(peakNo);
            //            return 0;
        }
        else
        {
            tempPeaks[peakNo].first = fitCurrentPeak;
            tempPeaks[peakNo].second = fitPeakAmplitude;
            tempPeaks[peakNo].sigma = f.GetParameter(2);
        }
    }

    // Process
    outPeaks.clear();
    for (int i = 0; i < tempPeaks.size(); i++)
    {
        if (tempPeaks[i].first > 0)
        {
            outPeaks.push_back(tempPeaks[i]);
        }
    }

    return outPeaks.size();
}

TH1D *hAPed[16];
TH1D *hBPed[16];
double gGainA[16];
double gPedA[16];
double gGainB[16];
double gPedB[16];

void UserDefine::SpCorrect()
{
    auto fileHist = new TFile("3/histo.root");
    auto fileHist2 = new TFile("../NoReflection/histo.root");
    TH1D *hA[16];
    TH1D *hB[16];
    TH1D *hANoRef[16];
    TH1D *hBNoRef[16];

    // TH1D *hAPed[16];
    // TH1D *hBPed[16];
    for (int i = 0; i < 16; i++)
    {
        hAPed[i] = new TH1D(Form("hAPed%d", i), Form("hAPed%d", i), 16384, 0, 65536);
        hBPed[i] = new TH1D(Form("hBPed%d", i), Form("hBPed%d", i), 16384, 0, 65536);
    }

    auto fileOrigin = new TFile("3/FEEBoard/Aligned.root");
    auto treeOrigin = (TTree *)(fileOrigin->Get("board"));
    double HGamp[32];
    double TDCTime[33];
    treeOrigin->SetBranchAddress("HGamp", HGamp);
    treeOrigin->SetBranchAddress("TDCTime", TDCTime);
    for (int jentry = 0; jentry < treeOrigin->GetEntries(); jentry++)
    {
        if (jentry % 10000 == 0)
            std::cout << jentry << std::endl;
        treeOrigin->GetEntry(jentry);
        if (TDCTime[32] != 0)
            continue;
        for (int i = 0; i < 16; i++)
        {
            int FiberA = GetFiberCode(i, SiPMAB::a);
            int FiberB = GetFiberCode(i, SiPMAB::b);
            hAPed[i]->Fill(HGamp[FiberA]);
            hBPed[i]->Fill(HGamp[FiberB]);
        }
    }
    std::ofstream foutA("GainA.txt");
    std::ofstream foutB("GainB.txt");

    for (int i = 0; i < 16; i++)
    {
        hA[i] = (TH1D *)(fileHist->Get(Form("hA%d", i)));
        hB[i] = (TH1D *)(fileHist->Get(Form("hB%d", i)));
        // hAPed[i] = (TH1D *)(fileHist->Get(Form("hAPed%d", i)));
        // hBPed[i] = (TH1D *)(fileHist->Get(Form("hBPed%d", i)));

        hANoRef[i] = (TH1D *)(fileHist2->Get(Form("hA%d", i)));
        hBNoRef[i] = (TH1D *)(fileHist2->Get(Form("hB%d", i)));

        std::vector<PeakFitValue> outPeaks;
        int nTemp1 = 0;
        double peakDev[20];
        double gain = 0;
        double ped = 0;

        FindPeaks(hAPed[i], outPeaks, 100);
        nTemp1 = (outPeaks.size() - 4) / 2 * 2;
        for (int j = 0; j < outPeaks.size(); j++)
        {
#ifdef DEBUG
            foutA << i << '\t' << outPeaks[j].first << '\t' << outPeaks[j].second << '\t' << outPeaks[j].sigma;
            if (j > 1)
                foutA << '\t' << outPeaks[j].first - outPeaks[j - 1].first;
            foutA << std::endl;
#endif

            if (j > 3)
                peakDev[j - 4] = outPeaks[j].first;
        }
        for (int j = 0; j < nTemp1 / 2; j++)
        {
            gain += peakDev[nTemp1 / 2 + j] - peakDev[j];
        }
        gain /= (nTemp1 / 2) * (nTemp1 / 2);
        ped = peakDev[0] - 4 * gain;
        foutA << i << '\t' << gain << '\t' << ped << std::endl;
        gGainA[i] = gain;
        gPedA[i] = ped;

        FindPeaks(hBPed[i], outPeaks, 100);
        nTemp1 = (outPeaks.size() - 4) / 2 * 2;
        for (int j = 0; j < outPeaks.size(); j++)
        {
#ifdef DEBUG
            foutB << i << '\t' << outPeaks[j].first << '\t' << outPeaks[j].second << '\t' << outPeaks[j].sigma;
            if (j > 1)
                foutB << '\t' << outPeaks[j].first - outPeaks[j - 1].first;
            foutB << std::endl;
#endif

            if (j > 3)
                peakDev[j - 4] = outPeaks[j].first;
        }
        gain = 0;
        for (int j = 0; j < nTemp1 / 2; j++)
        {
            gain += peakDev[nTemp1 / 2 + j] - peakDev[j];
        }
        gain /= (nTemp1 / 2) * (nTemp1 / 2);
        ped = peakDev[0] - 4 * gain;
        gGainB[i] = gain;
        gPedB[i] = ped;

        foutB << i << '\t' << gain << '\t' << ped << std::endl;

        // std::cout << hA[i] << '\t' << hB[i] << std::endl;
    }
    TH1D *hRefHist = new TH1D("hRefHist", "hRefHist", 1000, 0, 60000);
    auto tg = new TGraph();
    for (int i = 0; i < hA[7]->GetNbinsX(); i++)
    {
        hRefHist->SetBinContent(i, hA[7]->GetBinContent(i));
    }

    hRefHist->Rebin(10);
    hRefHist->Smooth();
    for (int i = 0; i < hRefHist->GetNbinsX(); i++)
        tg->SetPoint(i, hRefHist->GetBinCenter(i), hRefHist->GetBinContent(i));
    auto lambda = [tg](double *x, double *par) -> double
    {
        return par[0] * tg->Eval((x[0] - par[1]) / par[2]);
    };

    auto fun = new TF1("fun", lambda, 0, 60000, 3);
    fun->SetParameters(1, 0, 1);
    fun->SetParName(0, "A");
    fun->SetParName(1, "x0");
    fun->SetParName(2, "eta");

    auto c = new TCanvas("c", "c", 1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    for (int i = 0; i < 16; i++)
    {
        hA[i]->Rebin(10);
        hANoRef[i]->Rebin(10);
        hANoRef[i]->SetLineColor(kRed);

        auto axiA = hANoRef[i]->GetXaxis();
        axiA->Set(axiA->GetNbins(), (axiA->GetXmin() - gPedA[i]) / gGainA[i], (axiA->GetXmax() - gPedA[i]) / gGainA[i]);
        axiA = hA[i]->GetXaxis();
        axiA->Set(axiA->GetNbins(), (axiA->GetXmin() - gPedA[i]) / gGainA[i], (axiA->GetXmax() - gPedA[i]) / gGainA[i]);
        hANoRef[i]->SetTitle(";PE number;Counts");

        hANoRef[i]->DrawNormalized();
        hA[i]->DrawNormalized("same");
        hANoRef[i]->ResetStats();
        hANoRef[i]->SaveAs(Form("pdf/Compare-hA%d-histNoRef.root", i));
        hA[i]->ResetStats();
        hA[i]->SaveAs(Form("pdf/Compare-hA%d-histRef.root", i));
        c->SaveAs(Form("pdf/Compare-hA%d.pdf", i));
        c->SaveAs(Form("pdf/Compare-hA%d.root", i));

        // hA[i]->Fit(fun, "", "", 20000, 60000);
        hA[i]->Draw();
        hA[i]->Fit(fun);
        c->SaveAs(Form("pdf/hAFit%d.pdf", i));

        hAPed[i]->GetXaxis()->SetRangeUser(2000, 10000);
        hAPed[i]->Draw();
        c->SaveAs(Form("pdf/hAPed%d.pdf", i));

        hANoRef[i]->Draw();
        hANoRef[i]->Fit(fun);
        c->SaveAs(Form("pdf/hANoRefFit%d.pdf", i));

        // hB

        hB[i]->Rebin(10);
        hBNoRef[i]->Rebin(10);
        hBNoRef[i]->SetLineColor(kRed);

        auto axiB = hBNoRef[i]->GetXaxis();
        axiB->Set(axiB->GetNbins(), (axiB->GetXmin() - gPedB[i]) / gGainB[i], (axiB->GetXmax() - gPedB[i]) / gGainB[i]);
        axiB = hB[i]->GetXaxis();
        axiB->Set(axiB->GetNbins(), (axiB->GetXmin() - gPedB[i]) / gGainB[i], (axiB->GetXmax() - gPedB[i]) / gGainB[i]);
        hBNoRef[i]->SetTitle(";PE number;Counts");

        hBNoRef[i]->DrawNormalized();
        hB[i]->DrawNormalized("same");
        c->SaveAs(Form("pdf/Compare-hB%d.pdf", i));

        // hB[i]->Fit(fun, "", "", 20000, 60000);
        hB[i]->Draw();
        hB[i]->Fit(fun);
        c->SaveAs(Form("pdf/hBFit%d.pdf", i));

        hBNoRef[i]->Draw();
        hBNoRef[i]->Fit(fun);
        c->SaveAs(Form("pdf/hBNoRefFit%d.pdf", i));

        hBPed[i]->GetXaxis()->SetRangeUser(2000, 10000);
        hBPed[i]->Draw();
        c->SaveAs(Form("pdf/hBPed%d.pdf", i));
    }
}