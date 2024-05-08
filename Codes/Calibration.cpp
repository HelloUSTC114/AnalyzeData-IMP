#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TGraph.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

// #include "board.h"
//! TODO: Change to board.h in the future, if using in real program
#include "board.C"
#include "Decode.cpp"

class BoardCalibration
{
public:
    ~BoardCalibration() { Clear(); };

    void Clear();
    void Init(int boardNo, std::string sDataFolder = "../", std::string sResultFolder = "Cali/");
    bool ReadResultFile(int boardNo, std::string sResultFolder = "Cali/");
    bool ReadDataFile(std::string sDataFolder = "../");

    int DoCalibration(int boardNo, std::string sResultFolder = "Cali/", std::string sDataFolder = "../");

    void DrawPlots(std::string sFolder = "Cali/Plots/");
    void GetStripCaliFactor(std::string sFolder = "Cali/Plots/", bool drawFlag = 0);
    void DumpStripCali();

    board *GetBoard() { return fBoard; }

    // Calculate effective sipm hg value / sipm factor
    double CalcSiPMEffValue(int sipm);
    void CalcStripEffValue(int strip, double &weighta, double &weightb);

private:
    // Basic information
    int fBoardNo;
    std::string fsResultFolder, fsResultFileName, fsDataFolder, fsDataFileName;
    // TFile *fResultFile = NULL, *fDataFile = NULL;
    TFile *fResultFile = NULL;
    board *fBoard = NULL;
    bool fGenerateFlag = 0, fReadFlag = 0;
    bool fResultFlag = 0, fDataFlag = 0;

    // Channel information, used for SiPM gain calculation
    TH1D *fhChannels[32]{NULL};
    TH1D *fhPed[32]{NULL};
    TH1D *fhSPE[32]{NULL};        // Final converted SPE spectrum for SiPM
    TGraph *fgChannels[32]{NULL}; // Fit Using LG ADC-HG ADC (Y-X)
    TF1 *ffChannels[32]{NULL};    // Convert LG adc value into HG adc value

    TGraph *fgSPEFit[32]{NULL}; // Fit using HG ADC-SPE Number (Y-X)
    TF1 *ffSPEFit[32]{NULL};    // Convert HG ADC value into SPE number

    // Strip information, used for light collection calibration
    TH1D *fhStripsA[64]{NULL};
    TH1D *fhStripsB[64]{NULL};

    // Strip Factor Calculate
    double fPedStripA[64]{0};
    double fPedStripB[64]{0};
    double fSignalStripA[64]{0};
    double fSignalStripB[64]{0};
    double fFactorStripA[64]{0};
    double fFactorStripB[64]{0};

    std::vector<double> fvSiPMPed[32];
    std::vector<double> fvSiPMSignal[32];
    std::vector<double> fvSiPMFactor[32];
    double fSiPMPed[32], fSiPMSignal[32], fSiPMFactor[32];
    // For generate calibration file
    bool ReadDataFile(std::string sDataFolder, std::string sDataFile);
    bool GenerateResultFile(std::string sResultFolder, std::string sResultFile);

    int GenerateSiPMSpectrum();
    int GenerateStripSpectrum();

    // For read calibration file
    bool ReadResultFile(std::string sResultFolder, std::string sResultFile);

    // For Inside calculation
    double CalculateSPE(int sipm);
    double CalculateHG(int sipm, double threshold = 30000);

    // For calibrate strip spectrum
    double FitStripSpectrum(TH1D *hspectrum);
};

void BoardCalibration::Clear()
{
    if (fResultFile)
    {
        fResultFile->cd();
        for (int sipm = 0; sipm < 32; sipm++)
        {
            if (fhChannels[sipm])
            {
                if (fResultFile->IsWritable())
                    fhChannels[sipm]->Write();
                fhChannels[sipm] = NULL;
            }

            if (fhPed[sipm])
            {
                if (fResultFile->IsWritable())
                    fhPed[sipm]->Write();
                fhPed[sipm] = NULL;
            }

            if (fhSPE[sipm])
            {
                if (fResultFile->IsWritable())
                    fhSPE[sipm]->Write();
                fhSPE[sipm] = NULL;
            }

            if (fgChannels[sipm])
            {
                if (fResultFile->IsWritable())
                    fgChannels[sipm]->Write(Form("gsipm%d", sipm));
                delete fgChannels[sipm];
                fgChannels[sipm] = NULL;
            }

            if (ffChannels[sipm])
            {
                if (fResultFile->IsWritable())
                    ffChannels[sipm]->Write(Form("fsipm%d", sipm));
                delete ffChannels[sipm];
                ffChannels[sipm] = NULL;
            }

            if (fgSPEFit[sipm])
            {
                if (fResultFile->IsWritable())
                    fgSPEFit[sipm]->Write(Form("gspe%d", sipm));
                delete fgSPEFit[sipm];
                fgSPEFit[sipm] = NULL;
            }

            if (ffSPEFit[sipm])
            {
                if (fResultFile->IsWritable())
                    ffSPEFit[sipm]->Write(Form("fspe%d", sipm));
                delete ffSPEFit[sipm];
                ffSPEFit[sipm] = NULL;
            }
        }
        //! TODO: Add strip information storage
        for (int strip = 0; strip < 64; strip++)
        // for (int strip = 0; strip < 0; strip++)
        {
            if (fhStripsA[strip])
            {
                if (fResultFile->IsWritable())
                    fhStripsA[strip]->Write();
                fhStripsA[strip] = NULL;
            }
            if (fhStripsB[strip])
            {
                if (fResultFile->IsWritable())
                    fhStripsB[strip]->Write();
                fhStripsB[strip] = NULL;
            }
        }
        delete fResultFile;
    }

    // if (fDataFile)
    //     delete fDataFile;
    if (fBoard)
        delete fBoard;
    // fResultFile = fDataFile = NULL;
    fResultFile = NULL;
    fBoard = NULL;
    fDataFlag = fResultFlag = 0;
    fReadFlag = fGenerateFlag = 0;
}

void BoardCalibration::Init(int boardNo, std::string sDataFolder, std::string sResultFolder)
{
    if (fReadFlag)
        return;

    fBoardNo = boardNo;
    ReadDataFile(sDataFolder, Form("Board%d-Aligned.root", fBoardNo));
    GenerateResultFile(sResultFolder, Form("Board%d-Cali.root", fBoardNo));

    fResultFile->cd();
    for (int sipm = 0; sipm < 32; sipm++)
    {
        fhChannels[sipm] = new TH1D(Form("hsipm%d", sipm), Form("Spectrum of SiPM %d in board %d;ADC value;Counts", sipm, fBoardNo), 16384, 0, 65536);
        fhPed[sipm] = new TH1D(Form("hped%d", sipm), Form("Pedestal of SiPM %d in board %d;ADC value;Counts", sipm, fBoardNo), 16384, 0, 65536);
        fhSPE[sipm] = new TH1D(Form("hspe%d", sipm), Form("SPE spectrum of SiPM %d in board %d;ADC value;Counts", sipm, fBoardNo), 16384, -10, 350);

        fgChannels[sipm] = new TGraph();
        fgChannels[sipm]->SetName(Form("gsipm%d", sipm));
        fgChannels[sipm]->SetTitle(Form("Correlation between HG and LG for SiPM %d in board %d;HG ADC;LG ADC", sipm, fBoardNo));
        ffChannels[sipm] = new TF1(Form("fsipm%d", sipm), "pol1", 0, 65536);

        fgSPEFit[sipm] = new TGraph();
        fgSPEFit[sipm]->SetName(Form("gspe%d", sipm));
        fgSPEFit[sipm]->SetTitle(Form("Fit SPE for SiPM %d in board %d;PE Count;HG ADC value", sipm, fBoardNo));
        // ffSPEFit[sipm] = new TF1(Form("fspe%d", sipm), "pol1", 0, 100);
        ffSPEFit[sipm] = new TF1(Form("fspe%d", sipm), "pol1", 0, 65536);
    }
    for (int strip = 0; strip < 64; strip++)
    {
        fhStripsA[strip] = new TH1D(Form("hstripa%d", strip), Form("Spectrum of strip %d A in board %d;ADC value;Counts", strip, fBoardNo), 32768, 0, 131072);
        fhStripsB[strip] = new TH1D(Form("hstripb%d", strip), Form("Spectrum of strip %d B in board %d;ADC value;Counts", strip, fBoardNo), 32768, 0, 131072);
    }
    fGenerateFlag = 1;
    fReadFlag = 0;
}

bool BoardCalibration::ReadDataFile(std::string sDataFolder, std::string sDataFile)
{
    fsDataFolder = sDataFolder;
    fsDataFileName = sDataFile;
    // fDataFile = new TFile((fsDataFolder + "/" + fsDataFileName).c_str());
    fBoard = new board(fsDataFolder + "/" + fsDataFileName);

    // if (!fDataFile->IsOpen())
    //     return false;
    if (!(fBoard->fChain))
        return false;

    fDataFlag = 1;
    return true;
}

bool BoardCalibration::GenerateResultFile(std::string sResultFolder, std::string sResultFile)
{
    fsResultFolder = sResultFolder;
    fsResultFileName = sResultFile;
    fResultFile = new TFile((fsResultFolder + "/" + fsResultFileName).c_str(), "recreate");
    if (!fResultFile->IsOpen())
        return false;
    return true;
}

bool BoardCalibration::ReadResultFile(std::string sResultFolder, std::string sResultFile)
{
    if (fGenerateFlag)
        return false;
    if (fResultFlag)
        return false;

    fsResultFolder = sResultFolder;
    fsResultFileName = sResultFile;
    fResultFile = new TFile((fsResultFolder + "/" + fsResultFileName).c_str());
    if (!fResultFile->IsOpen())
        return false;

    for (int sipm = 0; sipm < 32; sipm++)
    {
        fhChannels[sipm] = (TH1D *)(fResultFile->Get(Form("hsipm%d", sipm)));
        if (!fhChannels[sipm])
            return false;

        fhPed[sipm] = (TH1D *)(fResultFile->Get(Form("hped%d", sipm)));
        if (!fhPed[sipm])
            return false;

        fhSPE[sipm] = (TH1D *)(fResultFile->Get(Form("hspe%d", sipm)));
        if (!fhSPE[sipm])
            return false;

        fgChannels[sipm] = (TGraph *)(fResultFile->Get(Form("gsipm%d", sipm)));
        if (!fgChannels[sipm])
            return false;

        ffChannels[sipm] = (TF1 *)(fResultFile->Get(Form("fsipm%d", sipm)));
        if (!ffChannels[sipm])
            return false;

        fgSPEFit[sipm] = (TGraph *)(fResultFile->Get(Form("gspe%d", sipm)));
        if (!fgSPEFit[sipm])
            return false;

        ffSPEFit[sipm] = (TF1 *)(fResultFile->Get(Form("fspe%d", sipm)));
        if (!ffSPEFit[sipm])
            return false;
    }
    //! TODO: Add strip information storage
    for (int strip = 0; strip < 64; strip++)
    // for (int strip = 0; strip < 0; strip++)
    {
        fhStripsA[strip] = (TH1D *)(fResultFile->Get(Form("hstripa%d", strip)));
        if (!fhStripsA[strip])
            return false;
        fhStripsB[strip] = (TH1D *)(fResultFile->Get(Form("hstripb%d", strip)));
        if (!fhStripsB[strip])
            return false;
    }
    fReadFlag = 1;
    fResultFlag = 1;

    fGenerateFlag = 0;

    return true;
}

bool BoardCalibration::ReadResultFile(int boardNo, std::string sResultFolder)
{
    auto flag = ReadResultFile(sResultFolder, Form("Board%d-Cali.root", boardNo));
    if (flag)
        fBoardNo = boardNo;
    return flag;
}

bool BoardCalibration::ReadDataFile(std::string sDataFolder)
{
    if (fDataFlag)
        return false;
    return ReadDataFile(sDataFolder, Form("Board%d-Aligned.root", fBoardNo));
}

auto gCanvas = new TCanvas("c", "c", 1);
void BoardCalibration::DrawPlots(std::string sFolder)
{
    if (!fResultFlag)
        return;
    // Draw sipm calibration
    for (int sipm = 0; sipm < 32; sipm++)
    {
        fhChannels[sipm]->Draw();
        gCanvas->SaveAs(Form("%s/B%dhsipm%d.jpg", sFolder.c_str(), fBoardNo, sipm));
        fhPed[sipm]->Draw();
        gCanvas->SaveAs(Form("%s/B%dhspe%d.jpg", sFolder.c_str(), fBoardNo, sipm));

        // LG HG
        TGraph tg;
        for (int point = 0; point < fgChannels[sipm]->GetN(); point++)
            tg.SetPoint(point, fgChannels[sipm]->GetPointY(point), fgChannels[sipm]->GetPointX(point));
        tg.SetTitle(fgChannels[sipm]->GetTitle());
        tg.GetXaxis()->SetTitle(fgChannels[sipm]->GetYaxis()->GetTitle());
        tg.GetYaxis()->SetTitle(fgChannels[sipm]->GetXaxis()->GetTitle());
        tg.Draw("AZP");
        ffChannels[sipm]->Draw("same");
        gCanvas->SaveAs(Form("%s/B%dgsipm%d.jpg", sFolder.c_str(), fBoardNo, sipm));

        // SPE
        tg.Set(0);
        for (int point = 0; point < fgSPEFit[sipm]->GetN(); point++)
            tg.SetPoint(point, fgSPEFit[sipm]->GetPointY(point), fgSPEFit[sipm]->GetPointX(point));
        tg.SetTitle(fgSPEFit[sipm]->GetTitle());
        tg.GetXaxis()->SetTitle(fgSPEFit[sipm]->GetYaxis()->GetTitle());
        tg.GetYaxis()->SetTitle(fgSPEFit[sipm]->GetXaxis()->GetTitle());
        tg.Draw("AZ*");
        ffSPEFit[sipm]->Draw("same");
        gCanvas->SaveAs(Form("%s/B%dgspe%d.jpg", sFolder.c_str(), fBoardNo, sipm));
    }

    // Draw strip calibration
    for (int strip = 0; strip < 64; strip++)
    {
        TH1D hA0(*fhStripsA[strip]);
        auto hA = &hA0;
        hA->Rebin(128);
        hA->Draw();
        gCanvas->SaveAs(Form("%s/B%dstripa%d.jpg", sFolder.c_str(), fBoardNo, strip));

        // auto hB = fhStripsB[strip];
        TH1D hB0(*fhStripsB[strip]);
        auto hB = &hB0;
        hB->Rebin(128);
        hB->Draw();
        gCanvas->SaveAs(Form("%s/B%dstripb%d.jpg", sFolder.c_str(), fBoardNo, strip));
    }
}

TF1 *gGaus = new TF1("ggaus", "gaus", 0, 131072);
TSpectrum *gsp = new TSpectrum();
#include <TStyle.h>
#include <TLatex.h>
auto gLatex = new TLatex();
void BoardCalibration::GetStripCaliFactor(std::string sFolder, bool drawFlag)
{
    gStyle->SetOptFit(111);
    if (!fResultFlag)
        return;
    for (int strip = 0; strip < 64; strip++)
    {
        gCanvas->cd();
        // Calculate pedestal
        int sipma = UserDefine::GetFiberCode(strip, UserDefine::a);
        TH1D ha0(*fhPed[sipma]);
        auto ha = &ha0;
        gsp->Search(ha);
        double peak0 = gsp->GetPositionX()[0];
        ha->Fit(gGaus, "Q", "", peak0 - 100, peak0 + 100);
        fPedStripA[strip] = gGaus->GetParameter(1);
        ha->Draw();
        if (drawFlag)
            gCanvas->SaveAs(Form("%s/B%dstripa-fit%d.jpg", sFolder.c_str(), fBoardNo, strip));
        // Calculate signal mean
        TH1D has0(*fhStripsA[strip]);
        auto has = &has0;
        has->Rebin(128);
        has->GetXaxis()->SetRangeUser(fPedStripA[strip] - 1000, 131072);
        has->ResetStats();
        auto x = FitStripSpectrum(has);
        has->Draw();
        fSignalStripA[strip] = x;
        fFactorStripA[strip] = fSignalStripA[strip] - fPedStripA[strip];
        // Related SiPM
        fvSiPMPed[sipma].push_back(fPedStripA[strip]);
        fvSiPMSignal[sipma].push_back(fSignalStripA[strip]);
        fvSiPMFactor[sipma].push_back(fFactorStripA[strip]);

        gLatex->DrawLatexNDC(0.2, 0.2, Form("%.1f,%.1f,%.1f", fPedStripA[strip], fSignalStripA[strip], fSignalStripA[strip] - fPedStripA[strip]));
        if (drawFlag)
            gCanvas->SaveAs(Form("%s/B%dstripa-signal%d.jpg", sFolder.c_str(), fBoardNo, strip));

        // Calculate pedestal B
        int sipmb = UserDefine::GetFiberCode(strip, UserDefine::b);
        TH1D hb0(*fhPed[sipma]);
        auto hb = &hb0;
        gsp->Search(hb);
        peak0 = gsp->GetPositionX()[0];
        hb->Fit(gGaus, "Q", "", peak0 - 100, peak0 + 100);
        fPedStripB[strip] = gGaus->GetParameter(1);
        hb->Draw();
        if (drawFlag)
            gCanvas->SaveAs(Form("%s/B%dstripb-fit%d.jpg", sFolder.c_str(), fBoardNo, strip));
        // Calculate signal mean B
        TH1D hbs0(*fhStripsB[strip]);
        auto hbs = &hbs0;
        hbs->Rebin(128);
        hbs->GetXaxis()->SetRangeUser(fPedStripA[strip] - 1000, 131072);
        hbs->ResetStats();
        x = FitStripSpectrum(hbs);
        hbs->Draw();
        fSignalStripB[strip] = x;
        fFactorStripB[strip] = fSignalStripB[strip] - fPedStripB[strip];
        // Related SiPM
        fvSiPMPed[sipmb].push_back(fPedStripB[strip]);
        fvSiPMSignal[sipmb].push_back(fSignalStripB[strip]);
        fvSiPMFactor[sipmb].push_back(fFactorStripB[strip]);

        gLatex->DrawLatexNDC(0.2, 0.2, Form("%.1f,%.1f,%.1f", fPedStripB[strip], fSignalStripB[strip], fSignalStripB[strip] - fPedStripB[strip]));
        if (drawFlag)
            gCanvas->SaveAs(Form("%s/B%dstripb-signal%d.jpg", sFolder.c_str(), fBoardNo, strip));
    }

    for (int sipm = 0; sipm < 32; sipm++)
    {
        fSiPMPed[sipm] = std::accumulate(fvSiPMPed[sipm].begin(), fvSiPMPed[sipm].end(), 0.0) / (fvSiPMPed[sipm].size());
        fSiPMSignal[sipm] = std::accumulate(fvSiPMSignal[sipm].begin(), fvSiPMSignal[sipm].end(), 0.0) / (fvSiPMSignal[sipm].size());
        fSiPMFactor[sipm] = std::accumulate(fvSiPMFactor[sipm].begin(), fvSiPMFactor[sipm].end(), 0.0) / (fvSiPMFactor[sipm].size());
    }
}

void BoardCalibration::DumpStripCali()
{
    std::cout << "Board " << fBoardNo << std::endl;
    // for (int strip = 0; strip < 64; strip++)
    // {
    //     std::cout << "strip " << strip << '\t' << " A: " << fPedStripA[strip] << '\t' << fSignalStripA[strip] << '\t' << fFactorStripA[strip] << std::endl;
    //     std::cout << "\t\t B: " << fPedStripB[strip] << '\t' << fSignalStripB[strip] << '\t' << fFactorStripB[strip] << std::endl;
    // }
    for (int sipm = 0; sipm < 32; sipm++)
    {
        std::cout << sipm << std::endl;
        for (int idx = 0; idx < fvSiPMPed[sipm].size(); idx++)
            std::cout << fvSiPMPed[sipm][idx] << '\t';
        std::cout << '\t' << fSiPMPed[sipm] << std::endl;
        for (int idx = 0; idx < fvSiPMSignal[sipm].size(); idx++)
            std::cout << fvSiPMSignal[sipm][idx] << '\t';
        std::cout << '\t' << fSiPMSignal[sipm] << std::endl;
        for (int idx = 0; idx < fvSiPMFactor[sipm].size(); idx++)
            std::cout << fvSiPMFactor[sipm][idx] << '\t';
        std::cout << '\t' << fSiPMFactor[sipm] << std::endl;
    }
    std::cout << std::endl;
}

double BoardCalibration::CalcSiPMEffValue(int sipm)
{
    double hg = CalculateHG(sipm);
    return (hg - fSiPMPed[sipm]) / (fSiPMFactor[sipm]);
}

void BoardCalibration::CalcStripEffValue(int strip, double &weighta, double &weightb)
{
    double hga = CalculateHG(UserDefine::GetFiberCode(strip, UserDefine::SiPMAB::a));
    double hgb = CalculateHG(UserDefine::GetFiberCode(strip, UserDefine::SiPMAB::b));
    weighta = (hga - fPedStripA[strip]) / fFactorStripA[strip];
    weightb = (hgb - fPedStripB[strip]) / fFactorStripB[strip];
}

int BoardCalibration::DoCalibration(int boardNo, std::string sResultFolder, std::string sDataFolder)
{
    if (!ReadResultFile(boardNo, sResultFolder))
    {
        Init(boardNo, sDataFolder, sResultFolder);
        if (!fGenerateFlag)
            return -1;
        GenerateSiPMSpectrum();
        GenerateStripSpectrum();
    }
    fResultFlag = 1;

    return 0;
}

int BoardCalibration::GenerateSiPMSpectrum()
{
    if (!fBoard || !(fBoard->fChain))
        return -1;

    int totalEntries = fBoard->fChain->GetEntries();
    int graphPointCounter[32]{0};
    for (int idx = 0; idx < totalEntries; idx++)
    // for (int idx = 0; idx < 100000; idx++)
    {
        if (idx % 1000 == 0)
            std::cout << "Processing: idx: " << idx << " in board: " << fBoardNo << std::endl;
        fBoard->GetEntry(idx);
        for (int sipm = 0; sipm < 32; sipm++)
        {
            fhPed[sipm]->Fill(fBoard->HGamp[sipm]);
            if (fBoard->TDCTime[sipm] > 0)
            {
                fhChannels[sipm]->Fill(fBoard->HGamp[sipm]);
                fgChannels[sipm]->SetPoint(graphPointCounter[sipm]++, fBoard->HGamp[sipm], fBoard->LGamp[sipm]);
            }
        }
    }

    for (int sipm = 0; sipm < 32; sipm++)
    {
        fgChannels[sipm]->Fit(ffChannels[sipm], "Q", "", 6000, 35000);
        double p0 = ffChannels[sipm]->GetParameter(0);
        double p1 = ffChannels[sipm]->GetParameter(1);
        ffChannels[sipm]->SetParameters(-p0 / p1, 1.0 / p1);

        std::vector<UserDefine::PeakFitValue> fitValue;
        UserDefine::FindPeaks(fhPed[sipm], fitValue, 2, 5, 0.001);
        for (int idx = 0; idx < fitValue.size(); idx++)
        {
            // std::cout << sipm << '\t' << idx << '\t' << fitValue[idx].first << '\t' << fitValue[idx].sigma << std::endl;
            fgSPEFit[sipm]->SetPoint(idx, idx, fitValue[idx].first);
        }
        // std::cout << std::endl;
        fgSPEFit[sipm]->Fit(ffSPEFit[sipm], "Q", "", 0, fitValue.size() - 1);
        p0 = ffSPEFit[sipm]->GetParameter(0);
        p1 = ffSPEFit[sipm]->GetParameter(1);
        ffSPEFit[sipm]->SetParameters(-p0 / p1, 1.0 / p1);
    }

    return 0;
}

int BoardCalibration::GenerateStripSpectrum()
{
    if (!fBoard || !(fBoard->fChain))
        return -1;

    int totalEntries = fBoard->fChain->GetEntries();
    for (int idx = 0; idx < totalEntries; idx++)
    {
        if (idx % 1000 == 0)
            std::cout << "Processing: idx: " << idx << " in board: " << fBoardNo << std::endl;
        fBoard->GetEntry(idx);

        std::vector<std::pair<int, double>> tempArray; // Array stores <channel, Effective HG value>
        std::map<int, double> mapArray;                // Array stores <channel, Effective HG value>
        for (int sipm = 0; sipm < 32; sipm++)
            if (fBoard->TDCTime[sipm] > 0)
            {
                double hgValid = CalculateHG(sipm, 0);
                fhSPE[sipm]->Fill(hgValid);
                if (hgValid > 10)
                {
                    tempArray.push_back({sipm, hgValid});
                    mapArray.insert({sipm, hgValid});
                }
            }
        //! TODO: Add decode for strip calibration
        std::sort(tempArray.begin(), tempArray.end(), [](std::pair<int, double> a, std::pair<int, double> b)
                  { return a.second > b.second; });
        if (tempArray.size() < 4 && tempArray.size() >= 2)
        {
            int strip;
            double weight;
            auto flag = UserDefine::Decode2SiPM(tempArray[0].first, tempArray[1].first, strip, weight, fBoard->HGamp);
            if (flag < 0)
                continue;

            fhStripsA[strip]->Fill(mapArray[UserDefine::GetFiberCode(strip, UserDefine::a)]);
            fhStripsB[strip]->Fill(mapArray[UserDefine::GetFiberCode(strip, UserDefine::b)]);
            // std::cout << tempArray[0].first << '\t' << tempArray[1].first << '\t' << UserDefine::GetFiberCode(strip, UserDefine::a) << '\t' << UserDefine::GetFiberCode(strip, UserDefine::b) << std::endl;
        }
        else if (tempArray.size() >= 4)
        {
            int strip1, strip2;
            double weight1, weight2;
            auto flag = UserDefine::Decode4SiPM(tempArray[0].first, tempArray[1].first, tempArray[2].first, tempArray[3].first, fBoard->HGamp, weight1, weight2, strip1, strip2);
            if (flag < 0)
                continue;

            fhStripsA[strip1]->Fill(mapArray[UserDefine::GetFiberCode(strip1, UserDefine::a)]);
            fhStripsB[strip1]->Fill(mapArray[UserDefine::GetFiberCode(strip1, UserDefine::b)]);

            fhStripsA[strip2]->Fill(mapArray[UserDefine::GetFiberCode(strip2, UserDefine::a)]);
            fhStripsB[strip2]->Fill(mapArray[UserDefine::GetFiberCode(strip2, UserDefine::b)]);
        }
    }
    return 0;
}

double BoardCalibration::CalculateSPE(int sipm)
{
    // Add judgement
    if (0)
        return 0;
    double spe = ffSPEFit[sipm]->Eval(CalculateHG(sipm));
    return spe;
}

double BoardCalibration::CalculateHG(int sipm, double threshold)
{
    // Add judgement
    if (0)
        return 0;

    double hgValid = 0;
    hgValid = fBoard->HGamp[sipm];
    if (hgValid > threshold)
        hgValid = ffChannels[sipm]->Eval(fBoard->LGamp[sipm]);
    return hgValid;
}

#include <Math/RootFinder.h>
#include <Math/RootFinderAlgorithms.h>
#include <Math/Functor.h>
double BoardCalibration::FitStripSpectrum(TH1D *hspectrum)
{
    auto h0 = hspectrum;
    static auto sGaus = new TF1("sgaus", "gaus", 0, 131072);
    static auto sPol0 = new TF1("spol0", "pol0", 0, 131072);
    static auto finder = new ROOT::Math::RootFinder();
    hspectrum->Fit(sPol0, "Q+", "", hspectrum->GetMean() * 0.7, hspectrum->GetMean() * 1.25);
    hspectrum->Fit(sGaus, "Q+", "", hspectrum->GetMean() * 1.25, 131072);

    TF1 funTemp("temp", "gaus(0)-[3]", 0, 131072);
    for (int i = 0; i < 3; i++)
        funTemp.SetParameter(i, sGaus->GetParameter(i));
    funTemp.SetParameter(3, sPol0->GetParameter(0) / 2.0);

    finder->SetMethod(ROOT::Math::RootFinder::kGSL_BISECTION);

    auto lambda = [&funTemp](double x)
    { return funTemp.Eval(x); };
    ROOT::Math::Functor1D f(lambda);

    finder->SetFunction(f, hspectrum->GetMean() * 1.25, 131072);
    finder->Solve();
    return finder->Root();
}

void Calibration(std::string sCaliFolder = "../Cali", std::string sDataFolder = "../Sum/")
{
    // for (int board : {0, 2, 4})
    gSystem->Exec(Form("mkdir %s", sCaliFolder.c_str()));
    gSystem->Exec(Form("mkdir %s/Plots", sCaliFolder.c_str()));
    // gSystem->Exec(Form("mkdir ../Cali/Plots", sCaliFolder));
    for (int board = 0; board < 8; board++)
    // for (int board = 0; board < 1; board++)
    {
        gSystem->Exec(Form("mkdir %s/Plots/B%d", sCaliFolder.c_str(), board));
        BoardCalibration a;
        // std::cout << a.ReadResultFile(0) << std::endl;
        a.DoCalibration(board, sCaliFolder, sDataFolder);
        a.DrawPlots(Form("%s/Plots/B%d", sCaliFolder.c_str(), board));
        a.GetStripCaliFactor(Form("%s/Plots/B%d", sCaliFolder.c_str(), board), 1);
        a.DumpStripCali();
        std::cout << a.ReadDataFile() << std::endl;
    }
    // a.DrawPlots();
    // a.Init(0);
    // a.GenerateSiPMSpectrum();
    // a.GenerateStripSpectrum();
}