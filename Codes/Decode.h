#ifndef DECODE_H
#define DECODE_H
class TH1D;
namespace UserDefine
{

    enum SiPMAB
    {
        a = 0,
        b,
    };

    int GetFiberCode(int stripNo, SiPMAB ab);

    int Decode2SiPM(int ch1NoSort, int ch2NoSort, int &firedStrip1, double &weight1, double *amp);

    int Decode4SiPM(int ch1NoSort, int ch2NoSort, int ch3NoSort, int ch4NoSort, double *amp, double &weight1, double &weight2, int &firedStrip1, int &firedStrip2);

    struct PeakFitValue;
    int FindPeaks(TH1D *hInput, std::vector<PeakFitValue> &outPeaks, int nRebins = 0, double maxPeakSigma = 10, double threshold = 0.05);

    void SpCorrect();

}

#endif