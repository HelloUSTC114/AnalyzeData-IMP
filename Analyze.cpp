#include <TCanvas.h>

auto c = new TCanvas("c", "c", 1);
void Analyze()
{
    // Open Sum/1stCorrection.root file
    TFile *f = new TFile("Sum/1stCorrection.root");
    if (f->IsZombie())
    {
        cout << "Error: Cannot open Sum/1stCorrection.root file" << endl;
        return;
    }
    // Get cor tree
    TTree *cor = (TTree *)f->Get("cor");
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
    if (1)
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

        // cor->Draw("TMath::ATan(TMath::Sqrt(TMath::Power(tanX,2)+TMath::Power(tanY,2)))/TMath::Pi()*180:TMath::ATan(TMath::Abs(tanX/tanY))/TMath::Pi()*180>>hr200(180,0,90,300,0,15)", "TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))<200&&TMath::Sqrt(TMath::Power((derivedPos[6]+derivedPos[0])/2.0-264.0,2)+TMath::Power((derivedPos[7]+derivedPos[1])/2.0-440.0,2))>150", "COLZ", 100e5);
        // TH2F *hr200 = (TH2F *)gDirectory->Get("hr200");
        // hr200->SetTitle("#theta vs #phi at 200 mm > r >150 mm ;#phi/degree;#theta/degree");
        // hr200->Draw("COLZ");
        // c->SaveAs("Sum/theta-phi-r200.jpg");
    }
}