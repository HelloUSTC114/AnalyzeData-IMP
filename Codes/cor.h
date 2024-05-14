//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  8 17:33:41 2024 by ROOT version 6.22/02
// from TTree cor/corrrection pos
// found on file: 1stCorrection.root
//////////////////////////////////////////////////////////

#ifndef cor_h
#define cor_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TVector3.h"

namespace john
{

    class cor
    {
    public:
        TTree *fChain;  //! pointer to the analyzed TTree or TChain
        Int_t fCurrent; //! current Tree number in a TChain

        // Fixed size dimensions of array or collections stored in the TTree if any.

        // Declaration of leaf types
        Double_t calcPos2[8];
        Int_t firedStrip1[8];
        Double_t weight1[8];
        Double_t weight1a[8];
        Double_t weight1b[8];
        Double_t weight2[8];
        Double_t weight2a[8];
        Double_t weight2b[8];
        Double_t calcPos3[8];
        Double_t eta[8];
        Int_t decodeFlag[8];
        Int_t decodeType[8];
        Double_t corDev[8];
        Double_t corEta[8];
        Int_t evenOdd[8];
        Double_t tanX;
        Double_t tanY;
        Double_t boardHG[8][32];
        Double_t boardLG[8][32];
        Double_t boardTDC[8][33];
        Double_t tdcTime[8];
        Double_t lastSeg[8];
        Double_t nextSeg[8];
        Double_t interval[8];
        Int_t processedCount[8];
        Int_t processedCh[8][32];
        Double_t processedAmp[8][32];
        Double_t processedTDC[8][33];
        Double_t derivedPos[8];
        Double_t fitSigma[8];
        TVector3 *v0;
        TVector3 *v1;
        TVector3 *v2;
        TVector3 *v3;
        TVector3 *vin;
        TVector3 *vout;
        TVector3 *vScat;
        Double_t thetaX;
        Double_t thetaY;
        Double_t theta;
        Double_t p0x;
        Double_t p1x;
        Double_t p0y;
        Double_t p1y;

        // List of branches
        TBranch *b_calcPos2;       //!
        TBranch *b_firedStrip1;    //!
        TBranch *b_weight1;        //!
        TBranch *b_weight1a;       //!
        TBranch *b_weight1b;       //!
        TBranch *b_weight2;        //!
        TBranch *b_weight2a;       //!
        TBranch *b_weight2b;       //!
        TBranch *b_calcPos3;       //!
        TBranch *b_eta;            //!
        TBranch *b_decodeFlag;     //!
        TBranch *b_decodeType;     //!
        TBranch *b_corDev;         //!
        TBranch *b_corEta;         //!
        TBranch *b_evenOdd;        //!
        TBranch *b_tanX;           //!
        TBranch *b_tanY;           //!
        TBranch *b_boardHG;        //!
        TBranch *b_boardLG;        //!
        TBranch *b_boardTDC;       //!
        TBranch *b_tdcTime;        //!
        TBranch *b_lastSeg;        //!
        TBranch *b_nextSeg;        //!
        TBranch *b_interval;       //!
        TBranch *b_processedCount; //!
        TBranch *b_processedCh;    //!
        TBranch *b_processedAmp;   //!
        TBranch *b_processedTDC;   //!
        TBranch *b_derivedPos;     //!
        TBranch *b_fitSigma;       //!
        TBranch *b_v0;             //!
        TBranch *b_v1;             //!
        TBranch *b_v2;             //!
        TBranch *b_v3;             //!
        TBranch *b_vin;            //!
        TBranch *b_vout;           //!
        TBranch *b_vScat;          //!
        TBranch *b_thetaX;         //!
        TBranch *b_thetaY;         //!
        TBranch *b_theta;          //!
        TBranch *b_p0x;            //!
        TBranch *b_p1x;            //!
        TBranch *b_p0y;            //!
        TBranch *b_p1y;            //!

        cor(std::string sFileName);
        cor(TTree *tree = 0);
        virtual ~cor();
        virtual Int_t Cut(Long64_t entry);
        virtual Int_t GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void Init(TTree *tree);
        virtual void Loop();
        virtual Bool_t Notify();
        virtual void Show(Long64_t entry = -1);
    };
}

#endif

#ifdef cor_cxx
namespace john
{

    inline cor::cor(std::string sFileName)
    {
        TTree *tree = NULL;
        TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(sFileName.c_str());
        if (!f || !f->IsOpen())
        {
            f = new TFile(sFileName.c_str());
        }
        f->GetObject("cor", tree);
        Init(tree);
    }
    cor::cor(TTree *tree) : fChain(0)
    {
        // if parameter tree is not specified (or zero), connect the file
        // used to generate this class and read the Tree.
        if (tree == 0)
        {
            TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("1stCorrection.root");
            if (!f || !f->IsOpen())
            {
                f = new TFile("1stCorrection.root");
            }
            f->GetObject("cor", tree);
        }
        Init(tree);
    }

    cor::~cor()
    {
        if (!fChain)
            return;
        delete fChain->GetCurrentFile();
    }

    Int_t cor::GetEntry(Long64_t entry)
    {
        // Read contents of entry.
        if (!fChain)
            return 0;
        return fChain->GetEntry(entry);
    }
    Long64_t cor::LoadTree(Long64_t entry)
    {
        // Set the environment to read one entry
        if (!fChain)
            return -5;
        Long64_t centry = fChain->LoadTree(entry);
        if (centry < 0)
            return centry;
        if (fChain->GetTreeNumber() != fCurrent)
        {
            fCurrent = fChain->GetTreeNumber();
            Notify();
        }
        return centry;
    }

    void cor::Init(TTree *tree)
    {
        // The Init() function is called when the selector needs to initialize
        // a new tree or chain. Typically here the branch addresses and branch
        // pointers of the tree will be set.
        // It is normally not necessary to make changes to the generated
        // code, but the routine can be extended by the user if needed.
        // Init() will be called many times when running on PROOF
        // (once per file to be processed).

        // Set object pointer
        v0 = 0;
        v1 = 0;
        v2 = 0;
        v3 = 0;
        vin = 0;
        vout = 0;
        vScat = 0;
        // Set branch addresses and branch pointers
        if (!tree)
            return;
        fChain = tree;
        fCurrent = -1;
        fChain->SetMakeClass(1);

        fChain->SetBranchAddress("calcPos2", calcPos2, &b_calcPos2);
        fChain->SetBranchAddress("firedStrip1", firedStrip1, &b_firedStrip1);
        fChain->SetBranchAddress("weight1", weight1, &b_weight1);
        fChain->SetBranchAddress("weight1a", weight1a, &b_weight1a);
        fChain->SetBranchAddress("weight1b", weight1b, &b_weight1b);
        fChain->SetBranchAddress("weight2", weight2, &b_weight2);
        fChain->SetBranchAddress("weight2a", weight2a, &b_weight2a);
        fChain->SetBranchAddress("weight2b", weight2b, &b_weight2b);
        fChain->SetBranchAddress("calcPos3", calcPos3, &b_calcPos3);
        fChain->SetBranchAddress("eta", eta, &b_eta);
        fChain->SetBranchAddress("decodeFlag", decodeFlag, &b_decodeFlag);
        fChain->SetBranchAddress("decodeType", decodeType, &b_decodeType);
        fChain->SetBranchAddress("corDev", corDev, &b_corDev);
        fChain->SetBranchAddress("corEta", corEta, &b_corEta);
        fChain->SetBranchAddress("evenOdd", evenOdd, &b_evenOdd);
        fChain->SetBranchAddress("tanX", &tanX, &b_tanX);
        fChain->SetBranchAddress("tanY", &tanY, &b_tanY);
        fChain->SetBranchAddress("boardHG", boardHG, &b_boardHG);
        fChain->SetBranchAddress("boardLG", boardLG, &b_boardLG);
        fChain->SetBranchAddress("boardTDC", boardTDC, &b_boardTDC);
        fChain->SetBranchAddress("tdcTime", tdcTime, &b_tdcTime);
        fChain->SetBranchAddress("lastSeg", lastSeg, &b_lastSeg);
        fChain->SetBranchAddress("nextSeg", nextSeg, &b_nextSeg);
        fChain->SetBranchAddress("interval", interval, &b_interval);
        fChain->SetBranchAddress("processedCount", processedCount, &b_processedCount);
        fChain->SetBranchAddress("processedCh", processedCh, &b_processedCh);
        fChain->SetBranchAddress("processedAmp", processedAmp, &b_processedAmp);
        fChain->SetBranchAddress("processedTDC", processedTDC, &b_processedTDC);
        fChain->SetBranchAddress("derivedPos", derivedPos, &b_derivedPos);
        fChain->SetBranchAddress("fitSigma", fitSigma, &b_fitSigma);
        fChain->SetBranchAddress("v0", &v0, &b_v0);
        fChain->SetBranchAddress("v1", &v1, &b_v1);
        fChain->SetBranchAddress("v2", &v2, &b_v2);
        fChain->SetBranchAddress("v3", &v3, &b_v3);
        fChain->SetBranchAddress("vin", &vin, &b_vin);
        fChain->SetBranchAddress("vout", &vout, &b_vout);
        fChain->SetBranchAddress("vScat", &vScat, &b_vScat);
        fChain->SetBranchAddress("thetaX", &thetaX, &b_thetaX);
        fChain->SetBranchAddress("thetaY", &thetaY, &b_thetaY);
        fChain->SetBranchAddress("theta", &theta, &b_theta);
        fChain->SetBranchAddress("p0x", &p0x, &b_p0x);
        fChain->SetBranchAddress("p1x", &p1x, &b_p1x);
        fChain->SetBranchAddress("p0y", &p0y, &b_p0y);
        fChain->SetBranchAddress("p1y", &p1y, &b_p1y);
        Notify();
    }

    Bool_t cor::Notify()
    {
        // The Notify() function is called when a new file is opened. This
        // can be either for a new TTree in a TChain or when when a new TTree
        // is started when using PROOF. It is normally not necessary to make changes
        // to the generated code, but the routine can be extended by the
        // user if needed. The return value is currently not used.

        return kTRUE;
    }

    void cor::Show(Long64_t entry)
    {
        // Print contents of entry.
        // If entry is not specified, print current entry
        if (!fChain)
            return;
        fChain->Show(entry);
    }
    Int_t cor::Cut(Long64_t entry)
    {
        // This function may be called from Loop.
        // returns  1 if entry is accepted.
        // returns -1 otherwise.
        return 1;
    }
}
#endif // #ifdef cor_cxx
