//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  6 11:07:52 2023 by ROOT version 6.22/02
// from TTree tsTree/Time stamp for all boards
// found on file: TS.root
//////////////////////////////////////////////////////////

#ifndef tsTree_h
#define tsTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class tsTree
{
public:
    TTree *fChain;  //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    UInt_t tsid;
    UChar_t tsBoardCount;
    UChar_t tsBoardStat[32]; //[tsBoardCount]
    Double_t ts[32];
    UChar_t tsFlag[32];

    // List of branches
    TBranch *b_tsid;         //!
    TBranch *b_tsBoardCount; //!
    TBranch *b_tsBoardStat;  //!
    TBranch *b_ts;           //!
    TBranch *b_tsFlag;       //!

    tsTree(std::string sFileName);
    tsTree(TTree *tree = 0);
    virtual ~tsTree();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop();
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef tsTree_cxx
inline tsTree::tsTree(std::string sFileName)
{
    TTree *tree = NULL;
    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(sFileName.c_str());
    if (!f || !f->IsOpen())
    {
        f = new TFile(sFileName.c_str());
    }
    f->GetObject("tsTree", tree);
    Init(tree);
}
tsTree::tsTree(TTree *tree) : fChain(0)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0)
    {
        TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("TS.root");
        if (!f || !f->IsOpen())
        {
            f = new TFile("TS.root");
        }
        f->GetObject("tsTree", tree);
    }
    Init(tree);
}

tsTree::~tsTree()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t tsTree::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t tsTree::LoadTree(Long64_t entry)
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

void tsTree::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("tsid", &tsid, &b_tsid);
    fChain->SetBranchAddress("tsBoardCount", &tsBoardCount, &b_tsBoardCount);
    fChain->SetBranchAddress("tsBoardStat", tsBoardStat, &b_tsBoardStat);
    fChain->SetBranchAddress("ts", ts, &b_ts);
    fChain->SetBranchAddress("tsFlag", tsFlag, &b_tsFlag);
    Notify();
}

Bool_t tsTree::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void tsTree::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}
Int_t tsTree::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef tsTree_cxx
