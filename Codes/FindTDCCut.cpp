void FindTDCCut()
{
    // Open file in ../Sum/GetPos.root
    TFile *file = new TFile("../Sum/GetPos.root");
    TTree *data = (TTree *)file->Get("data");

    auto c = new TCanvas("c", "c", 1);
    for (int board = 0; board < 8; board++)
    {
        data->Draw(Form("boardLG[%d][Iteration$]:boardTDC[%d][Iteration$]-tdcTime[%d]>>(100,-20,80,100,0,60000)", board, board, board), Form("boardTDC[%d]>0&&boardTDC[%d]!=tdcTime[%d]", board, board, board), "colz");
        auto sPNG = "TDCCut" + to_string(board) + ".png";
        c->SaveAs(sPNG.c_str());
    }
}