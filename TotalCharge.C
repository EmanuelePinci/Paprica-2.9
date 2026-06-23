#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>

using namespace std;

void TotalCharge(const vector<TString>& myfile) {

    TString mypath = "/Users/Utente/Desktop/PAPRICA/run_RAGGICOSMICI/Risultati/";

    if(gSystem->AccessPathName(mypath.Data())){
        cerr << "Error: Output folder not found: " << mypath << endl;
        cerr << "Run AllEventDisplay first to generate the .txt files." << endl;
        return;
    }

    // Histogram in PE units
    TH1D *htc = new TH1D("hTotalCharge",
                          "Total Charge Distribution;Total Charge [PE];Occurrence",
                          100, 0, 0);

    double ChargeTot_PE = 0.0;

    for(size_t f = 0; f < myfile.size(); f++){
        TString fileName = mypath + myfile[f];
        ifstream fileInput(fileName.Data());
        if (!fileInput.is_open()) {
            cerr << "Error: Cannot open " << fileName << endl;
            continue; // skip missing file, don't abort all runs
        }
        // Skip header line
        string header;
        getline(fileInput, header);

        while(fileInput >> ChargeTot_PE){
            htc->Fill(ChargeTot_PE);
        }
        fileInput.close();
        cout << "Loaded: " << fileName << endl;
    }

    TCanvas* ctc = new TCanvas("ctc", "Distribution of total charge [PE]", 800, 600);
    gStyle->SetOptStat(1111);
    htc->Draw();
    ctc->SaveAs(mypath + "Total_Charge_PE.png");
    cout << "Saved: " << mypath << "Total_Charge_PE.png" << endl;
}

// --- Wrapper: loads all cosmic ray runs in one call ---
void RunTotalCharge() {
    vector<TString> files = {
        "TotalCharge_run_127.txt",
        "TotalCharge_run_130.txt",
        "TotalCharge_run_131.txt",
        "TotalCharge_run_132.txt",
        "TotalCharge_run_133.txt",
        "TotalCharge_run_134.txt",
        "TotalCharge_run_135.txt",
        "TotalCharge_run_136.txt",
        "TotalCharge_run_137.txt"
    };
    TotalCharge(files);
}
