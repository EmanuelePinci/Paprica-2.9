#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <string>
#include <TROOT.h>

using namespace std;

// Aggiungiamo la costante per riscalare la relazione
//La soglia la lasciamo in adc counts xche applicata esternamente
const double GAIN_ADC_PER_PE = 300.0;

void AllEventDisplay(TString myfile, double soglia_ADC) {

    gROOT->SetBatch(kTRUE);

    // 0. Paths
    TString mypath  = "/Users/Utente/Desktop/PAPRICA/run_RAGGICOSMICI/";
    TString outpath = "/Users/Utente/Desktop/PAPRICA/run_RAGGICOSMICI/Risultati/";

    // Creo la cartella perchè non la ho
    if (gSystem->AccessPathName(outpath.Data())) {
        gSystem->MakeDirectory(outpath.Data());
        cout << "Created output folder: " << outpath << endl;
    }

    TString fileName = mypath + myfile;
    TString outfile  = outpath + myfile;

    // .txt file for total charge (in PE) — used by TotalCharge
    TString nomeFileTxt = outpath + Form("TotalCharge_%s.txt", myfile.ReplaceAll(".root", "").Data());
    ofstream outFileTxt(nomeFileTxt.Data());

    // 1. Channel mapping (X = ROW, Y = COL)
    int fisX[128], fisY[128];
    for(int i = 0; i < 128; i++) { fisX[i] = -1; fisY[i] = -1; }

    // 1st ROW
    fisX[94]=6; fisY[94]=3;   fisX[95]=7; fisY[95]=3;   fisX[64]=0; fisY[64]=0;   fisX[65]=1; fisY[65]=0;
    fisX[126]=6; fisY[126]=7; fisX[127]=7; fisY[127]=7; fisX[96]=0; fisY[96]=4;   fisX[97]=1; fisY[97]=4;
    // 2nd ROW
    fisX[93]=5; fisY[93]=3;   fisX[92]=4; fisY[92]=3;   fisX[67]=3; fisY[67]=0;   fisX[66]=2; fisY[66]=0;
    fisX[125]=5; fisY[125]=7; fisX[124]=4; fisY[124]=7; fisX[99]=3; fisY[99]=4;   fisX[98]=2; fisY[98]=4;
    // 3rd ROW
    fisX[90]=2; fisY[90]=3;   fisX[91]=3; fisY[91]=3;   fisX[68]=4; fisY[68]=0;   fisX[69]=5; fisY[69]=0;
    fisX[122]=2; fisY[122]=7; fisX[123]=3; fisY[123]=7; fisX[100]=4; fisY[100]=4; fisX[101]=5; fisY[101]=4;
    // 4th ROW
    fisX[89]=1; fisY[89]=3;   fisX[88]=0; fisY[88]=3;   fisX[71]=7; fisY[71]=0;   fisX[70]=6; fisY[70]=0;
    fisX[121]=1; fisY[121]=7; fisX[120]=0; fisY[120]=7; fisX[103]=7; fisY[103]=4; fisX[102]=6; fisY[102]=4;
    // 5th ROW
    fisX[86]=6; fisY[86]=2;   fisX[87]=7; fisY[87]=2;   fisX[72]=0; fisY[72]=1;   fisX[73]=1; fisY[73]=1;
    fisX[118]=6; fisY[118]=6; fisX[119]=7; fisY[119]=6; fisX[104]=0; fisY[104]=5; fisX[105]=1; fisY[105]=5;
    // 6th ROW
    fisX[85]=5; fisY[85]=2;   fisX[84]=4; fisY[84]=2;   fisX[75]=3; fisY[75]=1;   fisX[74]=2; fisY[74]=1;
    fisX[117]=5; fisY[117]=6; fisX[116]=4; fisY[116]=6; fisX[107]=3; fisY[107]=5; fisX[106]=2; fisY[106]=5;
    // 7th ROW
    fisX[82]=2; fisY[82]=2;   fisX[83]=3; fisY[83]=2;   fisX[76]=4; fisY[76]=1;   fisX[77]=5; fisY[77]=1;
    fisX[114]=2; fisY[114]=6; fisX[115]=3; fisY[115]=6; fisX[108]=4; fisY[108]=5; fisX[109]=5; fisY[109]=5;
    // 8th ROW
    fisX[81]=1; fisY[81]=2;   fisX[80]=0; fisY[80]=2;   fisX[79]=7; fisY[79]=1;   fisX[78]=6; fisY[78]=1;
    fisX[113]=1; fisY[113]=6; fisX[112]=0; fisY[112]=6; fisX[111]=7; fisY[111]=5; fisX[110]=6; fisY[110]=5;

    // 2. Open input ROOT file
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cout << "Error: Impossible to open file " << fileName << endl;
        return;
    }
    TTree *tree = (TTree*)file->Get("Tuple");
    if (!tree) {
        cout << "Error: TTree 'Tuple' not found!" << endl;
        return;
    }

    TFile* outFile = new TFile(outfile.Data(), "RECREATE");

    // 3. Branch initialization
    vector<int>    *iSiPM   = nullptr;
    vector<double> *ChargeHG = nullptr;
    double total_charge_PE   = 0.0;
    int nhit = 0, eventNumber = 0, bytes = 0;

    tree->SetBranchAddress("iSiPM",    &iSiPM);
    tree->SetBranchAddress("ChargeHG", &ChargeHG);
    tree->SetBranchAddress("nhit",     &nhit);

    // 4. Event loop
    gStyle->SetPaintTextFormat(".1f");
    outFileTxt << "//Total Charge [PE]" << endl;

    for(int i = 0; i < tree->GetEntries(); i++){

        bytes = tree->GetEntry(i);
        if (bytes <= 0) {
            cout << "Error: Empty or not existing event!" << endl;
            continue;
        } else if (iSiPM->size() != ChargeHG->size()) {
            cout << "Warning: ChargeHG and iSiPM size mismatch!" << endl;
            continue;
        }

        eventNumber = i;

        if(nhit > 6 && nhit < 50){

            TH2F *hMap = new TH2F(Form("hMap%d", eventNumber),
                                  Form("Event Display - Event %d [PE]", eventNumber),
                                  8, -0.5, 7.5, 8, -0.5, 7.5);

            for (size_t j = 0; j < iSiPM->size(); j++) {
                int id = iSiPM->at(j);
                double carica_hg = ChargeHG->at(j);

                if (id >= 0 && id < 128 && fisX[id] != -1) {
                    if (carica_hg > soglia_ADC) {
                        // ADC -> PE conversion
                        double carica_PE = carica_hg / GAIN_ADC_PER_PE;
                        total_charge_PE += carica_PE;
                        hMap->Fill(fisX[id], fisY[id], carica_PE);
                    }
                }
            }

            outFileTxt << total_charge_PE << endl;
            total_charge_PE = 0.0;

            hMap->SetMarkerSize(1.5);
            hMap->SetOption("colz text");
            hMap->GetZaxis()->SetTitle("Charge [PE]");
            hMap->Draw("colz text");
            hMap->Write();
            delete hMap;
        }
    }

    outFile->Close();
    outFileTxt.close();

    cout << "Processed " << tree->GetEntries() << " events. Output saved in " << outpath << endl;
    cout << "Gain used for ADC->PE conversion: " << GAIN_ADC_PER_PE << " ADC/PE" << endl;

    gROOT->SetBatch(kFALSE);
}

// --- Wrapper: processes all cosmic ray runs in one call ---
void RunAllEventDisplay(double soglia_ADC = 50.0) {
    vector<TString> runs = {
        "run_127.root",
        "run_130.root",
        "run_131.root",
        "run_132.root",
        "run_133.root",
        "run_134.root",
        "run_135.root",
        "run_136.root",
        "run_137.root"
    };
    for (const auto& run : runs) {
        cout << "\n=== Processing " << run << " ===" << endl;
        AllEventDisplay(run, soglia_ADC);
    }
    cout << "\nAll runs processed." << endl;
}
