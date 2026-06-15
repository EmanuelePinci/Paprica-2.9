#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <string>

using namespace std;

void EventDisplay(TString myfile, int eventNumber, double soglia_ADC) {
    
    // 0. We set the path where the file is

    TString mypath = "Dataset/Run_Cosmici/";
    TString fileName = mypath + myfile;

    // 1. First of all we initialize the mapping. 
    //    The mapping is initialized as electronics views the channels

    int fisX[128];
    int fisY[128];
    for(int i=0; i<128; i++) { 
        fisX[i] = -1; 
        fisY[i] = -1; 
    }

   // --- Trascription of the map (X = ROW, Y = COL) ---
    
    // 1st ROW on top
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


    // 2. We open the file and the TTree
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cout << "Errore: Impossibile aprire il file " << fileName << endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Tuple"); // TTree is named "Tuple"
    if (!tree) {
        cout << "Errore: TTree 'Tuple' non trovato!" << endl;
        return;
    }

    // 3. Each branch is a vector containig the values, so we initialize them 
    vector<int> *iSiPM = nullptr;
    vector<int> *ChargeHG = nullptr;
    
    tree->SetBranchAddress("iSiPM", &iSiPM);
    tree->SetBranchAddress("ChargeHG", &ChargeHG);

    // 4. Booking of the TH2F histo (8x8 pixel)
    // For centering the bin we set the range from -0.5 to 7.5
    TH2F *hMap = new TH2F("hMap", Form("Event Display (Spazio Fisico) - Evento %d; Colonna (X) ; Riga (Y)", eventNumber), 
                          8, -0.5, 7.5, 8, -0.5, 7.5);

// 5. Data extraction
    int bytes = tree->GetEntry(eventNumber); // We set the event from the tree
    if (bytes <= 0) {
        cout << "Errore: L'evento " << eventNumber << " non esiste nel file o e' vuoto!" << endl;
        return;
    }

    // We check if the size of iSiPM and ChargeHG is the same. We check only for HG (High Gain)
    // because we don't use LG (Low Gain)
    if (iSiPM->size() != ChargeHG->size()) {
        cout << "Attenzione: I vettori di ChargeHG e iSiPM hanno dimensioni diverse!" << endl;
        return;
    }

    // We fill the histogram 
    for (size_t i = 0; i < iSiPM->size(); i++) {
        int id = iSiPM->at(i);
        float carica_hg = ChargeHG->at(i); 
        
        // We check if the channel is in the matrix (so its index will not be -1)
        if (id >= 0 && id < 128 && fisX[id] != -1) {
            // We want to see only the channel over the threshold
            if (carica_hg > soglia_ADC) {
                int x_reale = fisX[id];
                int y_reale = fisY[id];
                
                hMap->Fill(x_reale, y_reale, carica_hg);
            }
        }
    }

    // 6. We send on screen the result
    TCanvas *c1 = new TCanvas("c1", "Analisi SiPM", 800, 600);
    c1->SetRightMargin(0.15); // Space for color bar
    gStyle->SetPaintTextFormat(".0f"); // Setting number format

    // "colz" --> to draw the color, "text" --> to write the value
    hMap->Draw("colz text");
    
    c1->Update();
}