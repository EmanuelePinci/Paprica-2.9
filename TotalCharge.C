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

void TotalCharge(const vector<TString>& myfile) {

    // 0. We set the path where the file is and we open the output file 

    TString mypath = "Dataset/Risultati_Cosmici/";
    
    // We make a control check to see if the folders exist

    if(gSystem->AccessPathName(mypath.Data())){
        cerr << "Error: The specified folder path does not exist: " << mypath << endl;
        return;
    }

    // Booking of the histo of total charge 

    TH1D *htc = new TH1D("hTotalCharge", "Total Charge Distribution;Total Charge [ADC counts];Occurrence", 100, 0, 0);

  // We want also to see the total charge (as the sum of each charge)collected by all SiPM in 
        // that event 
        
    double ChargeTot = 0.0; 
    
    for(size_t file=0; file<myfile.size(); file++){
        
        TString fileName = mypath + myfile[file];

        // 2. We open the data files for histo

        ifstream fileInput(fileName.Data());
        if (!fileInput.is_open()) {
            cerr << "Error: Impossible to open the file " << fileName << endl;
            return;
        }
        
        // We read always the first line because it's for intestation
        
        string header;
        getline(fileInput, header);

        // 4. Data extraction and filling of the histo

        // Now we cycle on all file's line and fill the total charge histo
        
        while(fileInput >> ChargeTot){
            
            // Filling histo of total charge
            htc->Fill(ChargeTot);
        
        }
        
        fileInput.close();
    } 

    TCanvas* ctc = new TCanvas("ctc", "Distribution of total charge", 800, 600);
    gStyle->SetOptStat(1111); // Mostra box statistiche completo
    htc->Draw();
    
    // We save the result as image
    ctc->SaveAs(mypath + "Total_Charge.png");
    
}