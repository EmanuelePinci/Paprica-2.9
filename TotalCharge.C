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

void TotalCharge(const vector<TString>& myfile, const string& Constraint = "") {

    // 0. We set the path where the file is and we open the output file 

    TString mypath = "Dataset/Risultati_Cosmici/";
    
    // We make a control check to see if the folders exist

    if(gSystem->AccessPathName(mypath.Data())){
        cerr << "Error: The specified folder path does not exist: " << mypath << endl;
        return;
    }

    // Booking of the histo of total charge (1D) and of RMS (2D) 

    TH1D *htc = new TH1D("hTotalCharge", "Total Charge Distribution;Total Charge [ADC counts];Occurrence", 200, 0, 150e3);
    TH1D *hthr = new TH1D("hTotalCharge_THR", "Total Charge Distribution with threshold;Total Charge [ADC counts];Occurrence", 200, 0, 150e3);

    TH2F *hRMS = new TH2F("hRMS", "Channel Correlation: RMS_{x} vs RMS_{y} ; RMS_{x}; RMS_{y}", 100, 0.0, 3.5, 100, 0.0, 3.5);
  
    // We want also to see the total charge (as the sum of each charge) collected by all SiPM in 
    // that event and the RMS of the event 
        
    int eventID = 0;

    double ChargeTot = 0.0;
    double charge_min = 0.0;
    double charge_max = 0.0;

    double RMS_x = .0, RMS_y = .0;

    // We inizialize the constraint formula for selecting histos later

    bool useConstraint = (!Constraint.empty());
    TFormula cutFormula;

    if (useConstraint) {
        cutFormula = TFormula("cutFormula", Constraint.c_str());
        
        if (cutFormula.Compile() != 0) {
            cerr << "Error: Non valid formula!" << endl;
            return;
        }
    }

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
        
        while(fileInput >> eventID >> ChargeTot >> RMS_x >> RMS_y){
            
            // Filling histo of total charge and of the correlation
            htc->Fill(ChargeTot);
            hRMS->Fill(RMS_x, RMS_y);

            // Filling histo of total charge with threshold
            // We also implement an automatic way to open the Heat Map 
            // of the charge deposited in each activated SiPM in a certain 
            // event. This is for event displaying

            if(useConstraint && cutFormula.Eval(RMS_x, RMS_y) != 0){
                hthr->Fill(ChargeTot);
                TString rootFileName = fileName;
                
                // ======================== OPENING ROOT FILE ==========================
                // We modify momentaneously the name of the myfile's element                
                
                rootFileName.ReplaceAll("TotalCharge_", "");
                rootFileName.ReplaceAll(".txt", ".root");

                TFile* fRoot = TFile::Open(rootFileName.Data());
                if(!fRoot || fRoot->IsZombie()){
                    cerr << "Error: Impossible to open ROOT file: " << rootFileName << endl;
                }else{
                    // ========================= OPENING HISTO ==========================
                    TH2F *h2D = (TH2F*)fRoot->Get(Form("hMap%d",eventID));
                    if(!h2D){
                        cerr << "Warning: Histo " << Form("hMap%d",eventID) << "not found or not existing!" << endl;
                    }else{
                        // ========================== DRAWING HISTO =======================
                        TCanvas* cDisplay = new TCanvas("cDisplay", Form("EventID:%d",eventID), 800, 600);
                        h2D->Draw("COLZ text");
                        cDisplay->Update();

                        cout << "Press [ENTER] to continue" << endl;
                        string stop;
                        getline(cin, stop);

                        delete cDisplay;
                    }
                }
                
                fRoot->Close();
                delete fRoot;
            }
            
            if(charge_min ==0 || charge_min > ChargeTot){
                charge_min = ChargeTot;
            } else if(charge_max <  ChargeTot){
                charge_max = ChargeTot;
            }    
        }
        
        fileInput.close();
    } 

    // Canvas creation and drawing for both the histos
    TCanvas* ctc = new TCanvas("ctc", "Distribution of total charge", 800, 600);
    TCanvas* cRMS = new TCanvas("cRMS", "Channel Correlation", 800, 600);
    /*
    // =================== LANDAU FIT =========================
    gStyle->SetOptStat(1111); // Entries, Mean, RMS
    gStyle->SetOptFit(1111);  // Fit parameters (Chi2, MPV, Sigma)

    // We execute the Landau fit 
    // "S": to save fit results into TFitResultPtr object
    // "Q": "Quiet" mode (optional)
    htc->Fit("landau", "S", "", charge_min, charge_max);

    // Changing fit line color
    TF1 *myfit = htc->GetFunction("landau");
    if(myfit) {
        myfit->SetLineColor(kBlue);
        myfit->SetLineWidth(3);
    }
*/
    // ====================== CANVAS DRAWING =====================================
    ctc->cd();
    htc->Draw();

    cRMS->cd();
    hRMS->Draw("COLZ");

    if(useConstraint){
        TCanvas* cthr = new TCanvas("cthr", "Distribution of total charge with threshold", 800, 600);
        cthr->cd();
        hthr->Draw();
        cthr->SaveAs(mypath + "Spettro_Carica_Totale_Threshold.png");
    }   
    
    // We save canvas as a picture
    ctc->SaveAs(mypath + "Spettro_Carica_Totale_FitLandau.png");
    cRMS->SaveAs(mypath + "Correlazione_Canali.png");
    
}