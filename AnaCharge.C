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

void AnaCharge(const vector<TString>& myfile, const string& Constraint = "") {

    // 0. We set the path where the file is and we open the output file 

    TString mypath = "Dataset/Risultati_Cosmici/";
    
    // We make a control check to see if the folders exist

    if(gSystem->AccessPathName(mypath.Data())){
        cerr << "Error: The specified folder path does not exist: " << mypath << endl;
        return;
    }

    // ---------------------------------------------------------------
    // 1. Theoretical estimate of MPV togliere diagonale
    // ---------------------------------------------------------------
    //mettere valore atteso e fare la banda più larga sfruttando radice di due
    double dEdx = 1.9;      // MeV/cm, muon MIP in plastic
    double Ly = 8000;       // photons/MeV, scintillating fiber light yield
    double QE = 0.35;       // SiPM quantum efficiency
    double eps_trap = 0.07/2.; // fiber trapping efficiency divided by two

    // Two geometrical cases: vertical (Δx = 1.1 cm) and diagonal (Δx * sqrt(2))
    double dx_vert = 1.0;                 // cm  (2.2cm / 2, average vertical path)

    double NPHE = dEdx * dx_vert * Ly * QE * eps_trap;

    // Relative uncertainties (added in quadrature)
    double rel_dEdx = 0.1 / dEdx;
    double rel_Ly = 0.10; // 10%
    double rel_QE = 0.015 / QE;
    double rel_trap = 0.50; // 10%

    double rel_err = sqrt(rel_dEdx * rel_dEdx + rel_Ly * rel_Ly + rel_QE * rel_QE 
                          + rel_trap * rel_trap);
    double NPHE_err = NPHE * rel_err;

    // Print theoretical calculation
    cout << "\n========================================" << endl;
    cout << "  Theoretical MPV estimate " << endl;
    cout << "========================================" << endl;
    cout << "  dE/dx        = " << dEdx << " MeV/cm" << endl;
    cout << "  Light Yield  = " << Ly << " photons/MeV" << endl;
    cout << "  QE           = " << QE << endl;
    cout << "  eps_trap     = " << eps_trap << endl;
    cout << "----------------------------------------" << endl;
    cout << "  Mean +/- err   :  " << NPHE << " +/- " << NPHE_err << " PHE" << endl;
    cout << "  Relative error :  " << rel_err * 100 << " %" << endl;
    cout << "========================================\n"
         << endl;

    // Booking of the histos of total charge (1D), of RMS (2D) and total charge over threshold (1D) 

    TH1D *htc = new TH1D("hTotalCharge", "Total Charge Distribution;Total Charge [PHE];Occurrence", 200, 0, 1000);
    TH1D *hthr = new TH1D("hTotalCharge_THR", "Total Charge Distribution with threshold;Total Charge [PHE];Occurrence", 200, 0, 1000);

    TH2F *hRMS = new TH2F("hRMS", "Channel Correlation: RMS_{x} vs RMS_{y} ; RMS_{x}; RMS_{y}", 350, 0.0, 3.5, 350, 0.0, 3.5);
  
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
                    }/*else{
                        // ========================== DRAWING HISTO =======================
                        TCanvas* cDisplay = new TCanvas("cDisplay", Form("EventID:%d",eventID), 800, 600);
                        h2D->Draw("COLZ text");
                        cDisplay->Update();

                        //cout << "Press [ENTER] to continue" << endl;
                        //string stop;
                        //getline(cin, stop);

                        delete cDisplay;
                    } IF YOU WANT VISUALIZE THE EVENT DECOMMENT THIS PART*/ 

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

    gStyle->SetOptStat(1111); // Entries, Mean, RMS
    gStyle->SetOptFit(1111);  // Fit parameters (Chi2, MPV, Sigma)

    double fit_min = htc->GetBinCenter(htc->GetMaximumBin()) - htc->GetRMS()*0.7;
        double fit_max = htc->GetBinCenter(htc->GetMaximumBin()) + 4.0 * htc->GetRMS();
        TF1 *fLandau_1 = new TF1("fLandau", "landau", fit_min, fit_max);
        fLandau_1->SetLineColor(kRed + 1);
        fLandau_1->SetLineWidth(2);
        htc->Fit(fLandau_1, "RQ");
        double MPV_measured = fLandau_1->GetParameter(1);
        double MPV_measured_err = fLandau_1->GetParError(1);

        cout << "----------------------------------------" << endl;
        cout << "  Landau fit result" << endl;
        cout << "  MPV measured = " << MPV_measured << " +/- " << MPV_measured_err << " PHE" << endl;
        cout << "----------------------------------------\n"
            << endl;

    // ====================== CANVAS DRAWING =====================================
    ctc->cd();
    htc->Draw("hist"); 
    htc->Draw("hist same"); 
    fLandau_1->Draw("same");
    ctc->Update();    

    cRMS->cd();
    hRMS->Draw("COLZ");


    

    if(useConstraint){
        TCanvas* cthr = new TCanvas("cthr", "Distribution of total charge with threshold", 800, 600);
        
        // =================== LANDAU FIT =========================
        gStyle->SetOptStat(1111); // Entries, Mean, RMS
        gStyle->SetOptFit(1111);  // Fit parameters (Chi2, MPV, Sigma)

        double fit_min = hthr->GetBinCenter(hthr->GetMaximumBin()) - hthr->GetRMS();
        double fit_max = hthr->GetBinCenter(hthr->GetMaximumBin()) + 4.5 * hthr->GetRMS();
        TF1 *fLandau = new TF1("fLandau", "landau", fit_min, fit_max);
        fLandau->SetLineColor(kRed + 1);
        fLandau->SetLineWidth(2);
        hthr->Fit(fLandau, "RQ");
        double MPV_measured = fLandau->GetParameter(1);
        double MPV_measured_err = fLandau->GetParError(1);

        cout << "----------------------------------------" << endl;
        cout << "  Landau fit result" << endl;
        cout << "  MPV measured = " << MPV_measured << " +/- " << MPV_measured_err << " PHE" << endl;
        cout << "----------------------------------------\n"
            << endl;
        
        // We draw a box which represents the range MPV +- resolution where for resolution we mean 
        // the Landau fit error
        double y_max = hthr->GetMaximum();
        TBox *band = new TBox(NPHE - NPHE_err, 0, NPHE + NPHE_err, y_max);
        band->SetFillStyle(1001); 
        band->SetFillColorAlpha(kGreen + 1, 0.20); 
        band->SetLineColor(kGreen + 3);
        band->SetLineWidth(1);

        // Legend
        TLegend *leg = new TLegend(0.58, 0.35, 0.98, 0.65);
        leg->SetBorderSize(1);
        leg->AddEntry(htc, Form("Data (%zu runs)", myfile.size()), "l"); // metti automatico
        leg->AddEntry(fLandau, "Landau fit", "l");
        leg->AddEntry(band, Form("Theory: [%.0f, %.0f] PHE", NPHE - NPHE_err, 
                                                             NPHE + NPHE_err), "f");
        

        cthr->cd();

        // 1. Disegna l'istogramma con l'opzione "hist" per creare assi e griglia
        hthr->Draw("hist"); 

        // 2. Forza lo stile pieno sul box e imposta l'alfa (0.15 = molto trasparente)
        band->SetFillStyle(1001); 
        band->SetFillColorAlpha(kGreen + 1, 0.15); 
        band->Draw("same");

        // 3. CRUCIALE: Ridisegna l'istogramma SOPRA la banda usando "same" 
        // per far riemergere le linee blu dei bin
        hthr->Draw("hist same"); 

        // 4. Ridisegna il fit per metterlo in primissimo piano
        fLandau->Draw("same");

        // 5. Spostiamo la legenda a sinistra per non sovrapporla allo stat bo
        leg->Draw("same");

        cthr->Update();      
        cthr->SaveAs(mypath + "Spettro_Carica_Totale_Threshold.png");
    }  
    
    // We save canvas as a picture
    ctc->SaveAs(mypath + "Spettro_Carica_Totale_FitLandau.png");
    cRMS->SaveAs(mypath + "Correlazione_Canali.png");
    
}