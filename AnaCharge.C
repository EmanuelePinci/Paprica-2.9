#include <TF1Convolution.h>
#include <TLatex.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <iostream>
#include <string>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>

using namespace std;

// Funzione di convoluzione Landau-Gauss (LanGau) ultra-stabile
Double_t langaufun(Double_t *x, Double_t *par) {
    // I parametri del fit sono:
    // par[0] = Landau sigma (larghezza intrinseca della Landau)
    // par[1] = Landau MPV (il picco reale della fisica)
    // par[2] = Area totale dell'istogramma (normalizzazione)
    // par[3] = Gauss sigma (lo smearing/risoluzione del rivelatore)
    
    Double_t InvSqrt2Pi = 0.39894228040143;
    Double_t mpshift  = -0.22278298;

    Double_t np = 100.0;      // Step di integrazione numerica
    Double_t sc =   5.0;      // Estensione dell'integrazione in unità di sigma gaussiana

    Double_t xx = x[0];
    Double_t mpc = par[1] - mpshift * par[0];

    Double_t xlow = xx - sc * par[3];
    Double_t xupp = xx + sc * par[3];
    Double_t step = (xupp - xlow) / np;

    Double_t sum = 0.0;
    for(Double_t i=1.0; i<=np/2; i++) {
        Double_t xxp = xlow + (i-0.5) * step;
        Double_t fland = TMath::Landau(xxp, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(xx, xxp, par[3]);

        xxp = xupp - (i-0.5) * step;
        fland = TMath::Landau(xxp, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(xx, xxp, par[3]);
    }

    return (par[2] * step * sum * InvSqrt2Pi / par[3]);
}

void AnaCharge(const vector<TString>& myfile, const string& Constraint = "") {

    // 0. We set the path where the file is and we open the output file 

    TString mypath = "/Users/simonenardi/Desktop/Paprica2026/run_cosmici/risultati_cosmici/";
    
    // We make a control check to see if the folders exist

    if(gSystem->AccessPathName(mypath.Data())){
        cerr << "Error: The specified folder path does not exist: " << mypath << endl;
        return;
    }

    // ---------------------------------------------------------------
    // 1. Theoretical estimate of MPV (Geometric path correction)
    // ---------------------------------------------------------------
    double dEdx = 1.9;         // MeV/cm, muon MIP in plastic
    double Ly = 8000;          // photons/MeV, scintillating fiber light yield
    double QE = 0.35;          // SiPM quantum efficiency
    double eps_trap = 0.07/2.; // fiber trapping efficiency divided by two

    // Path: vertical average is 1.0 cm. Maximum diagonal is 1.0 * sqrt(2)
    double dx_vert = 1.0;                 
    double NPHE = dEdx * dx_vert * Ly * QE * eps_trap;

    // Relative uncertainties
    double rel_dEdx = 0.1 / dEdx;
    double rel_Ly = 0.10; 
    double rel_QE = 0.015 / QE;
    double rel_trap = 0.10; 

    // Aggiungiamo l'effetto geometrico del cammino (fino a sqrt(2) ~ 41% in quadratura)
    double rel_geom = 0.10; 
    
    double rel_err = sqrt(rel_dEdx * rel_dEdx + rel_Ly * rel_Ly + rel_QE * rel_QE 
                          + rel_trap * rel_trap + rel_geom * rel_geom);
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

    TH1D *htc = new TH1D("hTotalCharge", "Total Charge Distribution;Total Charge [PHE];Occurrence", 200, 0, 800);
    TH1D *hthr = new TH1D("hTotalCharge_THR", "Total Charge Distribution with threshold;Total Charge [PHE];Occurrence", 120, 0, 600);

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

    double max_bin_center_1 = hthr->GetBinCenter(hthr->GetMaximumBin()); // Trova il picco (~190 PHE)                                           
    double integral_1 = hthr->Integral() * hthr->GetBinWidth(1);   // Area reale sotto la curva                                           

    // Impostiamo il range del fit visivo (dal grafico si vede che la struttura è tra 130 e 350)                                               
    double fit_min_1 = 130.0;
    double fit_max_1 = 350.0;

    // Creiamo il TF1 puntando alla funzione definita sopra (4 parametri)
    TF1 *fLangau_1 = new TF1("fLangau_1", langaufun, fit_min_1, fit_max_1, 4);
    fLangau_1->SetParNames("LandauSigma", "LandauMPV", "Area", "GaussSigma");

    // VALORI INIZIALI GUIDATI (Fondamentali per MINUIT)
    fLangau_1->SetParameter(0, 10.0);           // Stima iniziale Landau Sigma                                                                   
    fLangau_1->SetParameter(1, max_bin_center_1);  // Posiziona il MPV esattamente sul picco
    fLangau_1->SetParameter(2, integral_1);        // Normalizzazione basata sull'area reale                                                       
    fLangau_1->SetParameter(3, 11.0);           // Risoluzione gaussiana stimata dallo screenshot

    // Limiti di sicurezza per impedire ai parametri di convergere a valori assurdi o negativi
    fLangau_1->SetParLimits(0, 2.0, 25.0);
    fLangau_1->SetParLimits(1, 175.0, 205.0);
    fLangau_1->SetParLimits(3, 3.0, 25.0);

    fLangau_1->SetLineColor(kRed + 1);
    fLangau_1->SetLineWidth(2);

    // Eseguiamo il fit con l'opzione "R" (Range) e "B" (usa i limiti impostati)
    htc->Fit(fLangau_1, "RB");

    // Ora gli errori saranno perfettamente calcolati!
    double MPV_measured_1 = fLangau_1->GetParameter(1);
    double MPV_measured_err_1 = fLangau_1->GetParError(1);

    cout << "----------------------------------------" << endl;
    cout << "  Langau Fit Concluso Con Successo!" << endl;
    cout << "  MPV misurato = " << MPV_measured_1 << " +/- " << MPV_measured_err_1 << " PHE" << endl;
    cout << "----------------------------------------\n" << endl;

    // =================== LINEA VERTICALE PER MPV FIT =========================
    double mpv_fit_1 = fLangau_1->GetParameter(1); 
    double y_max_1 = htc->GetMaximum() * 1.05;
    double mpv_visivo_1 = fLangau_1->GetMaximumX();   // Il picco visivo (Effetto strumento)

    

    // Creiamo la linea: TLine(x1, y1, x2, y2)
    TLine *lineMPV_1 = new TLine(mpv_fit_1, 0, mpv_fit_1, y_max_1);

    // Stile della linea
    lineMPV_1->SetLineColor(kRed + 2); // Rosso scuro, coerente con il fit ma distinguibile
    lineMPV_1->SetLineStyle(2);        // Codice 2 = Tratteggiato (dashed)
    lineMPV_1->SetLineWidth(2);        // Spessore della linea

    // 2. Linea per il Picco Apparente (Strumentale) - Blu puntinata
    TLine *lineVisivo = new TLine(mpv_visivo_1, 0, mpv_visivo_1, y_max_1);
    lineVisivo->SetLineColor(kBlue + 1);
    lineVisivo->SetLineStyle(3); // 3 = Puntinata
    lineVisivo->SetLineWidth(2);



    // ====================== CANVAS DRAWING =====================================
    ctc->cd();
    htc->Draw("hist");
    // Aggiunta testo Bin e Bin Width
    TLatex lat1;
    lat1.SetNDC(); 
    lat1.SetTextSize(0.035);
    lat1.DrawLatex(0.15, 0.85, Form("N Bins: %d", htc->GetNbinsX()));
    lat1.DrawLatex(0.15, 0.80, Form("Bin Width: %.2f PHE", htc->GetBinWidth(1)));
    htc->Draw("hist same"); 
    fLangau_1->Draw("same");
    lineVisivo->Draw("same");   // Retta strumentale
    lineMPV_1->Draw("same");      // Linea verticale tratteggiata dell'MPV
    ctc->Update();    

    cRMS->cd();
    hRMS->Draw("COLZ");
    // Aggiunta testo Bin e Bin Width (2D)
    TLatex lat2;
    lat2.SetNDC();
    lat2.SetTextSize(0.035);
    lat2.DrawLatex(0.15, 0.85, Form("N Bins (X,Y): %d, %d", hRMS->GetNbinsX(), hRMS->GetNbinsY()));
    lat2.DrawLatex(0.15, 0.80, Form("Bin Width (X,Y): %.2f, %.2f", hRMS->GetXaxis()->GetBinWidth(1), hRMS->GetYaxis()->GetBinWidth(1)));
    
    cRMS->Update();


    

    if(useConstraint){
        TCanvas* cthr = new TCanvas("cthr", "Distribution of total charge with threshold", 800, 600);

	// =================== FIT LANGAU ULTRA-STABILE =========================
        gStyle->SetOptStat(1111); 
        gStyle->SetOptFit(1111);  

        double max_bin_center  = hthr->GetBinCenter(hthr->GetMaximumBin()); // Trova il picco (~190 PHE)
        double integral        = hthr->Integral() * hthr->GetBinWidth(1);   // Area reale sotto la curva

        // Impostiamo il range del fit visivo (dal grafico si vede che la struttura è tra 130 e 350)
        double fit_min = 130.0;
        double fit_max = 350.0;

        // Creiamo il TF1 puntando alla funzione definita sopra (4 parametri)
        TF1 *fLangau = new TF1("fLangau", langaufun, fit_min, fit_max, 4);
        fLangau->SetParNames("LandauSigma", "LandauMPV", "Area", "GaussSigma");

        // VALORI INIZIALI GUIDATI (Fondamentali per MINUIT)
        fLangau->SetParameter(0, 10.0);           // Stima iniziale Landau Sigma
        fLangau->SetParameter(1, max_bin_center);  // Posiziona il MPV esattamente sul picco
        fLangau->SetParameter(2, integral);        // Normalizzazione basata sull'area reale
        fLangau->SetParameter(3, 11.0);           // Risoluzione gaussiana stimata dallo screenshot

        // Limiti di sicurezza per impedire ai parametri di convergere a valori assurdi o negativi
        fLangau->SetParLimits(0, 2.0, 25.0);
        fLangau->SetParLimits(1, 175.0, 205.0);
        fLangau->SetParLimits(3, 3.0, 25.0);

        fLangau->SetLineColor(kRed + 1);
        fLangau->SetLineWidth(2);

        // Eseguiamo il fit con l'opzione "R" (Range) e "B" (usa i limiti impostati)
        hthr->Fit(fLangau, "RB");

        // Ora gli errori saranno perfettamente calcolati!
        double MPV_measured = fLangau->GetParameter(1);
        double MPV_measured_err = fLangau->GetParError(1);
        
        cout << "----------------------------------------" << endl;
        cout << "  Langau Fit Concluso Con Successo!" << endl;
        cout << "  MPV misurato = " << MPV_measured << " +/- " << MPV_measured_err << " PHE" << endl;
        cout << "----------------------------------------\n" << endl;

        
        // =================== PUNTO + ERRORE TEORICO =========================
        double y_max = hthr->GetMaximum();
        double y_pos = y_max * 0.5; // Posizioniamo il punto a metà altezza del grafico per visibilità
        
        // Creiamo un grafico con 1 punto: X = NPHE, Y = y_pos, Errore X = NPHE_err, Errore Y = 0
        TGraphErrors *gTheory = new TGraphErrors(1);
        gTheory->SetPoint(0, NPHE, y_pos);
        gTheory->SetPointError(0, NPHE_err, 0.0);
        
        // Stile del punto e della barra di errore
        gTheory->SetMarkerStyle(20);        // Cerchio pieno
        gTheory->SetMarkerSize(1.2);       // Dimensione punto
        gTheory->SetMarkerColor(kGreen+2);  // Colore del punto
        gTheory->SetLineColor(kGreen+2);    // Colore della barra d'errore
        gTheory->SetLineWidth(2);           // Spessore barra

	// =================== LINEA VERTICALE PER MPV FIT =========================
        double mpv_fit = fLangau->GetParameter(1); // Estrae il valore del picco trovato (184.2)
        double y_max_2 = hthr->GetMaximum() * 1.05;  // Prende l'altezza massima dell'istogramma + 5% di margine

        // Creiamo la linea: TLine(x1, y1, x2, y2)
        TLine *lineMPV = new TLine(mpv_fit, 0, mpv_fit, y_max_2);
        
        // Stile della linea
        lineMPV->SetLineColor(kRed + 2); // Rosso scuro, coerente con il fit ma distinguibile
        lineMPV->SetLineStyle(2);        // Codice 2 = Tratteggiato (dashed)
        lineMPV->SetLineWidth(2);        // Spessore della linea

	// =================== LEGENDA OTTIMIZZATA =========================
        // Coordinate ottimizzate per posizionarla bene sotto il box statistico
	TLegend *leg = new TLegend(0.60, 0.20, 0.90, 0.35);
        
        leg->SetBorderSize(0);            // Rimuove il bordo rigido per un look più pulito
        leg->SetFillStyle(0);             // Sfondo trasparente (non copre le linee del grafico)
        leg->SetTextFont(42);             // Font standard moderno di ROOT (Helvetica)
        leg->SetTextSize(0.035);          // Dimensione del testo in armonia con lo stat box

        // Aggiungiamo l'entry con i valori aggiornati, usando "#pm" per il simbolo ±
        leg->AddEntry(gTheory, Form("Theory MPV: %.1f #pm %.1f PHE", NPHE, NPHE_err), "pe");

        // Consiglio: se in futuro vorrai aggiungere anche i dati e il fit, ti basterà decommentare queste:
        // leg->AddEntry(hthr, Form("Data (%zu runs)", myfile.size()), "l"); 
        leg->AddEntry(fLangau, "Langau Fit", "l");
	// Aggiungiamo anche la linea nella legenda per renderla perfetta
        leg->AddEntry(lineMPV, Form("Fit MPV: %.1f PHE", mpv_fit), "l");
	
        cthr->cd();

        // 1. Disegna l'istogramma dati
        hthr->Draw("hist"); 

        // 2. Disegna il fit di Landau
        fLangau->Draw("same");

        // 3. Disegna il punto con la barra d'errore (opzione "P" obbligatoria per i TGraph)
        gTheory->Draw("P same");

	// 4. Linea verticale tratteggiata dell'MPV   
	lineMPV->Draw("same");
	
        // 5. Disegna la legenda
	leg->Draw("same");

	// Aggiunta testo Bin e Bin Width
        TLatex lat3;
        lat3.SetNDC();
        lat3.SetTextSize(0.035);
        lat3.DrawLatex(0.15, 0.85, Form("N Bins: %d", hthr->GetNbinsX()));
        lat3.DrawLatex(0.15, 0.80, Form("Bin Width: %.2f PHE", hthr->GetBinWidth(1)));

        cthr->Update();      
        cthr->SaveAs(mypath + "Spettro_Carica_Totale_Threshold.png");
    }  
    
    // We save canvas as a picture
    ctc->SaveAs(mypath + "Spettro_Carica_Totale_FitLandau.png");
    cRMS->SaveAs(mypath + "Correlazione_Canali.png");
    
}
