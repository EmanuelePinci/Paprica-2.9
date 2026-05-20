#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

void FitAllHistos() {

    // 1. Apriamo il file di input
    TFile *fileIn = new TFile("run_020.root", "READ");
    if (!fileIn || fileIn->IsZombie()) {
        std::cerr << "ERRORE: Impossibile aprire il file di input run_020.root" << std::endl;
        return; 
    }

    // 2. Creiamo il file di output per salvare i risultati
    TFile *fileOut = new TFile("SiPM_FitResults.root", "RECREATE");

    // 3. Ciclo su tutti gli SiPM da 64 a 127
    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm) {
        
        // Generiamo il nome dell'istogramma dinamicamente
        TString histoName = Form("charge_isipm%d", i_sipm);
        TH1D *h = (TH1D*) fileIn->Get(histoName);

        if (!h) {
            std::cerr << "Warning: Istogramma '" << histoName << "' non trovato. Salto al prossimo." << std::endl;
            continue;
        }

        h->Rebin(4);
        
        // Cerca i picchi
        TSpectrum *spec = new TSpectrum(10); 
        int nPeaks = spec->Search(h, 15, "goff", 0.001); // "goff" evita che TSpectrum disegni i marker rossi

        if (nPeaks < 2) {
            std::cerr << "Warning: TSpectrum non ha trovato abbastanza picchi in " << histoName << std::endl;
            delete spec;
            continue;
        }

        // Estraiamo le posizioni dei picchi e LE ORDINIAMO IN MODO CRESCENTE sull'asse X
        double *xPosRaw = spec->GetPositionX();
        std::vector<double> xPos(xPosRaw, xPosRaw + nPeaks);
        std::sort(xPos.begin(), xPos.end());

        // Dichiariamo i vettori per salvare i risultati del fit
        std::vector<double> means, meanErrors, sigmas, sigmaErrors, peakNumbers;

        // Limita a 5 picchi se ne ha trovati di più, oppure fitta quelli che ha trovato
        int peaksToFit = (nPeaks > 5) ? 5 : nPeaks;

        for(int p = 0; p < peaksToFit; p++) {
            double xc = xPos[p];
            double range = xc * 0.10;

            if (range < 5.0) range = 20.0;
            if (xc > 1000) range = xc * 0.10; 

            // Nome univoco della funzione per questo SiPM e questo picco
            TString funcName = Form("fit_sipm%d_peak%d", i_sipm, p);
            TF1 *fitFunc = new TF1(funcName, "gaus", xc - range, xc + range);

            // Parametri iniziali corretti e dinamici
            int binCenter = h->FindBin(xc);
            double peakHeight = h->GetBinContent(binCenter);
            double sigma_guess = h->GetStdDev() / 5.0; // Stima iniziale più piccola

            fitFunc->SetParameter(0, peakHeight);
            fitFunc->SetParameter(1, xc);
            fitFunc->SetParameter(2, sigma_guess);
            fitFunc->SetLineColor(kGreen);
            
            // Fit senza stampare a schermo il canvas ("0"), silenzioso ("Q"), range ("R")
            // e aggiunge la funzione all'istogramma ("+")
            h->Fit(fitFunc, "RQ0+");

            means.push_back(fitFunc->GetParameter(1));
            meanErrors.push_back(fitFunc->GetParError(1));
            sigmas.push_back(fitFunc->GetParameter(2));
            sigmaErrors.push_back(fitFunc->GetParError(2));
            peakNumbers.push_back(p);
        }

        // 4. Scriviamo l'istogramma con i fit gaussiani salvati al suo interno nel file di output
        fileOut->cd();
        h->Write();

        // 5. Creiamo e fittiamo il grafico Media vs Numero Picco
        if (means.size() > 1) {
            TString graphName = Form("gain_graph_isipm%d", i_sipm);
            TGraphErrors *graph = new TGraphErrors(means.size(), peakNumbers.data(), means.data(), 0, meanErrors.data());
            
            graph->SetName(graphName);
            graph->SetTitle(Form("Media vs Numero Picco - SiPM %d;Numero Picco;Media Carica [Units]", i_sipm));
            graph->SetMarkerStyle(21);
            graph->SetMarkerColor(kBlue);

            TF1 *lineFit = new TF1(Form("lineFit_sipm%d", i_sipm), "pol1", 0, means.size() - 1);
            lineFit->SetLineColor(kRed);
            
            // Fit lineare, anche qui silenzioso
            graph->Fit(lineFit, "RQ0");

            // Scriviamo il grafico nel file di output
            graph->Write();
        }

        delete spec;
    } // Fine ciclo SiPM

    // Chiusura file
    fileOut->Close();
    fileIn->Close();

    std::cout << "\n===== ANALISI COMPLETATA =====" << std::endl;
    std::cout << "Tutti i fit e i grafici sono stati salvati nel file: SiPM_FitResults.root" << std::endl;
}
