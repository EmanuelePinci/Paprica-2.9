void FitHisto(){

    // Usa un percorso relativo. ROOT cercherà il file nella stessa cartella da cui lanci lo script.
    TFile *file = new TFile("run_020.root", "READ");

    // Controllo di sicurezza: il file si è aperto correttamente?
    if (!file || file->IsZombie()) {
        std::cerr << "ERRORE: Impossibile aprire il file root. Controlla il percorso." << std::endl;
        return; 
    }

    // 2. Estraiamo l'istogramma direttamente dal file
    TH1D *h66 = (TH1D *) file->Get("charge_isipm81");
    TCanvas *c = new TCanvas("c", "Fit Histo", 800, 600);
    c->cd();

    // Controllo di sicurezza: l'istogramma esiste con questo nome?
    if (!h66) {
        std::cerr << "ERRORE: Istogramma 'charge_isipm66' non trovato nel file." << std::endl;
        file->Close();
        return;
    }

    // Da qui in poi, siamo sicuri che h66 non è nullo e il warning sparirà.

    
    // 3. Metodo per trovare e gestire i picchi
    
    h66->Rebin(4);
    // Prende il limite destro originale dell'istogramma (così non lo tagli)
    //    double maxX = h66->GetXaxis()->GetXmax();
    //double minX = h66->GetXaxis()->GetXmin();

    
    // Forza il grafico a essere disegnato (e analizzato) solo tra -50 e il massimo
    //h66->GetXaxis()->SetRangeUser(-50, 3000);

    TSpectrum *spec = new TSpectrum(5); // Puoi anche mettere 10 o 20, sarà la soglia a scartare i falsi positivi

    int nPeaks = spec->Search(h66, 15, "", 0.001);

    
    if (nPeaks == 0) {
      std::cerr << "ERRORE: TSpectrum non ha trovato nessun picco con l'attuale soglia." << std::endl;
      return;
    }

    // Estraiamo i risultati e li mettiamo in un vettore C++ standard per eleganza e sicurezza
    double *xPosRaw = spec->GetPositionX();

    // Ordiniamo i picchi da sinistra a destra ---
    std::vector<double> xPos(xPosRaw, xPosRaw + nPeaks);
    std::sort(xPos.begin(), xPos.end());

    // 4. Disegniamo il risultato
    
    h66->SetTitle("Fit of the peaks");
    h66->GetXaxis()->SetTitle("Charge [Units]");
    h66->GetYaxis()->SetTitle("Counts [Units]");
    h66->SetLineColor(kBlue);
    h66->Draw();

    //Dichiaro gli array per salvare i risulati
    std::vector<double> means, meanErrors, sigmas, sigmaErrors;
    std::vector<double> peakNumbers = {0, 1, 2, 3, 4}; // I numeri dei picchi (asse X)

    // Protezione: fittiamo fino a 5 picchi, o meno se ne ha trovati di meno
    int maxPeaks = std::min(nPeaks, 5);
    
    for(int i=0;i<5;i++)
      {
	double xc = xPos[i];

	double range = xc *0.10;
	if (range < 5.0) range = 20.0;
	if (xc > 1000) range = xc *0.10; 

	// 1. Creiamo un nome univoco per ogni funzione usando Form()
	TString funcName = Form("fit_%d", i);

	TF1 *fitFunc = new TF1(funcName.Data(), "gaus", xc - range, xc + range);



	fitFunc->SetParameters(0, h66->GetBinContent(h66->GetMaximumBin()));
	fitFunc->SetParameters(1, h66->GetBinCenter(h66->GetMaximumBin()));
	fitFunc->SetParameters(2, h66->GetStdDev());
	
    
	h66->Fit(fitFunc, "R+");
	fitFunc->SetLineColor(kRed);

	fitFunc->Draw("same");

	// 2. Salviamo i parametri estratti dal fit
	means.push_back(fitFunc->GetParameter(1));      // Media
	meanErrors.push_back(fitFunc->GetParError(1));  // Errore sulla media
	sigmas.push_back(fitFunc->GetParameter(2));     // Sigma
	sigmaErrors.push_back(fitFunc->GetParError(2)); // Errore sulla sigma
      }

    // Usiamo i vettori per popolare il grafico
    TCanvas *c2 = new TCanvas("c2", "Analisi Picchi", 800, 600);
    c2->cd();
    TGraphErrors *graph = new TGraphErrors(5, peakNumbers.data(), means.data(), 0, meanErrors.data());

    graph->SetTitle("Media vs Numero Picco;Numero Picco;Media Carica [Units]");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->Draw("AP"); // AP = Axis e Points

    // Il range è da 0 a 4 (i tuoi indici dei picchi)
    TF1 *lineFit = new TF1("lineFit", "pol1", 0, 4);

    graph->Fit(lineFit, "RQ");

    // 3. Estraiamo i parametri per visualizzarli (opzionale)
    double intercept = lineFit->GetParameter(0); // Offset (pedestal)
    double slope = lineFit->GetParameter(1);     // Guadagno (Guadagno del SiPM)

    std::cout << "Risultati del Fit Lineare:" << std::endl;
    std::cout << "Pedestal (p0): " << intercept << std::endl;
    std::cout << "Guadagno (p1): " << slope << " ADC channels/p.e." << std::endl;
    
    lineFit->Draw("same");
}
