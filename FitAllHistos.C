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

void FitAllHistos()
//first of all we load the file and we create the output file where we will save all the results of the fits and the graphs
{
    TFile *fileIn = new TFile("/Users/Utente/Desktop/run_020.root", "READ");
    if (!fileIn || fileIn->IsZombie())
    {
        std::cerr << "ERRORE: Impossibile aprire il file di input run_020.root" << std::endl;
        return;
    }

    TFile *fileOut = new TFile("/Users/Utente/Desktop/SiPM_FitResults.root", "RECREATE");

    // In this part we define the Multigraph, the idea is to create a multigraph to visualize all the gains of the differents simps in a single graph, this will be very useful later for recalibration of the electronics.

    TMultiGraph *mg = new TMultiGraph(); // Here we create the multigrph following the root instructions
    mg->SetTitle("Media vs Numero Picco - Tutti i SiPM;Numero Picco;Media Carica [Units]"); // we set the title
    TLegend *leg = new TLegend(0.1, 0.1, 0.3, 0.9); // legend for the graph

    // Since we want to read all the histograms from the 64 channels we create a loop to read and fit all of the data.
    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm)
    {
        TString histoName = Form("charge_isipm%d", i_sipm); // we create the name of the histogram to read, in this way we can read all the histograms in a loop without writing 64 times the name of the histogram
        TH1D *h = (TH1D *)fileIn->Get(histoName); // we read the histogram, if the histogram is not found we print a warning and we skip to the next one

        if (!h)
        {
            std::cerr << "Warning: Istogramma '" << histoName << "' non trovato. Salto al prossimo." << std::endl;
            continue;
        }

        h->Rebin(4);// We rebin the histogram to have less bins and a smoother distribution for the fit, this is important because the fit can be very sensitive to the noise in the histogram and with too many bins we can have a lot of noise that can affect the fit
// now 
        TSpectrum *spec = new TSpectrum(50);
        int nPeaks = spec->Search(h, 15, "goff", 0.001);

        if (nPeaks < 2)
        {
            std::cerr << "Warning: TSpectrum non ha trovato abbastanza picchi in " << histoName << std::endl;
            delete spec;
            continue;
        }

        double *xPosRaw = spec->GetPositionX();
        std::vector<double> xPos(xPosRaw, xPosRaw + nPeaks);
        std::sort(xPos.begin(), xPos.end());

        std::vector<double> means, meanErrors, sigmas, sigmaErrors, peakNumbers;

        int peaksToFit = (nPeaks > 5) ? 5 : nPeaks;

        for (int p = 0; p < peaksToFit; p++)
        {
            double xc = xPos[p];
            double range = xc * 0.10;

            if (range < 5.0)  range = 20.0;
            if (xc > 1000)    range = xc * 0.10;

            TString funcName = Form("fit_sipm%d_peak%d", i_sipm, p);
            TF1 *fitFunc = new TF1(funcName, "gaus", xc - range, xc + range);

            int binCenter = h->FindBin(xc);
            double peakHeight = h->GetBinContent(binCenter);
            double sigma_guess = h->GetStdDev() / 5.0;

            fitFunc->SetParameter(0, peakHeight);
            fitFunc->SetParameter(1, xc);
            fitFunc->SetParameter(2, sigma_guess);
            fitFunc->SetLineColor(kOrange + 1);

            h->Fit(fitFunc, "RQ+");

            means.push_back(fitFunc->GetParameter(1));
            meanErrors.push_back(fitFunc->GetParError(1));
            sigmas.push_back(fitFunc->GetParameter(2));
            sigmaErrors.push_back(fitFunc->GetParError(2));
            peakNumbers.push_back(p);
        }

        fileOut->cd();
        h->Write();

        if (means.size() > 1)
        {
            TString graphName = Form("gain_graph_isipm%d", i_sipm);
            TGraphErrors *graph = new TGraphErrors(means.size(), peakNumbers.data(), means.data(), 0, meanErrors.data());

            graph->SetName(graphName);
            graph->SetTitle(Form("SiPM %d", i_sipm));
            graph->SetMarkerStyle(21);

            int color = (i_sipm - 64) % 64 + 1;
            graph->SetMarkerColor(color);
            graph->SetLineColor(color);

            TF1 *lineFit = new TF1(Form("lineFit_sipm%d", i_sipm), "pol1", 0, means.size() - 1);
            lineFit->SetLineColor(color);
            graph->Fit(lineFit, "RQ0");

            mg->Add(graph);
            leg->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");

            fileOut->cd();
            graph->Write();
        }

        delete spec; // ← era mancante

    } // Fine ciclo SiPM ← parentesi graffa del for

    // Canvas e MultiGraph FUORI dal ciclo
    TCanvas *c = new TCanvas("c_all", "Tutti i SiPM", 1400, 800);
    mg->Draw("AP");
    leg->Draw();
    fileOut->cd();
    mg->Write("MultiGraph_tutti_sipm");
    c->Write("Canvas_tutti_sipm");

    // Chiusura file FUORI dal ciclo
    fileOut->Close();
    fileIn->Close();

    std::cout << "\n===== ANALISI COMPLETATA =====" << std::endl;
    std::cout << "Tutti i fit e i grafici sono stati salvati nel file: SiPM_FitResults.root" << std::endl;
}