#include <TFile.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

void FitPeaks(std::string numrun) {

    // We verify for security that the number of run is correct printing in on console
    
    cout << numrun << endl;

    /*  First of all we load the file and we create the output file where we will save all
        the results of the fits and the graphs
        PAY ATTENTION: file path refers to my own laptop you must to change it */ 
    
    TString fileInName = Form("/home/emanuele/Universita/Paprica 2.9/Dataset/Run_Data/run_%s.root", numrun.c_str());
    TFile *fileIn = new TFile(fileInName, "READ");
    if (!fileIn || fileIn->IsZombie())
    {
        std::cerr << "ERROR: It's impossible to open the input file run_" << numrun << ".root" << std::endl;
        return;
    }

    TString fileOutName = Form("/home/emanuele/Universita/Paprica 2.9/Dataset/Fit_Data/SiPM_Graphs_run_%s.root", numrun.c_str());
    TFile *fileOut = new TFile(fileOutName, "RECREATE");

    /*  Since we want to read all the histograms from the 64 channels we create a loop 
        to read and fit all of the data. */

    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm)
    {
        /* We create the name of the histogram to read, in this way we can read all the 
           histograms in a loop without writing 64 times the name of the histogram*/

        TString histoName = Form("charge_isipm%d", i_sipm);

        /* We read the histogram, if the histogram is not found we print a warning 
           and we skip to the next one*/

        TH1D *h = (TH1D *)fileIn->Get(histoName); 

        if (!h)
        {
            std::cerr << "Warning: Histogram '" << histoName << "' not found. Skipping to the next one." << std::endl;
            continue;
        }

        /* We rebin the histogram to have less bins and a smoother distribution for the fit, 
           this is important because the fit can be very sensitive to the noise in the histogram
           and with too many bins we can have a lot of noise that can affect the fit now */

        h->Rebin(4);
        TSpectrum *spec = new TSpectrum(50);
        int nPeaks = spec->Search(h, 15, "", 0.007);

        if (nPeaks < 2)
        {
            std::cerr << "Warning: TSpectrum doesn't find enough peaks in " << histoName << std::endl;
            delete spec;
            continue;
        }

        double *xPosRaw = spec->GetPositionX();
        std::vector<double> xPos(xPosRaw, xPosRaw + nPeaks);
        std::sort(xPos.begin(), xPos.end());

        std::vector<double> means, meanErrors, peakNumbers;
        int peaksToFit = (nPeaks > 5) ? 5 : nPeaks;

        for (int p = 0; p < peaksToFit; p++)
        {
            double xc = xPos[p];
            double range = xc * 0.10;
            if (range < 5.0)  range = 20.0;
            if (xc > 1000)    range = xc * 0.10;

            TF1 *fitFunc = new TF1(Form("fit_sipm%d_peak%d", i_sipm, p), "gaus", xc - range, xc + range);
            fitFunc->SetParameters(h->GetBinContent(h->FindBin(xc)), xc, h->GetStdDev() / 5.0);
            
            h->Fit(fitFunc, "RQ+");

            means.push_back(fitFunc->GetParameter(1));
            meanErrors.push_back(fitFunc->GetParError(1));
            peakNumbers.push_back(p);
            fileOut->cd();
            h->Write("", TObject::kOverwrite);
        }

        if (means.size() > 1)
        {
            TGraphErrors *graph = new TGraphErrors(means.size(), peakNumbers.data(), means.data(), 0, meanErrors.data());
            graph->SetName(Form("gain_graph_isipm%d", i_sipm));
            graph->SetTitle(Form("SiPM %d", i_sipm));
            
            //int color = (i_sipm - 64) % 64 + 1;
            // We set different colors from HSV map

            int idx = (i_sipm -64);
            float hue = (idx * 360.0) / 64.0; // Map index to hue (0-360)
            float saturation = 1.0;
            float value = 0.85;
            float r, g, b;
            TColor::HSV2RGB(hue, saturation, value, r, g, b);
            int color = TColor::GetColor(r, g, b);
            graph->SetMarkerColor(color);
            graph->SetLineColor(color);
            graph->SetMarkerStyle(21);

            fileOut->cd();
            graph->Write();
        }
        delete spec;
    }

    fileOut->Close();
    fileIn->Close();
}