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
// first of all we load the file and we create the output file where we will save all the results of the fits and the graphs
{
    TFile *fileIn = new TFile("/Users/Utente/Desktop/run_031.root", "READ");
    if (!fileIn || fileIn->IsZombie())
    {
        std::cerr << "ERRORE: Impossibile aprire il file di input run_031.root" << std::endl;
        return;
    }

    TFile *fileOut = new TFile("/Users/Utente/Desktop/SiPM_FitResults.root", "RECREATE");

    // In this part we define the Multigraph, the idea is to create a multigraph to visualize all the gains of the differents simps in a single graph, this will be very useful later for recalibration of the electronics.

    // Here we create the multigrph following the root instructions, since the data can be divided in two categories, the one that follow a good linear trend and the one that do not follow a good linear trend we create two different multigraphs to avoid graphical corruption of the graph with all the data, this is important for visualization because if we put all the data in a single graph it can be very difficult to visualize the results because of the outliers that can affect the scale of the graph and make it difficult to see the details of the good data.
    TMultiGraph *mg_good = new TMultiGraph();
    TMultiGraph *mg_bad = new TMultiGraph();

    TLegend *leg_good = new TLegend(0.1, 0.1, 0.3, 0.9);
    TLegend *leg_bad = new TLegend(0.1, 0.1, 0.3, 0.9);

    std::vector<double> slope; // <-- vector to store the slope(p1) values

    TH1D *h_P1 = new TH1D("h_P1", "Distribution of Gain Slopes (p1)", 100, 100, 400);

    double slopeMin = 220.0; // this values depends on the slopes so when we have a final range we change them, i've fixed the values from a poateriori analysis of thje graph
    double slopeMax = 280.0;

    // Since we want to read all the histograms from the 64 channels we create a loop to read and fit all of the data.
    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm)
    {
        TString histoName = Form("charge_isipm%d", i_sipm); // we create the name of the histogram to read, in this way we can read all the histograms in a loop without writing 64 times the name of the histogram
        TH1D *h = (TH1D *)fileIn->Get(histoName);           // we read the histogram, if the histogram is not found we print a warning and we skip to the next one
        // slope = 0;
        if (!h)
        {
            std::cerr << "Warning: Istogramma '" << histoName << "' non trovato. Salto al prossimo." << std::endl;
            continue;
        }

        h->Rebin(4); // We rebin the histogram to have less bins and a smoother distribution for the fit, this is important because the fit can be very sensitive to the noise in the histogram and with too many bins we can have a lot of noise that can affect the fit
                     // now we define TSpectrum to find the peaks, we set a high number of peaks to find because we want to be sure to find all the peaks, then we will select only the first 5 peaks for the fit because they are the most important for the gain calculation, we need to rember that each peak corresponds to a photon electron
        TSpectrum *spec = new TSpectrum(50);
        int nPeaks = spec->Search(h, 15, " ", 0.007); // the search function has for parameteres, h is the number of histogram, 15 is the sigma of the peaks, "goff" is the option to not draw the peaks on the histogram and 0.001 is the threshold for finding the peaks, this means that we will find only the peaks that are higher than 0.1% of the maximum value of the histogram

        if (nPeaks < 1) // if we find less than 2 peaks we print a warning and we skip to the next histogram because we need at least 2 peaks to calculate the gain
        {
            std::cerr << "Warning: TSpectrum non ha trovato abbastanza picchi in " << histoName << std::endl;
            delete spec;
            continue;
        }
        // double float create an array to store the x positions of the peaks, we need to sort the peaks in ascending order

        double *xPosRaw = spec->GetPositionX();
        std::vector<double> xPos(xPosRaw, xPosRaw + nPeaks);
        std::sort(xPos.begin(), xPos.end());

        std::vector<double> means, meanErrors, sigmas, sigmaErrors, peakNumbers; // we create vectors to store the results of the fits, we will use these vectors to create the graph later

        int peaksToFit = (nPeaks > 5) ? 5 : nPeaks; // if Tsprctrum finds more than 5 peaks we fit only the first 5 because they are the most important for the gain calculation, if we find less than 5 peaks we fit all the peaks found

        for (int p = 0; p < peaksToFit; p++)
        {
            double xc = xPos[p];      // takes the posiion x of the peak
            double range = xc * 0.10; // computes the range of the fit as the 10% of the peak position

            if (range < 5.0)
                range = 20.0;
            if (xc > 1000)
                range = xc * 0.10;

            TString funcName = Form("fit_sipm%d_peak%d", i_sipm, p);          // Tstrings creats a gaussian in Root with a unique name to avoid conflicts between differents fits
            TF1 *fitFunc = new TF1(funcName, "gaus", xc - range, xc + range); // create the gaussian function for the fit, the parameters are the name of the function, the formula of the function, the range of x for the fit.
                                                                              // these three lines estimates the first 3 parameters of the gaussian fit.
            int binCenter = h->FindBin(xc);
            double peakHeight = h->GetBinContent(binCenter);
            double sigma_guess = h->GetStdDev() / 5.0;
            // these are the limits for the fit parameter
            fitFunc->SetParameter(0, peakHeight);
            fitFunc->SetParameter(1, xc);
            fitFunc->SetParameter(2, sigma_guess);
            fitFunc->SetLineColor(kPink + 10); // Set the color of the fit function to a light pink for better visibility

            h->Fit(fitFunc, "RQ0"); // R stands for Range of the defined fit, Q stands for quiet mode to avoid printing the fit results on the console and + stands for to add the fit function to the list of functions of the histogram without deleting the previous ones, this is important because we want to fit all the peaks and we want to keep all the fit functions on the histogram for visualization

            means.push_back(fitFunc->GetParameter(1)); // we save the mean value of the fit
            meanErrors.push_back(fitFunc->GetParError(1));
            sigmas.push_back(fitFunc->GetParameter(2));
            sigmaErrors.push_back(fitFunc->GetParError(2));
            peakNumbers.push_back(p);
        }

        fileOut->cd();
        h->Write(); // we save the histogram with the fit functions in the output file for later visualization

        if (means.size() > 1)
        {
            TString graphName = Form("gain_graph_isipm%d", i_sipm);
            TGraphErrors *graph = new TGraphErrors(means.size(), peakNumbers.data(), means.data(), 0, meanErrors.data());

            graph->SetName(graphName);
            graph->SetTitle(Form("SiPM %d", i_sipm));
            graph->SetMarkerStyle(21); // the 21 is for a marker with a small full star
                                       // this int assignes a different color to each graph based on the SIMP number, this is important for visualization in the multigraph

            int idx = (i_sipm - 64);
            float hue = (idx * 360.0) / 64.0; // we map the SIMP number to a hue value between 0 and 360
            float saturation = 1.0;
            float value = 0.85;
            float r, g, b;
            TColor::HSV2RGB(hue, saturation, value, r, g, b);
            int color = TColor::GetColor(r, g, b); // we set a high saturation for better visibility
            graph->SetMarkerColor(color);
            graph->SetLineColor(color);

            // this part creats the linear fit
            TF1 *lineFit = new TF1(Form("lineFit_sipm%d", i_sipm), "pol1", 0, means.size() - 1);
            lineFit->SetLineColor(color);
            graph->Fit(lineFit, "RQ");
            //

            slope.push_back(lineFit->GetParameter(1)); // we save the slope of the linear fit in the vector to have an idea of the distribution of the gains, this will help us when our data will be taken in a low gain mode, from this modality we expect the majority of the gains function to not follow the linear pattarn resulting in more complex situations once we starte the tracking with the fibers

            std::cout << "SiPM " << i_sipm << " slope = " << slope.back() << std::endl;

            h_P1->Fill(slope.back());

            if (slope.back() >= slopeMin && slope.back() <= slopeMax)
            {
                mg_good->Add(graph);
                leg_good->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");
            }
            else
            {
                mg_bad->Add(graph);
                leg_bad->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");
            }

            fileOut->cd();
            graph->Write();
        }

        delete spec; // cleans tspectrum to avoid memory leakks

    } // end of the cicle for all the histograms

    // Canvas: the fisrt canvas is for the histogram of the gain distributions, while the other two are for the linear fit good vs bad

    TCanvas *c_dist = new TCanvas("c_dist", "Distribution of Gain Slopes (p1)", 800, 600);
    c_dist->cd();
    h_P1->SetLineColor(kBlue);
    h_P1->SetFillColor(kCyan);
    h_P1->Draw("HIST");
    fileOut->cd();
    h_P1->Write("h_P1");

    TCanvas *c_good = new TCanvas("c_good", "SiPM Buoni", 1400, 800);
    mg_good->Draw("AP"); // HIST stands for histogram, this is important to draw the graph with the histogram style, this will help us to visualize the distribution of the points and the fit line better, without this option the graph will be drawn with points and lines only and it can be difficult to visualize the results
    leg_good->Draw();
    fileOut->cd();
    mg_good->Write("MultiGraph_good");
    c_good->Write("Canvas_good");

    TCanvas *c_bad = new TCanvas("c_bad", "SiPM Fuori Range", 1400, 800);
    mg_bad->Draw("AP");
    leg_bad->Draw();
    fileOut->cd();
    mg_bad->Write("MultiGraph_bad");
    c_bad->Write("Canvas_bad");

    // we close the files in a clean way

    fileOut->Close();
    fileIn->Close();

    std::cout << "\n===== ANALISI COMPLETATA =====" << std::endl;
    std::cout << "Tutti i fit e i grafici sono stati salvati nel file: SiPM_FitResults.root" << std::endl;
}
// at this point we want to put all the slopes in a histogram to have an idea of the distribution of the gains, this will help us when our data will be taken in a low gain mode, from this modality we expect the majority of the gains function to not follow the linear pattarn resulting in more complex situations once we  starte the tracking with the fibers