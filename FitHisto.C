#include "TFile.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TDirectory.h"
#include <algorithm>
#include <vector>
#include <stdio.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
void FitHisto()
{

    // First of all, we need to load the file

    TFile *file = new TFile("/Users/Utente/Desktop/run_020.root", "READ");

    // Now we entry to the file

    file->cd();

    TCanvas *c = new TCanvas("c", "Fit Histo", 800, 600);
    c->cd();

    // Now we try to get the histogram 66 for testing fit

    TH1D *h66 = (TH1D *)gDirectory->Get("charge_isipm64;1");
    h66->Rebin(4); // We rebin the histogram to have less bins and a smoother distribution for the fit

    // Now we implement the fit method

    // First of all we want to select first 5 peaks linked to photoelectrons

    TSpectrum *spec = new TSpectrum(5); // 5 refers to maximum number of peaks attended by function but we use now only the first 5

    // We use Search method to find the peaks, Search method has 4 parameters: the histogram, the sigma of the peaks, the option for searching and the threshold for searching.
    // We set sigma to 15, graphic option  "" (to draw a triangle where the peak is) and threshold for detecting peaks to 0.001 (0.1% of the maximum value of the histogram)

    spec->Search(h66, 15, "", 0.001);

    // Now we obtain the peaks and then er order them in ascending order to select the first 5 peaks

    double *xPeaks = spec->GetPositionX(); // Get the x positions of the peaks

    std::sort(xPeaks, xPeaks + spec->GetNPeaks()); // Sort the peaks in ascending order so after this xPeaks[0] < xPeaks[1]

    // now we want to write a cicle for to select the first 5

    h66->SetTitle("Fit of the peaks");           // We set the title of the histogram and the axes
    h66->GetXaxis()->SetTitle("Charge [Units]"); // We set the title of the x axis
    h66->GetYaxis()->SetTitle("Counts [Units]"); // We set the title of the y axis
    h66->SetLineColor(kBlack);                   // We set the color of the histogram to black
    h66->Draw();                                 // We draw the histogram

    // now we declare the arrays to save the results
    std::vector<double> nPE = {0, 1, 2, 3, 4}; // Number of photoelectrons corresponding to the peaks
    std::vector<double> meanValues;            // To store the mean values of the fits
    std::vector<double> sigmaValues;           // To store the sigma values of the fits
    std::vector<double> errMeanValues;         // To store the
    std::vector<double> errsigma;              // To store the errors of the sigma values of the fits
    for (int i = 0; i < 5; i++)
    {
        double xc = xPeaks[i];   // We select the i-th peak as center for the fit
        double range = xc * 0.1; // We set the range of the fit to be 20% of the peak position, this is an arbitrary choice but it should be enough to include the peak and exclude the neighboring peaks
        if (range < 10)
            range = 20; // We set a minimum range of 20 to avoid problems with very low peaks
        if (xc > 1000)
            range = xc * 0.1; // We set a smaller range for very high peaks to avoid including too much background
        // Now we create the fit function, we use a Gaussian function for the fit, the parameters of the function are: the name of the function, the formula of the function, the range\
        of x for the fit and the number of parameters (3 for Gaussian)
        TString fitName = Form("fitFunc%d", i); // We create a unique name for each fit function to avoid conflicts

        TF1 *fitFunc = new TF1(fitName, "gaus", xc - range, xc + range); // We set the range of the fit to be around the first peak

        fitFunc->SetParameter(0, h66->GetBinContent(h66->FindBin(xc))); // We set the initial parameter of the fit function for the amplitude
        fitFunc->SetParameter(1, xc);
        // We set the initial parameter of the fit function for the mean value mu
        fitFunc->SetParameter(2, 5); // We set the initial parameter of the fit function for the std dev

        // Now we perform the fit

        h66->Fit(fitFunc, "R0+"); // We use the option "R0" to fit only in the range specified in the fit function and open only one fit

        // Now we can draw the histogram and the fit function to visualize the result
        fitFunc->SetLineColor(kGreen);                    // We set the color of the fit function to green
        fitFunc->Draw("same");                            // We draw the fit function on the same canvas
        meanValues.push_back(fitFunc->GetParameter(1));   // We save the mean value of the fit in the vector
        sigmaValues.push_back(fitFunc->GetParameter(2));  // We save the sigma value of the fit in the vector
        errsigma.push_back(fitFunc->GetParError(2));      // We save the error of the sigma value of the fit in the vector
        errMeanValues.push_back(fitFunc->GetParError(1)); // We save the error of the mean value of the fit in the vector
    }
    TCanvas *c2 = new TCanvas("c2", "Mean Values of the Fits", 800, 600);
    c2->cd();

    // Now we create a graph to plot the mean values of the fits as a function of the number of photoelectrons

    TGraphErrors *graph = new TGraphErrors(nPE.size(), nPE.data(), meanValues.data(), 0, errMeanValues.data()); // We create a graph with error bars, the x values are the number of photoelectrons, the y values are the mean values of the fits and the error bars are the errors of the mean values

    graph->SetTitle("Mean Values of the Fits;Number of Photoelectrons;Mean Value [Units]"); // We set the title of the graph and the axes
    graph->SetMarkerStyle(21);                                                              // We set the marker style to be a circle
    graph->SetMarkerSize(1);                                                                // We set the marker size to be 1
    graph->SetMarkerColor(kRed);                                                            // We set the marker color to red
    graph->Draw("AP");                                                                      // We draw the graph with axes and points

    // Now we compute the linear fit for the mean values as a function of the number of photoelectrons, this will give us the gain of the SiPM, the gain is the slope of the linear fit
    TF1 *linFit = new TF1("linFit", "pol1", -0.5, 4.5); // pol1 = a + b*x
    linFit->SetLineColor(kBlue);
    linFit->SetLineWidth(2);
    // since the root resolution is not very good we set a high number of points for the fit to have a smoother line
    linFit->SetNpx(100000000);
    graph->Fit(linFit, "R");

    printf("Fit lineare: intercetta = %.2f  pendenza (gain) = %.2f\n",
           linFit->GetParameter(0),  // intercetta
           linFit->GetParameter(1)); // pendenza = gain
}
