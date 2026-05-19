void FitHisto(){
    
    // First of all, we need to load the file
    
    TFile *file = new TFile("C:\Users\plele\OneDrive\Desktop\Università\Magistrale\Paprica 2.9\Dataset\run_020.root", "READ");

    // Now we entry to the file 

    file->cd();

    // Now we try to get the histogram 66 for testing fit 
    
    TH1D *h66 = (TH1D *) gDirectory->Get("charge_isipm66;1");

    // Now we implement the fit method

    // First of all we want to select first 5 peaks linked to photoelectrons

    TSpectrum *spec = new TSpectrum(6); // 6 refers to maximum number of peaks attended by function but we use now only the first 5

    // We use Search method to find the peaks, Search method has 4 parameters: the histogram, the sigma of the peaks, the option for searching and the threshold for searching.
    // We set sigma to 2, graphic option  "" (to draw a triangle where the peak is) and threshold for detecting peaks to 0.05 (5% of the maximum value of the histogram)
    
    spec->Search(h66, 2, "", 0.05); 

    // Now we select only the first peak for trying the fit

    double *xPeaks = spec->GetPositionX(); // Get the x positions of the peaks

    std::sort(xPeaks, xPeaks + spec->GetNPeaks()); // Sort the peaks in ascending order so after this xPeaks[0] < xPeaks[1] 

    double xc = xPeaks[0]; // We select the first peak as center for the fit
    double sigma = h66->GetStdDev() / 5; // We set the sigma for the fit as a fraction of the standard deviation of the histogram because having 5 peaks we have a std dev that is higher than the std dev of each peak
    
    // Now we create the fit function, we use a Gaussian function for the fit, the parameters of the function are: the name of the function, the formula of the function, the range of x for the fit and the number of parameters (3 for Gaussian)

    TF1 *fitFunc = new TF1("fitFunc", "gaus", xc - 3*sigma, xc + 3*sigma); // We set the range of the fit to be around the first peak, we choose 3sigma to obtain the 90% range of confidence of the Gaussian distribution
    
    *fitFunc->SetParameters(0, h66->GetBinContent(h66->GetMaximumBin())); // We set the initial parameter of the fit function for the amplitude
    *fitFunc->SetParameters(1, h66->GetBinCenter(h66->GetMaximumBin())); // We set the initial parameter of the fit function for the mean value mu
    *fitFunc->SetParameters(2, h66->GetStdDev()); // We set the initial parameter of the fit function for the std dev

    // Now we impose limits

    fitFunc->SetParLimits(0, -1*h66->GetBinContent(h66->GetMaximumBin()), h66->GetBinContent(h66->GetMaximumBin())); // We set the limit for the amplitude parameter to be between in range -FWHM and FWHM of the histogram
    fitFunc->SetParLimits(1, xc - sigma, xc + sigma); // We set the limit for the mean value parameter to be around the first peak
    fitFunc->SetParLimits(2, 0, 2 * sigma); // We set the limit for the std dev parameter to be between 0 and 2sigma to avoid unphysical values

    // Now we perform the fit

    h66->Fit(fitFunc, "R"); // We use the option "R" to fit only in the range specified in the fit function

    // Now we can draw the histogram and the fit function to visualize the result

    TCanvas *c = new TCanvas("c", "Fit Histo", 800, 600);
    h66 ->SetTitle("Fit of the peaks"); // We set the title of the histogram and the axes
    h66 ->GetXaxis()->SetTitle("Charge [Units]"); // We set the title of the x axis
    h66 ->GetYaxis()->SetTitle("Counts [Units]"); // We set the title of the y axis
    h66 ->SetLineColor(kBlack); // We set the color of the histogram to black
    fitFunc->SetLineColor(kGreen); // We set the color of the fit function to green

    h66->Draw(); // We draw the histogram
    fitFunc->Draw("same"); // We draw the fit function on the same canvas
      
    
    
    /*
    // Now we need to get the histograms from 64 up to 127 from the file

    for(int i = 64; i < 128; i++){
        TH1D *h = (TH1D *) gDirectory->Get()
    } */
}