#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>
#include <iostream>
#include <fstream> // <-- fstream is a c++ library that allows to read and write files

void ExtractGains(int mydac, std::string numrun){ // <-- mydac = DAC value for the run numrun

    /* We print the value of the DAC and the run number to check that we are analyzing 
       the correct file */

    cout << mydac << "\t" << numrun << endl; 

    // PAY ATTENTION: file path refers to my own laptop you must to change it 
    TString fileInName = Form("/home/emanuele/Universita/Paprica 2.9/Dataset/Fit_Data/SiPM_Graphs_run_%s.root", numrun.c_str());
    TFile *fileIn = new TFile(fileInName, "UPDATE");
    
    if (!fileIn || fileIn->IsZombie()) {
        std::cerr << "ERRORE: Impossibile aprire " << fileInName << std::endl;
        return;
    }

    /* We create a file in write mode where we will save all the results of the fits for 
       the parameter p1 (the gain). In particular in file there will be three columns:
       1. The name of the SiPM (from 64 to 127)
       2. The value of the DAC 
       3. The value of the gain (p1) for each SiPM */

    // We create a string with the name of the file 

    string myfilename = "/home/emanuele/Universita/Paprica 2.9/Dataset/Gain_Data/iSiPM_run_" 
                         + numrun + "_" + std::to_string(mydac) + "_Gain.txt";

    /* We open the file in write mode
       Next line defines the file where we will save the results of the fits for the gain (p1)
       for each SiPM */

    std::ofstream fileFitResult(myfilename); 
    
    if (!fileFitResult.is_open()){
        std::cerr << "ERROR: It's impossible to open the output file iSiPM_run_" + numrun + "_" 
                     + std::to_string(mydac) + "_Gain.txt" << std::endl;
        return;
    } 

    // The first line of data file for a more easy reading

    fileFitResult << "//SiPM index" << "\t" << "DAC" << "\t" << "p1 values" << "\t" 
                  << "p1 errors" << std::endl;

    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm) {
        TString graphName = Form("gain_graph_isipm%d", i_sipm);
        TGraphErrors *graph = (TGraphErrors*)fileIn->Get(graphName);

        if (!graph) {
            // If we don't find any peak we save in file a default value of 0 for gain and error

            fileFitResult << i_sipm << "\t" << mydac << "\t" << 0.0 << "\t" << 0.0 << std::endl;

            continue;
        }

        TF1 *lineFit = new TF1(Form("lineFit_sipm%d", i_sipm), "pol1", 0, graph->GetN() - 1);
        lineFit->SetLineColor(graph->GetLineColor());
        graph->Fit(lineFit, "RQ");

        /* Now we save in a variable the gain and the pedestal for each SiPM
           this will be useful later for the recalibration of the electronics, 
           we save these values in a root file. */
        
        double p0 = lineFit->GetParameter(0); // <-- Pedestal
        double p1 = lineFit->GetParameter(1); // <-- Gain
        double p1_error = lineFit->GetParError(1); // <-- Error on the gain 

        /* Now we put the values in file iSiPM_DAC_Gain.txt, we save the name of the SiPM, 
           the value of the DAC and the value of the gain (p1) for each SiPM */

        fileFitResult << i_sipm << "\t" << mydac << "\t" << p1 << "\t" << p1_error << std::endl;
        
        fileIn->cd();
        graph->Write("", TObject::kOverwrite);
    }

    fileFitResult.close();
    fileIn->Close();
}