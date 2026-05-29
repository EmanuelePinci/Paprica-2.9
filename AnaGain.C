#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <vector>
#include <sstream>

using namespace std;

void AnaGain(const vector<string>& skip_files) {
    
    // 1. Local directory configuration
    // PAY ATTENTION: Update this path to the location of your Gain_Data folder in your laptop
    
    string folder_path = "/home/emanuele/Universita/Paprica 2.9/Dataset/Gain_Data";

    // We make a control check to see if the folder exists
    if (gSystem->AccessPathName(folder_path.c_str())){
        cerr << "Error: The specified folder path does not exist: " << folder_path << endl;
        return;
    }

    // 2. Creation of the output .root file and initialization of TTree

    string out_file_path = folder_path + "/GainResults.root"; 
    TFile* outFile = new TFile(out_file_path.c_str(), "RECREATE");
    TTree* tree = new TTree("GainTree", "Tree containing DAC and Gain data");
    int channel;
    double dac, gain, gain_err;

    // In following lines char I and D indicate the type of the variable (Integer and Double respectively)
    tree->Branch("channel", &channel, "channel/I");
    tree->Branch("dac", &dac, "dac/D");
    tree->Branch("gain", &gain, "gain/D");
    tree->Branch("gain_err", &gain_err, "gain_err/D");

    // 3. Initialization of a C++ vector of TGraphErrors
    vector<TGraphErrors*> gr;
    TGraphErrors* gr_temp; // We also define an ausiliary TGraphErrors pointer to be used in the loop for filling the vector
    for(int i_sipm = 64; i_sipm <= 127; i_sipm++){
        gr_temp = new TGraphErrors();
        gr_temp->SetName(Form("Graph_ch%d", i_sipm));
        gr_temp->SetTitle(Form("Gain vs DAC - Channel %d; DAC Value; Gain", i_sipm));
        gr_temp->SetMarkerStyle(21);     
        gr_temp->SetMarkerSize(1.2);
        gr_temp->SetMarkerColor(kBlue);
        gr_temp->SetLineColor(kBlack);
        gr.push_back(gr_temp);
    } 

    cout << "\n--- Starting Data Analysis ---\n" << endl;

    // 4. Open the directory and loop over all .txt files

    void* dirp = gSystem->OpenDirectory(folder_path.c_str()); // 
    if(!dirp){
        cerr << "Error: Unable to open directory: " << folder_path << endl;
        return;
    }

    const char* entry;
    int file_counter = 0; // To keep track of how many files we processed
    while((entry = gSystem->GetDirEntry(dirp))){ // <-- This allows us to investigate all files in the folder
        TString fileName = entry;

        // We check if the file has a .txt extension file before processing it
        if(fileName.EndsWith(".txt")){
            string current_file = string(fileName.Data());
            if(find(skip_files.begin(), skip_files.end(), current_file) != skip_files.end()){
                continue;
            }else{
                file_counter++; // We increment file counter if file exists
                string full_path = folder_path + "/" + fileName.Data(); // We create the full path to the file 
                ifstream infile(full_path);
                if(infile.is_open()){
                    cout << "Processing file: " << full_path << endl;
                    string line;
                    int line_number = 0; // To keep track of line numbers for error reporting
                    // We skip first line of the file since it contains the header
                    if(getline(infile, line)){
                        line_number++;
                    }
                    // We analyze the remaining lines of the file 
                    while(getline(infile, line)){
                        line_number++;
                        if(line.empty()) continue; // Skip empty lines
                        istringstream iss(line); // This is used to parse the line into its components and to handle potential formatting issues
                        int channel_num;
                        double dac_value, gain_value, gain_error;
                        if(!(iss >> channel_num >> dac_value >> gain_value >> gain_error)){ // <-- This checks if the line was parsed correctly into four values
                            cerr << "Error parsing line " << line_number << " in file " << full_path << ": " << line << endl;
                            continue; // Skip malformed lines
                        }
                        // Fill the TTree
                        channel = channel_num;
                        dac = dac_value;
                        gain = gain_value;
                        gain_err = gain_error;
                        tree->Fill();
                    
                        // Fill the corresponding TGraphErrors
                        int graph_index = channel_num - 64; // Since channels start at 64
                        if(graph_index >= 0 && graph_index < gr.size()){
                            gr[graph_index]->SetPoint(gr[graph_index]->GetN(), dac_value, gain_value);
                            gr[graph_index]->SetPointError(gr[graph_index]->GetN() - 1, 0, gain_error); // Assuming no error in DAC value
                        } else {
                            cerr << "Warning: Channel number " << channel_num << " is out of expected range (64-127)." << endl;
                        }
                    }
                    infile.close();
                } else {
                    cerr << "Error: Unable to open file: " << full_path << endl;
                }
            }

        }
            
    }
    
    gSystem->FreeDirectory(dirp);
    
    // 5. Fit dei grafici e scrittura sul file .root
    cout << "\n Processed " << file_counter << " files. Performing linear fits..." << endl;
    
    outFile->cd(); // Ci posizioniamo nel file ROOT di output

    // We initialize a vector containing p1 parameter values
    vector<double> p1_values;
    
    // We loop over TGraphErrors to perform linear fit and save the graphs in the .root file
    for (auto* graph : gr) {
        if (graph->GetN() > 1) {
            graph->Sort(); // Sort the points by DAC value to ensure proper fitting
            graph->Fit("pol1", "Q"); // Fit lineare in modalità "Quiet" (silenziosa)
            p1_values.push_back(graph->GetFunction("pol1")->GetParameter(1)); // We save the p1 parameter value of the fit in the vector
        }
        graph->Write(); // Salva il grafico nel file ROOT
    }

    // 6. Scrittura del Tree e chiusura definitiva
    tree->Write();
    outFile->Close();
    
    cout << " Analysis complete! Data successfully saved in 'GainResults.root'." << endl;

    /* 7. Print in a canvas the p1 parameter values of the fits 
          in an histogram to check their distribution */
 
    TCanvas* c1 = new TCanvas("c1", "Distribution of p1 parameter values", 1200, 1000);
    c1->cd();
    TH1D* h_p1 = new TH1D("h_p1", "Distribution of SiPM's gain in all runs; Gain [units]; Occurrencies", 160, -1, 3); // Adjust the number of bins and range as needed
    h_p1->SetDirectory(0); // To prevent the histogram from being associated with the current directory (which is now closed)
    for (double val : p1_values) {
        h_p1->Fill(val);
    }
    h_p1->Draw();
   
}