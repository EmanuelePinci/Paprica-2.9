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

void Gain() {
    // 1. GitHub Repository configuration
    std::string repo_url = "https://github.com/EmanuelePinci/Paprica-2.9.git";
    std::string branch_name = "Gain_Data";
    std::string folder_path = "Paprica-2.9"; 
    
    // 2. Download the folder if it doesn't already exist
    if (gSystem->AccessPathName(folder_path.c_str())) {
        std::cout << "Data folder not found. Downloading from GitHub..." << std::endl;
        std::string command = "git clone -b " + branch_name + " --depth 1 " + repo_url;
        int exit_code = system(command.c_str());
        
        if (exit_code != 0) {
            std::cerr << "❌ Error during download. Exiting." << std::endl;
            return;
        }
    } else {
        std::cout << "✅ Data folder already exists. Skipping download." << std::endl;
    }

    // 3. Create the output .root file
    TFile* outFile = new TFile("GainResults.root", "RECREATE");

    // 4. Initialize the TTree
    TTree* tree = new TTree("GainTree", "Tree containing DAC and Gain data");
    int channel;
    double dac;
    double gain;
    double gain_err;
    
    tree->Branch("channel", &channel, "channel/I");
    tree->Branch("dac", &dac, "dac/D");
    tree->Branch("gain", &gain, "gain/D");
    tree->Branch("gain_err", &gain_err, "gain_err/D");

    // 5. Initialize an array of 64 TGraphErrors (for channels 64 to 127)
    const int num_channels = 64;
    const int ch_offset = 64; // Channels start at 64
    TGraphErrors* gr[num_channels];
    
    for(int i = 0; i < num_channels; i++) {
        gr[i] = new TGraphErrors();
        // Name the graph "Graph_ch64", "Graph_ch65", etc.
        gr[i]->SetName(Form("Graph_ch%d", i + ch_offset));
        gr[i]->SetTitle(Form("Gain vs DAC - Channel %d;DAC Value;Gain", i + ch_offset));
        gr[i]->SetMarkerStyle(20); // Solid circle marker
        gr[i]->SetMarkerSize(0.8);
        gr[i]->SetMarkerColor(kBlue);
        gr[i]->SetLineColor(kBlue);
    }

    std::cout << "\n--- Starting Data Analysis ---\n" << std::endl;

    // 6. Open the directory
    void* dirp = gSystem->OpenDirectory(folder_path.c_str());
    if (!dirp) {
        std::cerr << "❌ Error: Could not open directory " << folder_path << std::endl;
        return;
    }

    const char* entry;
    int file_counter = 0;
    
    // 7. Loop over all entries in the directory
    while ((entry = gSystem->GetDirEntry(dirp))) {
        TString fileName = entry;

        if (fileName.EndsWith(".txt")) {
            file_counter++;
            std::string full_path = folder_path + "/" + fileName.Data();
            
            std::ifstream infile(full_path);
            if (infile.is_open()) {
                
                std::string header_line;
                // Read and discard the first line (the comment "\\ valori")
                std::getline(infile, header_line); 
                
                // Read the remaining lines: 4 columns
                while (infile >> channel >> dac >> gain >> gain_err) {
                    
                    // Fill the TTree
                    tree->Fill();
                    
                    // Fill the corresponding TGraphErrors
                    if (channel >= 64 && channel <= 127) {
                        int idx = channel - ch_offset; // Map ch 64->0, ch 65->1, etc.
                        int n_pts = gr[idx]->GetN();   // Current number of points in the graph
                        
                        // Set X=dac, Y=gain
                        gr[idx]->SetPoint(n_pts, dac, gain);
                        // Set X_error=0, Y_error=gain_err
                        gr[idx]->SetPointError(n_pts, 0, gain_err);
                    }
                }
                infile.close();
            }
        }
    }

    gSystem->FreeDirectory(dirp);
    
    // 8. Fit the graphs with a linear function and write to the .root file
    std::cout << "✅ Processed " << file_counter << " files. Performing linear fits..." << std::endl;
    
    outFile->cd(); // Ensure we are writing to the file
    
    for(int i = 0; i < num_channels; i++) {
        // Only fit if the graph has at least 2 points
        if (gr[i]->GetN() > 1) {
            // "pol1" is a linear fit (y = p0 + p1*x)
            // "Q" means Quiet mode (don't print fit results for all 64 channels to the terminal)
            gr[i]->Fit("pol1", "Q"); 
        }
        // Write the graph (along with its fitted function) to the file
        gr[i]->Write(); 
    }

    // 9. Write the Tree and close the file
    tree->Write();
    outFile->Close();
    
    std::cout << "✅ Analysis complete! All graphs and the TTree are saved in 'GainResults.root'." << std::endl;
    std::cout << "   You can view them by typing: new TBrowser(); in your ROOT session." << std::endl;
}
