#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <vector>
#include <sstream>

using namespace std;

void ExtractDAC(string file_name, double target_gain){

    /* We open the .txt file containing Gain vs DAC fit results in order to invert fit function
       and obtain the DAC value corresponding to a specific gain value for each SIPM */
       
       string myfile = "/home/emanuele/Universita/Paprica 2.9/Dataset/Gain_Data/" + file_name;
       ifstream infile(myfile); // We open the file containing the fit results

    // We check if the file was opened succesfully

    if(!infile.is_open()){
        cerr << "Error: Unable to open file: " << myfile << endl;
        return;
    }        
    cout << "Processing file: " << myfile << " for a specific gain value: " << target_gain << endl;
        
    string line;

    // We made a check on file to check if it is empty

    if(!getline(infile, line)){
        cerr << "ERROR: Given data file is empty!" << endl;
        return;
    }
    
    /* We initialize two different vectors that identify the two CITIROC and we read the channel as 
       CITIROCs do, i.e. from 0 up to 31. We define VDAC as int because the CITIROC normally read
       values with 8 bit from 0 to 255 so double can make issues. So to avoid this we define these as
       int. */
    
       vector<int> CITIROC0(32, 0);
       vector<int> CITIROC1(32,0);
    
    // We define variables to read correctly the file

    int i_sipm = 0;
    double p0 = 0.0, p1 = 0.0, VDAC = 0.0;

    // We start to read file line by line 

    while(getline(infile, line)){
        
        if(line.empty()) continue;

        stringstream ss(line);

        // We fill the vectors with correct value
        if(ss >> i_sipm >> p0 >> p1){ // <-- We read the values ordinately

            /* Initial fit formula was GAIN = p1 * VDAC + p0, so we invert the formula in order
               to obtain VDAC. So we define VDAC */

               VDAC = (target_gain - p0) / p1;
            
            // We define the CITIROC mapping:
            // For CITIROC0 we have: 64  --> 0
            //                       95  --> 31
            // For CITIROC1 we have: 96  --> 0
            //                       127 --> 31

            if(i_sipm >= 64 && i_sipm <= 95){
                 int iglbsipm = i_sipm - 64;
                CITIROC0[iglbsipm] = static_cast<int>(round(VDAC));
            }else if(i_sipm >= 96 && i_sipm <= 127){
                int iglbsipm = i_sipm - 96;
                CITIROC1[iglbsipm] = static_cast<int>(round(VDAC));
            }

        }
    }

    // We close the input file

    infile.close();

    // Now we must to save the data in output file. We create first of all the output file
    
    ofstream outfile("iSipm_iglbSipm_DAC.txt");

    // Now we fill the file, the first line specifies the specific gain value kept into account

    if(outfile.is_open()){
        outfile << "// Gain value selected: " << target_gain << endl;
        outfile << "//i_SiPM" << "\t" << "iglbSiPM" << "\t" << "DAC" << endl;
        for(int i_SiPM = 64; i_SiPM < 128; i_SiPM++){
            if(i_SiPM < 96){
                outfile << i_SiPM << "\t" << i_SiPM - 64 << "\t" << CITIROC0[i_SiPM - 64] << endl;
            }else{
                outfile << i_SiPM << "\t" << i_SiPM - 96 << "\t" << CITIROC1[i_SiPM - 96] << endl;
            }
        }
        outfile.close();
        cout << "[INFO] Risultati salvati in 'iSipm_iglbSipm_DAC.txt'" << endl;
    }
}