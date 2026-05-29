#include <TFile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <string>
#include <TColor.h>
#include <TF1.h>
#include <TStyle.h>

using namespace std;

void GainMultiGraph(string folder_path, int graphx, int graphy, const vector<int>& i_sipm, bool save = false){
    
    if(i_sipm.empty()) {
        cout << "You select multigraph mode" << endl;
        cout << "Taking data from" << folder_path + "/GainResults.root" << endl;

        TString fileInName = Form("%s/GainResults.root", folder_path.c_str());
    
        TFile *fileIn = new TFile(fileInName, "UPDATE");
    
        if (!fileIn || fileIn->IsZombie()) {
            std::cerr << "ERRORE: Impossibile aprire " << fileInName << std::endl;
            return;
        }

        /* We initialize the multigraphs */

        TMultiGraph *g_g = new TMultiGraph();
    
        g_g->SetTitle("Good Gain; DAC [V]; Gain [Units]");

        TLegend *leg_good = new TLegend(0.1, 0.1, 0.3, 0.9);
        leg_good->SetNColumns(2); // We set the legend to have 2 columns to save space

        /* We fill the multigraphs with the proper graphs*/
        
        for(int i_graph = 64; i_graph < 128; i_graph++){
            // We set the graph name to get the graph from the file 
            TString graphName = Form("Graph_ch%d", i_graph);
            // We want also error bars
            TGraphErrors *graph = (TGraphErrors*)fileIn->Get(graphName);

            if (graph){
                graph->SetMarkerStyle(21); // We set the marker style to squares
                graph->SetMarkerSize(1.0);
                gStyle->SetOptFit(1111); // Show fit parameters on the graph

                // We set color from HSV map 

                int idx = (i_graph -64);
                float hue = (idx * 360.0) / 64.0; // Map index to hue (0-360)
                float saturation = 1.0;
                float value = 0.85;
                float r, g, b;
                TColor::HSV2RGB(hue, saturation, value, r, g, b);
                int color = TColor::GetColor(r, g, b);

                // We set the color of the graph

                graph->SetLineColor(color);
                graph->SetMarkerColor(color);
                
                TF1 *fitFunc = graph->GetFunction("pol1");
                if(fitFunc) {
                    fitFunc->SetLineColor(color);
                }

                // We write the graph on the multigraph
                g_g->Add(graph);
                leg_good->AddEntry(graph, Form("SiPM %d", i_graph), "lp");

            }else{
            
                std::cerr << "Warning: Impossible to find graph for SiPM " << i_graph << std::endl;
            
            } 

        }
        // We save and draw the multigraphs on canvas
        TCanvas *c_multig = new TCanvas("c_multig", "Multigraph of Gains", 1400, 800);
        c_multig->cd();
        g_g->Draw("APE"); // <-- A = axis, P = points, E = error bars
        leg_good->Draw();
        c_multig->Update(); // We update the canvas to show the multigraph

    }else if(i_sipm.size() == graphx * graphy) {
        cout << "Taking data from" << folder_path + "/GainResults.root" << endl;
        cout << "Putting" << graphx << "on x-axis and" << graphy << " on y-axis" << endl;
        cout << "Analyzing SiPMs with indices: ";
        for (const auto& idx : i_sipm) {
            cout << idx << " ";
        }cout << endl;
    
        /* We open the file with the gain results, we check if it is open correctly and if not we print
           an error message and we exit the function.
           PAY ATTENTION: file path refers to my own laptop you must to change it */
    
        TString fileInName = Form("%s/GainResults.root", folder_path.c_str());
    
        TFile *fileIn = new TFile(fileInName, "UPDATE");
    
        if (!fileIn || fileIn->IsZombie()) {
            std::cerr << "ERRORE: Impossibile aprire " << fileInName << std::endl;
            return;
        }

        /* We initialize the multigraphs */
        
        TMultiGraph *g_gs = new TMultiGraph();
    
        g_gs->SetTitle("Good Gain; DAC [V]; Gain [Units]");

        TLegend *leg_goods = new TLegend(0.1, 0.1, 0.3, 0.9);
        leg_goods->SetNColumns(2); // We set the legend to have 2 columns to save space

        /* We fill the multigraphs with the proper graphs*/
        
        for(int i_graph = 0; i_graph < graphx * graphy; i_graph++){

            // We define the index of the SiPM to analyze 
            int index = i_sipm[i_graph];
            // We set the graph name to get the graph from the file 
            TString graphName = Form("Graph_ch%d", index);
            // We want also error bars
            TGraphErrors *graph = (TGraphErrors*)fileIn->Get(graphName);

            if (graph){
                graph->SetMarkerStyle(21); // We set the marker style to squares
                graph->SetMarkerSize(1.0);
                gStyle->SetOptFit(1111); // Show fit parameters on the graph

                // We set color from HSV map 

                int idx = (i_graph -64);
                float hue = (idx * 360.0) / 64.0; // Map index to hue (0-360)
                float saturation = 1.0;
                float value = 0.85;
                float r, g, b;
                TColor::HSV2RGB(hue, saturation, value, r, g, b);
                int color = TColor::GetColor(r, g, b);

                // We set the color of the graph

                graph->SetLineColor(color);
                graph->SetMarkerColor(color);
                
                TF1 *fitFunc = graph->GetFunction("pol1");
                if(fitFunc) {
                    fitFunc->SetLineColor(color);
                }

                // We write the graph on the multigraph
                g_gs->Add(graph);
                leg_goods->AddEntry(graph, Form("SiPM %d", index), "lp");

            }else{
            
                std::cerr << "Warning: Impossible to find graph for SiPM " << index << std::endl;
            
            } 

        }
        // We save and draw the multigraphs on canvas
        TCanvas *c_multisel = new TCanvas("c_multisel", "Multigraph of Gains", 1400, 800);
        c_multisel->cd();
        g_gs->Draw("APE"); // <-- A = axis, P = points, E = error bars
        leg_goods->Draw();
        c_multisel->Update(); // We update the canvas to show the multigraph

        /* We create a canvas where we will draw graph matrix and the multigraph for selected channels*/

        TCanvas *c_matrix = new TCanvas("c_matrix", "Gain Graphs", 1000, 1000);
        c_matrix->Divide(graphx, graphy); // We divide the canvas in a grid of graphx by graphy

        // We fill the graph matrix

        for(int i_graph = 0; i_graph < graphx * graphy; i_graph++) {
            int index = i_sipm[i_graph];
            TString graphName = Form("Graph_ch%d", index);
            TGraphErrors *graph = (TGraphErrors*)fileIn->Get(graphName);

            if (graph){
                c_matrix->cd(i_graph + 1); // We select the pad corresponding to the current graph
                graph->SetTitle(Form("SiPM %d; DAC; Gain [Units]", index));
                graph->SetMarkerStyle(21);
                graph->SetMarkerColor(kBlue);
                graph->SetMarkerSize(1.0);
                gStyle->SetOptFit(1111); // Show fit parameters on the graph
                graph->Draw("APE"); // <-- A = axis, P = points, E = error bars
            }else{
                std::cerr << "Warning: Impossible to find graph for SiPM " << index << std::endl;
            }
        }

        // We save and draw the canvas with the graph matrix
        c_matrix->Update();
        c_matrix->Modified();
        if(save == true){
            c_matrix->SaveAs(Form("%s/GainGraphMatrix.png", folder_path.c_str()));
        }
    }else{
        cout << "Error: The number of SiPMs selected does not match the grid size (graphx * graphy)." << endl;
        cout << "Please check your inputs and try again." << endl;
        return;
    }

    std::cout << "\n===== ANALISI COMPLETATA =====" << std::endl;
}