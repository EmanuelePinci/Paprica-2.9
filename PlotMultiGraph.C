#include <TFile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

void PlotMultiGraph(std::string numrun){

    // We verify for security that the number of run is correct printing in on console
    
    cout << numrun << endl;

    // PAY ATTENTION: file path refers to my own laptop you must to change it 
    
    TString fileInName = Form("/home/emanuele/Universita/Paprica 2.9/Dataset/Fit_Data/SiPM_Graphs_run_%s.root", numrun.c_str());
    
    TFile *fileIn = new TFile(fileInName, "READ");
    
    if (!fileIn || fileIn->IsZombie()) {
        std::cerr << "ERRORE: Impossibile aprire " << fileInName << std::endl;
        return;
    }

    /* Here we create the multigraph following the root instructions, since the data can be divided
       in two categories, the one that follow a good linear trend and the one that do not follow a 
       good linear trend we create two different multigraphs to avoid graphical corruption of the 
       graph with all the data.*/

    TMultiGraph *mg_good = new TMultiGraph();
    
    mg_good->SetTitle("Good SiPMs; Peak Number; Mean Charge [Units]");
    
    TMultiGraph *mg_bad = new TMultiGraph();
    
    mg_bad->SetTitle("Bad SiPMs; Peak Number; Mean Charge [Units]");

    TLegend *leg_good = new TLegend(0.1, 0.1, 0.3, 0.9);
    
    TLegend *leg_bad = new TLegend(0.1, 0.1, 0.3, 0.9);

    /* At this point we want to put all the slopes in a histogram to have an idea of the 
       distribution of the gains */

    TH1D *h_P1 = new TH1D("h_P1", "Distribution of Gain Slopes (p1)", 100, 100, 400);

    // This values depends on the slopes so when we have a final range we change them
    
    double slopeMin = 220.0; 
    double slopeMax = 280.0;

    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm) {
        
        TString graphName = Form("gain_graph_isipm%d", i_sipm);
        
        TGraphErrors *graph = (TGraphErrors*)fileIn->Get(graphName);

        if (graph) {
            
            /* We now extract the slope (gain) from the fit line associated to the graph, 
               we will use this value to fill the histogram and to decide if the SiPM is good 
               or bad */

            TF1 *lineFit = graph->GetFunction(Form("lineFit_sipm%d", i_sipm));
            
            double slope = 0;

            if (lineFit) {
                slope = lineFit->GetParameter(1); // We extract the slope (gain) from the linear fit 
            } else {
                std::cerr << "Warning: Impossible to find linear fit function for SiPM " 
                          << i_sipm << std::endl;
                continue;
            }

            std::cout << "SiPM " << i_sipm << " slope = " << slope << std::endl;

            h_P1->Fill(slope);

            if (slope >= slopeMin && slope <= slopeMax)
            {
                mg_good->Add(graph);
                leg_good->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");
            }
            else
            {
                mg_bad->Add(graph);
                leg_bad->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");
            }
        }
    }

    TString fileOutName = Form("/home/emanuele/Universita/Paprica 2.9/Dataset/Fit_Data/SiPM_MultiGraph_run_%s.root", numrun.c_str());
    TFile *fileOut = new TFile(fileOutName, "RECREATE");
    
    // Canvas 1: Histogram of Gain Slopes

    TCanvas *c_dist = new TCanvas("c_dist", "Distribution of Gain Slopes (p1)", 800, 600);
    c_dist->cd();
    h_P1->SetLineColor(kBlue);
    h_P1->SetFillColor(kCyan);
    h_P1->Draw("HIST");
    fileOut->cd();
    h_P1->Write("h_P1");

    // Canvas 2: Good SiPMs
    
    TCanvas *c_good = new TCanvas("c_good", "SiPM Buoni", 1400, 800);
    
    // Check to avoid crash if mg_good is empty

    if (mg_good->GetListOfGraphs() && mg_good->GetListOfGraphs()->GetSize() > 0) {
        mg_good->Draw("AP"); 
        leg_good->Draw();
    }
    fileOut->cd();
    mg_good->Write("MultiGraph_good");
    c_good->Write("Canvas_good");

    // Canvas 3: Bad SiPMs

    TCanvas *c_bad = new TCanvas("c_bad", "SiPM Fuori Range", 1400, 800);
    
    // Check to avoid crash if mg_bad is empty
    
    if (mg_bad->GetListOfGraphs() && mg_bad->GetListOfGraphs()->GetSize() > 0) {
        mg_bad->Draw("AP");
        leg_bad->Draw();
    }
    fileOut->cd();
    mg_bad->Write("MultiGraph_bad");
    c_bad->Write("Canvas_bad");

    // We close the files in a clean way
    
    fileOut->Close();
    fileIn->Close();

    std::cout << "\n===== ANALISI COMPLETATA =====" << std::endl;
    std::cout << "MultiGraphs (Good/Bad) e Istogramma Slopes salvati in: " << fileOutName << std::endl;
}