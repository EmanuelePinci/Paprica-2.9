#include <TFile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

/* We define a struct to save data in first step. This is mandatory because we want analyze
   a series of run and not a single one. The problem is that varying DAC values we vary the 
   gain and so we vary the slopes, so we need to adapt the acceptance range accordingly.
   A struct is a feature of C and C++ that allows us to group together different variables. */

struct SiPMData{
    int i_sipm;
    TGraphErrors* graph;
    double slopeValue;
};

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
    
    mg_good->SetTitle(Form("Good SiPMs for run %s; Peak Number; Mean Charge [Units]", numrun.c_str()));
    
    TMultiGraph *mg_bad = new TMultiGraph();
    
    mg_bad->SetTitle(Form("Bad SiPMs for run %s; Peak Number; Mean Charge [Units]", numrun.c_str()));

    TLegend *leg_good = new TLegend(0.1, 0.1, 0.3, 0.9);
    
    TLegend *leg_bad = new TLegend(0.1, 0.1, 0.3, 0.9);

    /* At this point we want to put all the slopes in a histogram to have an idea of the 
       distribution of the gains. We set also the directory to visualize histogram. */

    TH1D *h_P1 = new TH1D("h_P1", Form("Distribution of Gain for run %s", numrun.c_str()), 100, 100, 400);
    h_P1->SetDirectory(0);

    // We define now a variable to store temporary the slope

    vector<SiPMData> temp_list;
    vector<double> slope;

    /*/ This values depends on the slopes so when we have a final range we change them
    
    double slopeMin = 220.0; 
    double slopeMax = 280.0;*/

    // Now we collect data and compute the slopes

    for (int i_sipm = 64; i_sipm <= 127; ++i_sipm) {
        
        TString graphName = Form("gain_graph_isipm%d", i_sipm);
        
        TGraphErrors *graph = (TGraphErrors*)fileIn->Get(graphName);

        if (graph) {
            
            /* We now extract the slope (gain) from the fit line associated to the graph, 
               we will use this value to fill the histogram and to decide if the SiPM is good 
               or bad */

            TF1 *lineFit = graph->GetFunction(Form("lineFit_sipm%d", i_sipm));

            if (lineFit) {
                slope.push_back(lineFit->GetParameter(1)); // We extract the slope (gain) from the linear fit
                h_P1->Fill(slope.back()); // We fill the histogram with the slope value
                temp_list.push_back({i_sipm, graph, slope.back()}); // We 

            } else {
                std::cerr << "Warning: Impossible to find linear fit function for SiPM " 
                          << i_sipm << std::endl;
                continue;
            }
        }
    }

    //Now we compute dinamically for each run the acceptance range 

    /* We consider as good SiPMs those with slope within 2 standard deviations from the mean,
      and as bad SiPMs those with slope outside this range */
            
    int maxBin = h_P1->GetMaximumBin();
    double PeakCenter = h_P1->GetBinCenter(maxBin);
    double PeakSigma = h_P1->GetRMS();
    double slopeMin = PeakCenter - 2 * PeakSigma;
    double slopeMax = PeakCenter + 2 * PeakSigma;

    // We now fill the multigraphs separating good and bad SiPMs
    for (const auto& data : temp_list) {
        int i_sipm = data.i_sipm;
        TGraphErrors* graph = data.graph;
        double slopeValue = data.slopeValue;
        TF1 *lineFit = graph->GetFunction(Form("lineFit_sipm%d", i_sipm));

        /* In next lines we define multigraph and graphics logic in order to have a better
            visualization of the data and of the results */

        graph->SetMarkerStyle(21); // <-- We set the marker style
            
        // We set different colors from HSV map

        int idx = (i_sipm -64);
        float hue = (idx * 360.0) / 64.0; // Map index to hue (0-360)
        float saturation = 1.0;
        float value = 0.85;
        float r, g, b;
        TColor::HSV2RGB(hue, saturation, value, r, g, b);
        int color = TColor::GetColor(r, g, b);

        // We set the color of the graph

        graph->SetLineColor(color);
        graph->SetMarkerColor(color);
        lineFit->SetLineColor(color);

        std::cout << "SiPM " << i_sipm << " slope = " << slopeValue << std::endl;

        //h_P1->Fill(slope.back());

        if (slopeValue >= slopeMin && slopeValue <= slopeMax)
        {
            mg_good->Add(graph);
            leg_good->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");
        }
        else if (slopeValue < slopeMin || slopeValue > slopeMax)
        {
            mg_bad->Add(graph);
            leg_bad->AddEntry(graph, Form("SiPM %d", i_sipm), "lp");
        }
    }

    TString fileOutName = Form("/home/emanuele/Universita/Paprica 2.9/Dataset/Fit_Data/SiPM_MultiGraph_run_%s.root", numrun.c_str());
    TFile *fileOut = new TFile(fileOutName, "RECREATE");
    
    // Canvas 1: Histogram of Gain Slopes

    TCanvas *c_dist = new TCanvas("c_dist", Form("Distribution of Gain for run %s", numrun.c_str()), 800, 600);
    c_dist->cd();
    h_P1->SetLineColor(kBlue);
    h_P1->Draw("HIST");
    fileOut->cd();
    h_P1->Write("h_P1");

    // Canvas 2: Good SiPMs
    
    TCanvas *c_good = new TCanvas("c_good", Form("Good SiPMs for run %s", numrun.c_str()), 1400, 800);
    
    // Check to avoid crash if mg_good is empty

    if (mg_good->GetListOfGraphs() && mg_good->GetListOfGraphs()->GetSize() > 0) {
        mg_good->Draw("AP"); 
        leg_good->Draw();
    }
    fileOut->cd();
    mg_good->Write("MultiGraph_good");
    c_good->Write("Canvas_good");

    // Canvas 3: Bad SiPMs

    TCanvas *c_bad = new TCanvas("c_bad", Form("Bad SiPMs for run %s", numrun.c_str()), 1400, 800);
    
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