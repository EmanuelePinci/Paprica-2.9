#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLine.h>
#include <TBox.h>
#include <TLegend.h>
#include <TLatex.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <TROOT.h>

using namespace std;

void TotalCharge(const vector<TString> &myfile)
{

    TString mypath = "/Users/Utente/Desktop/PAPRICA/run_RAGGICOSMICI/Risultati/";

    if (gSystem->AccessPathName(mypath.Data()))
    {
        cerr << "Error: Output folder not found: " << mypath << endl;
        cerr << "Run AllEventDisplay first to generate the .txt files." << endl;
        return;
    }

    // ---------------------------------------------------------------
    // 1. Theoretical estimate of MPV togliere diagonale
    // ---------------------------------------------------------------
    //mettere valore atteso e fare la banda più larga sfruttando radice di due
    double dEdx = 1.9;      // MeV/cm, muon MIP in plastic
    double Ly = 8000;       // photons/MeV, scintillating fiber light yield
    double QE = 0.35;       // SiPM quantum efficiency
    double eps_trap = 0.035; // fiber trapping efficiency divided by two

    // Two geometrical cases: vertical (Δx = 1.1 cm) and diagonal (Δx * sqrt(2))
    double dx_vert = 1.0;                 // cm  (2.2cm / 2, average vertical path)
    double dx_diag = dx_vert * sqrt(2.0); // cm  diagonal path

    double NPE_vert = dEdx * dx_vert * Ly * QE * eps_trap;
    double NPE_diag = dEdx * dx_diag * Ly * QE * eps_trap;
    double NPE_mean = 0.5 * (NPE_vert + NPE_diag);

    // Relative uncertainties (added in quadrature)
    double rel_dEdx = 0.1 / dEdx;
    double rel_Ly = 0.10; // 10%
    double rel_QE = 0.03 / QE;
    double rel_trap = 0.10; // 10%
    // Δx uncertainty: half the interval between vertical and diagonal
    double rel_dx = (dx_diag - dx_vert) / (2.0 * dx_vert);

    double rel_err = sqrt(rel_dEdx * rel_dEdx + rel_dx * rel_dx +
                          rel_Ly * rel_Ly + rel_QE * rel_QE + rel_trap * rel_trap);
    double NPE_err = NPE_mean * rel_err;

    // Print theoretical calculation
    cout << "\n========================================" << endl;
    cout << "  Theoretical MPV estimate" << endl;
    cout << "========================================" << endl;
    cout << "  dE/dx        = " << dEdx << " MeV/cm" << endl;
    cout << "  dx (vert)    = " << dx_vert << " cm  -> N_PE = " << NPE_vert << " PE" << endl;
   // cout << "  dx (diag)    = " << dx_diag << " cm  -> N_PE = " << NPE_diag << " PE" << endl;
    cout << "  Light Yield  = " << Ly << " photons/MeV" << endl;
    cout << "  QE           = " << QE << endl;
    cout << "  eps_trap     = " << eps_trap << endl;
    cout << "----------------------------------------" << endl;
    cout << "  Expected range : [" << NPE_vert << ", " << NPE_diag << "] PE" << endl;
    cout << "  Mean +/- err   :  " << NPE_mean << " +/- " << NPE_err << " PE" << endl;
    cout << "  Relative error :  " << rel_err * 100 << " %" << endl;
    cout << "========================================\n"
         << endl;

    // ---------------------------------------------------------------
    // 2. Fill histogram
    // ---------------------------------------------------------------
    TH1D *htc = new TH1D("hTotalCharge",
                         "Total Charge Distribution;Total Charge [PE];Occurrence",
                         100, 0, 0);

    double ChargeTot_PE = 0.0;
    for (size_t f = 0; f < myfile.size(); f++)
    {
        TString fileName = mypath + myfile[f];
        ifstream fileInput(fileName.Data());
        if (!fileInput.is_open())
        {
            cerr << "Error: Cannot open " << fileName << endl;
            continue;
        }
        string header;
        getline(fileInput, header);
        while (fileInput >> ChargeTot_PE)
        {
            htc->Fill(ChargeTot_PE);
        }
        fileInput.close();
        cout << "Loaded: " << fileName << endl;
    }

    // ---------------------------------------------------------------
    // 3. Landau fit
    // ---------------------------------------------------------------
    double fit_min = 170; // htc->GetMean() - htc->GetRMS();
    double fit_max = 500; // htc->GetMean() + 2.0 * htc->GetRMS();
    TF1 *fLandau = new TF1("fLandau", "landau", fit_min, fit_max);
    fLandau->SetLineColor(kRed + 1);
    fLandau->SetLineWidth(2);
    htc->Fit(fLandau, "RQ");

    double MPV_measured = fLandau->GetParameter(1);
    double MPV_measured_err = fLandau->GetParError(1);

    cout << "----------------------------------------" << endl;
    cout << "  Landau fit result" << endl;
    cout << "  MPV measured = " << MPV_measured << " +/- " << MPV_measured_err << " PE" << endl;
    cout << "  Expected range [" << NPE_vert << ", " << NPE_diag << "] PE" << endl;
    cout << "----------------------------------------\n"
         << endl;

    // ---------------------------------------------------------------
    // 4. Draw
    // ---------------------------------------------------------------
    TCanvas *ctc = new TCanvas("ctc", "Distribution of total charge [PE]", 900, 650);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1);

    htc->SetLineColor(kBlue + 1);
    htc->SetLineWidth(2);
    htc->Draw();
    fLandau->Draw("same");

    // Theoretical band [NPE_vert, NPE_diag]
    double ymax = htc->GetMaximum();
    TBox *band = new TBox(NPE_vert, 0, NPE_diag, ymax * 1.05);
    band->SetFillColorAlpha(kGreen + 1, 0.25);
    band->SetLineColor(kGreen + 3);
    band->SetLineWidth(1);
    band->Draw("same");

    // Measured MPV line
    TLine *lineMPV = new TLine(MPV_measured, 0, MPV_measured, ymax * 1.05);
    lineMPV->SetLineColor(kRed + 1);
    lineMPV->SetLineWidth(2);
    lineMPV->SetLineStyle(2); // dashed
    lineMPV->Draw("same");

    // Redraw histogram on top of band
    htc->Draw("same");
    fLandau->Draw("same");

    // Legend
    TLegend *leg = new TLegend(0.55, 0.60, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->AddEntry(htc, "Data (12 runs)", "l"); // metti automatico
    leg->AddEntry(fLandau, "Landau fit", "l");
    leg->AddEntry(lineMPV, Form("MPV = %.1f #pm %.1f PE", MPV_measured, MPV_measured_err), "l");
    leg->AddEntry(band, Form("Theory: [%.0f, %.0f] PE", NPE_vert, NPE_diag), "f");
    leg->Draw();

    ctc->SaveAs(mypath + "Total_Charge_PE.png");
    cout << "Saved: " << mypath << "Total_Charge_PE.png" << endl;
}

// --- Wrapper ---
void RunTotalCharge()
{
    vector<TString> files = {
        "TotalCharge_run_127.txt",
        "TotalCharge_run_130.txt",
        "TotalCharge_run_131.txt",
        "TotalCharge_run_132.txt",
        "TotalCharge_run_133.txt",
        "TotalCharge_run_134.txt",
        "TotalCharge_run_135.txt",
        "TotalCharge_run_136.txt",
        "TotalCharge_run_137.txt",
        "TotalCharge_run_138.txt",
        "TotalCharge_run_139.txt",
        "TotalCharge_run_140.txt",
        "TotalCharge_run_141.txt"};
    TotalCharge(files);
}
