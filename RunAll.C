#include "FitPeaks.C"
#include "ExtractGains.C"
#include "PlotMultiGraph.C"

void RunAll(const vector<int>& mydac, const vector<string>& numrun) {
    /* We check first of all if the size of mydac and numrun are equal*/

    if (mydac.size() != numrun.size()) {
        std::cerr << "ERROR: The size of mydac and numrun must be the same." << std::endl;
        return;
    }

    // Now we instanziate a cicle for analyzing all the runs together
    
    for(size_t i = 0; i < mydac.size(); ++i) {
        std::cout << "=== STARTING ANALYSIS FOR RUN " << numrun[i] << " ===" << std::endl;
        
        FitPeaks(numrun[i]);
        ExtractGains(mydac[i], numrun[i]);
        PlotMultiGraph(numrun[i]);
        
        std::cout << "=== ANALYSIS FOR RUN " << numrun[i] << " COMPLETED ===" << std::endl;
    }
}