#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>

{
    std::vector<std::string> filenames = {
        "u1_1.30-u3_1.69_T_164.txt",
        "u1_1.00-u3_1.69_T_164.txt",
        "u1_1.15-u3_1.69_T_164.txt",
        "u1_1.15-u3_1.3_T_164.txt",
        "u1_1.15-u3_1.4_T_164.txt",
        "u1_1.15-u3_1.5_T_164.txt",
        "u1_1.15-u3_1.5_T_170.txt",
        "u1_1.15-u3_1.5_T_175.txt",
        "u1_1.15-u3_1.5_T_181.txt", 
        "u1_1.15-u3_1.5_T_193.txt",
    };

    std::vector<std::vector<std::pair<double, double>>> allSetRanges = {
        {{9.3, 10.6}, {14, 15.6}, {19.2, 20.7}, {24.2, 25.95}, {29.8, 31}},
        {{9.5, 10.94}, {14.5, 15.8}, {19.45, 20.9}, {24.4, 25.98}, {29.8, 31.5}},
        {{9.4, 10.8}, {14.3, 15.7}, {19.3, 20.8}, {24.4, 25.98}, {29.8, 31.2}},
        {{9.4, 10.6}, {14.35, 15.5}, {19.3, 20.7}, {24.4, 25.9}, {29.7, 31}},
        {{9.4, 10.7}, {14.3, 15.6}, {19.35, 20.7}, {24.3, 25.9}, {30, 31}},
        {{9.4, 10.7}, {14.3, 15.6}, {19.35, 20.7}, {24.3, 25.9}, {30, 31}},
        {{9.3, 10.6}, {14.2, 15.3}, {19.1, 20.5}, {24.2, 25.5}, {29.4, 30.7}},
        {{9.3, 10.5}, {14.2, 15.3}, {19.2, 20.55}, {24.3, 25.55}, {29.4, 30.9}},
        {{9.3, 10.7}, {14, 15.5}, {19.2, 20.6}, {24.1, 25.6}, {29.3, 30.8}},
        {{9.3, 10.7}, {14, 15.5}, {19.2, 20.4}, {24.1, 25.4}, {29.2, 30.4}}
    };

    std::vector<std::vector<double>> allMeans;
    std::vector<std::vector<double>> allSigmas;
    std::vector<std::vector<double>> allDeltaMeans;
    std::vector<std::vector<double>> allUncertainties;
    std::vector<std::vector<double>> allStds;
    std::vector<double> finalmeans;
    std::vector<double> results;
    std::vector<double> rsigmas;
    
    TF1 *fnew = new TF1("fnew", "gaus");
    TF1 *poisson= new TF1("poissonFunc", "TMath::Poisson(x, [0])");

    for (size_t fileIndex = 0; fileIndex < filenames.size(); ++fileIndex) {
        
        std::string filename = filenames[fileIndex];
        std::stringstream titleStream;

        size_t u1Pos = filename.find("u1_");
        size_t u3Pos = filename.find("-u3_");
        size_t tPos = filename.find("_T_");

    if (u1Pos != std::string::npos && u3Pos != std::string::npos && tPos != std::string::npos) {
        double u1, u3, t;
        std::istringstream(filename.substr(u1Pos + 3, u3Pos - (u1Pos + 3))) >> u1;
        std::istringstream(filename.substr(u3Pos + 4, tPos - (u3Pos + 4))) >> u3;
        std::istringstream(filename.substr(tPos + 3)) >> t;
    
        titleStream << "U1 = " << std::fixed << std::setprecision(2) << u1 << " V, U3 = " << u3 << " V, T = " << std::fixed << std::setprecision(0) << t << " Celcius";} else {
        titleStream << "Upper Title Here"; 
    }

    std::string upperTitle = titleStream.str();

    TGraphErrors gr(filenames[fileIndex].c_str());
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.48);
    gStyle->SetStatY(0.9);

    gr.GetXaxis()->SetTitle("Acceleration Voltage [V]");
    gr.GetYaxis()->SetTitle("Current [nA]");
    gr.SetTitle(upperTitle.c_str());
    gr.SetLineColor(kBlack);
    gr.SetMarkerStyle(20); 
            
        

    for (int i = 0; i < gr.GetN(); ++i) {
        gr.SetPointError(i, 0.15, 0.005);
    }

    std::vector<double> means;
    std::vector<double> sigmas;
    std::vector<double> deltaMeans;
    std::vector<double> uncertainties;
    double stdmean[5] = {0, 0, 0, 0, 0};       
    double stdsigma[5] = {0, 0, 0, 0, 0};       

    for (int i = 0; i < 5; i++) {
        double rangeFirst = allSetRanges[fileIndex][i].first;
        double rangeSecond = allSetRanges[fileIndex][i].second;
        fnew->SetRange(rangeFirst, rangeSecond);
        gr.Fit("fnew", "QSR+");
        means.push_back(fnew->GetParameter(1));
        sigmas.push_back(fnew->GetParameter(2));
        stdmean[i] = fnew->GetParError(1);
        stdsigma[i] = fnew->GetParError(1);
    }

    //TCanvas *c1 = new TCanvas();
    //gr.Draw();
    std::string finalfilename = std::to_string(fileIndex) + "fitted.pdf";
    //c1->Print(finalfilename.c_str());


    for (size_t i = 1; i < means.size(); ++i) {
        double deltaMean = means[i] - means[i - 1];
        double uncertainty = std::sqrt(sigmas[i] * sigmas[i] + sigmas[i - 1] * sigmas[i - 1]);
        deltaMeans.push_back(deltaMean);
        uncertainties.push_back(uncertainty);
    }

    double wmean, wsigma;
    double totmean = 0;
    double totsigma = 0;
    double totweight = 0;
    for (int i = 0; i < 4; i++) {
        wmean = 1. / (stdmean[i] * stdmean[i]);
        wsigma = 1. / (stdsigma[i] * stdsigma[i]);
        totweight += wmean;
        totmean += deltaMeans[i] * wmean;
        totsigma += uncertainties[i] * wsigma;
    }

    totmean /= totweight;
    double finalsigma = 1. / std::sqrt(totsigma);

    allMeans.push_back(means);
    allSigmas.push_back(sigmas);
    allDeltaMeans.push_back(deltaMeans);
    allUncertainties.push_back(uncertainties);    
    results.push_back(totmean);
    rsigmas.push_back(finalsigma);
    
    }

    TCanvas *c2 = new TCanvas();
    TH1D *h1 = new TH1D("Parameters","Average Voltages - Gauss Fit",10,4.5,5.5); 
    for (int i=0; i<10; i++) h1->Fill(results[i] ,rsigmas[i]);
    h1->SetMarkerStyle(kFullCircle);
    h1->GetXaxis()->SetTitle("Average Delta Peak Voltage (V)");
    h1->GetYaxis()->SetTitle("Weighted Occurence");
    fnew->SetRange(4.9, 5.2);
    h1->Fit("fnew", "Q");
    //h1->Draw("HIST same E");

    
    double result = fnew->GetParameter(1);
    double sd_result = fnew->GetParameter(2);
        for(int i = 0; i < 10; i++){
        cout << results[i] << " +- " << rsigmas[i] << endl;
    } 

    cout << "First Excitation Energy = " << result << " +- " 
    << sd_result << " eV" <<endl;
    c2->Print("Franck-Hertz.pdf")    


}
