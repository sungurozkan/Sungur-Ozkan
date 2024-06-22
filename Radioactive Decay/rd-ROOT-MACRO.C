#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <sstream>


{
const float sigmat = 0.3;
const float ln2 = 0.6931;
const float truevalue = 55.6;
std::vector<float> results;
std::vector<float> resultsigmas;

TF1 *f1 = new TF1("f1","[0]*exp([1]*x)");
TF1 *f2 = new TF1("f2", "[0] + [1]*x");

gStyle->SetOptFit(1);
gStyle->SetStatX(0.9);
gStyle->SetStatY(0.9);
    std::vector<std::string> filenames = {
        "converted_data1.csv",
        "converted_data2.csv",
        "converted_data3.csv",
        "converted_data4.csv",
        "converted_data5.csv",
        "converted_data6.csv",
        "converted_data7.csv",
        "converted_data8.csv",
        "converted_data9.csv",
        "converted_data10.csv",
        "converted_data11.csv",
        "converted_data12.csv",
        "converted_data13.csv",
        "converted_data14.csv",
        "converted_data15.csv",
        "converted_data16.csv"
    };

for (size_t fileIndex = 0; fileIndex < filenames.size(); ++fileIndex) {
    TTree *t = new TTree("t", "t");
    t->ReadFile(filenames[fileIndex].c_str());
    int n = t->GetEntries();
    float * x, * y, * sx, * sy ;
    x = new float[n-1]; 
    y = new float[n-1]; 
    sx = new float[n-1]; 
    sy = new float[n-1];
    float time, interval;
    t->SetBranchAddress("T", &time);
    t->SetBranchAddress("I", &interval);

    for (int i = 0; i < n - 1; ++i) {
        t->GetEntry(i);
        float time_i = time;
        t->GetEntry(i+1);
        float time_ip1 = time;
        x[i] = (time_i + time_ip1) / 2;
        y[i] = 1 / interval;
        sx[i] = sigmat;
        sy[i] = sigmat / (interval * interval);
    }
    std::string voltage = std::to_string((fileIndex / 4) * 0.5 + 3); // Voltage starts from 3V and increments by 0.5V
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << std::stof(voltage); // Convert to float and set precision to 1 decimal place
    std::string voltageString = oss.str();    
    std::string squeezes = std::to_string((fileIndex % 4) + 2); // Squeezes start from 2 and increments by 1
    std::string finalfilename = voltageString + "V_" + squeezes + "Squeezes_fitted.pdf";
    std::string finalgraphname = voltageString + "V and " + squeezes + " Squeezes";

    TGraphErrors gr(n-1, x, y, sx, sy);
    gr.GetXaxis()->SetTitle("Time [s]");
    gr.GetYaxis()->SetTitle("1/Intervals [1/s]");
    gr.SetTitle(finalgraphname.c_str());
    f1->SetRange(0, 500);
    f1->SetParameter(0, 0.25);
    gr.Fit(f1, "S");
    TCanvas *c1 = new TCanvas();
    gr.Draw("A*");
    // Constructing the filename

    c1->Print(finalfilename.c_str());
    float lambda = f1->GetParameter(1);
    float lambdasigma = f1->GetParError(1);
    results.push_back(lambda);
    resultsigmas.push_back(lambdasigma);
}

float uppersum;
float lowersum;

for (int i=0; i<16; i++){

uppersum += results[i]/pow(resultsigmas[i], 2);
lowersum += 1/pow(resultsigmas[i], 2);

}

float finallambda = uppersum/lowersum;
float finallambdasigma = sqrt(1/lowersum);

float halflife = -ln2/finallambda;
float halflifesigma = (ln2*finallambdasigma)/(finallambda*finallambda);

std::cout << "Half-life: " << halflife << "s +- " << halflifesigma << "s" << std::endl;
std::cout << "Our result is " << (halflife-truevalue)/halflifesigma << "sigma away from true value." << std::endl;
float * u, * v, *z, *q;
u = new float[4]; 
v = new float[4]; 
z = new float[4]; 
q = new float[4]; 

std::vector<float> squeezeconst;
std::vector<float> voltageconst;

for (int i=0; i<4; i++){
squeezeconst.push_back(results[i]);
voltageconst.push_back(results[i*4]);
}
for (int i=0; i<4; i++){
v[i] = squeezeconst[i];
z[i] = voltageconst[i];
u[i] = i+2;
q[i] = 0.5*i+3.0;
}


TCanvas *c11 = new TCanvas();
TGraphErrors gr1(4, u, v);
gr1.GetXaxis()->SetTitle("Squeeze Amount");
gr1.GetYaxis()->SetTitle("Decay Constant");
gr1.SetTitle("Decay Constant As A Function Of Squeezes (3.0V)");
gr1.Draw("A*");
c11->Print("squeeze.pdf");


TCanvas *c12 = new TCanvas();
TGraphErrors gr2(4, q, z);
gr2.GetXaxis()->SetTitle("Voltage (V)");
gr2.GetYaxis()->SetTitle("Decay Constant");
gr2.SetTitle("Decay Constant As A Function Of Voltage (2 Squeezes)");


gr2.SetMarkerStyle(20); 
gr2.Draw("A*");
c12->Print("voltage.pdf");

}