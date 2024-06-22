#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>

{

double const rad = 0.2;
int const ncoil = 154;
float const mu0 = 1.257*pow(10, -6);
double const c = pow(0.8, 1.5);
float bc = (c*mu0*ncoil)/rad;
float qmtrue = 1.75*pow(10,11);

float rsigma = 0.002;
float vsigma = 1;
float isigma = 0.01;
float bsigma = 0.00000692;

TF1 *f1 = new TF1("f1", "[0] + [1]*x");
gStyle->SetOptFit(1);
gStyle->SetStatX(0.48);
gStyle->SetStatY(0.9);

TTree *t = new TTree("t", "t");
t->ReadFile("data.csv");
int n = t->GetEntries();
float * x, * y, * sx, * sy ;
x = new float[4]; 
y = new float[4]; 
sx = new float[4]; 
sy = new float[4];
float r,v,i;
t->SetBranchAddress("I", &i);
t->SetBranchAddress("V", &v);
std::vector<double> radius = {0.02, 0.03, 0.04, 0.05};   

std::vector<float> results;
std::vector<float> resultsigmas;


for (int j=0; j<2; j++){

for (int k = 0; k < 4; ++k) {
    t->GetEntry(k+j*4);
    x[k] = bc*bc*i*i;
    y[k] = (2*v)/(radius[k]*radius[k]);
    sx[k] = 2*bc*i*bsigma;
    sy[k] = sqrt(pow(((2*vsigma)/radius[k]*radius[k]),2) + pow((4*v*rsigma)/pow(radius[k],3),2));
}

TCanvas *c1 = new TCanvas();
TGraphErrors gr(4, x, y, sx, sy);
f1->SetParameter(1, 1.7*pow(10,11));
gr.Fit(f1,"S");
results.push_back(f1->GetParameter(1));
resultsigmas.push_back(f1->GetParError(1));
gr.Draw("A*");
f1->SetLineColor(kBlue);
gr.GetXaxis()->SetTitle("B^2 (T^2)");
gr.GetYaxis()->SetTitle("2V/r^2 (V/m^2)");
gr.SetTitle("Constant Electric Field");
std::string filename = "constE_" + std::to_string(j+1) + ".pdf";
c1->Print(filename.c_str());
}

for (int j=2; j<4; j++){

for (int k = 0; k < 4; ++k) {
    t->GetEntry(k+j*4);
    x[k] = bc*bc*i*i*(radius[k]*radius[k]);
    y[k] = (2*v);
    sx[k] = 2*bc*i*radius[k]*sqrt(pow((radius[k]*bsigma),2) + pow((bc*i*rsigma),2));
    sy[k] = 2*vsigma;
}

TCanvas *c2 = new TCanvas();
gStyle->SetOptFit(1);
TGraphErrors gr(4, x, y, sx, sy);
f1->SetParameter(1, 1.7*pow(10,11));
gr.Fit(f1,"S");
results.push_back(f1->GetParameter(1));
resultsigmas.push_back(f1->GetParError(1));
gr.Draw("A*");
f1->SetLineColor(kBlue);
gr.GetXaxis()->SetTitle("B^2*r^2 (T^2*m^2)");
gr.GetYaxis()->SetTitle("2V (V)");
gr.SetTitle("Constant Current");   
std::string filename = "constI_" + std::to_string(j-1) + ".pdf";
c2->Print(filename.c_str());
}

float uppersum1, uppersum2;
float lowersum1, lowersum2;

for (int i=0; i<2; i++){

uppersum1 += results[i]/pow(resultsigmas[i], 2);
lowersum1 += 1/pow(resultsigmas[i], 2);

}

for (int i=2; i<4; i++){

uppersum2 += results[i]/pow(resultsigmas[i], 2);
lowersum2 += 1/pow(resultsigmas[i], 2);

}

float finalqm1 = uppersum1/lowersum1;
float finalsigma1 = sqrt(1/lowersum1);

float finalqm2 = uppersum2/lowersum2;
float finalsigma2 = sqrt(1/lowersum2);
std::cout << "The Charge to Mass ratio with constant electric field is found to be: " << finalqm1 << " C/kg +- " << finalsigma1 << " C/kg" <<std::endl;
std::cout << "It is " << (qmtrue-finalqm1)/finalsigma1 << " away from true value." << std::endl;

std::cout << "The Charge to Mass ratio with constant current is found to be: " << finalqm2 << " C/kg +- " << finalsigma2 << " C/kg" <<std::endl;
std::cout << "It is " << (qmtrue-finalqm2)/finalsigma2 << " sigma away from true value." << std::endl;


}
