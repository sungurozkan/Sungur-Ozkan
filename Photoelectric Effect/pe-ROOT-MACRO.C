#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>

{

std::vector<std::string> filenames = {
    "blue.csv",
    "green.csv",
    "yellow.csv",
    "violet.csv",
    "turquoise.csv",
    };

std::vector<float> xintercepts;
std::vector<float> errors;
float trueplanck = 6.626*pow(10,-34);

gStyle->SetOptFit(1);
gStyle->SetStatX(0.9);
gStyle->SetStatY(0.9);
TF1 *f1 = new TF1("line", "[0] + [1]*x");


for (int fileIndex = 0; fileIndex < filenames.size(); ++fileIndex) {
TTree *t = new TTree("t", "t");
t->ReadFile(filenames[fileIndex].c_str());
int n = t->GetEntries();
float * x, * y, * sx, * sy ;
x = new float[n]; 
y = new float[n]; 
sx = new float[n]; 
sy = new float[n];
float a,v;
t->SetBranchAddress("A", &a);
t->SetBranchAddress("V", &v);
string color = filenames[fileIndex].substr(0, filenames[fileIndex].find_last_of('.'));
string finalfilename = color + "fitted.pdf"; 

for (int i = 0; i < n; ++i) {
    t->GetEntry(i);
    x[i] = v * pow(10, -3);
    y[i] = a * pow(10, -14);
    sx[i] = 0.1 * pow(10, -3);
    sy[i] = 0.1 * pow(10, -14);
}

TGraphErrors gr(n, x, y, sx, sy);

gr.SetTitle(("Voltage vs Current for " + color).c_str());
gr.GetXaxis()->SetTitle("Voltage(V)");
gr.GetYaxis()->SetTitle("Current(A)");

f1->SetRange(0,0.8);
gr.Fit(f1,"R");

float yintercept1 = f1->GetParameter(0);
float slope1 = f1->GetParameter(1);
float sigmayintercept1 = f1->GetParError(0);
float sigmaslope1 = f1->GetParError(1);

f1->SetRange(1,3.2);
gr.Fit(f1,"R+");

float yintercept2 = f1->GetParameter(0);
float slope2 = f1->GetParameter(1);
float sigmayintercept2 = f1->GetParError(0);
float sigmaslope2 = f1->GetParError(1);

float xintercept = (yintercept2 - yintercept1)/(slope1 - slope2);
xintercepts.push_back(xintercept);

float term1 = pow((yintercept2 - yintercept1) / pow((slope1 - slope2), 2), 2) * pow(sigmaslope1, 2);
float term2 = pow((yintercept2 - yintercept1) / pow((slope1 - slope2), 2), 2) * pow(sigmaslope2, 2);
float term3 = pow(1 / (slope1 - slope2), 2) * pow(sigmayintercept1, 2);
float term4 = pow(1 / (slope1 - slope2), 2) * pow(sigmayintercept2, 2);

float sigma_Vs = sqrt(term1 + term2 + term3 + term4);
errors.push_back(sigma_Vs);

cout << "Stopping potential for " << color << " is equal to " << xintercept << "V +- " << sigma_Vs << "V"<< endl;

float freqs[5] = {5.19, 5.49, 6.08, 6.88, 7.41};

TCanvas *c1 = new TCanvas();
gr.Draw("A*");
c1->Print(finalfilename.c_str());




}

const int ndata = 5;
double y2[ndata] = { 0.786, 0.407, 0.347, 0.404, 0.288};
double x2[ndata] = {5.19, 5.49, 6.08, 6.88, 7.41};
for (int i=0; i<5; i++){
    x2[i] = x2[i]*pow(10, 14);
}
double sy2[ndata] = {0.0026358, 0.0122895, 0.00611123, 0.0034491, 0.00693817};
double sx2[ndata] = {0,0,0,0,0};


TGraphErrors gr2(ndata, x2, y2, sx2, sy2);
gr2.SetTitle("Stopping Voltage vs Frequency");
gr2.GetYaxis()->SetTitle("Voltage(V)");
gr2.GetXaxis()->SetTitle("Frequency(Hz)");
f1->SetRange(1,1*pow(10,15));
gr2.Fit(f1,"R");

float wf = f1->GetParameter(0);
float sigmawf = f1->GetParError(0);
float hq = f1->GetParameter(1);
float sigmahq = f1->GetParError(1);

float q = (-1)*1.602*pow(10,-19);
float h = hq*q;
float error = (trueplanck-h)/(sigmahq*q);

cout << "Planck Constant: " << h << " J/Hz +- " << sigmahq*q << " J/Hz" << endl;
cout << "Work Function: " << wf << "V +- " << sigmawf << "V" << endl;
cout << "Planck Constant is" <<  error << " sigma away from true value" << endl;
TCanvas *c2 = new TCanvas();
gr2.Draw("A*");
c2->Print("freq.pdf");

}