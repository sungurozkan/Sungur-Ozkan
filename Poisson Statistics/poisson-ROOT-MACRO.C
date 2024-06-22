#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>

{
gStyle->SetOptFit(1111);
TF1 *fnew = new TF1("fnew", "gaus");
TF1 *f2 = new TF1("f2", "[0]*TMath::Poisson(x,[1])",0,30);
f2->SetLineColor(kBlue);

TCanvas *c1 = new TCanvas();
TGraphErrors cal("data_Sheet1.txt");
cal.GetXaxis()->SetTitle("High Voltage [V]");
cal.GetYaxis()->SetTitle("Number Of Counts");
cal.SetTitle("Number Of Counts In 100 Second Intervals");
cal.SetLineColor(kBlack);
cal.SetMarkerStyle(20); 
cal.Draw();
c1->Print("calibration.pdf");

f2->SetParName(0, "Amplitude Scale");
f2->SetParName(1, "Alpha");
fnew->SetParName(1, "Gauss Fit Mean");
fnew->SetParName(2, "Gauss Fit Sigma");

TCanvas *c2 = new TCanvas();
TH1F *h1 = new TH1F("h1","Sezium - 440V - 1/sn",10, 5, 30);  
std::ifstream file1("data_Sheet2.txt");
h1->GetXaxis()->SetTitle("Number of Counts");
h1->GetYaxis()->SetTitle("Weight Of Counts");

float datum1;
while (file1>>datum1) h1->Fill(datum1);
double meanSezyum1 = h1->GetMean();
double stdSezyum1 = h1->GetStdDev();
fnew->SetRange(5,30);
h1->Fit("fnew", "S");
f2->SetParameters(h1->GetMaximum(), h1->GetMean());
f2->SetRange(5,30);
h1->Fit("f2", "S+");
h1->Draw();
fnew->Draw("same");
f2->Draw("same");
double chiSquare11 = fnew->GetChisquare();
int ndf11 = fnew->GetNDF();
double chindf11 = chiSquare11 / ndf11;
double alpha1 = f2->GetParameter(1);

double chiSquare12 = f2->GetChisquare();
int ndf12 = f2->GetNDF();
double chindf12 = chiSquare12 / ndf12;
c2->Print("Seizum-1.pdf");



TCanvas *c3 = new TCanvas();
TH1F *h2 = new TH1F("h2","Sezium - 440V - 10/sn",10, 100, 180);
std::ifstream file2("data_Sheet3.txt");
h2->GetXaxis()->SetTitle("Number of Counts");
h2->GetYaxis()->SetTitle("Weight Of Counts");
float datum2;
while (file2>>datum2) h2->Fill(datum2);
double meanSezyum10 = h2->GetMean();
double stdSezyum10 = h2->GetStdDev();
fnew->SetRange(100,180);
h2->Fit("fnew", "S");
f2->SetParameters(h2->GetMaximum(), h2->GetMean());
f2->SetRange(100,180);
h2->Fit("f2", "S+");
h2->Draw();
fnew->Draw("same");
f2->Draw("same");
double chiSquare21 = fnew->GetChisquare();
int ndf21 = fnew->GetNDF();
double chindf21 = chiSquare21 / ndf21;

double chiSquare22 = f2->GetChisquare();
int ndf22 = f2->GetNDF();
double chindf22 = chiSquare22 / ndf22;
c3->Print("Seizum-10.pdf");
double alpha2 = f2->GetParameter(1);




TCanvas *c4 = new TCanvas();
TH1F *h3 = new TH1F("h3","Barium - 440V - 1/sn",10, 0, 10);
h3->GetXaxis()->SetTitle("Number of Counts");
h3->GetYaxis()->SetTitle("Weight Of Counts");
std::ifstream file3("data_Sheet4.txt");
float datum3;
while (file3>>datum3) h3->Fill(datum3);
h3->Rebin(10);
double meanBarium1 = h3->GetMean();
double stdBarium1 = h3->GetStdDev();
fnew->SetRange(0,10);
h3->Fit("fnew", "S");
f2->SetParameters(h3->GetMaximum(), h3->GetMean());
f2->SetRange(0,10);
h3->Fit("f2", "S+");
h3->Draw();
fnew->Draw("same");
//f2->Draw("same");
double chiSquare3 = fnew->GetChisquare();
int ndf3 = fnew->GetNDF();
double chindf3 = chiSquare3 / ndf3;
double chiSquare32 = f2->GetChisquare();
int ndf32 = f2->GetNDF();
double chindf32 = chiSquare32 / ndf32;
c4->Print("Barium-1.pdf");
double alpha3 = f2->GetParameter(1);



TCanvas *c5 = new TCanvas();
TH1F *h4 = new TH1F("h4","Barium - 440V - 10/sn",10, 10, 40);
std::ifstream file4("data_Sheet5.txt");
h4->GetXaxis()->SetTitle("Number of Counts");
h4->GetYaxis()->SetTitle("Weight Of Counts");
float datum4;
while (file4>>datum4) h4->Fill(datum4);

double meanBarium10 = h4->GetMean();
double stdBarium10 = h4->GetStdDev();
fnew->SetRange(10,40);
h4->Fit("fnew", "S");
f2->SetParameters(h4->GetMaximum(), h4->GetMean());
f2->SetRange(10,40);
h4->Fit("f2", "S+");
h4->Draw();
fnew->Draw("same");
f2->Draw("same");
double chiSquare4 = fnew->GetChisquare();
int ndf4 = fnew->GetNDF();
double chindf4 = chiSquare4 / ndf4;
double chiSquare42 = f2->GetChisquare();
int ndf42 = f2->GetNDF();
double chindf42 = chiSquare42 / ndf42;
c5->Print("Barium-10.pdf");
double alpha4 = f2->GetParameter(1);


double means[] = {meanSezyum1, meanSezyum10, meanBarium1, meanBarium10};
double sigmas[] = {stdSezyum1, stdSezyum10, stdBarium1, stdBarium10};
double sqr[4];
   for (int i = 0; i < 4; i++) {
        double x = std::sqrt(means[i]) / sigmas[i];
        sqr[i] = x;
    }

TCanvas *c6 = new TCanvas();
TGraphErrors gr(4, sigmas, sqr);
gr.SetTitle("Sqrt(mu) / Sigma as a function of Sigma");
gr.SetLineColor(kBlack);
gr.SetMarkerColor(kBlue);
gr.SetMarkerStyle(30); 
gr.GetXaxis()->SetTitle("sigma");
gr.GetYaxis()->SetTitle("sqrt(mu) / sigma");
TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x",0,12); 
fitFunc->SetParameter(0, initialGuessForParameter0);
fitFunc->SetParameter(1, initialGuessForParameter1);
gr.Fit("fitFunc", "R");
gr.Draw("A*");
fitFunc->Draw("same");
c6->Print("line.pdf");


// PART2


TF1 *f3 = new TF1("f3", "[1]*[0]*TMath::Exp(-[0]*x)");
TF1 *f4 = new TF1("f4", "[1]*[0]*[0]*x*TMath::Exp(-[0]*x)");
f3->SetParName(0, "Alpha");
f4->SetParName(0, "Alpha");
TCanvas *c7 = new TCanvas();
TH1F *h5 = new TH1F("h5","n=0",10, 0.0, 10.0);
h5->SetTitle("Histogram for n=0");

h5->GetXaxis()->SetTitle("Time Difference Between Successive Pulses");
h5->GetYaxis()->SetTitle("Frequency of Time Differences");
std::ifstream file5("data_Sheet9.txt");
float datum5;
while (file5>>datum5) h5->Fill(datum5);
f3->SetParameters(h5->GetMean(),h5->GetMaximum());
f3->SetRange(0.1 ,10);
h5->Fit("f3");
h5->Draw();
c7->Print("n0.pdf");

TCanvas *c8 = new TCanvas();
TH1F *h6 = new TH1F("h6","n=1",10, 0.0, 10.0);
h6->GetXaxis()->SetTitle("Time Difference Between Successive Pulses");
h6->GetYaxis()->SetTitle("Frequency of Time Differences");
h6->SetTitle("Histogram for n=1");
std::ifstream file6("data_Sheet10.txt");
float datum6;
while (file6>>datum6) h6->Fill(datum6);
f4->SetParameters(h6->GetMean(), h6->GetMaximum());
f4->SetRange(0,10);
h6->Fit("f4");
h6->Draw();
c8->Print("n1.pdf");






}