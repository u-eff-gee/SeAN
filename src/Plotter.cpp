#include "Plotter.h"
#include "Config.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TStyle.h"

#include <sstream>

using std::stringstream;

void Plotter::plot1DHistogram(const vector<double> &bins, const vector<double> &histogram, const string name, const string xaxis_label, const string yaxis_label){
	stringstream filename;
	filename << OUTPUT_DIR << name << ".pdf";
	stringstream canvasname;
	canvasname << name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), name.c_str(), 0, 0, 800, 500);
	
	TGraph *graph = new TGraph((int) bins.size(), &bins[0], &histogram[0]);
	graph->SetTitle(name.c_str());
	graph->GetHistogram()->GetXaxis()->SetTitle(xaxis_label.c_str());
	graph->GetHistogram()->GetYaxis()->SetTitle(yaxis_label.c_str());
	graph->Draw();

	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Plotter::plotMultiple1DHistograms(const vector<double> &bins, const vector< vector<double> > &histograms, const string name, const string xaxis_label, const string yaxis_label){
	stringstream filename;
	filename << OUTPUT_DIR << name << ".pdf";
	stringstream canvasname;
	canvasname << name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), name.c_str(), 0, 0, 800, 500);
	
	for(unsigned int i = 0; i < histograms.size(); ++i){
		TGraph *graph = new TGraph((int) bins.size(), &bins[0], &histograms[i][0]);
		if(i == 0){
			graph->SetTitle(name.c_str());
			graph->GetHistogram()->GetXaxis()->SetTitle(xaxis_label.c_str());
			graph->GetHistogram()->GetYaxis()->SetTitle(yaxis_label.c_str());

			graph->Draw();
		} else{
			graph->Draw("same");
		}
	}

	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Plotter::plotMultiple1DHistogramsAndSum(const vector<double> &bins, const vector< vector<double> > &histograms, const string name, const string xaxis_label, const string yaxis_label){
	vector<double> sum(bins.size(), 0.);

	for(unsigned int i = 0; i < histograms.size(); ++i){
		for(unsigned int j = 0; j < bins.size(); ++j){
			sum[j] += histograms[i][j];
		}	
	}	

	stringstream filename;
	filename << OUTPUT_DIR << name << ".pdf";
	stringstream canvasname;
	canvasname << name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), name.c_str(), 0, 0, 800, 500);

	TGraph *graph = new TGraph((int) bins.size(), &bins[0], &sum[0]);
	graph->SetTitle(name.c_str());
	graph->GetHistogram()->GetXaxis()->SetTitle(xaxis_label.c_str());
	graph->GetHistogram()->GetYaxis()->SetTitle(yaxis_label.c_str());

	graph->Draw();
	
	for(unsigned int i = 0; i < histograms.size(); ++i){
		TGraph *graph = new TGraph((int) bins.size(), &bins[0], &histograms[i][0]);
		graph->SetLineStyle(2);
		graph->Draw("same");
	}

	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}

void Plotter::plot2DHistogram(const vector<double> &bins1, const vector<double> &bins2, const vector< vector<double> > &histogram, const string name, const string xaxis_label, const string yaxis_label, const string zaxis_label){

	stringstream filename;
	filename << OUTPUT_DIR << name << ".pdf";
	stringstream canvasname;
	canvasname << name << "_canvas";

	unsigned int n1 = (unsigned int) bins1.size();
	unsigned int n2 = (unsigned int) bins2.size();

	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), name.c_str(), 0, 0, 800, 500);

	TGraph2D *graph = new TGraph2D(((int) n1)*((int) n2));	
	int npoint = 0;
	for(unsigned int i = 0; i < n1; ++i){
		for(unsigned int j = 0; j < n2; ++j){
			graph->SetPoint(npoint, bins1[i], bins2[j], histogram[i][j]);
			++npoint;
		}
	}

	graph->SetTitle(name.c_str());
	graph->GetHistogram()->GetXaxis()->SetTitle(xaxis_label.c_str());
	graph->GetHistogram()->GetYaxis()->SetTitle(yaxis_label.c_str());
	graph->GetHistogram()->GetZaxis()->SetTitle(zaxis_label.c_str());
	gStyle->SetPalette(55);
	graph->Draw("surf1");

	canvas->SetTheta(45.);
	canvas->SetPhi(250.);

	canvas->SaveAs(filename.str().c_str());
	delete canvas;
	delete graph;
}
