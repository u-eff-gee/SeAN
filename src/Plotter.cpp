#include "Plotter.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

#include <sstream>

using std::stringstream;

void Plotter::plot1DHistogram(vector<double> &bins, vector<double> &histogram, string name){
	stringstream filename;
	filename << PLOT_OUTPUT_DIR << name << ".pdf";
	stringstream canvasname;
	canvasname << name << "_canvas";
	TCanvas *canvas = new TCanvas(canvasname.str().c_str(), name.c_str(), 0, 0, 800, 500);
	//TLegend *legend = new TLegend(CS_PLOT_LEGEND_X1, CS_PLOT_LEGEND_Y1, CS_PLOT_LEGEND_X2, CS_PLOT_LEGEND_Y2);
	
	TGraph *graph = new TGraph(bins.size(), &bins[0], &histogram[0]);
	graph->Draw();

	//crossSection->plot_crosssection(energy_bins, crosssection_bins, target_name, canvas, legend, "Cross section");

	// legend->Draw();
	canvas->SaveAs(filename.str().c_str());
	delete canvas;
}
