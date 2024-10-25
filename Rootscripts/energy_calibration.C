#include <iostream>
#include <vector>
#include <string>


#include "include/strtools.h"


using namespace std;

// Get Hists from a file given a list of hist names
vector<TH1F*> getHists(TFile* f, vector<string>& vnames)
{
	vector<TH1F*> vh;
	for( auto it = vnames.begin(); it != vnames.end(); it++ ) {
		vh.push_back( (TH1F*)f->Get( (*it).c_str() ) );
	}
	return vh;
}

// Combine histograms (nomarlized by duration)
TH1F* add_histograms(vector<TH1F*>& vh)
{
	TH1F* h = NULL;
	for(int i = 0; i < vh.size(); i++){
		if(h == NULL) {
			h = (TH1F*)(vh[i])->Clone();
		}
		else {
			h->Add(vh[i]);
		}
	}
	return h;
}


// Takes in the root files that contain the energy calibration runs,
// fits the histograms, and extracts the mean and counts.
// This is done for both gam1 and gam2
// Cs137 has one peaks, (+one compton edge)
// Na22 has two peaks   (+two compton edge)
void energy_calibration(string sfNa, string sNaDur, string sfCs, string sCsDur, string sfout)
{
	// Open Save Folder
	TFile* fsave = new TFile(sfout.c_str(), "RECREATE");

	// Parse file names
	vector<string> vNa = parse_str(sfNa, '\n');
	vector<string> vCs = parse_str(sfCs, '\n');

	// Parse Durations
	vector<string> vNaDur = parse_str(sNaDur); vector<int> vNaDurI; for(auto it = vNaDur.begin(); it!=vNaDur.end(); it++) { vNaDurI.push_back(stoi(*it)); }
	vector<string> vCsDur = parse_str(sCsDur); vector<int> vCsDurI; for(auto it = vCsDur.begin(); it!=vCsDur.end(); it++) { vCsDurI.push_back(stoi(*it)); }



	// Open Up Na Files
	vector<TFile*> vfNa(vNa.size());
	for(int file = 0; file < vNa.size(); file++ ){
		vfNa[file] = new TFile(vNa[file].c_str(), "READ");
	}

	// Open Up Cs Files
	vector<TFile*> vfCs(vCs.size());
	for(int file = 0; file < vCs.size(); file++ ){
		vfCs[file] = new TFile(vCs[file].c_str(), "READ");
	}

	// If we have several Na Files open
	// Add gam1 and gam2 spectra together
	vector<string> vhnames = {"gam1", "gam2"};
	vector<TH1F*> hNa_gam1;
	vector<TH1F*> hNa_gam2;
	// Loop over files
	for(auto it = vfNa.begin(); it != vfNa.end(); it++)
	{
		// Get Histograms of interest
		vector<TH1F*> vh = getHists(*it, vhnames);
		hNa_gam1.push_back( vh[0] );
		hNa_gam2.push_back( vh[1] );
	}


	// If we have several Cs Files open
	// Add gam1 and gam2 spectra together
	vector<TH1F*> hCs_gam1;
	vector<TH1F*> hCs_gam2;
	// Loop over files
	for(auto it = vfCs.begin(); it != vfCs.end(); it++)
	{
		// Get Histograms of interest
		vector<TH1F*> vh = getHists(*it, vhnames);
		hCs_gam1.push_back( vh[0] );
		hCs_gam2.push_back( vh[1] );
	}

	// Get Added Hists
	TH1F* hNaGam1 = add_histograms(hNa_gam1);
	TH1F* hNaGam2 = add_histograms(hNa_gam2);
	TH1F* hCsGam1 = add_histograms(hCs_gam1);
	TH1F* hCsGam2 = add_histograms(hCs_gam2);
	
	// Save Hists
	fsave->cd();
	hNaGam1->Write("Na_Gam1_Added");
	hNaGam2->Write("Na_Gam2_Added");
	hCsGam1->Write("Cs_Gam1_Added");
	hCsGam2->Write("Cs_Gam2_Added");

	// Lets Fit the first one
	// hNaGam1->

}
