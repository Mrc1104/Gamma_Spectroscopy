#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>

#include "include/strtools.h"


using namespace std;

// Attenuation Coefficients
#define AL50  0.2279
#define AL125 0.1483
#define CU50  0.7492
#define CU125 0.4714



// Perform the background reduction
TH1F* background_sub(TH1F* h, TH1F* hbkg, double hdur = 10.0, double bkgdur=1.0, int xlow = 0, int xhigh=1023)
{
	TH1F* hret = (TH1F*)h->Clone();
	hret->Add(hbkg, -1*hdur/bkgdur);
	// Do proper error propogation
	for(int i = xlow; i < xhigh; i++) {
		double herr   = h->GetBinContent(i);
		double bkgerr = hbkg->GetBinContent(i);
		double scale  = hdur/bkgdur ;
		double err    = TMath::Sqrt(herr+scale*bkgerr);
		hret->SetBinError(i, err);
	}
	return hret;
}

// Fit peaks
void fit_peaks(TFile* file, TTree* t, TCanvas* c, const int pad, TH1F* h0, const int xlow, const int xhigh, const int ylow, const int yhigh, string name)
{
	c->cd(pad);
	TH1F* h = (TH1F*)h0->Clone();
	TLine *llow = new TLine(xlow, ylow, xlow, yhigh); llow->SetLineWidth(2); llow->SetLineColor(kBlack);
	TLine *lhigh = new TLine(xhigh, ylow, xhigh, yhigh); lhigh->SetLineWidth(2); lhigh->SetLineColor(kBlack);
	h->SetTitle(name.c_str());
	h->SetLineColor(kBlue);
	h->GetXaxis()->SetRange(xlow,xhigh+1);
	h->Draw();
	// llow->Draw("SAME"); lhigh->Draw("SAME");

	// Do a trial fit
	TF1* fitgaus = new TF1("fitgaus", "gaus(0)");
	h->Fit(fitgaus, "q");
	double* gptmp = fitgaus->GetParameters();
	double gp[5] = { gptmp[0],gptmp[1],gptmp[2], 0, 0 };

	// Fit with a gaus + pol1 using trial fit
	TF1* fit = new TF1("fit", "gaus(0)");
	fit->SetParameters(gp);
	TFitResultPtr pfr = h->Fit(fit,"qS");
	double *fp = fit->GetParameters();

	// Get Counts
	double peak_counts = h->Integral(xlow, xhigh);
	double peak_error = 0;
	for(int bin = xlow; bin <= xhigh; bin++ ) {
		peak_error += pow(h->GetBinError(bin), 2);
	}
	cout << "PEAK_ERROR: " << TMath::Sqrt(peak_error) << endl;


	double adjusted_counts = peak_counts;
	// if(adjusted_counts < 0) adjusted_counts = 0;
	double adjusted_counts_error = TMath::Sqrt( peak_error);

	if( name.find("GAM1_") != string::npos)
		cout << name << "\t|\tPeak Counts: " << peak_counts << "\tAdjusted: " << adjusted_counts << "\tAdjusted Error: " << adjusted_counts_error << endl;

	// Fill tree with peak info
	double bvals[] = { (double)xlow, (double)xhigh,fp[1], fp[2], adjusted_counts,adjusted_counts_error};
	auto branch = t->Branch(name.c_str(), &bvals, Form("%s/D:%s/D:%s/D:%s/D:%s/D:%s/D", "xlow", "xhigh", "mean", "sigma", "counts", "cts_error" ) );
	branch->Fill();
	t->Fill();
	
	// save
	file->cd();
	pfr->Write(Form("%s_fit", name.c_str()));
	h->Write(Form("%s_hist", name.c_str()));


	return;
}



vector<TH1F*> scale_hists(const vector<TH1F*> &vh, const vector<double> &vt, int xlow = 0, int xhigh=1023)
{
	vector<TH1F*> vret;
	if( vh.size() != vt.size() ){
		cout << "ERROR: SCALE_HISTS() EXPECTS EQUAL SIZED VECTORS!\nEXITING..." << endl;
		exit(-1);
	}
	for( int i = 0; i < vh.size(); i++ ) {
		vret.push_back( (TH1F*)vh[i]->Clone() );
		vret.back()->Scale(10.0 / vt[i] );
		for(int bin = xlow; bin < xhigh; bin++) {
			vret.back()->SetBinError(bin, TMath::Sqrt( vret.back()->GetBinContent(bin) * 10.0 / vt[i] ) );
		}
	}
	return vret;
}

// Get Hists from a file given a list of hist names
vector<TH1F*> getHists(TFile* f, vector<string>& vnames)
{
	vector<TH1F*> vh;
	for( auto it = vnames.begin(); it != vnames.end(); it++ ) {
		vh.push_back( (TH1F*)f->Get( (*it).c_str() ) );
	}
	return vh;
}



// string sdatafile -- path to data root file
// string spdata    -- path to txt files with thickness and run duration
// string sbkg      -- path to background rootfile
// string soutfile  -- path to save rootfile
// int elem         -- 0 (Al) 1 (Cu)
void attenuation(string sdatafile, string spdata, string sbkg, string soutfile, int elem)
{
	// Parse Elem
	if( elem != 0 && elem != 1 ) {
		cout << "ERROR: Not given a valid element value!\nEnter 0 (Al) or 1 (Cu)\nExiting..." << endl;
		exit(-1);
	}
	double coeff[2];
	if( elem == 0 ) { coeff[0] = AL50; coeff[1] = AL125; }
	if( elem == 1 ) { coeff[0] = CU50; coeff[1] = CU125; }

	// Parse data files and duration
	vector<string> vsdatafile = parse_str(sdatafile,'\n');

	// Open Up Files
	vector<TFile*> vfiles;
	for( auto it = vsdatafile.begin(); it != vsdatafile.end(); it++ )
	{
		vfiles.push_back( new TFile( (*it).c_str(), "READ") );
		vfiles.back()->Print();
	}

	// Get Histograms
	vector<string> vhnames = {"gam1", "gam2"};
	vector<TH1F*> vgam1;
	vector<TH1F*> vgam2;
	for( auto it = vfiles.begin(); it != vfiles.end(); it++ ) {
		vector<TH1F*>  vhists  = getHists( *it, vhnames );
		vgam1.push_back(vhists[0]);
		vgam2.push_back(vhists[1]);
	}

	// Read in file
	double pdata_time, pdata_th, pdata_therr;
	vector<double> vtime;
	// vector<pair<double,double>> vthickness;
	double vthick[255];
	double vtherr[255]; int kused = 0;
	ifstream fpdata(spdata);
	while(fpdata >> pdata_time >> pdata_th >> pdata_therr) {
		// cout << pdata_time << " " << pdata_th << " " << pdata_therr << endl;
		vtime.push_back(pdata_time);
		vthick[kused] = pdata_th;
		vtherr[kused] = pdata_therr;
		kused++;
	}


	// Scale Histograms
	vector<TH1F*> vgam1_scaled = scale_hists(vgam1, vtime);
	vector<TH1F*> vgam2_scaled = scale_hists(vgam2, vtime);

	// Remove Background
	TFile* fbkg = new TFile(sbkg.c_str(), "READ"); fbkg->Print();
	vector<TH1F*> vbkg = getHists( fbkg, vhnames );
	vector<TH1F*> vgam1_scaled_bkg;
	vector<TH1F*> vgam2_scaled_bkg;
	for(int i = 0; i < vgam1_scaled.size(); i++) {
		vgam1_scaled_bkg.push_back( background_sub(vgam1_scaled[i], vbkg[0], 10.0, 5.0) );
		vgam2_scaled_bkg.push_back( background_sub(vgam2_scaled[i], vbkg[1], 10.0, 5.0) );
	}

	// Fit the 511 & 1275 peaks
	int offset = 0;
	if(vgam1_scaled_bkg.size() % 2) offset = 1;
	TFile* fsave   = new TFile(soutfile.c_str(), "RECREATE");
	TTree* tgam1   = new TTree("Gam1", "Gam1 Data"  ); TCanvas* cgam1  = new TCanvas("cgam1"); cgam1->Divide(vgam1_scaled_bkg.size()/2, vgam1_scaled_bkg.size()/2+offset);
	TTree* tgam2   = new TTree("Gam2", "Gam2 Data"  ); TCanvas* cgam2  = new TCanvas("cgam2"); cgam2->Divide(vgam2_scaled_bkg.size()/2, vgam2_scaled_bkg.size()/2+offset);
	TTree* tgam1_2 = new TTree("Gam1_2", "Gam1 Data"); TCanvas* cgam1_2= new TCanvas("cgam1_2"); cgam1_2->Divide(vgam1_scaled_bkg.size()/2, vgam1_scaled_bkg.size()/2+offset);
	pair<int,int> gam1_range  = { 120, 170 };
	pair<int,int> gam1_range_2= { 300, 390 };
	pair<int,int> gam2_range  = { 120, 180 };
	vector<string> vgam1_names;
	vector<string> vgam2_names;
	vector<string> vgam1_2_names;
	for(int i = 0; i < vgam1_scaled_bkg.size(); i++) {
		string tmp1 = "GAM1_" + to_string(i);
		string tmp2 = "GAM2_" + to_string(i);
		string tmp3 = tmp1    + "_1275";
		vgam1_names.push_back(tmp1);
		vgam2_names.push_back(tmp2);
		vgam1_2_names.push_back(tmp3);
		fit_peaks(fsave, tgam1, cgam1, i+1, vgam1_scaled_bkg[i], gam1_range.first, gam1_range.second, 0,0, tmp1);
		fit_peaks(fsave, tgam2, cgam2, i+1, vgam2_scaled_bkg[i], gam2_range.first, gam2_range.second, 0,0, tmp2);
		fit_peaks(fsave, tgam1_2, cgam1_2, i+1, vgam1_scaled_bkg[i], gam1_range_2.first, gam1_range_2.second, 0,0, tmp3);
	// fit_peaks(fsave, tdata, cNa_Gam1, 2, hNa_Gam1_BKG, 300, 390, 0, 1600, "NA_GAM1_PEAK2");
	}
	tgam1->Write();
	tgam2->Write();
	tgam1_2->Write();

	// Obtain Counts (ONLY FOR GAM1)
	double cts_511[255];
	double err_511[255];
	int npoints = 0;
	for( auto it = vgam1_names.begin(); it != vgam1_names.end(); it++ ) {
		const char *bname = (*it).c_str();
		TBranch* b = tgam1->GetBranch(bname);
		double dvar[6];
		b->SetAddress(&dvar);
		b->GetEntry();
		cts_511[npoints] = dvar[4];
		err_511[npoints] = dvar[5];
		// cout << "\n\n" << vthick[npoints] << " " << vtherr_511[npoints] << endl;
		// cout << cts_511[npoints] << " " << err_511[npoints] << endl;
		npoints++;
	}
	
	// Get Counts of 1275 peak (ONLY FOR GAM1)
	double cts_1275[255];
	double err_1275[255];
	int npoints_2 = 0;
	for( auto it = vgam1_2_names.begin(); it != vgam1_2_names.end(); it++ ) {
		const char *bname = (*it).c_str();
		TBranch* b = tgam1_2->GetBranch(bname);
		double dvar[6];
		b->SetAddress(&dvar);
		b->GetEntry();
		cts_1275[npoints_2] = dvar[4];
		err_1275[npoints_2] = dvar[5];
		// cout << "\n\n" << vthick[npoints_2] << " " << vtherr_1275[npoints_2] << endl;
		// cout << cts_1275[npoints_2] << " " << err_1275[npoints_2] << endl;
		npoints_2++;
	}

	// Plot the sucker
	string element;
	if( elem == 0 ) element = "Al";
	else            element = "Cu";

	TCanvas* c = new TCanvas("Attenuation"); c->Divide(1,2); c->cd(1);
	// 511
	TGraphErrors* g = new TGraphErrors(npoints, vthick, cts_511, vtherr, err_511);
	g->SetTitle(Form("%s Attenuation (511 keV); Thickness [cm]; Counts", element.c_str()));
	g->SetMarkerStyle(8);
	g->Draw("AP");
	TFitResultPtr fexpo1 = g->Fit("expo","S");

	// 1275
	c->cd(2);
	TGraphErrors* g_2 = new TGraphErrors(npoints_2, vthick, cts_1275, vtherr, err_1275);
	g_2->SetTitle(Form("%s Attenuation (1.275 MeV); Thickness [cm]; Counts", element.c_str()));
	g_2->SetMarkerStyle(8);
	g_2->Draw("AP");
	TFitResultPtr fexpo2 = g_2->Fit("expo","S");

	g->Write(Form("%s_Attenuation_511",element.c_str()));
	fexpo1->Write(Form("%s_Attenuation_511_Fit",element.c_str()));
	g_2->Write(Form("%s_Attenuation_1275",element.c_str()));
	fexpo2->Write(Form("%s_Attenuation_1275_Fit",element.c_str()));

	return;
}
