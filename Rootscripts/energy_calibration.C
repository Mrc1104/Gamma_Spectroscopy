#include <iostream>
#include <vector>
#include <string>
#include <cmath>


#include "include/strtools.h"


using namespace std;

struct info_t
{
	double val[3];
	double err[3];
};


// Calculate Residuals
void c_residuals(info_t *arrres, const info_t *arrtrue, const info_t *arrmean, const TFitResultPtr pfit, const int n)
{
	TMatrixDSym cov  = pfit->GetCovarianceMatrix();
	double p[2] = { pfit->Parameter(0), pfit->Parameter(1) }; // [0] -- intercept; [1] -- slope
	double pe[2]= { pfit->ParError(0) , pfit->ParError(1)  };

	double cov_ab = cov[0][1];
	for( int i = 0; i < n; i++ ) {
		double fiterr = pow(arrmean->val[i]*pe[1], 2) + pow( pe[0], 2) + 2*arrmean->val[i]*cov_ab;
		double fitenr = arrmean->val[i]*p[1]+p[0];
		arrres->val[i] = fitenr - arrtrue->val[i];
		arrres->err[i] = TMath::Sqrt( fiterr*fiterr + arrtrue->err[i]*arrtrue->err[i]   );
		cout << "Residuals: " << arrres->val[i] << "\t+-\t" << arrres->err[i] << endl;
	}

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
	h->Fit(fitgaus);
	double* gptmp = fitgaus->GetParameters();
	double gp[5] = { gptmp[0],gptmp[1],gptmp[2], 0, 0 };

	// Fit with a gaus + pol1 using trial fit
	TF1* fit = new TF1("fit", "gaus(0)+pol1(3)");
	fit->SetParameters(gp);
	TFitResultPtr pfr = h->Fit(fit,"S");

	// Get Counts
	double* fp = fit->GetParameters();
	double intercept = fp[3];
	double slope     = fp[4];
	double peak_counts = h->Integral(xlow, xhigh);
	double line_counts = 0.5*slope*( xhigh*xhigh - xlow*xlow) + intercept*( xhigh - xlow ); // This is the integral of mx+b
	double adjusted_counts = peak_counts-line_counts;
	cout << name << "\t|\tPeak Counts: " << peak_counts << "\tLine Counts: " << line_counts << "\tAdjusted: " << adjusted_counts << endl;

	// Fill tree with peak info
	double bvals[] = { (double)xlow, (double)xhigh,fp[1], fp[2], adjusted_counts };
	auto branch = t->Branch(name.c_str(), &bvals, Form("%s/D:%s/D:%s/D:%s/D:%s/D", "xlow", "xhigh", "mean", "sigma", "counts" ) );
	branch->Fill();
	t->Fill();
	
	// save
	file->cd();
	pfr->Write(Form("%s_fit", name.c_str()));
	h->Write(Form("%s_hist", name.c_str()));


	return;
}

// Populate TCanvas
void add_to_canvas(TCanvas* c, const int pad, vector<TH1F*> vh, vector<int> vc, vector<string> vname, string name)
{
	c->cd(pad);
	TLegend* leg = new TLegend(0.7, 0.4, 0.89, 0.6);
	if( (vh.size() != vc.size()) && vh.size() == vname.size() ) {
		cout << "ERROR IN ADD_TO_CANVAS!\nSIZE MISMATCH\n";
		return;
	}
	for(int i = 0; i < vh.size(); i++) {
		vh[i]->SetLineColor(vc[i]);
		if(i == 0){
			vh[i]->Draw();
			vh[i]->SetTitle(name.c_str());
		}
		else
			vh[i]->Draw("SAME");
		leg->AddEntry(vh[i], vname[i].c_str(), "l");
	}
	leg->Draw();

}

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
void energy_calibration(string sfNa, string sNaDur, string sfCs, string sCsDur, string sbkg, string sbkgDur, string sfout)
{
	// Open Save Folder
	TFile* fsave = new TFile(sfout.c_str(), "RECREATE");

	// Parse file names
	vector<string> vNa = parse_str(sfNa, '\n');
	vector<string> vCs = parse_str(sfCs, '\n');
	vector<string> vbkg = parse_str(sbkg, '\n');

	// Parse Durations
	vector<string> vNaDur = parse_str(sNaDur); vector<int> vNaDurI; for(auto it = vNaDur.begin(); it!=vNaDur.end(); it++) { vNaDurI.push_back(stoi(*it)); }
	vector<string> vCsDur = parse_str(sCsDur); vector<int> vCsDurI; for(auto it = vCsDur.begin(); it!=vCsDur.end(); it++) { vCsDurI.push_back(stoi(*it)); }
	vector<string> vbkgDur = parse_str(sbkgDur); vector<int> vbkgDurI; for(auto it = vbkgDur.begin(); it!=vbkgDur.end(); it++) { vbkgDurI.push_back(stoi(*it)); }



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

	// Open Up background files
	vector<TFile*> vfbkg(vbkg.size());
	for(int file = 0; file < vbkg.size(); file++ ){
		vfbkg[file] = new TFile(vbkg[file].c_str(), "READ");
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

	// If we have several bkg Files open
	// Add gam1 and gam2 spectra together
	vector<TH1F*> hbkg_gam1;
	vector<TH1F*> hbkg_gam2;
	// Loop over files
	for(auto it = vfbkg.begin(); it != vfbkg.end(); it++)
	{
		// Get Histograms of interest
		vector<TH1F*> vh = getHists(*it, vhnames);
		hbkg_gam1.push_back( vh[0] );
		hbkg_gam2.push_back( vh[1] );
	}

	// Get Added Hists
	TH1F* hNaGam1 = add_histograms(hNa_gam1);
	TH1F* hNaGam2 = add_histograms(hNa_gam2);
	TH1F* hCsGam1 = add_histograms(hCs_gam1);
	TH1F* hCsGam2 = add_histograms(hCs_gam2);
	TH1F* hbkgGam1 = add_histograms(hbkg_gam1);
	TH1F* hbkgGam2 = add_histograms(hbkg_gam2);

	// Get Background subbed histograms
	TH1F* hNa_Gam1_BKG = background_sub(hNaGam1, hbkgGam1, 10.0, 5.0);
	TH1F* hNa_Gam2_BKG = background_sub(hNaGam2, hbkgGam2, 10.0, 5.0);
	TH1F* hCs_Gam1_BKG = background_sub(hCsGam1, hbkgGam1, 10.0, 5.0);
	TH1F* hCs_Gam2_BKG = background_sub(hCsGam2, hbkgGam2, 10.0, 5.0);

	
	// Save Hists
	fsave->cd();
	hNaGam1->Write("Na_Gam1_Added");
	hNaGam2->Write("Na_Gam2_Added");
	hCsGam1->Write("Cs_Gam1_Added");
	hCsGam2->Write("Cs_Gam2_Added");
	hbkgGam1->Write("bkg_Gam1_Added");
	hbkgGam2->Write("bkg_Gam2_Added");
	hNa_Gam1_BKG->Write("Na_Gam1_BKG");
	hNa_Gam2_BKG->Write("Na_Gam2_BKG");
	hCs_Gam1_BKG->Write("Cs_Gam1_BKG");
	hCs_Gam2_BKG->Write("Cs_Gam2_BKG");
	
	// Populate canvas
	TCanvas* cbkg = new TCanvas();cbkg->Divide(2,2);
	add_to_canvas(cbkg, 1, {hNaGam1, hbkgGam1, hNa_Gam1_BKG}, {4, 1, 2}, {"Na", "BKG", "No BKG"}, "NA_GAM1");
	add_to_canvas(cbkg, 2, {hNaGam2, hbkgGam2, hNa_Gam2_BKG}, {4, 1, 2}, {"Na", "BKG", "No BKG"}, "NA_GAM2");
	add_to_canvas(cbkg, 3, {hCsGam1, hbkgGam1, hCs_Gam1_BKG}, {4, 1, 2}, {"Cs", "BKG", "No BKG"}, "CS_GAM1");
	add_to_canvas(cbkg, 4, {hCsGam2, hbkgGam2, hCs_Gam2_BKG}, {4, 1, 2}, {"Cs", "BKG", "No BKG"}, "CS_GAM2");
	fsave->cd();
	cbkg->Write("BKG_SUBTRACTION");


	// Lets start Fitting
	TTree* tdata = new TTree("Data", "Peak Data Tree");

	TCanvas *cNa_Gam1 = new TCanvas("Na_GAM1_FIT"); cNa_Gam1->Divide(1,2);
	fit_peaks(fsave, tdata, cNa_Gam1, 1, hNa_Gam1_BKG, 120, 170, 0, 1600, "NA_GAM1_PEAK1");
	fit_peaks(fsave, tdata, cNa_Gam1, 2, hNa_Gam1_BKG, 300, 390, 0, 1600, "NA_GAM1_PEAK2");
	TCanvas *cNa_Gam2 = new TCanvas("Na_GAM2_FIT"); cNa_Gam2->Divide(1,2);
	fit_peaks(fsave, tdata, cNa_Gam2, 1, hNa_Gam2_BKG, 115, 180, 0, 1600, "NA_GAM2_PEAK1");
	fit_peaks(fsave, tdata, cNa_Gam2, 2, hNa_Gam2_BKG, 300, 390, 0, 1600, "NA_GAM2_PEAK2");
	TCanvas *cCs_Gam = new TCanvas("Cs_GAM_FIT"); cCs_Gam->Divide(1,2);
	fit_peaks(fsave, tdata, cCs_Gam, 1, hCs_Gam1_BKG, 150, 230, 0, 1600, "CS_GAM1_PEAK1");
	fit_peaks(fsave, tdata, cCs_Gam, 2, hCs_Gam2_BKG, 150, 230, 0, 1600, "CS_GAM2_PEAK1");
	tdata->Write();


	// Construct energy calibration
	double energies[] = {511.0, 1274.537, 661.657 }; // [0] - NA | [1] - NA |  [3] - Cs  => KeV
	double errors[]   = {0.0  ,    0.007,   0.003 };
	info_t energy  = { {energies[0], energies[1], energies[2]} , { errors[0], errors[1], errors[2]} };
	info_t measuredGAM1;
	info_t measuredGAM2;
	// vector<string> vbn = {"NA_GAM1_PEAK1","NA_GAM1_PEAK2","NA_GAM2_PEAK1","NA_GAM2_PEAK2","CS_GAM1_PEAK1","CS_GAM2_PEAK1"};
	vector<string> vbn = {"NA_GAM1_PEAK1","NA_GAM1_PEAK2","CS_GAM1_PEAK1"};
	double cts[vbn.size()];
	double mean[vbn.size()];
	double sigma[vbn.size()];
	double mean_error[vbn.size()];
	int npoint = 0;
	for ( auto it = vbn.begin(); it != vbn.end(); it++ ) {
		TBranch* b = tdata->GetBranch((*it).c_str());
		double dvar[5];
		b->SetAddress(&dvar);
		b->GetEntry();
		cts[npoint]  = dvar[4];
		mean[npoint] = dvar[2];
		sigma[npoint]= dvar[3];
		mean_error[npoint] = sigma[npoint]/TMath::Sqrt(cts[npoint]);
		//measuredGAM1[npoint]= {mean[npoint], mean_error[npoint]};
		measuredGAM1.val[npoint]= mean[npoint];
		measuredGAM1.err[npoint]= mean_error[npoint];
		npoint++;
	}
	TCanvas* cEC1 = new TCanvas("GAM 1 EC"); cEC1->cd();
	TGraphErrors* gCalib1 = new TGraphErrors(3, mean, energies, mean_error, errors);
	gCalib1->SetMarkerStyle(22);
	gCalib1->SetTitle("Energy Calibration: GAM1");
	gCalib1->GetXaxis()->SetTitle("ADC Bin");
	gCalib1->GetYaxis()->SetTitle("Energy [keV]");
	gCalib1->Draw("AP");
	TFitResultPtr pfr1 = gCalib1->Fit("pol1","S");
	vbn = {"NA_GAM2_PEAK1","NA_GAM2_PEAK2","CS_GAM2_PEAK1"};
	npoint = 0;
	for ( auto it = vbn.begin(); it != vbn.end(); it++ ) {
		TBranch* b = tdata->GetBranch((*it).c_str());
		double dvar[5];
		b->SetAddress(&dvar);
		b->GetEntry();
		cts[npoint]  = dvar[4];
		mean[npoint] = dvar[2];
		sigma[npoint]= dvar[3];
		mean_error[npoint] = sigma[npoint]/TMath::Sqrt(cts[npoint]);
		// measuredGAM2[npoint]= {mean[npoint], mean_error[npoint]};
		measuredGAM2.val[npoint]= mean[npoint];
		measuredGAM2.err[npoint]= mean_error[npoint];
		npoint++;
	}
	TCanvas* cEC2 = new TCanvas("GAM 2 EC"); cEC2->cd();
	TGraphErrors* gCalib2 = new TGraphErrors(3, mean, energies, mean_error, errors);
	gCalib2->SetMarkerStyle(23);
	gCalib2->SetMarkerColor(kBlue);
	gCalib2->SetTitle("Energy Calibration: GAM2");
	gCalib2->GetXaxis()->SetTitle("ADC Bin");
	gCalib2->GetYaxis()->SetTitle("Energy [keV]");
	gCalib2->Draw("AP");
	TF1* fec2 = new TF1("fec2", "pol1(0)");
	fec2->SetLineColor(kBlue);
	TFitResultPtr pfr2 = gCalib2->Fit("pol1","S");
	double pec2[] = {pfr2->Parameter(0),pfr2->Parameter(1)};
	fec2->SetParameters(pec2);
	pfr2 = gCalib2->Fit(fec2,"S");

	TCanvas* cCalib = new TCanvas("Energy Calibration"); cCalib->cd();
	TMultiGraph* mgEC = new TMultiGraph("mgEC", "Gam1 and Gam2 Energy Calibration");
	mgEC->Add(gCalib1);
	mgEC->Add(gCalib2);
	mgEC->GetXaxis()->SetTitle("ADC Bin");
	mgEC->GetYaxis()->SetTitle("Energy [keV]");
	mgEC->Draw("AP");
	TLegend* legEC1 = new TLegend(0.7, 0.4, 0.89, 0.6);
	legEC1->AddEntry(gCalib1, "Gam1", "p");
	legEC1->AddEntry(gCalib2, "Gam2", "p");
	legEC1->Draw();

	pfr1->Write("GAM1_EC_FIT");
	pfr2->Write("GAM2_EC_FIT");
	gCalib1->Write("GAM1_EC");
	gCalib2->Write("GAM2_EC");
	mgEC->Write("GAM1_GAM2_EC");


	// Calulate the residual for Gam1
	info_t gam1_res;
	info_t gam2_res;
	c_residuals(&gam1_res, &energy, &measuredGAM1, pfr1,3);
	TGraphErrors* gRes_GAM1 = new TGraphErrors(3, energy.val,gam1_res.val,energy.err, gam1_res.err);
	c_residuals(&gam2_res, &energy, &measuredGAM2, pfr2,3);
	TGraphErrors* gRes_GAM2 = new TGraphErrors(3, energy.val,gam2_res.val,energy.err, gam2_res.err);

	TCanvas* cRes = new TCanvas("EC Residuals"); cRes->cd();
	gRes_GAM1->SetMarkerStyle(22);
	gRes_GAM2->SetMarkerStyle(23);
	gRes_GAM2->SetMarkerColor(kBlue);
	TMultiGraph* mgRes = new TMultiGraph("mgRes", "Energy Calibration Residuals");
	TLine* lRes = new TLine(511, 0, 1300,0);
	mgRes->Add(gRes_GAM1);
	mgRes->Add(gRes_GAM2);
	mgRes->GetYaxis()->SetTitle("Residual [keV]");
	mgRes->GetXaxis()->SetTitle("Energy [keV]");
	mgRes->Draw("AP");
	lRes->Draw();
	TLegend* legRes= new TLegend(0.5, 0.6, 0.69, 0.8);
	legRes->AddEntry(gRes_GAM1, "Gam1", "p");
	legRes->AddEntry(gRes_GAM2, "Gam2", "p");
	legRes->Draw();

	cRes->Write("EC_Residual");

}
