#include <iostream>
#include <fstream>
#include <vector>
using namespace std;





double sum_hist(int xlow, int xhigh, string hname, TFile* file, double scale = 10.0)
{
	if(file == NULL){
		cout << "ERROR IN SUM_HIST()!\nTFILE CANNOT BE NULL\nEXITING...";
		exit(-1);
	}

	// Get Histogram that wants needs to be summed
	TH1F* h = (TH1F*)file->Get(hname.c_str());
	double sum = h->Integral(xlow, xhigh);
	sum *= 10.0/scale;
	// cout << "sum: " << sum << endl;
	return sum;
}



// parse string with vars seperated by a delim
vector<string> parse_str(string &str, char delim = ',')
{
	if(str.back() != delim)
		str += delim;
		
	std::vector<std::string> sfiles;
	char trash[755]; int ti = 0;
	memset(trash, 0, 255);
	for(int i = 0; i < str.size(); i++) {
		if( str[i] == delim ){
			sfiles.push_back(trash);
			// std::cout << sfiles.back() << std::endl;
			memset(trash, 0, 255);
			ti = 0;
		}
		else
			trash[ti++] = str[i];
	}
	return sfiles;
}





/***************************************************************/
// This will make four TGraphs
// 0) gam1 & g11 at 511
// 1) gam1 & g11 at 1.275
// 2) gam1 & g12 at 511
// 3) gam1 & g12 at 1.275
void angular_study(string sinfile, string sdeg, string scale, string soutfile)
{
	// Might need to parse files
	// parsed files is vs
	vector<string> vs = parse_str(sinfile, '\n');
	
	// Parse the degree list
	vector<string> vsdeg =  parse_str(sdeg);
	vector<int> vdeg;
	for(auto it=vsdeg.begin(); it!=vsdeg.end(); vdeg.push_back( stoi(*it) ), it++);

	// Parse the duration list
	vector<string> vsdur = parse_str(scale);
	vector<int> vdur;
	for(auto it=vsdur.begin(); it!=vsdur.end(); vdur.push_back( stoi(*it) ), it++);

	// Make histograms
	vector<TGraphErrors*> vg(4);
	for(auto it=vg.begin(); it!=vg.end(); *it=new TGraphErrors(), it++);

	// Open file(s)
	vector<TFile*> vfiles;
	int counter  = 0;
	for(auto it = vs.begin(); it != vs.end(); it++) {
		vfiles.push_back(  new TFile(it->c_str(), "READ" )  );
		vector<double> vy(4);
		vy[0] = ( sum_hist(126,168, "g11", vfiles.back(), vdur[counter]) );
		vy[1] = ( sum_hist(283,376, "g11", vfiles.back(), vdur[counter]) );
		vy[2] = ( sum_hist(126,168, "g12", vfiles.back(), vdur[counter]) );
		vy[3] = ( sum_hist(283,376, "g12", vfiles.back(), vdur[counter]) );
	
		vector<double> vey(4);
		for(int index = 0; index < 4; index++)
		{
			if(vy.at(index) == 0)
				vey.at(index) = 1.148;
			else
				vey.at(index) = TMath::Sqrt(vy[index]);
		}

		vg[0]->AddPoint(vdeg[counter], vy[0]);
		vg[1]->AddPoint(vdeg[counter], vy[1]);
		vg[2]->AddPoint(vdeg[counter], vy[2]);
		vg[3]->AddPoint(vdeg[counter], vy[3]);

		// Set Error (Counting statistic)
		vg[0]->SetPointError(counter, 2.5, vey[0] );
		vg[1]->SetPointError(counter, 2.5, vey[1] );
		vg[2]->SetPointError(counter, 2.5, vey[2] );
		vg[3]->SetPointError(counter, 2.5, vey[3] );

		counter++;
	}

	// Plot
	TCanvas* c = new TCanvas();
	c->Divide(2,2);
	c->cd(1); vg[0]->SetTitle("gam1 & g11 (511 keV)"  ); vg[0]->SetMarkerStyle(20); vg[0]->Draw("AP");
	c->cd(2); vg[1]->SetTitle("gam1 & g11 (1.275 keV)"); vg[1]->SetMarkerStyle(20); vg[1]->Draw("AP");
	c->cd(3); vg[2]->SetTitle("gam1 & g12 (511 keV)"  ); vg[2]->SetMarkerStyle(20); vg[2]->Draw("AP");
	c->cd(4); vg[3]->SetTitle("gam1 & g12 (1.275 keV)"); vg[3]->SetMarkerStyle(20); vg[3]->Draw("AP");


}
