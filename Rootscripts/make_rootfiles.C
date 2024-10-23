#include <iostream>
#include <string>
#include <fstream>
#include <vector>


using namespace std;

enum subtag { gam1 = 0, gam2 = 1, gtac = 2, g1_t = 3, g2_t = 4, g1_g2 = 5, g11 = 6, g21 = 7, g12 = 8, g22 = 9 };



/*
 * Func: xy_ints_to_tree
 * Brief:
 *	Opens a file that contains columnar xy data and writes it to a ROOT file
 *	stored as a TTree. It also gets the lowest and highest x values 
 *  (usefule for when making histograms out of the columnar data)
 * Comments:
 *	Assumes we will be reading in int data values
 */
void xy_ints_to_tree(TString& dfile, TTree* tree = NULL, TString bname = "XY_Data", TString xname = "Chan", TString yname = "Counts")
{
	// Check for tree	
	if(tree == NULL)
	{
		cout << "NEED TO HAVE TREE DEFINED BEFORE CALLING xy_ints_to_root(...)\n EXITING...\n";
		exit(-1);
	}
	// Create a TTree
	int vals[2], min=kMaxInt, max=kMinInt;
	tree->Branch(bname, &vals, Form("%s/I:%s/I", xname.Data(), yname.Data()) );
	
	// Open File
	std::ifstream inp; 
	inp.open(dfile, std::ios::in);
	while(inp >>vals[0]  >> vals[1] )
	{
		if(vals[0] > max) max = vals[0];
		if(vals[0] < min) min = vals[0];
		tree->Fill();
	}
	TString maxbname = bname+"MaxBin";
	TString minbname = bname+"MinBin";
	auto maxBranch = tree->Branch(maxbname,&max, "Maximum/I");
	auto minBranch = tree->Branch(minbname,&min, "Minimum/I");
	maxBranch->Fill();
	minBranch->Fill();
}

TString subtag_to_str(int n)
{
	string ret;
	switch(n)
	{
		case subtag::gam1:
			ret = "gam1";
			break;
		case subtag::gam2:
			ret = "gam2";
			break;
		case subtag::gtac:
			ret = "gtac";
			break;
		case subtag::g1_t:
			ret = "g1_t";
			break;
		case subtag::g2_t:
			ret = "g2_t";
			break;
		case subtag::g1_g2:
			ret = "g1_g2";
			break;
		case subtag::g11:
			ret = "g11";
			break;
		case subtag::g21:
			ret = "g21";
			break;
		case subtag::g12:
			ret = "g12";
			break;
		case subtag::g22:
			ret = "g22";
			break;
		default:
			cout << "Unknown subtag!\n Defaulting to 'UNKNOWN'\n";
			ret = "UNKNOWN";
	}
	return ret;
}

/*
 * Func: c_CountsVsChan_histo
 * Brief:
 *	Takes in x,y data stored inside a Tree (as a array[2]) and creates
 *	a simple Chan vs Counts histogram. Saves the histogram in the rootfile
 * Comments:
 *	Assumes we will be reading in int data values
 *	Assumes the TFile is already open
 */
TH1F* c_CountsVsChan_histo(TFile* file = NULL, TTree* T = NULL, TString bname = "BNAME", TString histoName = "HISTO") 
{
	if(T == NULL || file == NULL ){
		cout << "TTree and TFile cannot be NULL when calling c_CountsVsChan_histo!\n Exiting...";
		exit(-1);
	}
	file->cd();
	// Get Min and Max Chan number from the Tree
	int min, max;
	TString maxbname = bname+"MaxBin";
	TString minbname = bname+"MinBin";

	TBranch* bmin = T->GetBranch(minbname);
	TBranch* bmax = T->GetBranch(maxbname);
	bmin->SetAddress(&min);
	bmax->SetAddress(&max);
	bmin->GetEntry(0);
	bmax->GetEntry(0);
	// Make Histogram
	TH1F* h = new TH1F(histoName.Data(), histoName,max-min+1, min-.5, max+.5);
	int xy[2];
	TBranch* bdata = T->GetBranch(bname);
	bdata->SetAddress(&xy);
	for(long long int i = 0; i < T->GetEntries(); i++) {
		bdata->GetEntry(i);
		h->Fill(xy[0], xy[1]);
		h->SetBinError(xy[0], pow(float(xy[1]), 0.5));
	}
	h->GetXaxis()->SetTitle("Channels");
	h->GetYaxis()->SetTitle("Counts");
	std::cout << "Creating Counts vs Channel Histogram: " << histoName << std::endl;
	// h->Write(histoName, TObject::kOverwrite);
	return h;
}


void make_rootfiles(string path, string stag, TString sfout)
{
	TFile* fsave = new TFile(sfout, "RECREATE");
	TTree* tdata = new TTree("Data", "Data Tree");
	vector<TH1F*> vh;


	// make file path name for each subtag in main tag dir
	for( int tag = 0; tag < 10; tag++ ) {
		string sfile;
		int tagcounter = stoi(stag) + tag;
		string strtag  = to_string(tagcounter);
		if(strtag.size() != 4)
			strtag = "0"+strtag;
		if(path.back() == '/') {
			sfile = path+"tag"+strtag+".xy";
		}
		else {
			sfile = path+"/tag"+strtag+".xy";
		}
		TString tsfile(sfile);
		xy_ints_to_tree(tsfile, tdata, subtag_to_str(tag) );
		vh.push_back( c_CountsVsChan_histo(fsave, tdata, subtag_to_str(tag), subtag_to_str(tag) ) );
	}
	fsave->Write();
}
