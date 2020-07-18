void prettyTGraphErrors(TGraphErrors * gP,int color,int marker);
void prettyTH2( TH2D * h2 );
void prettyTH1( TH1D * h1 );
// ==========================================================================================================================
void jet_e_res(){
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	bool write_table = true;
	TString extra_label = "_smallb";//"_1C";
	TString jet_R = "10";
	// -----------------------------------------------------------------------------------
	float eta_bin[] = {0.0,0.5,1.0,2.0,3.0,4.0};	const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	int color[] = {1,2,62,8,92};
	//float eta_bin[] = {0.0,0.5,1.0,2.0,3.0,3.2,3.4,3.6,3.8,4.0};	const int size_eta_bin = sizeof(eta_bin)/sizeof(*eta_bin);
	//int color[] = {1,2,62,8,92,50,55,60,65,70,75,80,85};
	// -----------------------------------------------------------------------------------
	ifstream in_tab;
	// ---
	bool loaded_table_ch = false;
	double approx_mean_ch[size_eta_bin-1][30] = {{0}};
	double approx_sigma_ch[size_eta_bin-1][30] = {{0}};
	in_tab.open("tables/sigma_ch_R"+jet_R+Form("_eta_%i",size_eta_bin-1)+extra_label+".txt");
	if(in_tab){
		loaded_table_ch = true;
		for(int p = 4 ; p < 29 ; p++){
			for(int et = 0 ; et < size_eta_bin-1 ; et++){
				in_tab >> approx_mean_ch[et][p];
				in_tab >> approx_sigma_ch[et][p];
				//if(approx_sigma_ch[et][p]>0.04) approx_sigma_ch[et][p]=0.04;
				approx_sigma_ch[et][p]*=1.1;
			}
		}
		in_tab.close();
	}
	// ---
	bool loaded_table = false;
	double approx_mean[size_eta_bin-1][30] = {{0}};
	double approx_sigma[size_eta_bin-1][30] = {{0}};
        in_tab.open("tables/sigma_R"+jet_R+Form("_eta_%i",size_eta_bin-1)+extra_label+".txt");
        if(in_tab){
                loaded_table = true;
                for(int p = 4 ; p < 29 ; p++){
                        for(int et = 0 ; et < size_eta_bin-1 ; et++){
                                in_tab >> approx_mean[et][p];
				in_tab >> approx_sigma[et][p];
                                //if(approx_sigma[et][p]>0.04) approx_sigma[et][p]=0.04;
                        	approx_sigma[et][p]*=1.1;
			}
                }
                in_tab.close();
        }
	// -----------------------------------------------------------------------------------	
	TFile * F = new TFile("output_R"+jet_R+Form("_eta_%i",size_eta_bin-1)+extra_label+".root");
	TH2D ** h2_chjet_dpp_p_eta = new TH2D*[size_eta_bin-1]; // although it's called dpp, we'll read and plot dee
	TH2D ** h2_jet_dpp_p_eta   = new TH2D*[size_eta_bin-1];
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		h2_chjet_dpp_p_eta[et] = (TH2D*) F -> Get(Form("h1_chjet_er_et_vs_et_LB_%i",et));
		h2_jet_dpp_p_eta  [et] = (TH2D*) F -> Get(Form("h1_jet_er_et_vs_et_LB_%i",et));
	}
	// -----------------------------------------------------------------------------------
	TH1D *** h1_dpp_ch = new TH1D**[size_eta_bin-1];
	TH1D *** h1_dpp    = new TH1D**[size_eta_bin-1];
	TF1 *** f_dpp_ch   = new TF1 **[size_eta_bin-1];
	TF1 *** f_dpp      = new TF1 **[size_eta_bin-1];	
	TCanvas ** c_eta_ch = new TCanvas*[size_eta_bin-1];
	TCanvas ** c_eta    = new TCanvas*[size_eta_bin-1];

	double sigma_ch[size_eta_bin-1][30] = {{0}};	double sigma[size_eta_bin-1][30] = {{0}};
	double sig_E_ch[size_eta_bin-1][30] = {{0}};	double sig_E[size_eta_bin-1][30] = {{0}};
	double mom_ch  [size_eta_bin-1][30] = {{0}};    double mom  [size_eta_bin-1][30] = {{0}};
	int ctr_ch     [size_eta_bin-1]     =  {0} ;    int ctr     [size_eta_bin-1]     =  {0} ;

	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		// Charged Jets
		c_eta_ch[et] = new TCanvas(Form("c_eta_ch_%i",et),Form("CHARGED %.1f < |eta| < %.1f",eta_bin[et],eta_bin[et+1]),1200,900);
		c_eta_ch[et] -> Divide(5,5);
		int nBinsX_ch = h2_chjet_dpp_p_eta[et] -> GetXaxis() -> GetNbins();
		h1_dpp_ch[et] = new TH1D*[nBinsX_ch];
		f_dpp_ch [et] = new TF1*[nBinsX_ch];
		
		for(int p = 4 ; p < 29 ; p++){	
			h1_dpp_ch[et][p] = (TH1D*) h2_chjet_dpp_p_eta[et]->ProjectionY(Form("p_dpp_ch_%i_%i",et,p), p+1 , p+2 );
			h1_dpp_ch[et][p] -> SetTitle(Form("CH, %.1f < |#eta| < %.1f, %i < E < %i GeV",eta_bin[et],eta_bin[et+1],p,p+1));
			prettyTH1( h1_dpp_ch[et][p] );
			c_eta_ch[et] -> cd(p-3); gPad -> SetBottomMargin(0.23); gPad -> SetRightMargin(0.015); gPad -> SetLeftMargin(0.15);
			h1_dpp_ch[et][p] -> Draw("same");
			if(loaded_table_ch){
				double low = approx_mean_ch[et][p]-0.7*approx_sigma_ch[et][p];
				double high = approx_mean_ch[et][p]+1.3*approx_sigma_ch[et][p];
				if (low < -0.05){
					low = -0.05;
				}
				if (high > 0.05){
					high = 0.05;
				}
				f_dpp_ch[et][p] = new TF1(Form("f_dpp_ch_%i_%i",et,p),"gaus", low, high);
			}
			else
				f_dpp_ch[et][p] = new TF1(Form("f_dpp_ch_%i_%i",et,p),"gaus", -0.03, 0.03);

			h1_dpp_ch[et][p] -> Fit(Form("f_dpp_ch_%i_%i",et,p),"R");

			if(h1_dpp_ch[et][p]->GetMaximum()>40. && f_dpp_ch[et][p] -> GetParameter(2) != 0){
				sigma_ch[et][ctr_ch[et]] = (f_dpp_ch[et][p] -> GetParameter(2))*100.;
				sig_E_ch[et][ctr_ch[et]] = (f_dpp_ch[et][p] -> GetParError(2))*100.;
				mom_ch  [et][ctr_ch[et]] = (double)(p+p+1)/2.;
				ctr_ch  [et]++;
			}
			else {
				h1_dpp_ch[et][p] -> SetLineColor(92);
			}
		}
		// -----------------------------------------------------
		// Whole Jets
		c_eta[et] = new TCanvas(Form("c_eta_%i",et),Form("WHOLE, %.1f < |eta| < %.1f",eta_bin[et],eta_bin[et+1]),1200,900);
                c_eta[et] -> Divide(5,5);
		int nBinsX = h2_jet_dpp_p_eta[et] -> GetXaxis() -> GetNbins();
		h1_dpp[et] = new TH1D*[nBinsX];
                f_dpp [et] = new TF1*[nBinsX];
		
		for(int p = 4 ; p < 29 ; p++){
			h1_dpp[et][p] = (TH1D*) h2_jet_dpp_p_eta[et]->ProjectionY(Form("p_dpp_%i_%i",et,p), p+1 , p+2 );
                        h1_dpp[et][p] -> SetTitle(Form("%.1f < |#eta| < %.1f, %i < E < %i GeV",eta_bin[et],eta_bin[et+1],p,p+1));
                        prettyTH1( h1_dpp[et][p] );
			c_eta[et] -> cd(p-3); gPad -> SetBottomMargin(0.23); gPad -> SetRightMargin(0.015); gPad -> SetLeftMargin(0.15);
                        h1_dpp[et][p] -> Draw("same");

                        if(loaded_table){
				double low = approx_mean[et][p]-0.6*approx_sigma[et][p];
				double high = approx_mean[et][p]+1.3*approx_sigma[et][p];
				if (low < -0.1){
					low = -0.1;
				}
				if (high > 0.1){
					high = 0.1;
				}
                                f_dpp[et][p] = new TF1(Form("f_dpp_%i_%i",et,p),"gaus", low, high);
                	}
			else{
                                f_dpp[et][p] = new TF1(Form("f_dpp_%i_%i",et,p),"gaus",-0.03,0.03);
			}
                        h1_dpp[et][p] -> Fit(Form("f_dpp_%i_%i",et,p),"R");
	
                        if(h1_dpp[et][p]->GetMaximum()>60. && h1_dpp[et][p]->GetMaximum() > 1.7 * (h1_dpp[et][p]->GetEntries() / h1_dpp[et][p]->GetNbinsX())){
				sigma[et][ctr[et]] = (f_dpp[et][p] -> GetParameter(2))*100.;
				sig_E[et][ctr[et]] = (f_dpp[et][p] -> GetParError(2))*100.;
                                mom  [et][ctr[et]] = (double)(p+p+1)/2.;
                                ctr  [et]++;
			}
			else{
				h1_dpp[et][p] -> SetLineColor(92);
			}		
		}
	}
	// -----------------------------------------------------------------------------------
	if(write_table){
		ofstream out_tab;
		// ---
		out_tab.open("tables/sigma_ch_R"+jet_R+Form("_eta_%i",size_eta_bin-1)+extra_label+".txt");
		for(int p = 4 ; p < 29; p++){
			for(int et = 0 ; et < size_eta_bin-1 ; et++){
				out_tab << f_dpp_ch[et][p] -> GetParameter(1);
				out_tab << f_dpp_ch[et][p] -> GetParameter(2);
				if(et == size_eta_bin-2) out_tab << endl;
				else out_tab << "\t";
			}
		}
		out_tab.close();	
		// ---
		out_tab.open("tables/sigma_R"+jet_R+Form("_eta_%i",size_eta_bin-1)+extra_label+".txt");
                for(int p = 4 ; p < 29 ; p++){
                        for(int et = 0 ; et < size_eta_bin-1 ; et++){
                                out_tab /*<< " parameter 1: "*/ << f_dpp[et][p] -> GetParameter(1);
				out_tab /*<< " parameter 2: "*/ << f_dpp[et][p] -> GetParameter(2);
                                if(et == size_eta_bin-2) out_tab << endl;
                                else out_tab << "\t";
                        }
                }
                out_tab.close();
	}
	// -----------------------------------------------------------------------------------
	double xDum[] = {4,30};	double yDum[] = {0,6};
	TGraphErrors * gDummy = new TGraphErrors(2,xDum,yDum);
	prettyTGraphErrors(gDummy, 1 , 1);
	// -----------------------------------------------------------------------------------
	TGraphErrors ** g_dpp_p_ch = new TGraphErrors*[size_eta_bin-1];
	TGraphErrors ** g_dpp_p    = new TGraphErrors*[size_eta_bin-1];
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		g_dpp_p_ch[et] = new TGraphErrors(ctr_ch[et],mom_ch[et],sigma_ch[et],0,sig_E_ch[et]);
		g_dpp_p   [et] = new TGraphErrors(ctr   [et],mom   [et],sigma   [et],0,sig_E   [et]);
		prettyTGraphErrors( g_dpp_p_ch[et] , color[et] , 20);
		prettyTGraphErrors( g_dpp_p   [et] , color[et] , 21);	g_dpp_p[et] -> SetLineStyle(7);
	}
	// -----------------------------------------------------------------------------------
	TLegend * leg1 = new TLegend(0.6,0.35,0.8,0.6);
	leg1 -> SetLineColor(0);
	for(int et = 0 ; et < size_eta_bin-1 ; et++){leg1 -> AddEntry(g_dpp_p[et],Form("%.1f < |#eta| < %.1f",eta_bin[et],eta_bin[et+1]));}
	// -----------------------------------------------------------------------------------
	TCanvas * c1 = new TCanvas("c1");
	gDummy -> Draw("AP");
	gDummy -> SetTitle("Charged jets");
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		g_dpp_p_ch[et] -> Draw("samePLE3");
	}
	leg1 -> Draw("same");
	c1 -> Print("e_res"+extra_label+"_1.pdf(");
	// -----------------------------------------------------------------------------------
	TCanvas * c2 = new TCanvas("c2");
	gDummy -> Draw("AP");
	gDummy -> SetTitle("Whole jets");
	for(int et = 0 ; et < size_eta_bin-1 ; et++){
		g_dpp_p[et] -> Draw("samePLE3");
	}
	leg1 -> Draw("same");
	// -----------------------------------------------------------------------------------
	for(int et = 0 ; et < size_eta_bin-1 ; et++) c_eta_ch[et] -> Print("e_res"+extra_label+"_1.pdf");
	for(int et = 0 ; et < size_eta_bin-1 ; et++) c_eta   [et] -> Print("e_res"+extra_label+"_1.pdf");
	c2 -> Print("e_res"+extra_label+"_1.pdf)");
}
// ==========================================================================================================================
void prettyTGraphErrors(TGraphErrors * gP,int color,int marker){
	gP -> SetMarkerStyle(marker);
	gP -> SetMarkerColor(color);
	gP -> SetLineColor(color);
	gP -> SetLineWidth(2);
	gP -> SetFillColorAlpha(color,.15);

	gP -> GetXaxis() -> SetTitle("E [GeV]");
	gP -> GetXaxis() -> CenterTitle();
	gP -> GetXaxis() -> SetNdivisions(107);
	gP -> GetYaxis() -> SetTitle("#sigma(dE/E) [%]");
	gP -> GetYaxis() -> CenterTitle();
	gP -> GetYaxis() -> SetNdivisions(107);

	gP -> SetTitle("");
}
// ==========================================================================================================================
void prettyTH2( TH2D * h2 ){
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetXaxis() -> SetNdivisions(107);
	//h2 -> GetXaxis() -> SetTitleSize(0.06);
	//h2 -> GetXaxis() -> SetLabelSize(0.06);
	h2 -> GetYaxis() -> CenterTitle();
	h2 -> GetYaxis() -> SetNdivisions(107);
	//h2 -> GetYaxis() -> SetTitleSize(0.06);
        //h2 -> GetYaxis() -> SetLabelSize(0.06);
}
// ==========================================================================================================================
void prettyTH1( TH1D * h1 ){
        h1 -> GetXaxis() -> CenterTitle();
        h1 -> GetXaxis() -> SetNdivisions(107);
        h1 -> GetXaxis() -> SetTitleSize(0.06);
        h1 -> GetXaxis() -> SetLabelSize(0.06);
	h1 -> GetXaxis() -> SetTitleOffset(1.45);
        h1 -> GetYaxis() -> CenterTitle();
        h1 -> GetYaxis() -> SetNdivisions(107);
        h1 -> GetYaxis() -> SetTitleSize(0.06);
        h1 -> GetYaxis() -> SetLabelSize(0.06);
}
