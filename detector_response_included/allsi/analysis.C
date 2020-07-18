void prettyTH1( TH1F * h1 , int color , bool fill );
void prettyTH2( TH2F * h2 , TString xtit , TString ytit , TString tit , float labelsize );
void prettyTH2( TH2D * h2 , TString xtit , TString ytit , TString tit , float labelsize );
// ======================================================================================================================
void analysis(){
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gStyle -> SetOptStat(0);
	TString jet_rad = "10";
	// -----------------------------------------------------------------------------------------
	// Loading the data from the root file
	TFile * F = new TFile("$HOME/EIC/output/allsi/710.root");
	TTree * T = (TTree*) F -> Get("T");

	int njets, ntruthjets, nAlltruthjets, id[20], nComponent[20], matched_truthID[20], matched_truthNComponent[20], electron_truthPID,
	    all_truthID[20], all_truthNComponent[20],matched_charged_truthNComponent[20];
	float eta[20],phi[20],e[20],pt[20],matched_truthEta[20],matched_truthPhi[20],matched_truthE[20],matched_truthPt[20],
	matched_charged_truthEta[20],matched_charged_truthPhi[20],matched_charged_truthE[20],matched_charged_truthPt[20],
	electron_truthEta,electron_truthPhi,electron_truthE,electron_truthPt,electron_truthpX,electron_truthpY,electron_truthpZ,
	all_truthEta[20],all_truthPhi[20],all_truthE[20],all_truthPt[20],
	matched_Constituent_truthPt[20][25], matched_Constituent_truthPID[20][25], matched_Constituent_truthCharge[20][25], matched_Constituent_truthEta[20][25], matched_Constituent_truthPhi[20][25], matched_Constituent_truthE[20][25];

	T -> SetBranchAddress("njets"                          , &njets                                         );
	T -> SetBranchAddress("ntruthjets"                     , &ntruthjets                                    );
	T -> SetBranchAddress("nAlltruthjets"                  , &nAlltruthjets                                 );
	T -> SetBranchAddress("id"                             , id/*[njets]*/                                  ); // Reconstructed jets
	T -> SetBranchAddress("nComponent"                     , nComponent/*[njets]*/                          );
	T -> SetBranchAddress("eta"                            , eta/*[njets]*/                                 );
	T -> SetBranchAddress("phi"                            , phi/*[njets]*/                                 );
	T -> SetBranchAddress("e"                              , e/*[njets]*/                                   );
	T -> SetBranchAddress("pt"                             , pt/*[njets]*/                                  );
	T -> SetBranchAddress("matched_truthID"                , matched_truthID/*[ntruthjets]*/                ); // Matched Truth jets
	T -> SetBranchAddress("matched_truthNComponent"        , matched_truthNComponent/*[ntruthjets]*/        );
	T -> SetBranchAddress("matched_truthEta"               , matched_truthEta/*[ntruthjets]*/               );
	T -> SetBranchAddress("matched_truthPhi"               , matched_truthPhi/*[ntruthjets]*/               );
	T -> SetBranchAddress("matched_truthE"                 , matched_truthE/*[ntruthjets]*/                 );
	T -> SetBranchAddress("matched_truthPt"                , matched_truthPt/*[ntruthjets]*/                );
	T -> SetBranchAddress("matched_charged_truthNComponent", matched_charged_truthNComponent/*[ntruthjets]*/); // Matched Truth charged jets
	T -> SetBranchAddress("matched_charged_truthEta"       , matched_charged_truthEta/*[ntruthjets]*/       );
	T -> SetBranchAddress("matched_charged_truthPhi"       , matched_charged_truthPhi/*[ntruthjets]*/       );
	T -> SetBranchAddress("matched_charged_truthE"         , matched_charged_truthE /*[ntruthjets]*/        );
	T -> SetBranchAddress("matched_charged_truthPt"        , matched_charged_truthPt/*[ntruthjets]*/        );
	T -> SetBranchAddress("all_truthID"                    , all_truthID/*[nAlltruthjets]*/                 ); // All Truth jets
	T -> SetBranchAddress("all_truthNComponent"            , all_truthNComponent/*[nAlltruthjets]*/         );
	T -> SetBranchAddress("all_truthEta"                   , all_truthEta/*[nAlltruthjets]*/                );
	T -> SetBranchAddress("all_truthPhi"                   , all_truthPhi/*[nAlltruthjets]*/                );
	T -> SetBranchAddress("all_truthE"                     , all_truthE/*[nAlltruthjets]*/                  );
	T -> SetBranchAddress("all_truthPt"                    , all_truthPt/*[nAlltruthjets]*/                 );
	T -> SetBranchAddress("matched_Constituent_truthPID"   , matched_Constituent_truthPID	                );
	T -> SetBranchAddress("matched_Constituent_truthCharge", matched_Constituent_truthCharge                );
	T -> SetBranchAddress("matched_Constituent_truthEta"   , matched_Constituent_truthEta	                );
	T -> SetBranchAddress("matched_Constituent_truthPhi"   , matched_Constituent_truthPhi       	        );
	T -> SetBranchAddress("matched_Constituent_truthPt"    , matched_Constituent_truthPt	       	        );
	T -> SetBranchAddress("matched_Constituent_truthE"     , matched_Constituent_truthE	                );
	T -> SetBranchAddress("electron_truthEta"              , &electron_truthEta                             ); // Electron variables
	T -> SetBranchAddress("electron_truthPhi"              , &electron_truthPhi                             );
	T -> SetBranchAddress("electron_truthE"                , &electron_truthE                               );
	T -> SetBranchAddress("electron_truthPt"               , &electron_truthPt                              );
	T -> SetBranchAddress("electron_truthpX"               , &electron_truthpX                              );
	T -> SetBranchAddress("electron_truthpY"               , &electron_truthpY                              );
	T -> SetBranchAddress("electron_truthpZ"               , &electron_truthpZ                              );
	T -> SetBranchAddress("electron_truthPID"              , &electron_truthPID                             );

	int nEntries = T -> GetEntries();

	// -----------------------------------------------------------------------------------------
	// Declaring histograms
	int color[] = {1,2,62};

	TH1F * h1_jets_pT_Truth_Reco_SC_mix = new TH1F("h1_jets_pT_Truth_Reco_SC_mix",";p_{T,Reco} / p_{T,Truth};Counts"    ,100,0,1.6);	// several constituents, mixed
	TH1F * h1_parts_pT_Truth = new TH1F("h1_parts_pT_Truth", ";p_{T, truth}^{constituent} [GeV];Counts", 100, 0, 4);
	TH1F * h1_parts_pT_Truth_low = new TH1F("h1_parts_pT_Truth_low", ";p_{T, truth}^{constituent} [GeV];Counts", 100, 0, 4);
	prettyTH1( h1_jets_pT_Truth_Reco_SC_mix , 2 , true );
	prettyTH1( h1_parts_pT_Truth, 1, true);
	prettyTH1( h1_parts_pT_Truth_low, 2, true);


	int count = 0;
	// -----------------------------------------------------------------------------------------
	// Loop over entries
	for(int ent = 0 ; ent < 110000; ent++){
		T -> GetEntry(ent);		
		if( ent%10000 == 0 ) cout << "Event " << ent << "\tof " << nEntries << endl;
		if( ntruthjets != njets ){ cout << "Something went wrong. Skipping event." << endl; continue;}
		for(int jet = 0 ; jet < ntruthjets ; jet++){
			TLorentzVector truth_jet, truth_charged_jet, reco_jet;
			truth_jet        .SetPtEtaPhiE(matched_truthPt        [jet],matched_truthEta        [jet],matched_truthPhi        [jet],matched_truthE        [jet]);
			truth_charged_jet.SetPtEtaPhiE(matched_charged_truthPt[jet],matched_charged_truthEta[jet],matched_charged_truthPhi[jet],matched_charged_truthE[jet]);
			reco_jet.SetPtEtaPhiE(pt                     [jet],eta                     [jet],phi                     [jet],e                     [jet]);

			if( matched_truthNComponent[jet] > 1 && matched_truthNComponent[jet] != matched_charged_truthNComponent[jet] && matched_charged_truthNComponent[jet] != 0 ){
				h1_jets_pT_Truth_Reco_SC_mix -> Fill(pt[jet]/matched_truthPt[jet]);

				if ( matched_truthNComponent[jet] < 26 ){
					for (int part = 0; part < matched_truthNComponent[jet]; part ++){
						h1_parts_pT_Truth -> Fill(matched_Constituent_truthE[jet][part]);
						if (pt[jet]/matched_truthPt[jet] < 0.6){
							h1_parts_pT_Truth_low -> Fill(matched_Constituent_truthE[jet][part]);
						}	
					}
				}
				else if ( matched_truthNComponent[jet] > 10000 ){
					count += 1;
				}
				
			}
		}
	}
	cout << count << endl;
	// ------------------------------------------------------
	TCanvas * c8 = new TCanvas("c8","c8",1300,900);
	gPad -> SetBottomMargin(0.15);
	gPad -> SetLeftMargin(0.15);
	h1_jets_pT_Truth_Reco_SC_mix -> Draw();
	c8 -> Print("images/results_R"+jet_rad+".pdf(");
	
	TLegend * leg1 = new TLegend(0.5, 0.8, 0.9, 0.86);
	leg1 -> AddEntry(h1_parts_pT_Truth, "all const from >1 const mixed");
	leg1 -> AddEntry(h1_parts_pT_Truth_low, "const from jets with >1 const mixed and p_{T,rec}^{jet}/p_{T,tru}^{jet} < 0.6");
	TCanvas * c1 = new TCanvas("c1", "c1", 1300, 900);
	gPad -> SetBottomMargin(0.15);
	gPad -> SetLeftMargin(0.15);
	h1_parts_pT_Truth -> Draw("hist");
	h1_parts_pT_Truth_low -> Draw("samehist");
	leg1 -> Draw("same");
	c1 -> Print("images/results_R"+jet_rad+".pdf)");
}
// ======================================================================================================================
void prettyTH1( TH1F * h1 , int color , bool fill ){
	h1 -> SetLineColor(color);
	h1 -> SetMarkerColor(color);
	if(fill) h1 -> SetFillColorAlpha(color,0.15);

	h1 -> SetLineWidth(2);

	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetTitleSize(0.05);
	h1 -> GetXaxis() -> SetLabelSize(0.05);

	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetYaxis() -> SetTitleSize(0.05);
	h1 -> GetYaxis() -> SetLabelSize(0.05);
}
// ======================================================================================================================
void prettyTH2( TH2F * h2 , TString xtit , TString ytit , TString tit , float labelsize ){
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetXaxis() -> SetTitleSize(labelsize);
	h2 -> GetXaxis() -> SetLabelSize(labelsize);
	//h2 -> GetXaxis() -> SetNdivisions(107);
	if(xtit!="") h2 -> GetXaxis() -> SetTitle(xtit);

	h2 -> GetYaxis() -> CenterTitle();
	h2 -> GetYaxis() -> SetTitleSize(labelsize);
	h2 -> GetYaxis() -> SetLabelSize(labelsize);
	//h2 -> GetYaxis() -> SetNdivisions(107);
	if(ytit!="") h2 -> GetYaxis() -> SetTitle(ytit);

	h2 -> SetTitle(tit);
}
// ======================================================================================================================
void prettyTH2( TH2D * h2 , TString xtit , TString ytit , TString tit , float labelsize ){
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetXaxis() -> SetTitleSize(labelsize);
	h2 -> GetXaxis() -> SetLabelSize(labelsize);
	//h2 -> GetXaxis() -> SetNdivisions(107);
	if(xtit!="") h2 -> GetXaxis() -> SetTitle(xtit);

	h2 -> GetYaxis() -> CenterTitle();
	h2 -> GetYaxis() -> SetTitleSize(labelsize);
	h2 -> GetYaxis() -> SetLabelSize(labelsize);
	//h2 -> GetYaxis() -> SetNdivisions(107);
	if(ytit!="") h2 -> GetYaxis() -> SetTitle(ytit);

	h2 -> SetTitle(tit);
}
