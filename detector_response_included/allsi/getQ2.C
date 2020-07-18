void getQ2(){
	
	TFile *F = new TFile("$HOME/EIC/output/allsi/710.root");
	TTree *T = (TTree *)F->Get("T");
	
	float electron_truthEta, electron_truthPhi, electron_truthE, electron_truthPt, Q2;
	T->SetBranchAddress("electron_truthEta", &electron_truthEta);
	T->SetBranchAddress("electron_truthPhi", &electron_truthPhi);
	T->SetBranchAddress("electron_truthE", &electron_truthE);
	T->SetBranchAddress("electron_truthPt", &electron_truthPt);

	int nEntries = T->GetEntries();
		
	TH1F *h_Q2 = new TH1F("Q2", "", 200, 16, 500);
	TH1F *h_Q2_large_eta = new TH1F("Q2 for large eta", "", 20, 0, 0.5);
	TH1::SetDefaultSumw2();
	
	vector <float> arr;

	for (int ent = 0; ent < nEntries; ent ++){
	
		T->GetEntry(ent);
		TLorentzVector e_vec;
		e_vec.SetPtEtaPhiE(electron_truthPt, electron_truthEta, electron_truthPhi, electron_truthE);
		Q2 = 2. * 20 * e_vec.E() * (1 - TMath::Abs(e_vec.CosTheta()));
		h_Q2->Fill(Q2);
		if (Q2 < 0.2 || Q2 == 16){
			cout << "event, Q2, eta, e, pt, phi: " << ent << ", " << Q2 << ", " << e_vec.Eta() << ", " << e_vec.E() << ", " << e_vec.Pt() << ", " << e_vec.Phi() << endl;
			cout << TMath::Abs(e_vec.CosTheta()) << endl;
		}
		if (Q2 > 16 && Q2 < 500){
			arr.push_back(Q2);
		}
		if (electron_truthEta > 3.5){
			h_Q2_large_eta->Fill(Q2);
		}
	}

	float mean = h_Q2->GetMean();
	cout << "mean: " << mean << endl;

	int entries = h_Q2->Integral();
	cout << "entries: " << entries << endl;

	float min = *min_element(arr.begin(), arr.end());
	cout << "min: " << min << endl;

	float max = *max_element(arr.begin(), arr.end());
	cout << "max: " << max << endl;

	TCanvas *c1 = new TCanvas();
	h_Q2->Draw();
	c1->Print("Q2.pdf(");
	TCanvas *c2 = new TCanvas();
	h_Q2_large_eta->Draw();
	c2->Print("Q2.pdf)");
}
