void Lambda_Decay(){
    TFile *f = new TFile("Lc_thCy0_T2cm_Ge5mrad7cm.root");
    
    TTree *tree = (TTree*)f->Get("MCDecayTree");
    
    Double_t p_lambda;
    tree->SetBranchAddress("Lambda_cplus_TRUEP_E", &p_lambda);

    TH1F *lambda_p   = new TH1F("lambda","p distribution of Lambda",100,400000,7000000);
    TH1F *p= new TH1F("p", "proton distribution",100, 0, 4000);
    TH1F *K= new TH1F("K", "K",100, 0, 4000);
    TH1F *pi= new TH1F("pi", "pi",100, 0, 3000);  
    Int_t nentries = (Int_t)tree->GetEntries();
    cout << nentries << endl;

    TLorentzVector lambda(0.,0.,0,2.286);
    Double_t masse[3]= {0.938, 0.494, 0.140};
    TGenPhaseSpace event;
    Double_t weight;

    for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
      lambda_p->Fill(p_lambda);
      
      lambda.SetXYZM(0.,0.,p_lambda*1E-3,2.286);   // esprimo tutto in GeV
      event.SetDecay(lambda, 3, masse);
      weight = event.Generate();
        TLorentzVector *pProton = event.GetDecay(0);
        TLorentzVector *pK    = event.GetDecay(1);
        TLorentzVector *pPi    = event.GetDecay(2);

        p->Fill(pProton->E(), weight);
        K->Fill(pK->E(), weight);
        pi->Fill(pPi->E(), weight);
   }
   lambda_p->Draw();
   TCanvas *c2 = new TCanvas("c2","c2");
   c2->cd();
   p->Draw();

   TCanvas *c3 = new TCanvas("c3","c3");
   c3->cd();
   K->SetLineColor(kRed);
   K->Draw("same");

   TCanvas *c4 = new TCanvas("c4","c4");
   c4->cd();
   pi->SetLineColor(kCyan);
   pi->Draw("same");
}