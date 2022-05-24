void Lambda_Decay(){
    TFile *f = new TFile("Lc_thCy0_T2cm_Ge5mrad7cm.root");
    
    TTree *tree = (TTree*)f->Get("MCDecayTree");
    
    Double_t E_lambda, px_lambda, py_lambda, pz_lambda;
    tree->SetBranchAddress("Lambda_cplus_TRUEP_E", &E_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_X", &px_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Y", &py_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Z", &pz_lambda);

    TH1F *lambda_p   = new TH1F("#Lambda_{c}^{+}","Momentum distribution of #Lambda_{c}^{+}",100,400,7000);
    TH1F *p= new TH1F("p", "Energy distribution of protons",100, 0, 4000);
    TH1F *K= new TH1F("K^{-}", "Energy distribution of K^{-}",100, 0, 4000);
    TH1F *pi= new TH1F("#pi^{+}", "Energy distribution of #pi^{+}",100, 0, 3000);  
    Int_t nentries = (Int_t)tree->GetEntries();
    cout << nentries << endl;

    TLorentzVector lambda(0.,0.,0,2.286);
    Double_t masse[3]= {0.938, 0.494, 0.140};
    TGenPhaseSpace event;
    Double_t weight;

    for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
      
      
      lambda.SetPxPyPzE(px_lambda*1E-3,py_lambda*1E-3,pz_lambda*1E-3,E_lambda*1E-3);   // esprimo tutto in GeV
      lambda_p->Fill(lambda.P());
      //cout << lambda.M() << endl;
      event.SetDecay(lambda, 3, masse);
      weight = event.Generate();
        TLorentzVector *pProton = event.GetDecay(0);
        TLorentzVector *pK    = event.GetDecay(1);
        TLorentzVector *pPi    = event.GetDecay(2);

        p->Fill(pProton->E(), weight);
        K->Fill(pK->E(), weight);
        pi->Fill(pPi->E(), weight);
   }
   TCanvas *c1 = new TCanvas("c1","c1");
   lambda_p->GetXaxis()->SetTitle("p[GeV/c]");
   lambda_p->GetYaxis()->SetTitle("#frac{dN}{dp}");
   lambda_p->Draw();
   c1->SaveAs("lambda_distrib.png");
   TCanvas *c2 = new TCanvas("c2","c2");
   c2->cd();
   p->GetXaxis()->SetTitle("E[GeV]");
   p->GetYaxis()->SetTitle("#frac{dN}{dE}");
   K->SetLineColor(kSpring);
   p->Draw();
   c2->SaveAs("p_distrib.png");

   TCanvas *c3 = new TCanvas("c3","c3");
   c3->cd();
   K->GetXaxis()->SetTitle("E[GeV]");
   K->GetYaxis()->SetTitle("#frac{dN}{dE}");
   K->SetLineColor(kRed);
   K->Draw("same");
   c3->SaveAs("K_distrib.png");

   TCanvas *c4 = new TCanvas("c4","c4");
   c4->cd();
   pi->GetXaxis()->SetTitle("E[GeV]");
   pi->GetYaxis()->SetTitle("#frac{dN}{dE}");
   pi->SetLineColor(kCyan);
   pi->Draw("same");
   c4->SaveAs("pi_distrib.png");
}
