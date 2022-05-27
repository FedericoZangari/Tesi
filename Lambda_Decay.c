void Lambda_Decay(){
    TFile *f = new TFile("Lc_thCy0_T2cm_Ge5mrad7cm.root");
    TTree *tree = (TTree*)f->Get("MCDecayTree");

    TFile *f1= new TFile("tree_lambda.root","recreate");
    TTree t1("t1","momentum of lambda and decay product");
    
    Double_t E_lambda, px_lambda, py_lambda, pz_lambda;
    Double_t E_p, px_p, py_p, pz_p;
    Double_t E_K, px_K, py_K, pz_K;
    Double_t E_pi, px_pi, py_pi, pz_pi;
    tree->SetBranchAddress("Lambda_cplus_TRUEP_E", &E_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_X", &px_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Y", &py_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Z", &pz_lambda);

    t1.Branch("lambda_E",&E_lambda,"lambda_E/D");
    t1.Branch("lambda_px",&px_lambda,"lambda_px/D");
    t1.Branch("lambda_py",&py_lambda,"lambda_py/D");
    t1.Branch("lambda_pz",&pz_lambda,"lambda_pz/D");
    t1.Branch("p_E",&E_p,"p_E/D");
    t1.Branch("p_px",&px_p,"p_px/D");
    t1.Branch("p_py",&py_p,"p_py/D");
    t1.Branch("p_pz",&pz_p,"p_pz/D");
    t1.Branch("K_E",&E_K,"K_E/D");
    t1.Branch("K_px",&px_K,"K_px/D");
    t1.Branch("K_py",&py_K,"K_py/D");
    t1.Branch("K_pz",&pz_K,"K_pz/D");
    t1.Branch("pi_E",&E_pi,"pi_E/D");
    t1.Branch("pi_px",&px_pi,"pi_px/D");
    t1.Branch("pi_py",&py_pi,"pi_py/D");
    t1.Branch("pi_pz",&pz_pi,"pi_pz/D");

    TH1F *lambda_p   = new TH1F("#Lambda_{c}^{+}","Momentum distribution of #Lambda_{c}^{+}",100,400,7000);
    TH1F *p= new TH1F("p", "Momentum distribution of protons",100, 0, 4000);
    TH1F *K= new TH1F("K^{-}", "Momentum distribution of K^{-}",100, 0, 4000);
    TH1F *pi= new TH1F("#pi^{+}", "Momentum distribution of #pi^{+}",100, 0, 3000);  
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
      E_lambda =  E_lambda*1E-3;
      px_lambda = px_lambda*1E-3;
      py_lambda = py_lambda*1E-3;
      pz_lambda = pz_lambda*1E-3;

      event.SetDecay(lambda, 3, masse);
      weight = event.Generate();
      TLorentzVector *pProton = event.GetDecay(0);
      TLorentzVector *pK    = event.GetDecay(1);
      TLorentzVector *pPi    = event.GetDecay(2);

      E_p = pProton->E();
      px_p = pProton->Px();
      py_p = pProton->Py();
      pz_p = pProton->Pz();

      E_K = pK->E();
      px_K = pK->Px();
      py_K = pK->Py();
      pz_K = pK->Pz();

      E_pi = pPi->E();
      px_pi = pPi->Px();
      py_pi = pPi->Py();
      pz_pi = pPi->Pz();

      p->Fill(pProton->P(), weight);
      K->Fill(pK->P(), weight);
      pi->Fill(pPi->P(), weight);

      t1.Fill();
   }
   t1.Write();

   TCanvas *c1 = new TCanvas("c1","c1");
   lambda_p->GetXaxis()->SetTitle("p[GeV/c]");
   lambda_p->GetYaxis()->SetTitle("#frac{dN}{dp}");
   lambda_p->Draw();
   c1->SaveAs("lambda_distrib.png");
   TCanvas *c2 = new TCanvas("c2","c2");
   c2->cd();
   p->GetXaxis()->SetTitle("p[GeV/c]");
   p->GetYaxis()->SetTitle("#frac{dN}{dp}");
   K->SetLineColor(kSpring);
   p->Draw();
   c2->SaveAs("p_distrib.png");

   TCanvas *c3 = new TCanvas("c3","c3");
   c3->cd();
   K->GetXaxis()->SetTitle("p[GeV/c]");
   K->GetYaxis()->SetTitle("#frac{dN}{dp}");
   K->SetLineColor(kRed);
   K->Draw("same");
   c3->SaveAs("K_distrib.png");

   TCanvas *c4 = new TCanvas("c4","c4");
   c4->cd();
   pi->GetXaxis()->SetTitle("p[GeV/c]");
   pi->GetYaxis()->SetTitle("#frac{dN}{dp}");
   pi->SetLineColor(kCyan);
   pi->Draw("same");
   c4->SaveAs("pi_distrib.png");

  
}
