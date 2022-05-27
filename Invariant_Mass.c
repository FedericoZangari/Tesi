void Invariant_Mass(){
    TFile *f = new TFile("Lc_thCy0_T2cm_Ge5mrad7cm.root");    
    TTree *tree = (TTree*)f->Get("MCDecayTree");
    
    Double_t E_lambda, px_lambda, py_lambda, pz_lambda;
    tree->SetBranchAddress("Lambda_cplus_TRUEP_E", &E_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_X", &px_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Y", &py_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Z", &pz_lambda);
         
    Int_t nentries = (Int_t)tree->GetEntries();
    cout << nentries << endl;

    TLorentzVector lambda(0.,0.,0,2.286);
    TLorentzVector lambda_gen(0.,0.,0,2.286);
    
    
    TH1F *lambda_M_gen   = new TH1F("#Lambda_{c}^{+}"," Invariant mass of generated #Lambda_{c}^{+}",100, 1.5,3);
    
    Double_t masse[3]= {0.938, 0.494, 0.140};
    TGenPhaseSpace event;
    Double_t weight;

    double B= 1.4, L=3.4, D= 0.4, sigma_trk = 10E-6;
    //double B= 1.1, L=1.7, D= 0.4, sigma_trk = 10E-6;
    Double_t sigma= 2*sigma_trk/0.3/B/L/D;
    //Double_t sigma= 3.5E-5;

    for(int j=0; j < 500; j++){    
        for (Int_t i=0; i<nentries; i++) {
            tree->GetEntry(i);
      
            lambda.SetPxPyPzE(px_lambda*1E-3,py_lambda*1E-3,pz_lambda*1E-3,E_lambda*1E-3);   // esprimo tutto in GeV
      
            event.SetDecay(lambda, 3, masse);
            weight = event.Generate();

            TLorentzVector *pProton = event.GetDecay(0);
            TLorentzVector *pK    = event.GetDecay(1);
            TLorentzVector *pPi    = event.GetDecay(2);
        
            TLorentzVector p_gen = *pProton;
            TLorentzVector K_gen = *pK;
            TLorentzVector pi_gen = *pPi;


            p_gen.SetRho(gRandom->Gaus(pProton->P(), sigma*pProton->P()*pProton->P()));
            p_gen.SetE(sqrt(p_gen.P()*p_gen.P()+masse[0]*masse[0]));
         
            K_gen.SetRho(gRandom->Gaus(pK->P(), sigma*pK->P()*pK->P()));
            K_gen.SetE(sqrt(K_gen.P()*K_gen.P()+masse[1]*masse[1]));
        
            pi_gen.SetRho(gRandom->Gaus(pPi->P(), sigma*pPi->P()*pPi->P()));
            pi_gen.SetE(sqrt(pi_gen.P()*pi_gen.P()+masse[2]*masse[2]));
        

            lambda_gen = p_gen + K_gen + pi_gen;
            lambda_M_gen->Fill(lambda_gen.M(),weight);
                
        }
    }
   TCanvas *c1 =new TCanvas("c1","c1");
   lambda_M_gen->GetXaxis()->SetTitle("m_{#Lambda}[GeV/c^{2}]");
   lambda_M_gen->GetYaxis()->SetTitle("#frac{dN}{dm}");
   
   lambda_M_gen->Draw();
   c1->SaveAs("invariant_mass.png");
   
}
