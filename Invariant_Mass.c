void Invariant_Mass(){
    //string title = "Lc_thCy0_T2cm_Ge5mrad7cm.root";
    //string title = "Lc_thCy0_T2cm_Ge7mrad7cm.root";
    //string title = "Lc_thCy0_T2cm_Si5mrad7cm.root";
    string title = "Lc_thCy0_T2cm_Si7mrad7cm.root";
    string material = title.substr (14,2);
    string bending = title.substr (16,5);
    TFile *f = new TFile(title.c_str());     
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
    
    string title2 = " Invariant mass of generated #Lambda_{c}^{+} " + material +" "+ bending;
    TH1F *lambda_M_gen   = new TH1F("#Lambda_{c}^{+}",title2.c_str(),40, 2,2.55);
    
    Double_t masse[3]= {0.938, 0.494, 0.140};
    TGenPhaseSpace event;
    Double_t weight;

    Double_t B= 1.4, L=3.4, D= 0.4, sigma_trk = 10E-6, N=4;
    //double B= 1.1, L=1.7, D= 0.4, sigma_trk = 10E-6;
    Double_t sigma_p= 2*sigma_trk/0.3/B/L/D;
    Double_t sigma_theta= sigma_trk/D*sqrt(12*(N-1)/N/(N+1));

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


            p_gen.SetRho(gRandom->Gaus(pProton->P(), sigma_p*pProton->P()*pProton->P()));
            p_gen.SetTheta(gRandom->Gaus(pProton->Theta(), sigma_theta));
            p_gen.SetPhi(gRandom->Gaus(pProton->Phi(), sigma_theta));
            p_gen.SetE(sqrt(p_gen.P()*p_gen.P()+masse[0]*masse[0]));
         
            K_gen.SetRho(gRandom->Gaus(pK->P(), sigma_p*pK->P()*pK->P()));
            K_gen.SetTheta(gRandom->Gaus(pK->Theta(), sigma_theta));
            K_gen.SetPhi(gRandom->Gaus(pK->Phi(), sigma_theta));
            K_gen.SetE(sqrt(K_gen.P()*K_gen.P()+masse[1]*masse[1]));
        
            pi_gen.SetRho(gRandom->Gaus(pPi->P(), sigma_p*pPi->P()*pPi->P()));
            pi_gen.SetTheta(gRandom->Gaus(pPi->Theta(), sigma_theta));
            pi_gen.SetPhi(gRandom->Gaus(pPi->Phi(), sigma_theta));
            pi_gen.SetE(sqrt(pi_gen.P()*pi_gen.P()+masse[2]*masse[2]));
        

            lambda_gen = p_gen + K_gen + pi_gen;
            lambda_M_gen->Fill(lambda_gen.M(),weight);
                
        }
    }

    TCanvas *c1 =new TCanvas("c1","c1");
    lambda_M_gen->GetXaxis()->SetTitle("m_{#Lambda}[GeV/c^{2}]");
    lambda_M_gen->GetYaxis()->SetTitle("#frac{dN}{dm}");
    lambda_M_gen->Fit("gaus");
    gStyle->SetOptFit(0101);
    lambda_M_gen->Draw();
    TLegend *l = new TLegend(0.2,0.7,0.35,0.89);
    l->AddEntry((TObject*)0, TString::Format("D = %g m" , D), "");
    l->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l->Draw("same");
    
    c1->SaveAs("invariant_mass.png");
   
}
