void dispersioni_xy(){
    string title;
    title = "Lc_thCy0_T2cm_Ge5mrad7cm.root";
    //title = "Lc_thCy0_T2cm_Ge7mrad7cm.root";
    //title = "Lc_thCy0_T2cm_Si5mrad7cm.root";
    //title = "Lc_thCy0_T2cm_Si7mrad7cm.root";
    string material = title.substr (14,2);
    string bending = title.substr (16,5);
    cout << endl << material << "\t" << bending << endl;
    Double_t theta_bend = stod(bending.substr(0,1));
    TFile *f = new TFile(title.c_str());    
    TTree *tree = (TTree*)f->Get("MCDecayTree");

    Double_t E_lambda, px_lambda, py_lambda, pz_lambda;
    tree->SetBranchAddress("Lambda_cplus_TRUEP_E", &E_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_X", &px_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Y", &py_lambda);
    tree->SetBranchAddress("Lambda_cplus_TRUEP_Z", &pz_lambda);
        
    Double_t nentries = (Double_t)tree->GetEntries();
    //cout << nentries << endl;

    TLorentzVector lambda(0.,0.,0,2.286);
               
    Double_t masse[3]= {0.938, 0.494, 0.140};
    TGenPhaseSpace event;
    Double_t weight;

    //Double_t B= 1.4, L=3.4;
    Double_t B= 1.1, L=1.7;
    vector<Double_t> counts;
    Double_t y;                               // tutte le dispersioni in y sono in cm
    Double_t x;

    TGraph *p_distrib = new TGraph();
    TGraph *k_distrib = new TGraph();
    TGraph *pi_distrib = new TGraph();
    int k=0, w=0, z=0;

    for(int j=0; j < 100; j++){    
        for (Int_t i=0; i<nentries; i++) {
            tree->GetEntry(i);
      
            lambda.SetPxPyPzE(px_lambda*1E-3,py_lambda*1E-3,pz_lambda*1E-3,E_lambda*1E-3);   // esprimo tutto in GeV
      
            event.SetDecay(lambda, 3, masse);
            weight = event.Generate();

            TLorentzVector *p = event.GetDecay(0);
            TLorentzVector *K    = event.GetDecay(1);
            TLorentzVector *Pi    = event.GetDecay(2);

            y = 0.5 + theta_bend *1.E-3*(L+1.)*1.E2;
            y += 0.3*L*L*B/2./p->P()*1.E2;
            y += p->Py() / p->Pz() * (L+1.)*1.E2;
            y += gRandom->Gaus(0,0.05);
            x = p->Px() / p->Pz() * (L+1)*1.E2;
            x += gRandom->Gaus(0,0.05);
            p_distrib->SetPoint(k, x, y);
            k++;

            y = 0.5 + theta_bend *1.E-3*(L+1.)*1.E2;
            y += -0.3*L*L*B/2./K->P()*1.E2;
            y += K->Py()/K->Pz() *(L+1.)*1.E2;
            y += gRandom->Gaus(0,0.05);
            x = K->Px() / K->Pz() * (L+1)*1E2;
            x += gRandom->Gaus(0,0.05);
            k_distrib->SetPoint(w, x, y);
            w++;
                    
            y = 0.5 + theta_bend *1.E-3*(L+1.)*1.E2;
            y += 0.3*L*L*B/2./Pi->P()*1.E2;
            y += Pi->Py() / Pi->Pz() * (L+1.)*1.E2;
            y += gRandom->Gaus(0,0.05);
            x = Pi->Px() / Pi->Pz() * (L+1)*1E2;
            x += gRandom->Gaus(0,0.05);
            pi_distrib->SetPoint(z, x, y);
            z++;  
        }
    }
    TCanvas *c = new TCanvas("c","c");
    
    p_distrib->SetMarkerStyle(8);
    p_distrib->SetMarkerColorAlpha(kRed,0.25);
	p_distrib->SetMarkerSize(1);

    k_distrib->SetMarkerStyle(8);
    k_distrib->SetMarkerColorAlpha(kOrange,0.25);
	k_distrib->SetMarkerSize(1);
    
    pi_distrib->SetMarkerStyle(8);
    pi_distrib->SetMarkerColorAlpha(kSpring,0.25);
	pi_distrib->SetMarkerSize(1);

    TMultiGraph *g = new TMultiGraph();
    g->Add(pi_distrib);
    g->Add(k_distrib);
    g->Add(p_distrib);
    g->SetTitle("Dispersione dopo il magnete");
    g->Draw("AP");
    g->GetXaxis()->SetTitle("x [cm]");
    g->GetYaxis()->SetTitle("y [cm]");
    
    TLegend *l = new TLegend(0.8,0.1,0.89,0.29);
    l->AddEntry(p_distrib, "p", "p");
    l->AddEntry(k_distrib, "K^{-}", "p");
    l->AddEntry(pi_distrib, "#pi^{+}", "p");
    l->Draw("same");
    
    c->SaveAs("dispersioni_xy.png");
}