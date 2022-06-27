void Momentum_Reconstruction(){
    string title;
    //title = "Lc_thCy0_T2cm_Ge5mrad7cmx5.root";
    title = "Lc_thCy0_T2cm_Ge7mrad7cmx5.root";
    //title = "Lc_thCy0_T2cm_Si5mrad7cmx5.root";
    //title = "Lc_thCy0_T2cm_Si7mrad7cmx5.root";
    string material = title.substr (14,2);
    string bending = title.substr (16,5);
    cout << endl << material << "\t" << bending << endl;
    
    TFile *f = new TFile(title.c_str());    
    TTree *tree = (TTree*)f->Get("MCDecayTree");

    Double_t E_p, px_p, py_p, pz_p;
    Double_t E_k, px_k, py_k, pz_k;
    Double_t E_pi, px_pi, py_pi, pz_pi;
    tree->SetBranchAddress("pplus_TRUEP_E", &E_p);
    tree->SetBranchAddress("pplus_TRUEP_X", &px_p);
    tree->SetBranchAddress("pplus_TRUEP_Y", &py_p);
    tree->SetBranchAddress("pplus_TRUEP_Z", &pz_p);

    tree->SetBranchAddress("Kminus_TRUEP_E", &E_k);
    tree->SetBranchAddress("Kminus_TRUEP_X", &px_k);
    tree->SetBranchAddress("Kminus_TRUEP_Y", &py_k);
    tree->SetBranchAddress("Kminus_TRUEP_Z", &pz_k);

    tree->SetBranchAddress("piplus_TRUEP_E", &E_pi);
    tree->SetBranchAddress("piplus_TRUEP_X", &px_pi);
    tree->SetBranchAddress("piplus_TRUEP_Y", &py_pi);
    tree->SetBranchAddress("piplus_TRUEP_Z", &pz_pi);
        
    Double_t nentries = (Double_t)tree->GetEntries();
    cout << nentries << endl;
    TLorentzVector p, K, pi;
    TLorentzVector p_gen, K_gen, pi_gen;
    TLorentzVector lambda_gen;
    

    //Double_t B= 1.4, L=3.4, R_B = 3;
    Double_t B= 1.1, L=1.7, R_B=3;
    
                                  
    
    Double_t theta_y, theta_x, theta_y_in, theta_y_out, theta_x_mis;            //tutti gli angoli sono espressi in radianti
    Double_t theta_bend_1 = 50.E-6;
    Double_t theta_bend_2 = stod(bending.substr(0,1))*1.E-3;

    Double_t y, y_real;                                    // tutte le dispersioni in y sono in cm
    vector<Double_t> track_y;
    Double_t x, x_real;
    vector<Double_t> track_x;
    vector<Double_t> z;
    Double_t p_mis;
    Double_t D= 50;
    Double_t sigma_p;
    Double_t masse[3]= {0.938272081, 0.493677, 0.13957039};
    vector<bool> index = {0,0,0};

    TGraph * sigma = new TGraph();
    TRandom3 *gRandom = new  TRandom3(1);
    gRandom->SetSeed(1);
    string title2 = " Invariant mass of #Lambda_{c}^{+} generated from tracks (" + material +" "+ bending+ ")";
    TH1F *lambda_M_gen   = new TH1F("#Lambda_{c}^{+}",title2.c_str(),40, 2,2.5);
      
    for (Int_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
      
        p.SetPxPyPzE(px_p*1E-3,py_p*1E-3,pz_p*1E-3,E_p*1E-3);   // esprimo tutto in GeV
        K.SetPxPyPzE(px_k*1E-3,py_k*1E-3,pz_k*1E-3,E_k*1E-3);
        pi.SetPxPyPzE(px_pi*1E-3,py_pi*1E-3,pz_pi*1E-3,E_pi*1E-3);

        TLorentzVector v= p;
        for(int j= 0; j < 3; j++){                  // eseguo un ciclo per ogni particella
            
            theta_y = v.Py() / v.Pz() + theta_bend_1 + theta_bend_2;        // angolo vero della particella lungo y, considerando l'apertura del decadimento e le deviazioni dei 2 cristalli
            theta_x = v.Px() / v.Pz();
            y =  theta_bend_1*100.*1.E2  +  gRandom->Gaus(0, 0.05); 
            x =  gRandom->Gaus(0, 0.05); 

            TGraph *trk_y = new TGraph();
            TGraph *trk_x = new TGraph();
            for(int k=0; k < 4; k++){                  // un ciclo per tracker
                z.push_back(1.*100. - D +  k * D / 3.);   // z misurata rispetto al secondo cristallo
                y_real= y + tan(theta_y) * z[k];
                x_real= x + tan(theta_x) * z[k];
                track_y.push_back( gRandom->Gaus( y_real, 10E-4 ));             // genero le tracce
                track_x.push_back( gRandom->Gaus( x_real, 10E-4 ));
                theta_y += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );          //aggiungo una deviazione dovuta allo scattering multiplo sul tracker
                theta_x += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );
                trk_y->SetPoint(k, z[k], track_y[k]);
                trk_x->SetPoint(k, z[k], track_x[k]);
            }

            trk_y->Fit("pol1","Q");
            TF1 *fit = trk_y->GetFunction("pol1");
            theta_y_in = atan( fit->GetParameter(1));
            //cout << theta_y << "\t" << theta_y_in << endl;
            z.clear();
            track_y.clear();

            if(j != 1) theta_y += 0.3 * L * B / v.P();                    // aggiungo la deflessione del magnete all'angolo y vero
            if(j == 1) theta_y -= 0.3 * L * B / v.P();

            for(int k=0; k < 4; k++){                  // un ciclo per tracker
                z.push_back(1*100 + L +  k * D / 3.);
                y_real= y + tan( theta_y ) * z[k];
                x_real= x + tan( theta_x ) * z[k];
                track_y.push_back( gRandom->Gaus( y_real, 10E-4 ));             // genero le tracce dopo il magnete
                track_x.push_back( gRandom->Gaus( x_real, 10E-4 )); 
                theta_y += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );       //scattering multiplo
                theta_x += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );
                trk_y->SetPoint(k, z[k], track_y[k]);
                trk_x->SetPoint(k+4, z[k], track_x[k+4]);
            }
            
            trk_y->Fit("pol1","Q");
            fit = trk_y->GetFunction("pol1");
            theta_y_out = atan( fit->GetParameter(1));

            trk_x->Fit("pol1","Q");
            fit = trk_x->GetFunction("pol1");
            theta_x_mis = atan( fit->GetParameter(1));
            
            //cout << theta_x << "\t" << theta_x_mis << endl;
            z.clear();
            track_y.clear();
            track_x.clear();

            p_mis = 0.3* L * B / abs(theta_y_out - theta_y_in);     //in GeV
            sigma_p = abs(p_mis - v.P());
            if(track_y[0] < R_B){
                index[j] = 1;
                sigma->SetPoint(i*3+j, v.P(), sigma_p );           // considero solo le particelle che riescono ad uscire dal magnete
            }
            if(j==0) v=K;
            if(j==1) v=pi;
            Double_t pz = p_mis / sqrt( pow( tan(theta_x_mis),2 ) + pow( tan(theta_y_in),2) + 1 );
            //cout << pz << endl;
            if(j==0) p_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
            if(j==1) K_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
            if(j==2) pi_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
        }
        lambda_gen = p_gen + K_gen + pi_gen;
        if(index[0]==1 && index[1]==1 &&  index[2]==1 )lambda_M_gen->Fill(lambda_gen.M());
        //cout << lambda_gen.M() << endl;
        
    }
    
    TCanvas *c = new TCanvas("c","c");
    TF1 * func = new TF1("func", "[0]*x^2", 0, 4000);
    func->SetParameter(0, 1./L/B);
    TF1 *f2 = new TF1("sigma","[0]/0.3*x^2",0, 4000);
    f2->SetParameter(0, 2* 10E-6/L/B/D*100);
    //f2->SetParameter(0, 3*sqrt(2./5.)* 10E-6/L/B/D*100);
    f2->SetLineColor(kSpring);
    TF1 *f3 = new TF1("sigma2","2*[0]/0.3*x^2",0, 4000);
    f3->SetParameter(0, 10E-6/L/B/L);
    f3->SetLineColor(kBlue);
    sigma->Fit(func);
    sigma->SetMarkerSize(1);
    TLegend *l = new TLegend(0.15,0.65,0.4,0.89);
    l->AddEntry((TObject*)0, TString::Format("D = %g m" , D/100.), "");
    l->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l->AddEntry((TObject*)0, TString::Format("R_{B} = %g cm", R_B), "");
    l->AddEntry(sigma,"Incertezze generate", "p");
    l->AddEntry(func,"Fit incertezze generate", "l");
    l->AddEntry(f2, "#sigma_{p}#propto 1/LD", "l");
    l->AddEntry(f3, "#sigma_{p}#propto 1/L^{2}", "l");
    sigma->SetTitle("Incertezza impulso");
    sigma->GetXaxis()->SetTitle("p[GeV/c]");
    sigma->GetYaxis()->SetTitle("#sigma_{p} [GeV/c]");
    sigma->Draw("AP");
    func->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    l->Draw("same");
    cout << "1/(LD)\t" << 2/0.3/L/B/D*10E-6*100 << endl;
    cout << "1/L^2\t" << 2/0.3/L/B/L*10E-6 << endl;


    TCanvas *c1 =new TCanvas("c1","c1", 900, 550);
    lambda_M_gen->GetXaxis()->SetTitle("m_{#Lambda}[GeV/c^{2}]");
    lambda_M_gen->GetYaxis()->SetTitle("#frac{dN}{dm}");
    lambda_M_gen->Draw();
    TLegend *l1 = new TLegend(0.15,0.65,0.4,0.89);
    l1->AddEntry((TObject*)0, TString::Format("D = %g m" , D/100.), "");
    l1->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l1->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l1->AddEntry((TObject*)0, TString::Format("R_{B} = %g cm", R_B), "");
    l1->Draw("same");
    
    c->SaveAs("Incertezze_p_gen.png");
}
