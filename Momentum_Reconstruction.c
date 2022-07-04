using namespace RooFit;

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
    TLorentzVector p, K, pi, v;
    TLorentzVector p_gen, K_gen, pi_gen;
    TLorentzVector lambda_gen;
    

    //Double_t B= 1.4, L=3.4, R_B = 4;
    Double_t B= 1.1, L=1.7, R_B=4;                     
    
    Double_t theta_y, theta_y_new, theta_x, theta_y_in, theta_y_out, theta_x_mis;            //tutti gli angoli sono espressi in radianti
    Double_t theta_bend_1 = 50.E-6;
    Double_t theta_bend_2 = stod(bending.substr(0,1))*1.E-3;

    Double_t y, y_B, y_real;                                    // tutte le dispersioni in y sono in cm
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
    //string title2 = " Invariant mass of #Lambda_{c}^{+}, only angles contribution  (" + material +" "+ bending+ ")";
    string title2 = " Invariant mass of #Lambda_{c}^{+}  (" + material +" "+ bending+ ")";
    TH1F *lambda_M_gen   = new TH1F("#Lambda_{c}^{+}",title2.c_str(),40, 2,2.55);
    lambda_M_gen->Sumw2();
    TGraph *trk_y = new TGraph();
    TGraph *trk_x = new TGraph();
    //TGraph *trk_y1 = new TGraph();
    //TGraph *trk_x1 = new TGraph();
    TF1 *fit1 = new TF1("fit1", "pol1", 0, 101);
    TF1 *fit2 = new TF1("fit2", "pol1", 99 + L*100, 100+L*100 +D+1);
    TCanvas *c2 = new TCanvas("c2","c2", 1000, 800);
    c2->Divide(1,2);
    
    Int_t rec = 0;
      
    for (Int_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
      
        p.SetPxPyPzE(px_p*1E-3,py_p*1E-3,pz_p*1E-3,E_p*1E-3);   // esprimo tutto in GeV
        K.SetPxPyPzE(px_k*1E-3,py_k*1E-3,pz_k*1E-3,E_k*1E-3);
        pi.SetPxPyPzE(px_pi*1E-3,py_pi*1E-3,pz_pi*1E-3,E_pi*1E-3);

        v = p;
        y =  theta_bend_1*100.*1.E2  +  gRandom->Gaus(0, 0.05); 
        x =  gRandom->Gaus(0, 0.05);

        for(int j= 0; j < 3; j++){                  // eseguo un ciclo per ogni particella
            
            theta_y = atan( v.Py() / v.Pz()) + theta_bend_1 + theta_bend_2;        // angolo vero della particella lungo y, considerando l'apertura del decadimento e le deviazioni dei 2 cristalli
            theta_x = atan( v.Px() / v.Pz());
            theta_y_new = theta_y;
             
            //cout << theta_y << "\t" << v.P() << "\t" << y << endl;
            
        
            for(int k=0; k < 4; k++){                  // un ciclo per tracker
                z.push_back(1.*100. - D +  k * D / 3.);   // z misurata rispetto al secondo cristallo
                y_real= y + tan(theta_y_new) * z[k];
                x_real= x + tan(theta_x) * z[k];
                
                track_y.push_back( gRandom->Gaus( y_real, 10E-4 ));             // genero le tracce
                track_x.push_back( gRandom->Gaus( x_real, 10E-4 ));
                theta_y_new += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );          //aggiungo una deviazione dovuta allo scattering multiplo sul tracker
                theta_x += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );
                trk_y->SetPoint(k, z[k], track_y[k]);
                trk_x->SetPoint(k, z[k], track_x[k]);
                //cout << z[k] << "\t" <<  x_real<< endl;
            }

            trk_y->Fit(fit1,"QR");
            theta_y_in = atan( fit1 ->GetParameter(1));
            //cout << theta_y << "\t" << theta_y_in << endl;
            
            
            y_real= y + tan(theta_y_new) * (1*100 + L*100);
            if(j != 1) {theta_y_new = theta_y_new + 0.3 * L * B / v.P();     y_B = 0.3*L*L*B/2./v.P()*1.E2 ; }                  // aggiungo la deflessione del magnete all'angolo y vero
            if(j == 1) {theta_y_new = theta_y_new - 0.3 * L * B / v.P();     y_B = - 0.3*L*L*B/2./v.P()*1.E2 ; }

            
            for(int k=0; k < 4; k++){                  // un ciclo per tracker
                z.push_back(1*100 + L*100 +  k * D / 3.);
                
                x_real= x + tan( theta_x ) * z[k+4];
                //y_real= y + tan( theta_y ) * z[k] + y_B;
                if(k==0) y_real += y_B ;
                else    y_real += tan( theta_y_new ) * D / 3.;
                //cout << x_real << endl;
                track_y.push_back( gRandom->Gaus( y_real, 10E-4 ));             // genero le tracce dopo il magnete
                track_x.push_back( gRandom->Gaus( x_real, 10E-4 )); 
                theta_y_new += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );       //scattering multiplo
                theta_x += gRandom->Gaus( 0,  13.6 / (v.Beta() * v.P()*1000.) * sqrt(10E-4/9.37) );
                trk_y->SetPoint(k+4, z[k+4], track_y[k+4]);
                trk_x->SetPoint(k+4, z[k+4], track_x[k+4]);
                //cout << z[k+4] << "\t" << track_y[k+4] << "\t\t" << y_real[k+4+1] << endl;
            }
            //cout << endl;
            trk_y->Fit(fit2,"QR");
            theta_y_out = atan( fit2->GetParameter(1));
            
            trk_x->Fit("pol1","Q");
            TF1 * fit = trk_x->GetFunction("pol1");
            theta_x_mis = atan( fit->GetParameter(1));
            
            //cout << theta_x << "\t" << theta_x_mis << endl;
            

            p_mis = 0.3* L * B / abs(theta_y_out - theta_y_in);     //in GeV
            sigma_p = abs(p_mis - v.P());
            if(track_y[4] < R_B){
                index[j] = 1;
                sigma->SetPoint(i*3+j, v.P(), sigma_p );           // considero solo le particelle che riescono ad uscire dal magnete
            }   //else{index[j]=0;}
            
            Double_t pz = p_mis / sqrt( pow( tan(theta_x_mis),2 ) + pow( tan(theta_y_in),2) + 1 );
            //cout << pz << endl;
            if(j==0) p_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
            if(j==1) K_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
            if(j==2) pi_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
            /*
            Double_t pz = p_mis / sqrt( pow( tan(theta_x),2 ) + pow( tan(theta_y),2) + 1 );
            if(j==0) p_gen.SetPxPyPzE(pz*tan(theta_x),pz*tan(theta_y),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));       // studio il contrib dell'impulso
            if(j==1) K_gen.SetPxPyPzE(pz*tan(theta_x),pz*tan(theta_y),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));
            if(j==2) pi_gen.SetPxPyPzE(pz*tan(theta_x),pz*tan(theta_y),pz,sqrt(p_mis*p_mis+masse[j]*masse[j]));

            Double_t pz = v.P() / sqrt( pow( tan(theta_x_mis),2 ) + pow( tan(theta_y_in),2) + 1 );
            if(j==0) p_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(v.P()*v.P()+masse[j]*masse[j]));        // studio il contrib degli angoli
            if(j==1) K_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(v.P()*v.P()+masse[j]*masse[j]));
            if(j==2) pi_gen.SetPxPyPzE(pz*tan(theta_x_mis),pz*tan(theta_y_in),pz,sqrt(v.P()*v.P()+masse[j]*masse[j]));
            */
            if(i==nentries-1 && j == 2) {
                TGraph * trk_x1 = trk_x;
                TGraph * trk_y1 = trk_y;
                trk_y1->SetMarkerSize(2);
                trk_y1->SetMarkerColor(kBlue);
                trk_y1->SetMarkerStyle(43);
                trk_y1->SetTitle("Track reconstruction of a #pi^{+} in y direction");
                trk_y1->GetXaxis()->SetTitle("z [cm]");
                trk_y1->GetYaxis()->SetTitle("y [cm]");
                trk_x1->SetTitle("Track reconstruction of a #pi^{+} in x direction");
                trk_x1->GetXaxis()->SetTitle("z [cm]");
                trk_x1->GetYaxis()->SetTitle("x [cm]");
                c2->cd(1);
                trk_y1->Draw("APC");
                fit1->Draw("same");
                fit2->Draw("same");
                TLegend *l2 = new TLegend(0.15,0.65,0.4,0.89);
                l2->AddEntry((TObject*)0, TString::Format("p = %g GeV/c" , pi_gen.P()), "");
                l2->AddEntry(trk_y1,"Generated tracks", "pl");
                l2->AddEntry(fit1,"Fit of tracks", "l");
                l2->Draw("same");
                //TCanvas *c3 = new TCanvas("c3","c3", 900, 550);
                c2->cd(2);
                trk_x1->SetMarkerSize(2);
                trk_x1->SetMarkerColor(kBlue);
                trk_x1->SetMarkerStyle(43);
                trk_x1->Draw("APL");
            }
            if(j==0) v=K;
            if(j==1) v=pi;
            z.clear();
            track_y.clear();
            track_x.clear();
        }
        lambda_gen = p_gen + K_gen + pi_gen;
        if(index[0]==1 && index[1]==1 &&  index[2]==1 ) {lambda_M_gen->Fill(lambda_gen.M()); rec++;}
        //cout << lambda_gen.M() << endl;
        index = {0,0,0};
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
    TLegend *l = new TLegend(0.15,0.6,0.4,0.89);
    l->AddEntry((TObject*)0, TString::Format("D = %g m" , D/100.), "");
    l->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l->AddEntry((TObject*)0, TString::Format("R_{B} = %g cm", R_B), "");
    l->AddEntry(sigma,"Generated errors", "p");
    l->AddEntry(func,"Fit generated errors", "l");
    l->AddEntry(f2, "#sigma_{p}#propto 1/LD", "l");
    l->AddEntry(f3, "#sigma_{p}#propto 1/L^{2}", "l");
    string title3 = "Momentum errors (" + material +" "+ bending+ ")";
    sigma->SetTitle(title3.c_str());
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

    RooRealVar m("m", "m", 0, 10);
    RooDataHist h("h", "h", RooArgSet(m), Import(*lambda_M_gen));
    SumW2Error(true);
    RooRealVar mean("mean", "mean", 2.286, 2, 2.55);
    RooRealVar sigma1("sigma1", "sigma1", 0.01, 0., 1);
    RooGaussian gauss1("gauss1", "gauss1", m, mean, sigma1);
    
    RooRealVar sigma2("sigma2", "sigma2", 0.1, 0., 1);
    RooGaussian gauss2("gauss2", "gauss2", m, mean, sigma2);
    RooRealVar frac("frac", "fraction of component 1 in signal", 0.5, 0., 1.);
    RooAddPdf model("model","2 gaus",RooArgList(gauss1,gauss2),frac);
   
    model.chi2FitTo(h);
    
    RooPlot *frame = m.frame(Title(title2.c_str()));
    h.plotOn(frame);
    model.plotOn(frame);
    
    frame->GetXaxis()->SetTitle("m(pK^{-}#pi^{+})[GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("#frac{dN}{dm}");
    frame->Draw();
    Double_t sigma_fit= sqrt(pow(frac.getValV()*sigma1.getValV(),2) + pow((1-frac.getValV())*sigma2.getValV(),2));

    TLegend *l1 = new TLegend(0.15,0.65,0.4,0.89);
    l1->AddEntry((TObject*)0, TString::Format("D = %g m" , D/100.), "");
    l1->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l1->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l1->AddEntry((TObject*)0, TString::Format("R_{B} = %g cm", R_B), "");
    l1->AddEntry((TObject*)0, TString::Format("#sigma_{m} = %g MeV/c^{2}", sigma_fit*1000.), "");
    l1->AddEntry("h","Generated mass", "P");
    l1->AddEntry("model","Fit of mass distribution", "L");
    l1->Draw("same");
    
    //c->SaveAs("Incertezze_p_gen.png");
    //c1->SaveAs("Invariant_mass_distribution.png");
    //c1->SaveAs("Invariant_mass_distribution_p_contrib.png");
    //c1->SaveAs("Invariant_mass_distribution_theta_contrib.png");
    //c2->SaveAs("generated_tracks_pi.png");

    cout << rec << endl << rec/nentries*100 << endl;
    //cout << v.P() << endl;
    /*    TCanvas *c2 = new TCanvas("c2","c2", 900, 550);
    trk_y1->SetMarkerSize(2);
    trk_y1->SetMarkerColor(kBlue);
    trk_y1->SetMarkerStyle(43);
    trk_y1->Draw("APC");
    //fit1->Draw("same");
    //fit2->Draw("same");
    TCanvas *c3 = new TCanvas("c3","c3", 900, 550);
    trk_x1->SetMarkerSize(2);
    trk_x1->SetMarkerColor(kBlue);
    trk_x1->SetMarkerStyle(43);
    trk_x1->Draw("APL");*/
}
