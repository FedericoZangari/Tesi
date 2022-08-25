using namespace RooFit;

void Invariant_Mass(){
    //string title = "Lc_thCy0_T2cm_Ge5mrad7cmx5.root";
    //string title = "Lc_thCy0_T2cm_Ge7mrad7cmx5.root";
    string title = "Lc_thCy0_T2cm_Si5mrad7cmx5.root";
    //string title = "Lc_thCy0_T2cm_Si7mrad7cmx5.root";
    string material = title.substr (14,2);
    string bending = title.substr (16,5);
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
         
    Int_t nentries = (Int_t)tree->GetEntries();
    cout << nentries << endl;

    TLorentzVector lambda(0.,0.,0,2.286);
    TLorentzVector lambda_gen(0.,0.,0,2.286);
    
    string title2 = " Invariant mass of generated #Lambda_{c}^{+} " + material +" "+ bending;
    TH1F *lambda_M_gen   = new TH1F("#Lambda_{c}^{+}",title2.c_str(),30, 2,2.55);
    lambda_M_gen->Sumw2();    
    
    Double_t masse[3]= {0.938272081, 0.493677, 0.13957039};
    //TGenPhaseSpace event;
    //Double_t weight;

    //Double_t B= 1.4, L=3.4, D= 0.8, sigma_trk = 10E-6, N=4;
    double B= 1.1, L=1.7, D= 0.8, sigma_trk = 10E-6,  N=4;
    //Double_t sigma_p= 2*sigma_trk/0.3/B/L/D;
    Double_t sigma_p= 2*sigma_trk/0.3/B/L/L;
    Double_t sigma_theta= sigma_trk/D*sqrt(12*(N-1)/N/(N+1));
    //Double_t sigma_theta= sigma_trk/D*sqrt(2);
    Double_t sigma;

    TLorentzVector *p = new TLorentzVector();
    TLorentzVector *K = new TLorentzVector();
    TLorentzVector *Pi = new TLorentzVector();

    gRandom->SetSeed(1);

    
        for (Int_t i=0; i<nentries; i++) {
            tree->GetEntry(i);
      
            p->SetPxPyPzE(px_p*1E-3,py_p*1E-3,pz_p*1E-3,E_p*1E-3);   // esprimo tutto in GeV
            K->SetPxPyPzE(px_k*1E-3,py_k*1E-3,pz_k*1E-3,E_k*1E-3);
            Pi->SetPxPyPzE(px_pi*1E-3,py_pi*1E-3,pz_pi*1E-3,E_pi*1E-3);
        
            TLorentzVector p_gen = *p;
            TLorentzVector K_gen = *K;
            TLorentzVector pi_gen = *Pi;


            p_gen.SetRho(gRandom->Gaus(p->P(), sigma_p*p->P()*p->P()));
            p_gen.SetTheta(gRandom->Gaus(p->Theta(), sigma_theta));
            p_gen.SetPhi(gRandom->Gaus(p->Phi(), sigma_theta));
            p_gen.SetE(sqrt(p_gen.P()*p_gen.P()+masse[0]*masse[0]));
         
            K_gen.SetRho(gRandom->Gaus(K->P(), sigma_p*K->P()*K->P()));
            K_gen.SetTheta(gRandom->Gaus(K->Theta(), sigma_theta));
            K_gen.SetPhi(gRandom->Gaus(K->Phi(), sigma_theta));
            K_gen.SetE(sqrt(K_gen.P()*K_gen.P()+masse[1]*masse[1]));
        
            pi_gen.SetRho(gRandom->Gaus(Pi->P(), sigma_p*Pi->P()*Pi->P()));
            pi_gen.SetTheta(gRandom->Gaus(Pi->Theta(), sigma_theta));
            pi_gen.SetPhi(gRandom->Gaus(Pi->Phi(), sigma_theta));
            pi_gen.SetE(sqrt(pi_gen.P()*pi_gen.P()+masse[2]*masse[2]));
        

            lambda_gen = p_gen + K_gen + pi_gen;
            lambda_M_gen->Fill(lambda_gen.M());
                
        }
    

    TCanvas *c1 =new TCanvas("c1","c1", 900, 550);


    RooRealVar x("x", "x", 0, 10);
    RooDataHist h("h", "h", RooArgSet(x), Import(*lambda_M_gen));
    SumW2Error(true);
    RooRealVar mean("mean", "mean", 2.286, 2, 2.55);
    RooRealVar sigma1("sigma1", "sigma1", 0.01, 0., 1);
    RooGaussian gauss1("gauss1", "gauss1", x, mean, sigma1);
    
    RooRealVar sigma2("sigma2", "sigma2", 0.1, 0., 1);
    RooGaussian gauss2("gauss2", "gauss2", x, mean, sigma2);
    RooRealVar frac("frac", "fraction of component 1 in signal", 0.5, 0., 1.);
    RooAddPdf model("model","2 gaus",RooArgList(gauss1,gauss2),frac);
    //model.fitTo(h);
    
    model.chi2FitTo(h);
    //gauss1.fitTo(h);
    RooPlot *frame = x.frame(Title(title2.c_str()));
    h.plotOn(frame);
    model.plotOn(frame);
    RooHist* hresid = frame->residHist() ;
    //RooChi2MCSModule chi();
    //frame->addPlotable(hresid,"P") ;
    //frame2->addPlotable(hresid,"P") ;*/
    //gauss1.plotOn(frame);
    frame->GetXaxis()->SetTitle("m_{#Lambda}[GeV/c^{2}]");
    frame->GetYaxis()->SetTitle("#frac{dN}{dm}");
    frame->Draw();
    
    /*lambda_M_gen->Fit("gaus");
    gStyle->SetOptFit(0101);
    lambda_M_gen->Draw();*/
    sigma= sqrt(pow(frac.getValV()*sigma1.getValV(),2) + pow((1-frac.getValV())*sigma2.getValV(),2));
    TLegend *l = new TLegend(0.2,0.7,0.35,0.89);
    l->AddEntry((TObject*)0, TString::Format("D = %g m" , D), "");
    l->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l->AddEntry((TObject*)0, TString::Format("#sigma_{m} = %g MeV/c^{2}", sigma*1000.), "");
    l->Draw("same");
    
    cout << sigma*1E3 << endl;
    cout << frac.getValV() << endl;
    cout << "chi^2 = " << frame->chiSquare() << endl ;
    //RooAbsReal * chi2_o_ndf = model.createChi2(h, Range("fullRange"), Extended(true), DataError(RooAbsData::Poisson));
    //cout << "chi2_o_ndf: " << chi2_o_ndf->getVal() << endl;
    
    
    //c1->SaveAs("invariant_mass.png");
    string filename = "incertezze_masse_old_tot_"+material+bending+"2.dat";
    ofstream out;
    out.open(filename.c_str(),ios::app);
    out << sigma*1E3 << "\t" << D << "\t" << B << "\t" << L  << endl;
   
}