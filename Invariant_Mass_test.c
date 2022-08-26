using namespace RooFit;

void Invariant_Mass_test(){
    //string title = "Lc_thCy0_T2cm_Ge5mrad7cmx5.root";
    //string title = "Lc_thCy0_T2cm_Ge7mrad7cmx5.root";
    string title = "Lc_thCy0_T2cm_Si5mrad7cm.root";
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
    TH1F *lambda_M_gen   = new TH1F("#Lambda_{c}^{+}",title2.c_str(),30, 2000,2.55*1000);
    lambda_M_gen->Sumw2();    
    
    Double_t masse[3]= {938, 493.7, 137.57};
    //TGenPhaseSpace event;
    //Double_t weight;

    //Double_t B= 1.4, L=3.4, D= 0.8, sigma_trk = 10E-6, N=4;
    double B= 1.1, L=1.7, D= 0.8, sigma_trk = 10E-6,  N=4;
    //Double_t sigma_p= 2*sigma_trk/0.3/B/L/D;
    Double_t sigma_p= 2*sigma_trk/0.3/B/L/L/1000;
    Double_t sigma_theta= sigma_trk/D*sqrt(2);
    //Double_t sigma_theta= sigma_trk/D*sqrt(2);
    Double_t sigma;

    TLorentzVector *p = new TLorentzVector();
    TLorentzVector *K = new TLorentzVector();
    TLorentzVector *Pi = new TLorentzVector();

    gRandom->SetSeed(1);

    double theta_x, theta_y, p_smear, px_smear, py_smear;

    
        for (Int_t i=0; i<nentries; i++) {
            tree->GetEntry(i);
      
            p->SetPxPyPzE(px_p,py_p,pz_p,E_p);   // tutto in MeV
            K->SetPxPyPzE(px_k,py_k,pz_k,E_k);
            Pi->SetPxPyPzE(px_pi,py_pi,pz_pi,E_pi);
        
            TLorentzVector p_gen = *p;
            TLorentzVector K_gen = *K;
            TLorentzVector pi_gen = *Pi;
            
            p_smear = gRandom->Gaus(p->P(), sigma_p * p->P()*p->P());
            px_smear = px_p + gRandom->Gaus(0, sigma_theta)*pz_p;
            py_smear = py_p + gRandom->Gaus(0, sigma_theta)*pz_p;
            p_gen.SetPxPyPzE(px_smear, py_smear, pz_p*p_smear/p->P(), sqrt(px_smear*px_smear + py_smear*py_smear + pz_p*p_smear/p->P()*pz_p*p_smear/p->P() + masse[0]*masse[0]));
            
            p_smear = gRandom->Gaus(K->P(), sigma_p * K->P()*K->P());
            px_smear = px_k + gRandom->Gaus(0, sigma_theta)*pz_k;
            py_smear = py_k + gRandom->Gaus(0, sigma_theta)*pz_k;
            K_gen.SetPxPyPzE(px_smear, py_smear, pz_k*p_smear/K->P(), sqrt(px_smear*px_smear + py_smear*py_smear + pz_k*p_smear/K->P()*pz_k*p_smear/K->P() + masse[1]*masse[1]));
            
            p_smear = gRandom->Gaus(Pi->P(), sigma_p * Pi->P()*Pi->P());
            px_smear = px_pi + gRandom->Gaus(0, sigma_theta)*pz_pi;
            py_smear = py_pi + gRandom->Gaus(0, sigma_theta)*pz_pi;
            pi_gen.SetPxPyPzE(px_smear, py_smear, pz_pi*p_smear/Pi->P(), sqrt(px_smear*px_smear + py_smear*py_smear + pz_pi*p_smear/Pi->P()*pz_pi*p_smear/Pi->P() + masse[2]*masse[2]));
         
            /*
            p_gen.SetRho(gRandom->Gaus(p->P(), sigma_p*p->P()*p->P()));
            
            p_gen.SetE(sqrt(p_gen.P()*p_gen.P()+masse[0]*masse[0]));
         
            K_gen.SetRho(gRandom->Gaus(K->P(), sigma_p*K->P()*K->P()));
           
            K_gen.SetE(sqrt(K_gen.P()*K_gen.P()+masse[1]*masse[1]));
        
            pi_gen.SetRho(gRandom->Gaus(Pi->P(), sigma_p*Pi->P()*Pi->P()));
            
            pi_gen.SetE(sqrt(pi_gen.P()*pi_gen.P()+masse[2]*masse[2]));

            /*
            p_smear = gRandom->Gaus(p->P(), sigma_p * p->P()*p->P());
            p_gen.SetPxPyPzE(px_p*1E-3* p_smear / p->P(),py_p*1E-3* p_smear / p->P(),pz_p*1E-3* p_smear / p->P(),sqrt(p_smear*p_smear+masse[0]*masse[0]));

            p_smear = gRandom->Gaus(K->P(), sigma_p * K->P()*K->P());
            K_gen.SetPxPyPzE(px_k*1E-3* p_smear / K->P(), py_k*1E-3* p_smear / K->P(), pz_k*1E-3* p_smear / K->P(),sqrt(p_smear*p_smear+masse[1]*masse[1]));

            p_smear = gRandom->Gaus(Pi->P(), sigma_p * Pi->P()*Pi->P());
            pi_gen.SetPxPyPzE(px_pi*1E-3* p_smear / Pi->P(), py_pi*1E-3* p_smear / Pi->P(), pz_pi*1E-3* p_smear / Pi->P(), sqrt(p_smear*p_smear+masse[2]*masse[2]));
         
            

            /*
            p_gen.SetRho(gRandom->Gaus(p->P(), sigma_p*p->P()*p->P()));
            theta_x = atan( tan(p->Phi())* cos(p->Theta()) );
            theta_x = gRandom->Gaus(theta_x, sigma_theta);
            theta_y = atan( tan(p->Phi())* sin(p->Theta()) );
            theta_y = gRandom->Gaus(theta_y, sigma_theta);
            p_gen.SetTheta( atan( tan(theta_y) / tan(theta_x)) );
            p_gen.SetPhi( atan( tan(theta_x) / cos(p_gen.Theta())) );
            p_gen.SetE(sqrt(p_gen.P()*p_gen.P()+masse[0]*masse[0]));
         
            K_gen.SetRho(gRandom->Gaus(K->P(), sigma_p*K->P()*K->P()));
            theta_x = atan( tan(K->Phi())* cos(K->Theta()) );
            theta_x = gRandom->Gaus(theta_x, sigma_theta);
            theta_y = atan( tan(K->Phi())* sin(K->Theta()) );
            theta_y = gRandom->Gaus(theta_y, sigma_theta);
            K_gen.SetTheta( atan( tan(theta_y) / tan(theta_x)) );
            K_gen.SetPhi( atan( tan(theta_x) / cos(K_gen.Theta())) );
            K_gen.SetE(sqrt(K_gen.P()*K_gen.P()+masse[1]*masse[1]));
        
            pi_gen.SetRho(gRandom->Gaus(Pi->P(), sigma_p*Pi->P()*Pi->P()));
            theta_x = atan( tan(Pi->Phi())* cos(Pi->Theta()) );
            theta_x = gRandom->Gaus(theta_x, sigma_theta);
            theta_y = atan( tan(Pi->Phi())* sin(Pi->Theta()) );
            theta_y = gRandom->Gaus(theta_y, sigma_theta);
            pi_gen.SetTheta( atan( tan(theta_y) / tan(theta_x)) );
            pi_gen.SetPhi( atan( tan(theta_x) / cos(pi_gen.Theta())) );
            pi_gen.SetE(sqrt(pi_gen.P()*pi_gen.P()+masse[2]*masse[2]));
        */

            lambda_gen = p_gen + K_gen + pi_gen;
            lambda_M_gen->Fill(lambda_gen.M());
                
        }
    

    TCanvas *c1 =new TCanvas("c1","c1", 900, 550);


    
    
    lambda_M_gen->Fit("gaus");
    gStyle->SetOptFit(0101);
    lambda_M_gen->Draw();
    TF1 * fit = lambda_M_gen->GetFunction("gaus");
    sigma = fit->GetParameter(2);
    TLegend *l = new TLegend(0.2,0.7,0.35,0.89);
    l->AddEntry((TObject*)0, TString::Format("D = %g m" , D), "");
    l->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    l->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    l->AddEntry((TObject*)0, TString::Format("#sigma_{m} = %g MeV/c^{2}", sigma), "");
    l->Draw("same");
    
    cout << sigma << endl;
    
    
    
    //c1->SaveAs("invariant_mass.png");
    string filename = "incertezze_masse_old_tot_"+material+bending+"2_test.dat";
    ofstream out;
    out.open(filename.c_str(),ios::app);
    out << sigma << "\t" << D << "\t" << B << "\t" << L  << endl;
   
}