void Accettanza(){

    for(int z=0; z<4; z++){

        string title;
        if (z==0) title = "Lc_thCy0_T2cm_Ge5mrad7cm.root";
        if (z==1) title = "Lc_thCy0_T2cm_Ge7mrad7cm.root";
        if (z==2) title = "Lc_thCy0_T2cm_Si5mrad7cm.root";
        if (z==3) title = "Lc_thCy0_T2cm_Si7mrad7cm.root";
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

        Double_t B= 1.4, L=3.4;
        vector<Double_t> counts;
        Double_t y1= 0.5;                               // tutte le dispersioni in y sono in cm
        Double_t y2=0, y3=0, y4=0, y5=0;
       
        for(int k = 0; k < 2; k++){
            counts = {0.,0.,0.,0.};
            y2= theta_bend *1.E-3*(L+1.)*1.E2;

            for(int j=0; j < 500; j++){    
                for (Int_t i=0; i<nentries; i++) {
                    tree->GetEntry(i);
      
                    lambda.SetPxPyPzE(px_lambda*1E-3,py_lambda*1E-3,pz_lambda*1E-3,E_lambda*1E-3);   // esprimo tutto in GeV
      
                    event.SetDecay(lambda, 3, masse);
                    weight = event.Generate();

                    TLorentzVector *p = event.GetDecay(0);
                    TLorentzVector *K    = event.GetDecay(1);
                    TLorentzVector *Pi    = event.GetDecay(2);

                    y3 = 0.3*L*L*B/2./p->P()*1.E2;
                    y4 = (p->Theta())*(L+1.)*1.E2;
                    y5= gRandom->Gaus(0,0.05);
                    if(y1+y2+y3+y4+y5 < 2) counts[0]++;
                    if(y1+y2+y3+y4+y5 < 2.5) counts[1]++;
                    if(y1+y2+y3+y4+y5 < 3) counts[2]++;
                    if(y1+y2+y3+y4+y5 < 4) counts[3]++;

                    y3 = 0.3*L*L*B/2./K->P()*1.E2;
                    y4 =(K->Theta())*(L+1.)*1.E2;
                    y5= gRandom->Gaus(0,0.05);
                    if(y1+y2+y3+y4+y5 < 2) counts[0]++;
                    if(y1+y2+y3+y4+y5 < 2.5) counts[1]++;
                    if(y1+y2+y3+y4+y5 < 3) counts[2]++;
                    if(y1+y2+y3+y4+y5 < 4) counts[3]++;

                    y3 = 0.3*L*L*B/2./Pi->P()*1.E2;
                    y4 = (Pi->Theta())*(L+1.)*1.E2;
                    y5= gRandom->Gaus(0,0.05);
                    if(y1+y2+y3+y4+y5 < 2) counts[0]++;
                    if(y1+y2+y3+y4+y5 < 2.5) counts[1]++;
                    if(y1+y2+y3+y4+y5 < 3) counts[2]++;
                    if(y1+y2+y3+y4+y5 < 4) counts[3]++;
               
                }
            }
            cout << "B= " << B << "T \t L= " << L << " m\t";
            cout << setprecision(2) << counts[0]/nentries/3./500.*100. <<  "\t" << counts[1]/nentries/3./500.*100. << "\t" << counts[2]/nentries/3./500.*100. << "\t" << counts[3]/nentries/3./500.*100. <<endl;
            B= 1.1; L=1.7;

        }
    }
}