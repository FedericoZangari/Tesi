void Rate(){
    TFile* file = new TFile("hist_mb_1000.root");
    //TFile* file = new TFile("hist_Lc_1000.root");
    int trk_id = 2;
    int trk_nmr = 1;
    int Event_number = 1000;
    string t = "Occupancy"+ to_string(trk_id) + "_" + to_string(trk_nmr) + "_" + to_string(Event_number);
    TH2F * rate =  (TH2F*)file->Get(t.c_str());
    TCanvas *c = new TCanvas("c", "c", 850, 600);
    c->SetRightMargin(0.2);
    

    //rate->Rebin2D( 10./0.055, 10./0.055);
    //rate->Scale(1./ 0.000025);
    rate->Scale(1000/100);      //in occupancy salvo il numero di volte il pix è trigged *100 / num_event
    double n_pix = 10;
    rate->Rebin2D(n_pix,n_pix);
    rate->Scale(1.E6*0.1822);   //scalo per il rate di eventi di minimum bias
    rate->Scale(1./1000.);   //scalo per il numero di eventi studiati
    rate->Scale(1/(n_pix*n_pix*5.5*5.5E-6));      // scalo per l'area di 9 pixel espressa in cm^2
    rate->Scale(1/1000000.);                    //  esprimo in MHz
    rate->GetXaxis()->SetTitle("x [mm]");
    rate->GetYaxis()->SetTitle("y [mm]");
    rate->GetZaxis()->SetTitleOffset(1.5);
    rate->GetZaxis()->SetTitle("rate of hits [MHz/cm^{2}]");
    string title = "Rate of events, tracker " + to_string(trk_id) + "_" + to_string(trk_nmr) + ", " + to_string(Event_number) +" background events";
    //string title = "Rate of events, tracker " + to_string(trk_id) + "_" + to_string(trk_nmr) + ", " + to_string(Event_number) +" #Lambda_{c}^{+} events";
    rate->SetTitle(title.c_str());
    rate->GetZaxis()->SetRangeUser(0, 3.700000);

    //rate->GetXaxis()->SetRangeUser(-25,25);
    //rate->GetYaxis()->SetRangeUser(-25,25);
    gStyle->SetPalette(kRainBow);
    rate->SetStats(0);
    rate->Draw("colz");
    
    //t = "rate_mb_" + to_string(trk_id) + "_" + to_string(trk_nmr) + "_" + to_string(Event_number) + "_zoom.png";
    t = "rate_mb_" + to_string(trk_id) + "_" + to_string(trk_nmr) + "_" + to_string(Event_number) + ".png";
    c->SaveAs(t.c_str());
    
}
