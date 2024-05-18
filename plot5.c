void plot5(){
    //string title = "incertezze_masse_track_Ge5mrad";
    //string title = "incertezze_masse_track_Ge7mrad";
    //string title = "incertezze_masse_track_Si5mrad";
    //string title = "incertezze_masse_track_Si7mrad";
    string title = "incertezze_masse_theta_contrib_Si7mrad";    
    string title11 = "incertezze_masse_p_contrib_Si7mrad.dat";
    string title12 = "incertezze_masse_track_Si7mrad.dat";
    string material = title.substr (31,2);
    string bending = title.substr (33,5);
    string title1 = title + ".dat";

    ifstream f;
    f.open(title1.c_str(),ios::in);
    ifstream f1;
    f1.open(title11.c_str(),ios::in);
    ifstream f2;
    f2.open(title12.c_str(),ios::in);
    Double_t D, B, L, chi;
    Double_t sigma;
    TGraph *c1= new TGraph();
    TGraph *c2= new TGraph();
    TGraph *t1= new TGraph();
    
    
    int i= 0,j=0, k=0;
    for (;;){
        f >> sigma >> D >> B >> L >> chi;
        c1->SetPoint(i, D, sigma);i++;
        f1 >> sigma >> D >> B >> L >> chi;
        c2->SetPoint(j, D, sigma);j++;
        f2 >> sigma >> D >> B >> L >> chi;
        t1->SetPoint(k, D, sigma);k++;
        
        if(f.eof()){
			cout << "End of file reached "<< endl;
			break;
		}		
	}
    
    
    TCanvas *c =new TCanvas("c1","c1");
    c1->Draw("APL");
    c2->Draw("samePL");
    t1->Draw("samePL");
    c1->SetMarkerStyle(43);
	c1->SetMarkerSize(1.5);
    c2->SetMarkerStyle(41);
	c2->SetMarkerSize(1.5);
    c2->SetMarkerColor(kBlue);
    c1->SetMarkerColor(kOrange);
    string title2 = "Invariant mass uncertainty from tracks " + material + " " + bending;
    c1->SetTitle(title2.c_str());
    c1->GetXaxis()->SetTitle("D [cm]");
	c1->GetYaxis()->SetTitle("#sigma_{M} [MeV]");
    c1->SetLineColor(kOrange);
    c2->SetLineColor(kBlue);
    t1->SetLineColor(kGreen+3);

    
    c1->GetYaxis()->SetRangeUser(9,85);

    TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
    leg->AddEntry(c1, "only angles contribution", "lp");
    leg->AddEntry(c2, "only momentum contribution", "lp");
    leg->AddEntry(t1, "total", "lp");
    leg->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    leg->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    leg->AddEntry((TObject*)0, TString::Format("R_{B} = %g cm", 4.), "");
    leg->Draw();

    c->SaveAs("incertezze_masse_track_contributions_Si7mrad.png");
}