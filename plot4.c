void plot4(){
    string title = "incertezze_masse_track_Ge5mrad";
    //string title = "incertezze_masse_track_Ge7mrad";
    //string title = "incertezze_masse_track_Si5mrad";
    //string title = "incertezze_masse_track_Si7mrad";
    string material = title.substr (23,2);
    string bending = title.substr (25,5);
    string title1 = title + ".dat";
    ifstream f;
    f.open(title1.c_str(),ios::in);
    Double_t D, B, L, chi;
    Double_t sigma;
    TGraph *f1= new TGraph();
    TGraph *g1= new TGraph();
    
    
    int i= 0,j=0;
    for (;;){
        f >> sigma >> D >> B >> L >> chi;
        if( B==1.4) {f1->SetPoint(i, D, sigma);i++;}
        if( B==1.1) {g1->SetPoint(j, D, sigma);j++;}
        
        if(f.eof()){
			cout << "End of file reached "<< endl;
			break;
		}		
	}
    
    
    TCanvas *c1 =new TCanvas("c1","c1");
    f1->Draw("APL");
    g1->Draw("samePL");
    f1->SetMarkerStyle(43);
	f1->SetMarkerSize(1.5);
    g1->SetMarkerStyle(41);
	g1->SetMarkerSize(1.5);
    f1->SetMarkerColor(kRed);
    g1->SetMarkerColor(kGreen+3);
    string title2 = "Invariant mass uncertainty from tracks " + material + " " + bending;
    f1->SetTitle(title2.c_str());
    f1->GetXaxis()->SetTitle("D [cm]");
	f1->GetYaxis()->SetTitle("#sigma_{M} [MeV]");
    f1->SetLineColor(kRed+2);
    g1->SetLineColor(kGreen+3);

    
    f1->GetYaxis()->SetRangeUser(9,100);

    TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
    leg->AddEntry(f1, "B=1.4T, L=3.4m", "lp");
    leg->AddEntry(g1, "B=1.1T, L=1.7m", "lp");
    leg->Draw();

    title = title +".png";
    c1->SaveAs(title.c_str());
}