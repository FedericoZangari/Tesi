void plot3(){
    string title = "incertezze_masse_Ge5mrad";
    //string title = "incertezze_masse_Ge7mrad";
    //string title = "incertezze_masse_Si5mrad";
    //string title = "incertezze_masse_Si7mrad";
    string material = title.substr (17,2);
    string bending = title.substr (19,5);
    string title1 = title + ".dat";
    ifstream f;
    f.open(title1.c_str(),ios::in);
    Double_t D, B, L;
    int index;
    Double_t sigma;
    TGraph *f1= new TGraph();
    TGraph *g1= new TGraph();
    TGraph *f2= new TGraph();
    TGraph *g2= new TGraph();
    
    int i= 0,j=0,k=0,l=0;
    for (;;){
        f >> sigma >> D >> B >> L >> index;
        if(index ==1 && B==1.4) {f1->SetPoint(i, D, sigma);i++;}
        if(index ==2 && B==1.4) {g1->SetPoint(j, D, sigma);j++;}
        if(index ==1 && B==1.1) {f2->SetPoint(k, D, sigma);k++;}
        if(index ==2 && B==1.1) {g2->SetPoint(l, D, sigma);l++;}
        
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
    g1->SetMarkerStyle(43);
	g1->SetMarkerSize(1.5);
    f1->SetMarkerColor(kRed);
    g1->SetMarkerColor(kGreen+3);
    string title2 = "Invariant mass uncertainty " + material + " " + bending;
    f1->SetTitle(title2.c_str());
    f1->GetXaxis()->SetTitle("D [m]");
	f1->GetYaxis()->SetTitle("#sigma_{M} [MeV]");
    f1->SetLineColor(kRed+2);
    g1->SetLineColor(kGreen+3);

    f2->Draw("samePL");
    g2->Draw("samePL");
    f2->SetMarkerStyle(41);
	f2->SetMarkerSize(1.5);
    g2->SetMarkerStyle(41);
	g2->SetMarkerSize(1.5);
    f2->SetMarkerColor(kRed);
    g2->SetMarkerColor(kGreen+3);
    f2->SetLineColor(kRed+2);
    g2->SetLineColor(kSpring-7);
    f1->GetYaxis()->SetRangeUser(10,105);

    TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
    leg->AddEntry(f1, "B=1.4T, L=3.4m, #sigma_{p}#propto 1/LD", "lp");
    leg->AddEntry(g1, "B=1.4T, L=3.4m, #sigma_{p}#propto 1/L^{2}", "lp");
    leg->AddEntry(f2, "B=1.1T, L=1.7m, #sigma_{p}#propto 1/LD", "lp");
    leg->AddEntry(g2, "B=1.1T, L=1.7m, #sigma_{p}#propto 1/L^{2}", "lp");
    leg->Draw();

    title = title +".png";
    c1->SaveAs(title.c_str());
}