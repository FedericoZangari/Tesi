void plot6(){
    string title1 = "incertezze_masse_old_theta_Si7mrad.dat";    
    string title2 = "incertezze_masse_old_p_Si7mrad.dat";
    string title3 = "incertezze_masse_old_tot_Si7mrad.dat";
    //string title1 = "incertezze_masse_old_theta_Si7mrad2.dat";    
    //string title2 = "incertezze_masse_old_p_Si7mrad2.dat";
    //string title3 = "incertezze_masse_old_tot_Si7mrad2.dat";
    string material = title1.substr (27,2);
    string bending = title1.substr (29,5);

    ifstream f1, f2, f3;
    f1.open(title1.c_str(),ios::in);
    f2.open(title2.c_str(),ios::in);
    f3.open(title3.c_str(),ios::in);
    Double_t D, B, L;
    Double_t sigma;


    TGraph *c1= new TGraph();
    TGraph *c2= new TGraph();
    TGraph *c3= new TGraph();
    
    
    int i= 0,j=0, k=0;
    for (;;){
        f1 >> sigma >> D >> B >> L;
        c1->SetPoint(i, D, sigma);i++;

        f2 >> sigma >> D >> B >> L;
        c2->SetPoint(j, D, sigma);j++;

        f3 >> sigma >> D >> B >> L;
        c3->SetPoint(k, D, sigma);k++;
        
        if(f1.eof()){
			cout << "End of file reached "<< endl;
			break;
		}		
	}
    
    
    TCanvas *c =new TCanvas("c1","c1");
    c1->Draw("APL");
    c2->Draw("samePL");
    c3->Draw("samePL");

    c1->SetMarkerStyle(43);
	c1->SetMarkerSize(1.5);
    c1->SetMarkerColor(kOrange);
    c1->SetLineColor(kOrange);

    c2->SetMarkerStyle(41);
	c2->SetMarkerSize(1.5);
    c2->SetMarkerColor(kBlue);
    c2->SetLineColor(kBlue);
    
    c3->SetLineColor(kGreen+3);

    string title = "Invariant mass uncertainty contributions " + material + " " + bending;
    c1->SetTitle(title.c_str());
    c1->GetXaxis()->SetTitle("D [m]");
	c1->GetYaxis()->SetTitle("#sigma_{M} [MeV]");    
    c1->GetYaxis()->SetRangeUser(9,80);

    TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
    leg->AddEntry(c1, "only angles contribution", "lp");
    leg->AddEntry(c2, "only momentum contribution", "lp");
    leg->AddEntry(c3, "total", "lp");
    leg->AddEntry((TObject*)0, TString::Format("B = %g T", B), "");
    leg->AddEntry((TObject*)0, TString::Format("L = %g m", L), "");
    leg->AddEntry((TObject*)0, "#sigma_{p}#propto 1/LD", "");
    //leg->AddEntry((TObject*)0, "#sigma_{p}#propto 1/L^{2}", "");
    leg->Draw();

    c->SaveAs("incertezze_masse_old_tot_Si7mrad.png");
    //c->SaveAs("incertezze_masse_old_tot_Si7mrad2.png");
}