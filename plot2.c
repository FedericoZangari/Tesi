void plot2(){
    TCanvas * plot = new TCanvas("plot","plot", 1100,600);
    plot->Divide(2,1);
    plot->cd(1);
    double B=1.1, L=1.7, D=0.4, sig =1E-5;
    TF1 *f = new TF1("sigma","2*[0]/0.3*x",100, 1000);
    f->SetParameter(0,sig/L/B/D*100);
    f->SetTitle("Momentum uncertainty as a function of p");
    f->GetXaxis()->SetTitle("p[GeV/c]");
    f->GetYaxis()->SetTitle("#sigma_{p}/p % [GeV/c]");
    TF1 *g = new TF1("sigma2","2*[1]*x/(0.3)",100, 1000);
    g->SetParameter(1,sig/L/L/B*100);
    g->SetLineColor(kSpring);
    
    TLegend *l = new TLegend(0.2,0.7,0.45,0.89);
    l->AddEntry((TObject*)0, "D=0.4m", "");
    l->AddEntry((TObject*)0, "B=1.1T", "");
    l->AddEntry((TObject*)0, "L=1.7m", "");
    l->AddEntry(f, "#sigma_{p}#propto 1/LD", "l");
    l->AddEntry(g, "#sigma_{p}#propto 1/L^{2}", "l");
    
    f->Draw();
    g->Draw("SAME");
    l->Draw("same");
    f->GetYaxis()->SetRangeUser(0,9.);
    
    plot->cd(2);
    B=1.4; L=3.4; 
    TF1 *f2 = new TF1("sigma","2*[0]/0.3*x",100, 1000);
    f2->SetParameter(0,sig/L/D/B*100);
    f2->SetTitle("Momentum uncertainty as a function of p");
    f2->GetXaxis()->SetTitle("p[GeV/c]");
    f2->GetYaxis()->SetTitle("#sigma_{p}/p % [GeV/c]");
    TF1 *g2 = new TF1("sigma2","2*[0]/0.3*x",100, 1000);
    g2->SetParameter(0,sig/L/L/B*100);
    g2->SetLineColor(kSpring);
    TLegend *l2 = new TLegend(0.2,0.7,0.45,0.89);
    l2->AddEntry((TObject*)0, "D=0.4m", "");
    l2->AddEntry((TObject*)0, "B=1.4T", "");
    l2->AddEntry((TObject*)0, "L=3.4m", "");
    l2->AddEntry(f2, "#sigma_{p}#propto 1/LD", "l");
    l2->AddEntry(g2, "#sigma_{p}#propto 1/L^{2}", "l");
    f2->GetYaxis()->SetRangeUser(0.,9);
    f2->Draw();
    l2->Draw("same");
    g2->Draw("SAME");
    

    plot->SaveAs("incertezza2 p.png");

    

}