void plot(){
    TCanvas * plot = new TCanvas("plot","plot", 1100,600);
    plot->Divide(2,1);
    plot->cd(1);
    double B=1.1, L=1.7, p=600, sig =1E-5;
    TF1 *f = new TF1("sigma","2*[0]/(0.3*x)",0.3, 0.7);
    //f->SetParameter(0,p*p*sig/L/B);
    f->SetParameter(0,p*sig/L/B*100);
    f->SetTitle("Momentum uncertainty as a function of D");
    f->GetXaxis()->SetTitle("D[m]");
    f->GetYaxis()->SetTitle("#sigma_{p}/p % [GeV/c]");
    TF1 *g = new TF1("sigma2","2*[1]/(0.3)",0.3, 0.7);
    g->SetParameter(1,p*sig/L/L/B*100);
    g->SetLineColor(kSpring);
    
    TLegend *l = new TLegend(0.5,1,"","brNDC");
    l->AddEntry((TObject*)0, "p=600GeV/c", "");
    l->AddEntry((TObject*)0, "B=1.1T", "");
    l->AddEntry((TObject*)0, "L=1.7m", "");
    l->AddEntry(f, "#sigma_{p}#propto 1/LD", "l");
    l->AddEntry(g, "#sigma_{p}#propto 1/L^{2}", "l");
    
    f->Draw();
    g->Draw("SAME");
    l->Draw("same");
    f->GetYaxis()->SetRangeUser(0.,10);
    
    plot->cd(2);
    B=1.4; L=3.4; 
    TF1 *f2 = new TF1("sigma","2*[0]/(0.3*x)",0.3, 0.7);
    f2->SetParameter(0,p*sig/L/B*100);
    f2->SetTitle("Momentum uncertainty as a function of D");
    f2->GetXaxis()->SetTitle("D[m]");
    f2->GetYaxis()->SetTitle("#sigma_{p}/p % [GeV/c]");
    TF1 *g2 = new TF1("sigma2","2*[1]/(0.3)",0.3, 0.7);
    g2->SetParameter(1,p*sig/L/L/B*100);
    cout << p*p*sig/L/L/B*2/0.3 << endl;
    g2->SetLineColor(kSpring);
    TLegend *l2 = new TLegend(0.2,1,"","brNDC");
    l2->AddEntry((TObject*)0, "p=600GeV/c", "");
    l2->AddEntry((TObject*)0, "B=1.4T", "");
    l2->AddEntry((TObject*)0, "L=3.4m", "");
    l2->AddEntry(f2, "#sigma_{p}#propto 1/LD", "l");
    l2->AddEntry(g2, "#sigma_{p}#propto 1/L^{2}", "l");
    f2->GetYaxis()->SetRangeUser(0.,10);
    f2->Draw();
    l2->Draw("same");
    g2->Draw("SAME");

    plot->SaveAs("incertezza p.png");

    

}