void plot_ode() {
    TCanvas *c2 = new TCanvas("c2", "Trasporto Primario Fotoni", 1000, 800);
    c2->Divide(2,1);

    TGraph *g_rk4 = new TGraph("fluenza_terma_z.dat", "%lg %lg");
    TGraph *g_ana = new TGraph("fluenza_terma_z.dat", "%lg %*lg %lg");
    TGraph *g_ter = new TGraph("fluenza_terma_z.dat", "%lg %*lg %*lg %lg");

    c2->cd(1);
    gPad->SetGrid();

    g_ana->SetTitle("Decadimento Fluenza Primaria;Profondita' z (cm);\\Phi(z) / \\Phi_{0}");
    g_ana->SetLineColor(kBlue);
    g_ana->SetLineWidth(4);
    
    g_rk4->SetLineColor(kRed);
    g_rk4->SetLineStyle(7);
    g_rk4->SetLineWidth(3);

    g_ana->Draw("AL");
    g_rk4->Draw("L SAME");

    TLegend *leg = new TLegend(0.50, 0.70, 0.88, 0.88);
    leg->AddEntry(g_ana, "Analitica", "l");
    leg->AddEntry(g_rk4, "Numerica (RK4)", "l");
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetBorderSize(1);
    leg->Draw();

    c2->cd(2);
    gPad->SetGrid();

    g_ter->SetTitle("Calcolo TERMA;Profondita' z (cm);TERMA (MeV)");
    g_ter->SetLineColor(kGreen+2);
    g_ter->SetLineWidth(3);
    g_ter->Draw("AL");
}
