void plot_mu() {
    TCanvas *c_mu = new TCanvas("c_mu", "Interpolazione mu e mu_en", 800, 600);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGrid();

    // Lettura Dati
    TGraph *g_mu_exact  = new TGraph("confronto_mu.dat", "%lg %lg %*lg");
    TGraph *g_mu_interp = new TGraph("confronto_mu.dat", "%lg %*lg %lg");
    TGraph *g_en_exact  = new TGraph("confronto_mu_en.dat", "%lg %lg %*lg");
    TGraph *g_en_interp = new TGraph("confronto_mu_en.dat", "%lg %*lg %lg");

    // Styling mu Totale (Blu)
    g_mu_exact->SetLineColor(kBlue+2);   g_mu_exact->SetLineWidth(5);
    g_mu_interp->SetLineColor(kCyan);    g_mu_interp->SetLineWidth(2); g_mu_interp->SetLineStyle(7);

    // Styling mu_en (Rosso/Arancio)
    g_en_exact->SetLineColor(kRed+2);    g_en_exact->SetLineWidth(5);
    g_en_interp->SetLineColor(kOrange);  g_en_interp->SetLineWidth(2); g_en_interp->SetLineStyle(7);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Coefficienti NIST e Interpolazione Chebyshev;Energia (MeV);Coefficiente (cm^{2}/g)");
    mg->Add(g_mu_exact, "L"); mg->Add(g_mu_interp, "L");
    mg->Add(g_en_exact, "L"); mg->Add(g_en_interp, "L");
    mg->Draw("A");

    TLegend *leg = new TLegend(0.50, 0.65, 0.88, 0.88);
    leg->AddEntry(g_mu_exact, "\\mu/\\rho Totale (NIST)", "l");
    leg->AddEntry(g_mu_interp, "\\mu/\\rho Totale (Interp)", "l");
    leg->AddEntry(g_en_exact, "\\mu_{en}/\\rho Assorbimento (NIST)", "l");
    leg->AddEntry(g_en_interp, "\\mu_{en}/\\rho Assorbimento (Interp)", "l");
    leg->Draw();
}
