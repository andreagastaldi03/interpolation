void plot_mu() {
    TCanvas *c_mu = new TCanvas("c_mu", "Coefficienti di Attenuazione NIST", 800, 600);
    
    // Attiviamo la scala logaritmica su entrambi gli assi (fondamentale per mu!)
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGrid();

    // ==========================================
    // 1. LETTURA DATI
    // ==========================================
    // File 1: Attenuazione Totale (mu)
    TGraph *g_mu_exact  = new TGraph("confronto_mu.dat", "%lg %lg %*lg");
    TGraph *g_mu_interp = new TGraph("confronto_mu.dat", "%lg %*lg %lg");

    // File 2: Assorbimento di Energia (mu_en)
    TGraph *g_mu_en_exact  = new TGraph("confronto_mu_en.dat", "%lg %lg %*lg");
    TGraph *g_mu_en_interp = new TGraph("confronto_mu_en.dat", "%lg %*lg %lg");

    // ==========================================
    // 2. STYLING DELLE CURVE
    // ==========================================
    // MU TOTALE (Blu)
    g_mu_exact->SetLineColor(kBlue+2);
    g_mu_exact->SetLineWidth(5);       // Linea spessa di background per il dato esatto

    g_mu_interp->SetLineColor(kCyan);
    g_mu_interp->SetLineStyle(7);      // Tratteggiata per mostrare l'interpolazione che vi si adagia
    g_mu_interp->SetLineWidth(2);

    // MU ENERGIA (Rosso/Arancio)
    g_mu_en_exact->SetLineColor(kRed+2);
    g_mu_en_exact->SetLineWidth(5);

    g_mu_en_interp->SetLineColor(kOrange+1);
    g_mu_en_interp->SetLineStyle(7);
    g_mu_en_interp->SetLineWidth(2);

    // ==========================================
    // 3. MULTIGRAPH E PLOT
    // ==========================================
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Confronto Dati NIST vs Interpolazione Chebyshev;Energia (MeV);Coefficiente Massico (cm^{2}/g)");
    
    mg->Add(g_mu_exact, "L");
    mg->Add(g_mu_interp, "L");
    mg->Add(g_mu_en_exact, "L");
    mg->Add(g_mu_en_interp, "L");

    mg->Draw("A");

    // ==========================================
    // 4. LEGENDA
    // ==========================================
    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->AddEntry(g_mu_exact,  "\\mu/\\rho Dato Esatto NIST", "l");
    leg->AddEntry(g_mu_interp, "\\mu/\\rho Interpolato", "l");
    leg->AddEntry(g_mu_en_exact,  "\\mu_{en}/\\rho Dato Esatto", "l");
    leg->AddEntry(g_mu_en_interp, "\\mu_{en}/\\rho Interpolato", "l");
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
    leg->SetBorderSize(1);
    leg->Draw();
}
