void plot_ode() {
    TCanvas *c2 = new TCanvas("c2", "Propagazione ODE e Rilascio Energia", 800, 800);
    c2->Divide(1, 2);

    // ==========================================
    // PAD 1: Validazione Fluenza (Colonna 2 e 3)
    // ==========================================
    c2->cd(1);
    gPad->SetGrid();
    
    // Leggiamo z (1) e Fluenza Esatta (3)
    TGraph *g_phi_exact = new TGraph("fluenza_terma_z.dat", "%lg %*lg %lg");
    // Leggiamo z (1) e Fluenza RK4 (2)
    TGraph *g_phi_rk4   = new TGraph("fluenza_terma_z.dat", "%lg %lg %*lg");

    g_phi_exact->SetTitle("Decadimento Fluenza Primaria;Profondita' z (cm);\\Phi(z) / \\Phi_{0}");
    g_phi_exact->SetLineColor(kBlue);
    g_phi_exact->SetLineWidth(4);

    g_phi_rk4->SetLineColor(kRed);
    g_phi_rk4->SetLineStyle(7);
    g_phi_rk4->SetLineWidth(3);

    g_phi_exact->Draw("AL");
    g_phi_rk4->Draw("L SAME");

    TLegend *leg1 = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg1->AddEntry(g_phi_exact, "Fluenza Analitica Esatta", "l");
    leg1->AddEntry(g_phi_rk4, "Fluenza Numerica (RK4)", "l");
    leg1->SetTextSize(0.045);
    leg1->SetBorderSize(1);
    leg1->Draw();

    // ==========================================
    // PAD 2: TERMA e KERMA calcolati dall'ODE (Colonna 4 e 5)
    // ==========================================
    c2->cd(2);
    gPad->SetGrid();
    
    // Leggiamo z (1) e TERMA (4)
    TGraph *g_terma = new TGraph("fluenza_terma_z.dat", "%lg %*lg %*lg %lg");
    // Leggiamo z (1) e KERMA (5)
    TGraph *g_kerma = new TGraph("fluenza_terma_z.dat", "%lg %*lg %*lg %*lg %lg");

    g_terma->SetLineColor(kBlue);
    g_terma->SetLineWidth(3);
    g_terma->SetLineStyle(7);

    g_kerma->SetLineColor(kGreen+2);
    g_kerma->SetLineWidth(3);
    g_kerma->SetLineStyle(7);

    // Usiamo TMultiGraph per auto-scalare perfettamente l'asse Y
    TMultiGraph *mg2 = new TMultiGraph();
    mg2->SetTitle("Rilascio e Assorbimento Energetico;Profondita' z (cm);Energia (a.u.)");
    mg2->Add(g_terma, "L");
    mg2->Add(g_kerma, "L");
    mg2->Draw("A");

    TLegend *leg2 = new TLegend(0.50, 0.70, 0.88, 0.88);
    leg2->AddEntry(g_terma, "TERMA", "l");
    leg2->AddEntry(g_kerma, "Collision KERMA", "l");
    leg2->SetTextSize(0.045);
    leg2->SetBorderSize(1);
    leg2->Draw();
}
