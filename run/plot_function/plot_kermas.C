void plot_kermas() {
    TCanvas *c_kerma = new TCanvas("c_kerma", "Analisi KERMA e ODE", 800, 600);
    c_kerma->SetGrid();

    // 1 ODE (Colonna 5)
    TGraph *g_kerma_1ode = new TGraph("fluenza_terma_z.dat", "%lg %*lg %*lg %*lg %lg");
    // 2 ODE: p (2), s (3), tot (4)
    TGraph *g_kerma_p   = new TGraph("fluenza_prim_scat.dat", "%lg %lg %*lg %*lg %*lg %*lg");
    TGraph *g_kerma_s   = new TGraph("fluenza_prim_scat.dat", "%lg %*lg %lg %*lg %*lg %*lg");
    TGraph *g_kerma_tot = new TGraph("fluenza_prim_scat.dat", "%lg %*lg %*lg %lg %*lg %*lg");

    // Styling
    g_kerma_1ode->SetLineColor(kBlack);   g_kerma_1ode->SetLineWidth(2); g_kerma_1ode->SetLineStyle(7);
    g_kerma_p->SetLineColor(kBlue);       g_kerma_p->SetLineWidth(3);
    g_kerma_s->SetLineColor(kOrange+1);   g_kerma_s->SetLineWidth(3);
    g_kerma_tot->SetLineColor(kRed+1);    g_kerma_tot->SetLineWidth(4);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Scomposizione del KERMA (Dinamica Fotoni);Profondita' z (cm);Energia (a.u.)");
    mg->Add(g_kerma_1ode, "L"); mg->Add(g_kerma_p, "L");
    mg->Add(g_kerma_s, "L");    mg->Add(g_kerma_tot, "L");

    mg->Draw("A");
    mg->GetXaxis()->SetRangeUser(0.0, 15.0);

    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->AddEntry(g_kerma_1ode, "KERMA (Modello senza scatter)", "l");
    leg->AddEntry(g_kerma_p, "KERMA Primario (K_{p})", "l");
    leg->AddEntry(g_kerma_s, "KERMA Scatter (K_{s})", "l");
    leg->AddEntry(g_kerma_tot, "KERMA Totale (K_{tot} = K_{p} + K_{s})", "l");
    leg->Draw();
}
