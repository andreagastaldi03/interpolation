void plot_confronto_dosi_smallStep() {
    TCanvas *c_dosi = new TCanvas("c_dosi_small", "Dose: Riemann vs Gauss-Kronrod", 1000, 800);
    c_dosi->SetGrid();

    // Da profilo_dose.dat (Riemann) -> z (1), Dose (3)
    TGraph *g_dose_riemann = new TGraph("profilo_dose_small.dat", "%lg %*lg %lg");
    // Da profilo_dose_gk.dat -> z (1), TERMA (2), Dose_K_Coll (5), Kerma_Tot (6)
    TGraph *g_terma        = new TGraph("profilo_dose_gk_small.dat", "%lg %lg %*lg %*lg %*lg %*lg %*lg");
    TGraph *g_kerma        = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %lg %*lg %*lg %*lg %*lg");
    TGraph *g_kerma_tot    = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %lg %*lg");
    TGraph *g_dose_k_coll  = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %lg %*lg %*lg");
    TGraph *g_dose_k_tot   = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %*lg %lg");
    TGraph *g_dose_t_tot   = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    TGraph *g_dose_t       = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %lg %*lg %*lg %*lg");
    TGraph *g_dose_k_tot_e = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    // Styling (Curve sorgente in background)
    g_terma->SetLineColor(kAzure-3); g_terma->SetLineWidth(4); g_terma->SetLineStyle(3);
    g_kerma_tot->SetLineColor(kGreen+2); g_kerma_tot->SetLineWidth(4); g_kerma_tot->SetLineStyle(3);
    g_kerma->SetLineColor(kP10Green); g_kerma->SetLineWidth(4); g_kerma->SetLineStyle(3);
   
    // Curve calcolate
    g_dose_riemann->SetLineColor(kMagenta); g_dose_riemann->SetLineWidth(2); g_dose_riemann->SetLineStyle(9);
    g_dose_k_coll->SetLineColor(kRed); g_dose_k_coll->SetLineWidth(2);
    g_dose_t_tot->SetLineColor(kP8Orange); g_dose_t_tot->SetLineWidth(2);
    g_dose_k_tot->SetLineColor(kP6Grape); g_dose_k_tot->SetLineWidth(2);
    g_dose_t->SetLineColor(kP6Gray); g_dose_t->SetLineWidth(2);
    g_dose_k_tot_e->SetLineColor(kP10Cyan); g_dose_k_tot_e->SetLineWidth(2);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Confronto Integratori (Riemann vs Gauss-Kronrod);Profondita' z (cm);Dose (a.u.)");
    mg->Add(g_terma, "L");        mg->Add(g_kerma, "L");
    mg->Add(g_kerma_tot, "L");    mg->Add(g_dose_k_tot, "L");
    mg->Add(g_dose_riemann, "L"); mg->Add(g_dose_k_coll, "L");
    mg->Add(g_dose_t_tot, "L");   mg->Add(g_dose_t, "L");
    mg->Add(g_dose_k_tot_e, "L");
    mg->Draw("A");
    mg->GetXaxis()->SetRangeUser(0.0, 15.0);

    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->AddEntry(g_terma, "TERMA", "l");
    leg->AddEntry(g_kerma_tot, "KERMA (ODE flusso primario e scatter)", "l");
    leg->AddEntry(g_kerma, "KERMA (ODE solo flusso primario)", "l");
    leg->AddEntry(g_dose_riemann, "Dose K_{coll} (Somma di Riemann)", "l");
    leg->AddEntry(g_dose_k_coll, "Dose K_{coll} (Gauss-Kronrod)", "l");
    leg->AddEntry(g_dose_t_tot, "Dose TERMA (Flusso primario e scatter)", "l");
    leg->AddEntry(g_dose_t, "Dose TERMA (Funzione continua)", "l");
    leg->AddEntry(g_dose_k_tot, "Dose KERMA (ODE flusso primario e scatter)", "l");
    leg->AddEntry(g_dose_k_tot_e, "Dose KERMA (ODE flusso primario e scatter, senza 3 gen)", "l");
    leg->Draw();
}
