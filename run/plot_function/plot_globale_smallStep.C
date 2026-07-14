void plot_globale_smallStep() {
    TCanvas *c_all = new TCanvas("c_all_small", "Panoramica Completa Radioterapia", 1000, 800);
    c_all->SetGrid();

    // Lettura da profilo_dose_gk.dat
    // Colonne: 1=z, 2=Terma, 3=Kerma1, 4=DoseT, 5=DoseK, 6=KermaTot, 7=DoseKTot, 8=DoseTTot
    TGraph *g_terma     = new TGraph("profilo_dose_gk_small.dat", "%lg %lg %*lg %*lg %*lg %*lg %*lg");
    TGraph *g_kerma1    = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %lg %*lg %*lg %*lg %*lg");
    TGraph *g_kerma_tot = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %lg %*lg");
    
    TGraph *g_dose_t    = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %lg %*lg %*lg %*lg");
    TGraph *g_dose_k    = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %lg %*lg %*lg");
    TGraph *g_dose_ktot = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %*lg %lg");
    TGraph *g_dose_ttot = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");
    TGraph *g_dose_kten = new TGraph("profilo_dose_gk_small.dat", "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg");

    // SORGENTI (Sfondi tratteggiati)
    g_terma->SetLineColor(kGray+2);      g_terma->SetLineWidth(3); g_terma->SetLineStyle(3);
    g_kerma1->SetLineColor(kTeal-6);     g_kerma1->SetLineWidth(3); g_kerma1->SetLineStyle(7);
    g_kerma_tot->SetLineColor(kGreen+3); g_kerma_tot->SetLineWidth(3); g_kerma_tot->SetLineStyle(7);

    // DOSI (Linee continue e colorate in primo piano)
    g_dose_t->SetLineColor(kBlue);       g_dose_t->SetLineWidth(3);
    g_dose_k->SetLineColor(kP10Yellow);  g_dose_k->SetLineWidth(3);
    g_dose_ktot->SetLineColor(kRed);     g_dose_ktot->SetLineWidth(3); 
    g_dose_ttot->SetLineColor(kP8Pink);  g_dose_ttot->SetLineWidth(3);
    g_dose_kten->SetLineColor(kP10Brown); g_dose_kten->SetLineWidth(3);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Integrazione Globale: Sorgenti KERMA/TERMA e Dose Assorbita;Profondita' z (cm);Energia/Dose (a.u.)");
    
    mg->Add(g_terma, "L");     mg->Add(g_kerma1, "L");    mg->Add(g_kerma_tot, "L");
    mg->Add(g_dose_t, "L");    mg->Add(g_dose_k, "L");    mg->Add(g_dose_ktot, "L");
    mg->Add(g_dose_ttot, "L"); mg->Add(g_dose_kten, "L");

    mg->Draw("A");
    mg->GetXaxis()->SetRangeUser(0.0, 15.0);

    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->AddEntry(g_terma, "TERMA", "l");
    leg->AddEntry(g_kerma1, "KERMA (Solo Primari)", "l");
    leg->AddEntry(g_kerma_tot, "KERMA (ODE, P+S)", "l"); 
    leg->AddEntry(g_dose_t, "Dose da TERMA continuo", "l");
    leg->AddEntry(g_dose_k, "Dose da KERMA (Solo Primari)", "l");
    leg->AddEntry(g_dose_ktot, "Dose da KERMA (flusso primario e scatter)", "l");
    leg->AddEntry(g_dose_ttot, "Dose da TERMA (flusso primario e scatter)", "l");
    leg->AddEntry(g_dose_kten, "Dose da KERMA (flusso primario e scatter, senza 3 gen)");
    leg->SetBorderSize(1);
    leg->Draw();
}
