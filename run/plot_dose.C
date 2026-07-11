void plot_dose() {
    TCanvas *c_phys = new TCanvas("c_phys", "Fisica della Dose (PDD)", 800, 600);
    c_phys->SetGrid();

    TGraph *g_terma = new TGraph("profilo_dose_gk.dat", "%lg %lg");
    
    TGraph *g_kerma = new TGraph("profilo_dose_gk.dat", "%lg %*lg %lg");
    
    TGraph *g_dose_k  = new TGraph("profilo_dose_gk.dat", "%lg %*lg %*lg %lg");
    
    TGraph *g_dose_t  = new TGraph("profilo_dose_gk.dat", "%lg %*lg %*lg %*lg %lg");

    g_terma->SetLineColor(kBlue);
    g_terma->SetLineStyle(7); // Tratteggiata
    g_terma->SetLineWidth(2);
    
    g_kerma->SetLineColor(kGreen+2);
    g_kerma->SetLineStyle(7); // Tratteggiata
    g_kerma->SetLineWidth(3);

    g_dose_k->SetLineColor(kRed);
    g_dose_k->SetLineWidth(4);

    g_dose_t->SetLineColor(kOrange);
    g_dose_t->SetLineStyle(7);
    g_dose_t->SetLineColor(3);


    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Profilo di Dose (PDD);Profondita' z (cm);Energia Rilasciata / Assorbita (a.u.)");
    
    // Aggiungiamo i tre grafici al contenitore
    mg->Add(g_terma, "L");
    mg->Add(g_kerma, "L");
    mg->Add(g_dose_k, "L");
    mg->Add(g_dose_t, "L");

    // Disegniamo con "A" per creare gli assi calibrati su tutti i dati
    mg->Draw("A");
    
    // Limitiamo la visualizzazione ai primi 15 cm per vedere bene il build-up
    mg->GetXaxis()->SetRangeUser(0.0, 15.0); 

    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->AddEntry(g_terma, "TERMA", "l");
    leg->AddEntry(g_kerma, "Collision KERMA", "l");
    leg->AddEntry(g_dose_k, "Dose (from KERMA)", "l");
    leg->AddEntry(g_dose_t, "Dose (from TERMA)", "l");
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->Draw();
}
