void plot_mu() {
    // 1. Crea la finestra grafica (Canvas)
    TCanvas *c1 = new TCanvas("c1", "Confronto Attenuazione NIST vs Chebyshev", 800, 600);
    
    // 2. Imposta la scala logaritmica su entrambi gli assi
    c1->SetLogx();
    c1->SetLogy();
    c1->SetGrid(); // Aggiunge la griglia per leggere meglio i valori

    // 3. Legge i dati dal file generato dal tuo programma
    // Il primo TGraph legge colonna 1 (E) e colonna 2 (NIST)
    TGraph *g_nist = new TGraph("confronto_mu.dat", "%lg %lg");
    
    // Il secondo TGraph legge colonna 1 (E) e colonna 3 (Chebyshev), saltando la 2 (%*lg)
    TGraph *g_cheb = new TGraph("confronto_mu.dat", "%lg %*lg %lg");

    // 4. Stile della linea NIST (Blu, solida)
    g_nist->SetTitle("Attenuazione Fotoni in Acqua;Energia E (MeV);#mu/\\rho (cm^{2}/g)");
    g_nist->SetLineColor(kBlue);
    g_nist->SetLineWidth(3);
    
    // 5. Stile della linea Chebyshev (Rossa, tratteggiata)
    g_cheb->SetLineColor(kRed);
    g_cheb->SetLineStyle(7);
    g_cheb->SetLineWidth(3);

    // 6. Disegna i grafici (AL = Assi e Linea, L SAME = Linea sovrapposta al precedente)
    g_nist->Draw("AL");
    g_cheb->Draw("L SAME");

    // 7. Aggiunge una legenda in alto a destra
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(g_nist, "Dati NIST (Ponte Log-Log)", "l");
    leg->AddEntry(g_cheb, "Interpolazione Chebyshev", "l");
    leg->SetBorderSize(0);
    leg->Draw();
}
