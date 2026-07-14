#!/bin/bash

echo "Avvio il plot dei dati dosimetrici e ODE."

# Verifica che il file dei dati esista prima di lanciare ROOT
if [ ! -f confronto_mu.dat ]; then
    echo "Errore: Il file confronto_mu.dat non esiste. Lancia prima il motore dosimetrico."
    exit 1
fi

if [ ! -f fluenza_terma_z.dat ]; then
    echo "Errore: Il file fluenza_terma_z.dat non esiste. Lancia prima il motore dosimetrico."
    exit 1
fi

# Lancia ROOT eseguendo direttamente la macro
# -l : toglie il banner testuale iniziale di ROOT
root -l plot_confronto_dosi_noSc.C plot_kermas_noSc.C plot_globale_noSc.C
