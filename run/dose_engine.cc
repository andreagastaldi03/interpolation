#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "Interpolation/interpolation.hh"
#include "runge_kutta.hh"


int main ()
{
    // 1. lettura file di testo contenente attenuazione
    std::vector<double> E_data;
    std::vector<double> mu_data;
    std::vector<double> mu_en_data;

    std::ifstream file("nist_water.txt");
    if (!file.is_open()) {
        std::cerr << "Errore: Impossibile trovare nist_water.txt" 
            << std::endl;
        return 1;
    }

    double E_val, mu_val, mu_en_val;
    while (file >> E_val >> mu_val >> mu_en_val) {
        E_data.push_back(E_val);
        mu_data.push_back(mu_val);
        mu_en_data.push_back(mu_en_val);
    }
    file.close();
    std::cout << "Letti " << E_data.size() << " punti dal database NIST."
        << std::endl;

    const double E_min = E_data.front(); // MeV
    const double E_max = E_data.back();  // MeV

    // 2. interpolatore loglog sui dati
    // questa funzione lambda prende un E qualsiasi e calcola il mu_rho 
    auto mu_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_data.begin(), E_data.end(), E_req);
        // trova il primo punto E_data >= E_req

        // gestione dei bordi
        if (it == E_data.begin()) return mu_data.front();
        if (it == E_data.end()) return mu_data.back();

        // indici dei due punti in cui cade E_req
        int i = std::distance(E_data.begin(), it);
        double E1 = E_data[i-1], E2 = E_data[i];
        double m1 = mu_data[i-1], m2 = mu_data[i];

        // interpolazione loglog
        // y=exp(log(y1)+(log(y2)-log(y1))*(log(x)-log(x1))/(log(x2)-log(x1))) 
        double log_E1 = std::log(E1);
        double log_E2 = std::log(E2);
        double log_E_req = std::log(E_req);
        double log_m1 = std::log(m1);
        double log_m2 = std::log(m2);

        double log_m_req = log_m1 + (log_m2 - log_m1) * (log_E_req - log_E1)
            / (log_E2 - log_E1);
        return std::exp(log_m_req);
    };

    auto mu_en_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_data.begin(), E_data.end(), E_req);
        // trova il primo punto E_data >= E_req

        // gestione dei bordi
        if (it == E_data.begin()) return mu_en_data.front();
        if (it == E_data.end()) return mu_en_data.back();

        // indici dei due punti in cui cade E_req
        int i = std::distance(E_data.begin(), it);
        double E1 = E_data[i-1], E2 = E_data[i];
        double m1 = mu_en_data[i-1], m2 = mu_en_data[i];

        // interpolazione loglog
        // y=exp(log(y1)+(log(y2)-log(y1))*(log(x)-log(x1))/(log(x2)-log(x1))) 
        double log_E1 = std::log(E1);
        double log_E2 = std::log(E2);
        double log_E_req = std::log(E_req);
        double log_m1 = std::log(m1);
        double log_m2 = std::log(m2);

        double log_m_req = log_m1 + (log_m2 - log_m1) * (log_E_req - log_E1)
            / (log_E2 - log_E1);
        return std::exp(log_m_req);
    };
    // 3. creazione grid sull'intervallo usando tassellazione
    std::vector<double> tasselli = {E_data.front(), 0.05, 0.1, 1, 5, E_data.back()};
    std::vector<size_t> gradi = {20, 20, 20, 20, 20};
    
    Interpolation::Grid1D grid(
            Interpolation::make_discretization_info<Interpolation::details::log_0_maps>(tasselli, gradi));

    auto function_for_chebyshev = [&](double E) {
        return std::log(mu_raw_loglog(E));
    };

    auto function_for_chebyshev_en = [&](double E) {
        return std::log(mu_en_raw_loglog(E));
    };

    std::vector<double> fj = Discretize<std::vector<double>, double>(
            grid,
            function_for_chebyshev,
            [](size_t n) {return std::vector<double>(n, 0.); }
    );

    std::vector<double> fj_en = Discretize<std::vector<double>, double>(
            grid,
            function_for_chebyshev_en,
            [](size_t n) {return std::vector<double>(n, 0.); }
    );

    std::cout << "Dati NIST discretizzati in una Grid1D a " << gradi.size() 
        << " tasselli." << std::endl;

    //double E_test = 1.75;
    //double exact_mu = mu_raw_loglog(E_test);

    //double log_interpolated_mu = grid.interpolate<double, std::vector<double>>(
    //    E_test, fj, []() -> double {return 0.;}
    //);

    //double interpolated_mu = std::exp(log_interpolated_mu);

    // 4. stampa di dati per plot confronto mu ricostruiti da dati NIST e mu interp
    std::ofstream out("confronto_mu.dat");
    for(int i=0; i<1000; ++i) {
        double E = E_min * std::exp(i * std::log(E_max/E_min) / (1000-1));
        
        double log_interp = grid.interpolate<double, std::vector<double>>(
            E, fj, []() -> double {return 0.;}    
        );

        double val_interp = std::exp(log_interp);
    out << E << " " << mu_raw_loglog(E) << " " << val_interp  << std::endl;
    }
    out.close();
    std::cout << "Dati salvati in confronto_mu.dat" << std::endl;

    std::ofstream out_en("confronto_mu_en.dat");
    for(int i=0; i<1000; ++i) {
        double E = E_min * std::exp(i * std::log(E_max/E_min) / (1000-1));
        
        double log_interp = grid.interpolate<double, std::vector<double>>(
            E, fj_en, []() -> double {return 0.;}    
        );

        double val_interp = std::exp(log_interp);
    out_en << E << " " << mu_en_raw_loglog(E) << " " << val_interp  << std::endl;
    }
    out_en.close();

    std::cout << "Dati salvati in confronto_mu.dat e confronto_mu_en.dat" 
        << std::endl;
    
    // 5. analisi globale degli errori su tutto l'intervallo
    int n_test = 10000;
    double max_abs_err = 0.;
    double max_rel_err = 0.;
    double min_abs_err = 100.;
    double min_rel_err = 1.;

    double sum_rel_err = 0.;
    double E_worst_abs = 0.;
    double E_worst_rel = 0.;
    double E_best_abs = 0.;
    double E_best_rel = 0.;

    double max_abs_err_en = 0.;
    double max_rel_err_en = 0.;
    double min_abs_err_en = 100.;
    double min_rel_err_en = 1.;

    double sum_rel_err_en = 0.;
    double E_worst_abs_en = 0.;
    double E_worst_rel_en = 0.;
    double E_best_abs_en = 0.;
    double E_best_rel_en = 0.;

    for(int i=0; i<n_test; i++) {
        double E = E_min * std::exp(i * std::log(E_max / E_min) / (n_test - 1));

        double exact = mu_raw_loglog(E);
        double log_interp = grid.interpolate<double, std::vector<double>>(
            E, fj, []() -> double {return 0.;}
        );
        double interp = std::exp(log_interp);

        double abs_err = std::abs(exact - interp);
        double rel_err = abs_err / exact;

        if (abs_err > max_abs_err) {
            max_abs_err = abs_err;
            E_worst_abs = E;
        }
        if (rel_err > max_rel_err) {
            max_rel_err = rel_err;
            E_worst_rel = E;
        }

        if (abs_err < min_abs_err) {
            min_abs_err = abs_err;
            E_best_abs = E;
        }
        if (rel_err < min_rel_err) {
            min_rel_err = rel_err;
            E_best_rel = E;
        }
        sum_rel_err += rel_err;

        double exact_en = mu_en_raw_loglog(E);
        double log_interp_en = grid.interpolate<double, std::vector<double>>(
            E, fj_en, []() -> double {return 0.;}
        );
        double interp_en = std::exp(log_interp_en);

        double abs_err_en = std::abs(exact_en - interp_en);
        double rel_err_en = abs_err_en / exact_en;

        if (abs_err_en > max_abs_err_en) {
            max_abs_err_en = abs_err_en;
            E_worst_abs_en = E;
        }
        if (rel_err_en > max_rel_err_en) {
            max_rel_err_en = rel_err_en;
            E_worst_rel_en = E;
        }

        if (abs_err_en < min_abs_err_en) {
            min_abs_err_en = abs_err_en;
            E_best_abs_en = E;
        }
        if (rel_err_en < min_rel_err_en) {
            min_rel_err_en = rel_err_en;
            E_best_rel_en = E;
        }
        sum_rel_err_en += rel_err_en;
    }
    double mean_rel_err = sum_rel_err / n_test;
    double mean_rel_err_en = sum_rel_err_en / n_test;

    std::cout << "\nRisultati del Test di Interpolazione, " << n_test 
        << " punti valutati." << std::endl;
    //std::cout << "Energia richiesta: " << E_test << " MeV" << std::endl;
    //std::cout << "Valore Esatto:      " << exact_mu << std::endl;
    //std::cout << "Valore Interpolato: " << interpolated_mu << std::endl;
    //std::cout << "Errore assoluto:    " << std::abs(exact_mu - interpolated_mu) 
    //    << std::endl;
    std::cout << "Errore Assoluto Massimo: " << max_abs_err << " (E = " 
        << E_worst_abs << " MeV)" << std::endl;
    std::cout << "Errore Relativo Massimo: " << max_rel_err * 100 << "% (E = " 
        << E_worst_rel << " MeV)" << std::endl;
    std::cout << "Errore Assoluto Minimo: " << min_abs_err << " (E = " 
        << E_best_abs << " MeV)" << std::endl;
    std::cout << "Errore Relativo Minimo: " << min_rel_err * 100 << "% (E = " 
        << E_best_rel << " MeV)" << std::endl;
    std::cout << "Errore Relativo Medio: " << mean_rel_err * 100 << "%" << std::endl;
    std::cout << "Errore Assoluto Massimo (mu_en): " << max_abs_err_en << " (E = " 
        << E_worst_abs_en << " MeV)" << std::endl;
    std::cout << "Errore Relativo Massimo (mu_en): " << max_rel_err_en * 100 << "% (E = " 
        << E_worst_rel_en << " MeV)" << std::endl;
    std::cout << "Errore Assoluto Minimo (mu_en): " << min_abs_err_en << " (E = " 
        << E_best_abs_en << " MeV)" << std::endl;
    std::cout << "Errore Relativo Minimo (mu_en): " << min_rel_err_en * 100 << "% (E = " 
        << E_best_rel_en << " MeV)" << std::endl;
    std::cout << "Errore Relativo Medio (mu_en): " << mean_rel_err_en * 100 << "%" << std::endl;

////////////////////////////////////////////////////////////////////////////////////
// FINE PARTE CALCOLO MU
// INIZIO CALCOLO ODE RK
////////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n" << std::endl;

    // 1. energia e mu del test
    double E_beam = 1.75;
    double log_mu = grid.interpolate<double, std::vector<double>>(
        E_beam, fj, []() -> double {return 0.;});
    double log_mu_en = grid.interpolate<double, std::vector<double>>(
        E_beam, fj_en, []() -> double {return 0.;});
    double mu_const = std::exp(log_mu);
    double mu_en_const = std::exp(log_mu_en);

    std::cout << "Simulazione fascio a " << E_beam << " MeV / mu = " << mu_const
        << " / mu_en = " << mu_en_const << std::endl;

    // 2. funzione rhs (right hand side)
    // solutore passa a un vettore state, ne leggiamo il valore e lo sovrasciviamo
    // col valore della derivata
    rk::rk_rhs_t<rk::vd<1>> rhs = [mu_const](double z, rk::vd<1>& state) {
        double current_phi = state(0); // leggiamo fluenza attuale
        
        state[0] = -mu_const * current_phi; // scriviamo la derivata
                                            // dPhi/dz = - mu Phi
    };

    // 3. condizioni iniziali
    rk::vd<1> phi_initial;
    phi_initial[0] = 1.0; // fluenza a z=0, sul bordo del fantoccio, 100%
    std::vector<double> z_points;
    for(int i=0; i <=50; i++) {
        z_points.push_back(i * 0.5); // passi di 0.5 cm
    }
    // TimeInfo gestisce la griglia, chiamato time ma per noi è spazio
    rk::TimeInfo spatial_info(z_points);

    // 4. inizializzazione solutore
    rk::RungeKutta<rk::vd<1>, 4> solver(
        rk::PreImplementedTableau::RKOriginal, // tab butcher per RK4
        phi_initial,                           // condizione iniziale a z=0
        spatial_info,                          // punti spaziali
        rhs                                    // equazione
    );

    // 5. callback, esportiamo i dati passo passo mentre risolve
    std::ofstream out_ode("fluenza_terma_z.dat");
    std::vector<double> z_history;
    std::vector<double> terma_history;
    std::vector<double> kerma_history;

    double max_abs_err_ode = 0.0;
    double max_rel_err_ode = 0.0;
    double worst_z_abs = 0.0;
    double worst_z_rel = 0.0;

    rk::rk_callback_t<rk::vd<1>> callback = [&](double z, const rk::TimeInfo&, 
            const rk::vd<1>& phi) {
        double fluenza_attuale = phi(0);
        double esatta = std::exp(-mu_const * z);
        double terma_attuale = mu_const * E_beam * fluenza_attuale;
        double kerma_attuale = mu_en_const * E_beam * fluenza_attuale;
        
        double abs_err = std::abs(fluenza_attuale-esatta);
        double rel_err = (esatta > 0.0) ? (abs_err / esatta) : 0.0;

        if (abs_err > max_abs_err_ode) {
            max_abs_err_ode = abs_err;
            worst_z_abs = z;
        }
        if (rel_err > max_rel_err_ode) {
            max_rel_err_ode = rel_err;
            worst_z_rel = z;
        }
        z_history.push_back(z);
        terma_history.push_back(terma_attuale);
        kerma_history.push_back(kerma_attuale);

        out_ode << z << " " << fluenza_attuale << " " << esatta << " "  
            << terma_attuale << " " << kerma_attuale << "\n";
    };
    solver.AddCallback(callback);
    solver.CallbackOnlyOnTimeStamp(); // stampa solo ai z_points definiti da noi

    // 6. risoluzione
    std::cout << "Avvio propagazione del fascio nel mezzo." << std::endl;
    solver();
    out_ode.close();
    std::cout << "Propagazione completata " << std::endl;
    std::cout << "Dati salvati fluenza_terma.dat" << std::endl;

    std::cout << "Errore Assoluto Massimo: " << max_abs_err_ode << " (z = " <<
        worst_z_abs << " cm)" << std::endl;
    std::cout << "Errore Relativo Massimo: " << max_rel_err_ode * 100.0 << " % (z = " 
        << worst_z_rel << " cm)" << std::endl;

////////////////////////////////////////////////////////////////////////////////////
// FINE CALCOLO ODE RK
// INIZIO CALCOLO DOSE
////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n" << std::endl;

    // 1. definizione del kernel 
    // a = controlla la pendenza / A fattore di normalizzazione
    double a_pr = 2.0; // elettroni primari, decadimento rapido
    double A_pr = a_pr / 2.0;

    // funzione kernel semplificata
    auto dose_kernel = [A_pr, a_pr](double distanza) {
        return A_pr * std::exp(-a_pr * std::abs(distanza));
    };

    // 2. esecuzione dell'integrale discreto
    std::ofstream out_dose_riemann("profilo_dose.dat");
    std::vector<double> dose_history_riemann;

    std::cout << "Calcolo dose con Somma di Riemann." << std::endl;

    for (size_t i=0; i<z_history.size(); i++) {
       double z = z_history[i];
       double dose_z = 0.0;

        // somma di riemann, integro su tutti gli z'
        for (size_t j=0; j<z_history.size()-1; j++) {
            double z_prime = z_history[j];
            double terma_prime = terma_history[j];
            // passo di integrazione dz'
            double dz_prime = z_history[j+1] - z_history[j];
            // distanza tra il punto di interazione e il punto di calcolo
            double distanza = z - z_prime;

            dose_z += terma_prime * dose_kernel(distanza) * dz_prime;
        }
        dose_history_riemann.push_back(dose_z);

        out_dose_riemann << z << " " << terma_history[i] << " " << dose_z << "\n";
    }
    out_dose_riemann.close();

    // 1. definizione kernel 
    double c_p = 0.7; // 70% componente elettroni primari
    double c_s = 0.3; // 30% componente scatter
    
    double a_p = 5.0; // elettroni primari, decadimento rapido
    double A_p = c_p * (a_p / 2.0);
    
    double a_s_fwd = 0.3; // scatter forward lento
    double a_s_bwd = 1.2; // veloce decadimento all'indietro

    double A_s = c_s / ((1.0 / a_s_fwd) + (1.0 / a_s_bwd));

    auto kernel = [A_p, a_p, A_s, a_s_fwd, a_s_bwd](double distanza) {
        double term_p = A_p * std::exp(-a_p * std::abs(distanza));
        double term_s = 0.0;

        if (distanza >= 0) {
            term_s = A_s * std::exp(-a_s_fwd * distanza);
        } else {
            term_s = A_s * std::exp(a_s_bwd * distanza);
        }
        return term_p + term_s;
    };

    double a_fwd = 2.0;
    double a_bwd = 8.0;
    double A = 1.0 / ((1.0 / a_fwd) + (1.0 / a_bwd));

    auto kernel_elettroni = [A, a_fwd, a_bwd](double distanza) {
        if (distanza >= 0.0) return A * std::exp(-a_fwd * distanza);
        else                 return A * std::exp(a_bwd * distanza);
    };

    // 2. definizione terma e kerma continuo - ponte spaziale (spezzata)
    auto terma_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return terma_history.front();
        if (it == z_history.end()) return terma_history.back();

        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double t0 = terma_history[i-1], t1 = terma_history[i];

        return t0 + (t1 - t0) * (z_req - z0) / (z1 - z0);
    };

    auto kerma_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return kerma_history.front();
        if (it == z_history.end()) return kerma_history.back();

        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double k0 = kerma_history[i-1], k1 = kerma_history[i];

        return k0 + (k1 - k0) * (z_req - z0) / (z1 - z0);
    };

    // 3. definizione terma continuo - grid tassellazione chebishev
    std::vector<double> tasselli_z = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0};
    std::vector<size_t> gradi_z = {15, 15, 15, 15, 15};

    Interpolation::Grid1D grid_terma(
        Interpolation::make_discretization_info<Interpolation::details::identity_maps>(tasselli_z, gradi_z)        
    );

    std::vector<double> fj_terma = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma,
        terma_raw_linear,
        [](size_t n) {return std::vector<double>(n, 0.);}        
    );
    
    std::vector<double> fj_kerma = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma,
        kerma_raw_linear,
        [](size_t n) {return std::vector<double>(n, 0.);}        
    );
    
    // 4. definizione terma continuo - funzione finale liscia
    auto terma_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(
            z_req, fj_terma, []() -> double {return 0.;}        
        );
    };
    std::cout << "TERMA discretizzato spazialmente su Grid1D." << std::endl;

    auto kerma_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(
            z_req, fj_kerma, []() -> double {return 0.;}        
        );
    };
    std::cout << "KERMA discretizzato spazialmente su Grid1D." << std::endl;

    // 5. integrale di gauss-kronrod
    std::ofstream out_dose("profilo_dose_gk.dat");
    double z_max = z_history.back(); // integriamo su tutto lo spessore
    double toll_rel = 1.0e-5;
    double toll_abs = 1.0e-8;

    std::cout << "Calcolo dose con quadratura Gauss-Kronrod." << std::endl;

    for (double z=0.0; z<=15.0; z+=0.1) { // integro su 15 cm ignorando padding 
                                          // finale
        // definisco integranda per lo specifico z
        auto integranda_k = [&](double z_prime) {
            return kerma_continuo(z_prime) * kernel_elettroni(z - z_prime);
        };

        // contributo da monte
        double dose_fwd_k = 0.0;

        if (z > 0.0) {
            dose_fwd_k = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_k, 0.0, z, toll_rel, toll_abs);
        }

        // contributo da valle
        double dose_bwd_k = 0.0;
        if (z < z_max) {
            dose_bwd_k = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_k, z, z_max, toll_rel, toll_abs);
        }

        double dose_z_k = dose_fwd_k + dose_bwd_k;

        auto integranda_t = [&](double z_prime) {
            return terma_continuo(z_prime) * kernel(z - z_prime);
        };

        // contributo da monte
        double dose_fwd_t = 0.0;

        if (z > 0.0) {
            dose_fwd_t = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_t, 0.0, z, toll_rel, toll_abs);
        }

        // contributo da valle
        double dose_bwd_t = 0.0;
        if (z < z_max) {
            dose_bwd_t = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_t, z, z_max, toll_rel, toll_abs);
        }

        double dose_z_t = dose_fwd_t + dose_bwd_t;

        out_dose << z << " " << terma_continuo(z) << " " << kerma_continuo(z) << 
           " " << dose_z_k << " " << dose_z_t << "\n";
    }
    out_dose.close();

    std::cout << "Calcolo dose completato." << std::endl;
    std::cout << "Dati salvati in profilo_dose.dat e profilo_dose_gk.dat" 
        << std::endl; 
  
    return 0;
}
