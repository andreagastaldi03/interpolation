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

    std::vector<double> E_csda_data;
    std::vector<double> Rcsda_data;
 
    std::ifstream file("nist_water.txt");
    if (!file.is_open()) {
        std::cerr << "Errore: Impossibile trovare nist_water.txt" 
            << std::endl;
        return 1;
    }

    std::ifstream file_csda("estar_water.txt");
    double Ec, Rc;
    while (file_csda >> Ec >> Rc) {
        E_csda_data.push_back(Ec);
        Rcsda_data.push_back(Rc);
    }
    file_csda.close();
    std::cout << "Letti " << E_csda_data.size() << " punti dal database NIST ESTAR."
        << std::endl;

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
/*    const double E_min_r = E_csda_data.front();
    const double E_max_r = E_csda_data.back();
*/
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

    auto Rcsda_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_csda_data.begin(), E_csda_data.end(), E_req);
        if (it == E_csda_data.begin()) return Rcsda_data.front();
        if (it == E_csda_data.end())   return Rcsda_data.back();

        int i = std::distance(E_csda_data.begin(), it);
        double E1 = E_csda_data[i-1], E2 = E_csda_data[i];
        double R1 = Rcsda_data[i-1],  R2 = Rcsda_data[i];

        double log_R = std::log(R1) + (std::log(R2) - std::log(R1))
            * (std::log(E_req) - std::log(E1)) / (std::log(E2) - std::log(E1));
        return std::exp(log_R);
    };

    // 3. creazione grid sull'intervallo usando tassellazione
    std::vector<double> tasselli = {E_data.front(), 0.05, 0.1, 1, 5, E_data.back()};
    std::vector<size_t> gradi = {20, 20, 20, 20, 20};
    
    /*std::vector<double> tasselli_csda = {
        E_csda_data.front(), 0.05, 0.1, 0.5, 1, 2, 5, 10, 100, E_csda_data.back()
    };
    std::vector<size_t> gradi_csda = {20, 20, 20, 20, 20, 20, 20, 20, 20};

    Interpolation::Grid1D grid_csda(
        Interpolation::make_discretization_info<Interpolation::details::log_0_maps>(tasselli_csda, gradi_csda));
*/
    Interpolation::Grid1D grid(
        Interpolation::make_discretization_info<Interpolation::details::log_0_maps>(tasselli, gradi));

    auto function_for_chebyshev = [&](double E) {
        return std::log(mu_raw_loglog(E));
    };

    auto function_for_chebyshev_en = [&](double E) {
        return std::log(mu_en_raw_loglog(E));
    };
/*
    auto function_for_chebyshev_r = [&](double E) {
        return std::log(Rcsda_raw_loglog(E));
    };
*/
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
/*
    std::vector<double> fj_r = Discretize<std::vector<double>, double>(
            grid_csda,
            function_for_chebyshev_r,
            [](size_t n) {return std::vector<double>(n, 0.); }
    );
*/
    std::cout << "Dati NIST discretizzati in una Grid1D a " << gradi.size() 
        << " tasselli." << std::endl;

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
/*
    std::ofstream out_r("confronto_rcsda.dat");
    for(int i=0; i<1000; ++i) {
        double E = E_min_r * std::exp(i * std::log(E_max_r/E_min_r) / (1000-1));
        
        double log_interp = grid_csda.interpolate<double, std::vector<double>>(
            E, fj_r, []() -> double {return 0.;}    
        );

        double val_interp = std::exp(log_interp);
    out_r << E << " " << Rcsda_raw_loglog(E) << " " << val_interp  << std::endl;
    }
    out_r.close();
*/
    auto mu_continuo = [&](double E_req) {
        double log_mu = grid.interpolate<double, std::vector<double>>(
            E_req, fj, []() -> double {return 0.;}        
        );
        return std::exp(log_mu);
    };

    auto mu_en_continuo = [&](double E_req) {
        double log_mu_en = grid.interpolate<double, std::vector<double>>(
            E_req, fj_en, []() -> double {return 0.;}        
        );
        return std::exp(log_mu_en);
    };
/*
    auto Rcsda_continuo = [&](double E_req) {
        double log_r = grid_csda.interpolate<double, std::vector<double>>(
            E_req, fj_r, []() -> double {return 0.;}        
        );
        return std::exp(log_r);
    };
*/
    std::cout << "Dati salvati in confronto_mu.dat, confronto_mu_en.dat e confronto_rcsda.dat" << std::endl;
    
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
/*
    double max_abs_err_r = 0.;
    double max_rel_err_r = 0.;
    double min_abs_err_r = 100.;
    double min_rel_err_r = 1.;

    double sum_rel_err_r = 0.;
    double E_worst_abs_r = 0.;
    double E_worst_rel_r = 0.;
    double E_best_abs_r = 0.;
    double E_best_rel_r = 0.;
*/
    for(int i=0; i<n_test; i++) {
        double E = E_min * std::exp(i * std::log(E_max / E_min) / (n_test - 1));
        //double E_r = E_min_r * std::exp(i * std::log(E_max_r / E_min_r) 
        //        / (n_test - 1));

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
/*
        double exact_r = Rcsda_raw_loglog(E_r);
        double log_interp_r = grid.interpolate<double, std::vector<double>>(
            E_r, fj_r, []() -> double {return 0.;}
        );
        double interp_r = std::exp(log_interp_r);

        double abs_err_r = std::abs(exact_r - interp_r);
        double rel_err_r = abs_err_r / exact_r;

        if (abs_err_r > max_abs_err_r) {
            max_abs_err_r = abs_err_r;
            E_worst_abs_r = E_r;
        }
        if (rel_err_r > max_rel_err_r) {
            max_rel_err_r = rel_err_r;
            E_worst_rel_r = E_r;
        }

        if (abs_err_r < min_abs_err_r) {
            min_abs_err_r = abs_err_r;
            E_best_abs_r = E_r;
        }
        if (rel_err_r < min_rel_err_r) {
            min_rel_err_r = rel_err_r;
            E_best_rel_r = E_r;
        }
        sum_rel_err_r += rel_err_r;
 */
    }
    double mean_rel_err = sum_rel_err / n_test;
    double mean_rel_err_en = sum_rel_err_en / n_test;
    //double mean_rel_err_r = sum_rel_err_r / n_test;

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
    std::cout << "\n" << std::endl;
    std::cout << "Errore Assoluto Massimo (mu_en): " << max_abs_err_en << " (E = " 
        << E_worst_abs_en << " MeV)" << std::endl;
    std::cout << "Errore Relativo Massimo (mu_en): " << max_rel_err_en * 100 
        << "% (E = " << E_worst_rel_en << " MeV)" << std::endl;
    std::cout << "Errore Assoluto Minimo (mu_en): " << min_abs_err_en << " (E = " 
        << E_best_abs_en << " MeV)" << std::endl;
    std::cout << "Errore Relativo Minimo (mu_en): " << min_rel_err_en * 100 
        << "% (E = " << E_best_rel_en << " MeV)" << std::endl;
    std::cout << "Errore Relativo Medio (mu_en): " << mean_rel_err_en * 100 
        << "%" << std::endl;
/*    std::cout << "\n" << std::endl;
    std::cout << "Errore Assoluto Massimo (r_csda): " << max_abs_err_r << " (E = " 
        << E_worst_abs_r << " MeV)" << std::endl;
    std::cout << "Errore Relativo Massimo (r_csda): " << max_rel_err_r * 100 
        << "% (E = " << E_worst_rel_r << " MeV)" << std::endl;
    std::cout << "Errore Assoluto Minimo (r_csda): " << min_abs_err_r << " (E = " 
        << E_best_abs_r << " MeV)" << std::endl;
    std::cout << "Errore Relativo Minimo (r_csda): " << min_rel_err_r * 100 
        << "% (E = " << E_best_rel_r << " MeV)" << std::endl;
    std::cout << "Errore Relativo Medio (r_csda): " << mean_rel_err_r * 100 
        << "%" << std::endl;
*/

////////////////////////////////////////////////////////////////////////////////////
// FINE PARTE CALCOLO MU
// INIZIO CALCOLO ODE RK
////////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n" << std::endl;

    // 1. energia e mu del test
    double E_beam = 1.75;
/*    double log_mu = grid.interpolate<double, std::vector<double>>(
        E_beam, fj, []() -> double {return 0.;});
    double log_mu_en = grid.interpolate<double, std::vector<double>>(
        E_beam, fj_en, []() -> double {return 0.;});
    double mu_const = std::exp(log_mu);
    double mu_en_const = std::exp(log_mu_en);*/
    double mu_const = mu_continuo(E_beam);
    double mu_en_const = mu_en_continuo(E_beam);

    double E_beam_sc = 0.60;
/*    double log_mu_sc = grid.interpolate<double, std::vector<double>>(
        E_beam_sc, fj, []() -> double {return 0.;});
    double log_mu_en_sc = grid.interpolate<double, std::vector<double>>(
        E_beam_sc, fj_en, []() -> double {return 0.;});
    double mu_const_sc = std::exp(log_mu_sc);
    double mu_en_const_sc = std::exp(log_mu_en_sc);*/
    double mu_const_sc = mu_continuo(E_beam_sc);
    double mu_en_const_sc = mu_en_continuo(E_beam_sc);

    double f_fwd = 0.80;
    double mu_scatt_create = mu_const - mu_en_const;

    std::cout << "Simulazione fascio. " << std::endl;
    std::cout << "Parametri fisici: " << E_beam << " MeV / mu = " << mu_const
        << " / mu_en = " << mu_en_const << " -> Parametri fascio primario" 
        << std::endl;
    std::cout << "Parametri fisici: " << E_beam_sc << " MeV / mu = " << mu_const_sc
        << " / mu_en = " << mu_en_const_sc << " -> Parametri fascio scatterato" 
        << std::endl;

    // 2. funzione rhs (right hand side)
    // solutore passa a un vettore state, ne leggiamo il valore e lo sovrasciviamo
    // col valore della derivata
    rk::vd<2> phi_init; // condizioni iniziali
    phi_init[0] = 1.0;  // fluenza primaria 1
    phi_init[1] = 0.0;  // fluenza scatterata 0 

    rk::rk_rhs_t<rk::vd<2>> rhs = [mu_const, mu_const_sc, mu_scatt_create, f_fwd,
        E_beam, E_beam_sc](double /*z*/, rk::vd<2>& state) {
        double phi_p = state(0); 
        double phi_s = state(1);

        double dp_dz = - mu_const * phi_p;
        double ds_dz = (mu_scatt_create * (E_beam / E_beam_sc) * f_fwd) 
            * phi_p - mu_const_sc * phi_s; 
        
        state[0] = dp_dz;
        state[1] = ds_dz;
    };
    
    // 3. condizioni iniziali
    std::vector<double> z_points;
    for(int i=0; i <=50; i++) {
        z_points.push_back(i * 0.5); // passi di 0.5 cm
    }
    // TimeInfo gestisce la griglia, chiamato time ma per noi è spazio
    rk::TimeInfo spatial_info(z_points);

    // 4. inizializzazione solutore
    rk::RungeKutta<rk::vd<2>, 4> solver(
        rk::PreImplementedTableau::RKOriginal,  // tab butcher per RK4
        phi_init,                               // condizione iniziale a z=0
        spatial_info,                           // punti spaziali
        rhs                                     // equazione
    );
   
    rk::rk_rhs_t<rk::vd<1>> rhs0 = [mu_const](double /*z*/, rk::vd<1>& state) {
        double current_phi = state(0);      // leggiamo fluenza attuale
        
        state[0] = -mu_const * current_phi; // scriviamo la derivata
                                            // dPhi/dz = - mu Phi
    };
    
    rk::vd<1> phi_initial;
    phi_initial[0] = 1.0; // fluenza a z=0, sul bordo del fantoccio, 100%

    rk::RungeKutta<rk::vd<1>, 4> solver0(
        rk::PreImplementedTableau::RKOriginal, // tab butcher per RK4
        phi_initial,                           // condizione iniziale a z=0
        spatial_info,                          // punti spaziali
        rhs0                                   // equazione
    );

    // 5. callback, esportiamo i dati passo passo mentre risolve
    std::ofstream out_ode0("fluenza_terma_z.dat");
    std::ofstream out_ode("fluenza_prim_scat.dat");
    
    std::vector<double> z_history;
    std::vector<double> terma_history;
    std::vector<double> kerma_history;
    
    std::vector<double> phi_p_history;
    std::vector<double> phi_s_history;
    std::vector<double> kerma_tot_history;

    double max_abs_err_ode = 0.0;
    double max_rel_err_ode = 0.0;
    double worst_z_abs = 0.0;
    double worst_z_rel = 0.0;

    rk::rk_callback_t<rk::vd<1>> callback0 = [&](double z, const rk::TimeInfo&, 
            const rk::vd<1>& phi0) {
        double fluenza_attuale = phi0(0);
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

        out_ode0 << z << " " << fluenza_attuale << " " << esatta << " "  
            << terma_attuale << " " << kerma_attuale << "\n";
    };
    solver0.AddCallback(callback0);
    solver0.CallbackOnlyOnTimeStamp(); // stampa solo ai z_points definiti da noi

    rk::rk_callback_t<rk::vd<2>> callback = [&](double z, const rk::TimeInfo&, 
            const rk::vd<2>& phi) {
        double phi_p = phi(0);
        double phi_s = phi(1);

        // KERMA locale depositato dai primari
        double kerma_p = mu_en_const * E_beam * phi_p;

        // KERMA locale depositato dagli scatterati (che viaggiano e poi assorbono)
        double kerma_s = mu_en_const_sc * E_beam_sc * phi_s;

        // Collision KERMA Totale
        double kerma_tot = kerma_p + kerma_s;

        //z_history.push_back(z);
        phi_p_history.push_back(phi_p);
        phi_s_history.push_back(phi_s);
        kerma_tot_history.push_back(kerma_tot);

        out_ode << z << " " << kerma_p << " " << kerma_s << " " << kerma_tot 
            << " " << phi_p << " " << phi_s << "\n";
    };

    solver.AddCallback(callback);
    solver.CallbackOnlyOnTimeStamp();

    // 6. risoluzione
    std::cout << "Avvio propagazione del fascio nel mezzo." << std::endl;
    solver0();
    out_ode0.close();
    solver();
    out_ode.close();
    std::cout << "Propagazione completata " << std::endl;
    std::cout << "Dati salvati fluenza_terma.dat" << std::endl;
    std::cout << "Dati salvati fluenza_prim_scat.dat" << std::endl;

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
/*    double a_pr = 2.0; // elettroni primari, decadimento rapido
    double A_pr = a_pr / 2.0;

    // funzione kernel semplificata
    auto dose_kernel = [A_pr, a_pr](double distanza) {
        return A_pr * std::exp(-a_pr * std::abs(distanza));
    };

    double w_el = mu_en_const / mu_const;
    double w_sc = 1.0 - w_el;
    
    double a_fwd_el = 4.8;
    double a_bwd_el = 8.0;
    double A_el = 1.0 / ((1.0 / a_fwd_el) + (1.0 / a_bwd_el));

    double a_s_fwd = 0.8;
    double a_s_bwd = 1.2;
    double A_sc = 1.0 / ((1.0 / a_s_fwd) + (1.0 / a_s_bwd));

    auto kernel_terma = [w_el, w_sc, A_el, a_fwd_el, a_bwd_el, A_sc, a_s_fwd, 
         a_s_bwd](double distanza) {
        double term_el = 0.0;
        if (distanza >= 0) {
            term_el = (w_el * A_el) * std::exp(- a_fwd_el * distanza);
        } else {
            term_el = (w_el * A_el) * std::exp(a_bwd_el * distanza);
        }

        double term_sc = 0.0;
        if (distanza >= 0) {
            term_sc = (w_sc * A_sc) * std::exp(- a_s_fwd * distanza);
        } else {
            term_sc = (w_sc * A_sc) * std::exp(a_s_bwd * distanza);
        }

        return term_el + term_sc;
    };

    auto kernel_elettroni = [A_el, a_fwd_el, a_bwd_el](double distanza) {
        if (distanza >= 0.0) return A_el * std::exp(-a_fwd_el * distanza);
        else                 return A_el * std::exp(a_bwd_el * distanza);
    };
*/
    double k_fwd_el = 2.5; // TODO calibrare
    double k_bwd_el = 4.0; // TODO calibrare

    double k_s_fwd = 1.0; // TODO calibrare
    double k_s_bwd = 1.8; // TODO calibrare

    double R_csda_primari = Rcsda_raw_loglog(E_beam);
    double R_csda_scatter = Rcsda_raw_loglog(E_beam_sc);

    // crea kernel elettronico normalizzato per una data energia, usando Rcsda
    struct kernelExp {
        double A, a_fwd, a_bwd;
        double operator()(double distanza) const {
            if (distanza>=0.0) return A * std::exp(-a_fwd * distanza);
            else               return A * std::exp( a_bwd * distanza);
        }
    };

    auto make_kernel_elettroni = [&](double R) {
        double a_fwd = k_fwd_el / R;
        double a_bwd = k_bwd_el / R;
        double A = 1.0 / ((1.0 / a_fwd) + (1.0 / a_bwd)); // normalizzazione
        return kernelExp{A, a_fwd, a_bwd};
    };

    kernelExp kernel_el_primari = make_kernel_elettroni(R_csda_primari);
    kernelExp kernel_el_scatter = make_kernel_elettroni(R_csda_scatter);

    std::cout << "Kernel elettronico primari:  a_fwd=" << kernel_el_primari.a_fwd
        << " a_bwd=" << kernel_el_primari.a_bwd
        << " (R_csda=" << R_csda_primari << " cm)" << std::endl;
    std::cout << "Kernel elettronico scatter:  a_fwd=" << kernel_el_scatter.a_fwd
        << " a_bwd=" << kernel_el_scatter.a_bwd
        << " (R_csda=" << R_csda_scatter << " cm)" << std::endl;

    // kernel scatter fotonico per TERMA, basato sul libero cammino medio 1/mu
    double lambda_sc = 1.0 / mu_const_sc;
    double a_s_fwd = k_s_fwd / lambda_sc;
    double a_s_bwd = k_s_bwd / lambda_sc;
    double A_sc = 1.0 / ((1.0 / a_s_fwd) + (1.0 / a_s_bwd));

    // pesi energetici
    double w_el = mu_en_const / mu_const;
    double w_sc = 1.0 - w_el;

    // kernel TERMA elettornico + scatter
    auto kernel_terma_fisico = [=](double distanza) {
        double term_el = w_el * kernel_el_primari(distanza);
        double term_sc = 0.0;
        if (distanza >= 0.0) term_sc = (w_sc * A_sc) * std::exp(-a_s_fwd * distanza);
        else                 term_sc = (w_sc * A_sc) * std::exp(a_s_bwd * distanza);
        return term_el + term_sc;
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

    auto kerma_tot_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return kerma_tot_history.front();
        if (it == z_history.end()) return kerma_tot_history.back();

        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double k0 = kerma_tot_history[i-1], k1 = kerma_tot_history[i];

        return k0 + (k1 - k0) * (z_req - z0) / (z1 - z0);
    };

    auto phi_p_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return phi_p_history.front();
        if (it == z_history.end()) return phi_p_history.back();

        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double k0 = phi_p_history[i-1], k1 = phi_p_history[i];

        return k0 + (k1 - k0) * (z_req - z0) / (z1 - z0);
    };

    auto phi_s_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return phi_s_history.front();
        if (it == z_history.end()) return phi_s_history.back();

        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double k0 = phi_s_history[i-1], k1 = phi_s_history[i];

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
 
    std::vector<double> fj_kerma_tot = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma,
        kerma_tot_raw_linear,
        [](size_t n) {return std::vector<double>(n, 0.);}        
    );  
 
    std::vector<double> fj_phi_p = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma,
        phi_p_raw_linear,
        [](size_t n) {return std::vector<double>(n, 0.);}        
    );
 
    std::vector<double> fj_phi_s = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma,
        phi_s_raw_linear,
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

    auto kerma_tot_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(
            z_req, fj_kerma_tot, []() -> double {return 0.;}        
        );
    };
    std::cout << "KERMA TOT discretizzato spazialmente su Grid1D." << std::endl;

    auto phi_p_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(
            z_req, fj_phi_p, []() -> double {return 0.;}        
        );
    };
    std::cout << "Flusso Primario discretizzato spazialmente su Grid1D." 
        << std::endl;

    auto phi_s_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(
            z_req, fj_phi_s, []() -> double {return 0.;}        
        );
    };
    std::cout << "Flusso Scatter discretizzato spazialmente su Grid1D." 
        << std::endl;

    auto kerma_p_continuo = [&](double z_req) {
        return mu_en_const * E_beam * phi_p_continuo(z_req);
    };
    auto kerma_s_continuo = [&](double z_req) {
        return mu_const_sc * E_beam_sc * phi_s_continuo(z_req);
        // uso mu_const anzichè mu_en perchè mu_en scarterebbe tutti i
        // fotoni terziari, la terza generazione, allora approx con tutta
        // energia assorbita
    };

    // 5. integrale di gauss-kronrod
    std::ofstream out_dose("profilo_dose_gk.dat");
    double toll_rel = 1.0e-5;
    double toll_abs = 1.0e-8;
    double z_max = z_history.back(); // integriamo su tutto lo spessore
 
    std::ofstream out_dose_riemann("profilo_dose.dat");
    std::vector<double> dose_history_riemann;

    std::cout << "Calcolo dose con Somma di Riemann." << std::endl;

    std::vector<double> z_riemann_check, dose_riemann_check;
    for (double z = 0.0; z <= 15.0; z += 0.1) {
       double dose_z = 0.0;
       double dz_fine = 0.02;

        // somma di riemann, integro su tutti gli z'
        for (double z_prime = 0.0; z_prime <= z_max; z_prime += dz_fine) {
            dose_z += kerma_continuo(z_prime) * kernel_el_primari(z - z_prime) * dz_fine;
        }
        z_riemann_check.push_back(z);
        dose_riemann_check.push_back(dose_z);

        out_dose_riemann << z << " " << kerma_continuo(z) << " " << dose_z << "\n";
    }
    out_dose_riemann.close();
 

    std::cout << "Calcolo dose con quadratura Gauss-Kronrod." << std::endl;

    std::vector<double> z_check, dose_k_check, dose_t_check, dose_kt_check;
    
    for (double z=0.0; z<=15.0; z+=0.1) { // integro su 15 cm ignorando padding 
                                          // finale
        // definisco integranda per lo specifico z
        auto integranda_k = [&](double z_prime) {
            return kerma_continuo(z_prime) * kernel_el_primari(z - z_prime);
        };

        double dose_fwd_k = (z > 0.0)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_k, 0.0, z, toll_rel, toll_abs)
            : 0.0;
        double dose_bwd_k = (z < z_max)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_k, z, z_max, toll_rel, toll_abs)
            : 0.0;
        double dose_z_k = dose_fwd_k + dose_bwd_k;

        auto integranda_t = [&](double z_prime) {
            return terma_continuo(z_prime) * kernel_terma_fisico(z - z_prime);
        };

        double dose_fwd_t = (z > 0.0)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_t, 0.0, z, toll_rel, toll_abs)
            : 0.0;
        double dose_bwd_t = (z < z_max)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_t, z, z_max, toll_rel, toll_abs)
            : 0.0;
        double dose_z_t = dose_fwd_t + dose_bwd_t;

        auto integranda_kt = [&](double z_prime) {
            return kerma_p_continuo(z_prime) * kernel_el_primari(z - z_prime)
                + kerma_s_continuo(z_prime) * kernel_el_scatter(z - z_prime);
        };

        double dose_fwd_kt = (z > 0.0)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_kt, 0.0, z, toll_rel, toll_abs)
            : 0.0;
        double dose_bwd_kt = (z < z_max)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_kt, z, z_max, toll_rel, toll_abs)
            : 0.0;
        double dose_z_kt = dose_fwd_kt + dose_bwd_kt;

        z_check.push_back(z);
        dose_k_check.push_back(dose_z_k);
        dose_t_check.push_back(dose_z_t);
        dose_kt_check.push_back(dose_z_kt);

        out_dose << z << " " << terma_continuo(z) << " " << kerma_continuo(z) << 
           " " << dose_z_t << " " << dose_z_k << " " << 
           (kerma_p_continuo(z) + kerma_s_continuo(z)) << " " << dose_z_kt << "\n";
    }
    out_dose.close();

    std::cout << "Calcolo dose completato." << std::endl;
    std::cout << "Dati salvati in profilo_dose.dat e profilo_dose_gk.dat" 
        << std::endl; 

    auto integra_trapezi = [](const std::vector<double>& z, 
            const std::vector<double>& f) {
        double s = 0.0;
        for (size_t i = 0; i + 1 < z.size(); i++) {
            s += 0.5 * (f[i] + f[i+1]) * (z[i+1] - z[i]);
        }
        return s;
    };

    std::vector<double> kcol_tot_check;
    for (double z : z_check) {
        kcol_tot_check.push_back(kerma_p_continuo(z)+kerma_s_continuo(z));
    }

    // energia "di terza generazione" mai contabilizzata dal modello a 2 stati
    auto terza_generazione = [&](double z) {
        return (mu_const_sc - mu_en_const_sc) * E_beam_sc * phi_s_continuo(z);
    };

    double integrale_terza_gen = 0.0;
    double dz = 0.05;
    for (double z = 0.0; z <= 15.0; z += dz) {
        integrale_terza_gen += terza_generazione(z) * dz;
    }
    std::cout << "Energia di terza generazione non contabilizzata: "
        << integrale_terza_gen << std::endl;

    std::cout << "Riferimento, K_col totale (primari+scatter): "
        << integra_trapezi(z_check, kcol_tot_check) << std::endl;
    std::cout << "Riemann (K_col primari, kernel elettronico):  "
        << integra_trapezi(z_riemann_check, dose_riemann_check) << std::endl;
    std::cout << "GK K_col (K_col primari, kernel elettronico): "
        << integra_trapezi(z_check, dose_k_check) << std::endl;
    std::cout << "GK TERMA (TERMA, kernel elettr.+scatter):     "
        << integra_trapezi(z_check, dose_t_check) << std::endl;
    std::cout << "GK K_tot (K_col p.+s., kernel separati):      "
        << integra_trapezi(z_check, dose_kt_check) << std::endl;
   
    return 0;
}
