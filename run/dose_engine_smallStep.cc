#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "Interpolation/interpolation.hh"
#include "runge_kutta.hh"


int main ()
{
    ////////////////////////////////////////////////////////////////////////////////
    /// CREAZIONE DATI COEFFICIENTI DI ATTENUAZIONE
    ////////////////////////////////////////////////////////////////////////////////

    // 1. lettura file di testo contenente attenuazione e range elettroni
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
    if (!file_csda.is_open()) {
        std::cerr << "Errore: Impossibile trovare estar_water.txt" 
            << std::endl;
        return 1;
    }

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


    // 2. interpolatore loglog sui dati
    // questa funzione lambda prende un E qualsiasi e calcola il mu_rho 
    auto mu_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_data.begin(), E_data.end(), E_req);

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

        if (it == E_data.begin()) return mu_en_data.front();
        if (it == E_data.end()) return mu_en_data.back();

        int i = std::distance(E_data.begin(), it);
        double E1 = E_data[i-1], E2 = E_data[i];
        double m1 = mu_en_data[i-1], m2 = mu_en_data[i];
 
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

    // 4. stampa di dati per plot confronto mu ricostruiti da dati NIST e mu interpolato
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

    std::cout << "Dati salvati in confronto_mu.dat, confronto_mu_en.dat" << std::endl;
    
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

////////////////////////////////////////////////////////////////////////////////////
/// INIZIO CALCOLO ODE RK
////////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n" << std::endl;

    // 1. def parametri in uso nel codice
    double E_beam = 1.75;
    double mu_const = mu_continuo(E_beam);
    double mu_en_const = mu_en_continuo(E_beam);

    // definizione f_fwd e en beam scatterato, usando sez urto klein nishina
    double m_e_c2 = 0.511; // MeV
    double alpha = E_beam / m_e_c2;

    auto E_prime = [&](double theta) {
        return E_beam / (1.0 + alpha * (1.0 - std::cos(theta)));
    }; // energia di fotone scatterato

    auto dsigma_dOmega = [&](double theta) {
        double Ep = E_prime(theta);
        double r = Ep / E_beam;
        double sin2 = std::sin(theta) * std::sin(theta);
        return r*r * (r + 1.0 / r - sin2);
    };

    auto integranda_energia = [&](double theta) {
        return E_prime(theta) * dsigma_dOmega(theta) * std::sin(theta);
    };

    auto integranda_norm = [&](double theta) {
        return dsigma_dOmega(theta) * std::sin(theta);
    };

    double tol = 1e-8;
    double num_E   = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
            integranda_energia, 0.0, M_PI, tol, tol);
    double den_E   = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
            integranda_norm,    0.0, M_PI, tol, tol);
    double num_fwd = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
            integranda_energia, 0.0, M_PI/2.0, tol, tol);
    
    double E_beam_sc = num_E / den_E; // media pesata su tutto angolo solido E'
    double f_fwd = num_fwd / num_E; // fattore di fotoni scatterati forward e non persi

    // parametri calibrati su energia del fascio
    double mu_const_sc = mu_continuo(E_beam_sc); // attenuazione totale
    double mu_en_const_sc = mu_en_continuo(E_beam_sc); // attenuazione energia assorbita

    double mu_scatt_create = mu_const - mu_en_const; // coefficiente energia scatterata

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

    // primo approccio al problema, decadimento exp di fotoni primari solo  
    rk::vd<1> phi_initial;
    phi_initial[0] = 1.0; // fluenza a z=0, sul bordo del fantoccio, 100%

    rk::rk_rhs_t<rk::vd<1>> rhs0 = [mu_const](double /*z*/, rk::vd<1>& state) {
        double current_phi = state(0);      // leggiamo fluenza attuale
        
        state[0] = -mu_const * current_phi; // scriviamo la derivata
                                            // dPhi/dz = - mu Phi
    };

    // approccio due stati, flusso fotoni primari e secondari
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
    
    // 3. specificazione problema
    std::vector<double> z_points;
    for(int i=0; i <=100; i++) {
        z_points.push_back(i * 0.25); // passi di 0.5 cm
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

    rk::RungeKutta<rk::vd<1>, 4> solver0(
        rk::PreImplementedTableau::RKOriginal, // tab butcher per RK4
        phi_initial,                           // condizione iniziale a z=0
        spatial_info,                          // punti spaziali
        rhs0                                   // equazione
    );

    // 5. callback, esportiamo i dati passo passo mentre risolve
    std::ofstream out_ode0("fluenza_terma_z_small.dat");
    std::ofstream out_ode("fluenza_prim_scat_small.dat");
    
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
    std::cout << "Propagazione completata con small step" << std::endl;
    std::cout << "Dati salvati fluenza_terma_small.dat" << std::endl;
    std::cout << "Dati salvati fluenza_prim_scat_small.dat" << std::endl;

    std::cout << "Errore Assoluto Massimo: " << max_abs_err_ode << " (z = " <<
        worst_z_abs << " cm)" << std::endl;
    std::cout << "Errore Relativo Massimo: " << max_rel_err_ode * 100.0 << " % (z = " 
        << worst_z_rel << " cm)" << std::endl;

////////////////////////////////////////////////////////////////////////////////////
/// INIZIO CALCOLO DOSE
////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n" << std::endl;

    // 1. definizione del kernel di convoluzione 
    double k_fwd_el = 4.46; 
    double k_bwd_el = 15.60; 

    double k_s_fwd = 0.7; 
    double k_s_bwd = 2.9; 

    double R_csda_primari = Rcsda_raw_loglog(E_beam);
    double R_csda_scatter = Rcsda_raw_loglog(E_beam_sc);

    // kernel elettronico normalizzato per una data energia, usando Rcsda
    // usato per dose direttamente depositata da elettroni
    struct kernelExp {
        double A, a_fwd, a_bwd;
        double operator()(double distanza) const {
            if (distanza>=0.0) return A * std::exp(-a_fwd * distanza);
            else               return A * std::exp( a_bwd * distanza);
        }
    };

    auto make_kernel_elettroni = [&](double R) {
        // parametro exp dato da reale cammino fisico della particella + costante di peso
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
    std::cout << "\n" << std::endl;

    // kernel scatter fotonico, basato sul libero cammino medio 1/mu
    double lambda_sc = 1.0 / mu_const_sc;

    // parametro exp dato da reale cammino fisico della particella + costante di peso
    double a_s_fwd = k_s_fwd / lambda_sc;
    double a_s_bwd = k_s_bwd / lambda_sc;
    double A_sc = 1.0 / ((1.0 / a_s_fwd) + (1.0 / a_s_bwd));

    // pesi energetici per componente scatter e primario
    double w_el = mu_en_const / mu_const;
    double w_sc = 1.0 - w_el;

    // kernel TERMA elettornico + scatter, unico kernel per il terma complessivo
    auto kernel_terma_fisico = [=](double distanza) {
        double term_el = w_el * kernel_el_primari(distanza);
        double term_sc = 0.0;
        if (distanza >= 0.0) term_sc = (w_sc * A_sc) * std::exp(-a_s_fwd * distanza);
        else                 term_sc = (w_sc * A_sc) * std::exp(a_s_bwd * distanza);
        return term_el + term_sc;
    };

    double k_s2_fwd = 2.8;
    double k_s2_bwd = 4.8;

    // kernel specifico per modellizzare dose fotoni scatterati di 
    // generazione successiva (doppio scatter o oltre)
    // viene modellata direttamente dose depositata da fotoni scatterati
    // e dose invece rilasciata lontano da ulteriore scatter
    auto build_kernel_terma_scatter_gen2 = [&]() {
        double w_el_s = mu_en_const_sc / mu_const_sc;
        double w_sc_s = 1.0 - w_el_s;
        double lambda_s2 = 1.0 / mu_const_sc; 
        double a_s2_fwd = k_s2_fwd / lambda_s2, a_s2_bwd = k_s2_bwd / lambda_s2;
        double A_s2 = 1.0 / (1.0/a_s2_fwd + 1.0/a_s2_bwd);

        return [=](double d) {
            double el = w_el_s * kernel_el_scatter(d); 
            double sc = (d>=0.0) ? (w_sc_s*A_s2)*std::exp(-a_s2_fwd*d) : (w_sc_s*A_s2)*std::exp(a_s2_bwd*d);
            return el + sc;
        };
    };
    auto kernel_terma_gen2 = build_kernel_terma_scatter_gen2();

   // 2. definizione terma e kerma continuo - creazione funzione raw
   // per kerma (da singola eq diff), terma (considerando unicamente decadimento fotoni),
   // kerma (doppia ode), flusso primario e scatter.
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

    // 4. definizione funzione continua - funzione finale liscia
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
    std::cout << "KERMA (from Primario e Scatter) discretizzato spazialmente su Grid1D." << std::endl;

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
        // fotoni terziari, la seconda generazione, allora approx con tutta
        // energia assorbita
    };

    auto kerma_s_continuo_en = [&](double z_req) {
        return mu_en_const_sc * E_beam_sc * phi_s_continuo(z_req);
    };

    std::cout << "\n" << std::endl;

    // 5. integrale di gauss-kronrod
    double toll_rel = 1.0e-5;
    double toll_abs = 1.0e-8;
    double z_max = z_history.back(); 
 
    // primo approccio base, per verificare correttezza calcolo
    std::ofstream out_dose_riemann("profilo_dose_small.dat");
    std::vector<double> dose_history_riemann;

    std::cout << "Calcolo dose con Somma di Riemann, small step." << std::endl;

    std::vector<double> z_riemann_check, dose_riemann_check;
    for (double z = 0.0; z <= 15.0; z += 0.1) {
       double dose_z = 0.0;
       double dz_fine = 0.01;

        // somma di riemann, integro su tutti gli z'
        for (double z_prime = 0.0; z_prime <= z_max; z_prime += dz_fine) {
            dose_z += kerma_continuo(z_prime) * kernel_el_primari(z - z_prime) * dz_fine;
        }
        z_riemann_check.push_back(z);
        dose_riemann_check.push_back(dose_z);

        out_dose_riemann << z << " " << kerma_continuo(z) << " " << dose_z << "\n";
    }
    out_dose_riemann.close();

    std::cout << "Calcolo dose con quadratura Gauss-Kronrod, small step." << std::endl;
    
    std::ofstream out_dose("profilo_dose_gk_small.dat");
    std::vector<double> z_check, dose_k_check, dose_t_check, dose_kt_check, dose_t_tot_check, dose_kt_en_check;
    
    for (double z=0.0; z<=15.0; z+=0.1) { // integro su 15 cm ignorando 
                                          // comportamenti anomali a fine dominio
        // definisco integranda per lo specifico z
        // Dose per Kerma calcolato da singola eq diff
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

        // Dose per Terma calcolato da flusso fotoni decadimento exp
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

        // Dose per Terma, usando flusso fotoni primari e secondari
        auto integranda_terma = [&](double z_prime) {
            double kerma_p_prime = mu_en_const * E_beam * phi_p_continuo(z_prime);
            // solo elettroni
            double terma_s_prime = mu_const_sc * E_beam_sc * phi_s_continuo(z_prime); 
            // TERMA (mu totale) dello scatter

            return kerma_p_prime * kernel_el_primari(z - z_prime)        
                   + terma_s_prime * kernel_terma_gen2(z - z_prime);     
           // dose da elettroni primari e fotoni scatterati calc separatamente
           // maggior realismo e adesione alla realtà 
        };
        double dose_fwd_terma = (z > 0.0)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_terma, 0.0, z, toll_rel, toll_abs)
            : 0.0;
        double dose_bwd_terma = (z < z_max)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_terma, z, z_max, toll_rel, toll_abs)
            : 0.0;
        double dose_z_terma = dose_fwd_terma + dose_bwd_terma;

        // calcolo della dose da Kerma, ma considerando per il contributo scatter solo 
        // energia rilasciata da elettroni prodotti da fotoni scatterati, ulteriore 
        // scatter eliminato e perdita conseguente di dose
        // dose persa calcolata come energia di terza generazione in fondo.
        auto integranda_kt_en = [&](double z_prime) {
            return kerma_p_continuo(z_prime) * kernel_el_primari(z - z_prime)
                + kerma_s_continuo_en(z_prime) * kernel_el_scatter(z - z_prime);
        };
        double dose_fwd_kt_en = (z > 0.0)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_kt_en, 0.0, z, toll_rel, toll_abs)
            : 0.0;
        double dose_bwd_kt_en = (z < z_max)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_kt_en, z, z_max, toll_rel, toll_abs)
            : 0.0;
        double dose_z_kt_en = dose_fwd_kt_en + dose_bwd_kt_en;

        // dose da Kerma considerando flusso fotoni primario e scatter, fotoni scatter
        // depositano tutta l'energia localmente, non possibilità ulteriore scatter
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

        // salvataggio valori calcolati per ogni z
        z_check.push_back(z);
        dose_k_check.push_back(dose_z_k);
        dose_t_check.push_back(dose_z_t);
        dose_kt_check.push_back(dose_z_kt);
        dose_t_tot_check.push_back(dose_z_terma);
        dose_kt_en_check.push_back(dose_z_kt_en);

        // salvataggio valori in file apposito per plot
        out_dose << z << " " << terma_continuo(z) << " " << kerma_continuo(z) << 
           " " << dose_z_t << " " << dose_z_k << " " << 
           (kerma_p_continuo(z) + kerma_s_continuo(z)) << " " << dose_z_kt << 
           " " << dose_z_terma << " " << dose_z_kt_en << "\n";
    }
    out_dose.close();

    std::cout << "Calcolo dose completato." << std::endl;
    std::cout << "Dati salvati in profilo_dose_small.dat e profilo_dose_gk_small.dat" 
        << std::endl; 

    std::cout << "\n" << std::endl;

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
    std::cout << "Riferimento, K_col totale (primari+scatter): "
        << integra_trapezi(z_check, kcol_tot_check) << std::endl;
    std::cout << "Riemann (K_col primari, kernel elettronico): "
        << integra_trapezi(z_riemann_check, dose_riemann_check) << std::endl;
    std::cout << "GK K_col (K_col primari, kernel elettronico): "
        << integra_trapezi(z_check, dose_k_check) << std::endl;
    std::cout << "GK TERMA (TERMA, w_sc su primario): "
        << integra_trapezi(z_check, dose_t_check) << std::endl;
    std::cout << "GK TERMA (TERMA, phi_p e phi_s separati): "
        << integra_trapezi(z_check, dose_t_tot_check) << std::endl;
    std::cout << "GK K_tot (K_col p.+s., kernel separati): "
        << integra_trapezi(z_check, dose_kt_check) << std::endl;
    std::cout << "GK K_tot (K_col p.+s., kernel separati, perdita 2 gen): "
        << integra_trapezi(z_check, dose_kt_en_check) << std::endl;

    std::vector<double> z_check_full, dose_kt_check_full, dose_terma_v3_check_full;
    for (double z = 0.0; z <= z_max; z += 0.1) {
    // stessa logica di calcolo di dose_z_kt e dose_z_terma, ma su dominio esteso
        auto integranda_terma = [&](double z_prime) {
            double kerma_p_prime = mu_en_const * E_beam * phi_p_continuo(z_prime);
            double terma_s_prime = mu_const_sc * E_beam_sc * phi_s_continuo(z_prime); 
            return kerma_p_prime * kernel_el_primari(z - z_prime)        
                   + terma_s_prime * kernel_terma_gen2(z - z_prime);        
        };
        double dose_fwd_terma = (z > 0.0)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_terma, 0.0, z, toll_rel, toll_abs)
            : 0.0;
        double dose_bwd_terma = (z < z_max)
            ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(
                    integranda_terma, z, z_max, toll_rel, toll_abs)
            : 0.0;
        double dose_z_terma = dose_fwd_terma + dose_bwd_terma;


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

        z_check_full.push_back(z);
        dose_kt_check_full.push_back(dose_z_kt);
        dose_terma_v3_check_full.push_back(dose_z_terma);
    }
    std::cout << "GK K_tot (dominio esteso): " 
        << integra_trapezi(z_check_full, dose_kt_check_full) << std::endl;
    std::cout << "GK TERMA (dominio esteso): " 
        << integra_trapezi(z_check_full, dose_terma_v3_check_full) << std::endl;



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

    auto Dose_kcol_z = [&](double z) {
        auto integ = [&](double zp) { return kerma_continuo(zp) * kernel_el_primari(z - zp); };
        double fwd = (z > 0.0)    ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ, 0.0, z, 1e-5, 1e-8) : 0.0;
        double bwd = (z < z_max)  ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ, z, z_max, 1e-5, 1e-8) : 0.0;
        return fwd + bwd;
    };

    double integrale_outer_GK = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(Dose_kcol_z, 0.0, 15.0, 1e-5, 1e-8);

    std::cout << "GK annidata (outer+inner): " << integrale_outer_GK << std::endl;
    std::cout << "Trapezi su griglia 0.1cm:  " << integra_trapezi(z_check, dose_k_check) 
        << std::endl;

////////////////////////////////////////////////////////////////////////////////////
/// FIT SU VALORI PARAMETRI LIBERI
////////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n" << std::endl;
    // fit sul parametro a_s_fwd e _bwd, dovuti alla propagazione dei fotoni scatter
    // l'idea è di associare stessa forma di dose dovuta a due calcoli e modelli 
    // simili ma diversi 
    auto build_kernel_terma = [&](double ksf, double ksb) {
        double a_sf = ksf / lambda_sc, a_sb = ksb / lambda_sc;
        double A_sc_loc = 1.0 / ((1.0 / a_sf) + (1.0 / a_sb));
        return [=](double d) {
            double el = w_el * kernel_el_primari(d);
            double sc = (d >= 0.0) ? (w_sc * A_sc_loc) * std::exp(-a_sf * d)
                                   : (w_sc * A_sc_loc) * std::exp( a_sb * d);
            return el + sc;
        };
    };

    auto costo_forma = [&](double ksf, double ksb) {
        auto ker = build_kernel_terma(ksf, ksb);
        double s2 = 0.0;
        for (size_t i=0.0; i < z_check.size(); i++) {
            double z = z_check[i];
            auto integ = [&](double zp) {return terma_continuo(zp) * ker(z-zp);};
            double d_fwd = (z>0.0) ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ,0.0,z,1e-4,1e-6) : 0.0;
            double d_bwd = (z<z_max) ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ,z,z_max,1e-4,1e-6) : 0.0;
            double scarto = (d_fwd + d_bwd) - dose_kt_check[i];
            s2 += scarto * scarto; 
        }
        return s2;
    };

    // grid search grezza
    double best_ksf = k_s_fwd, best_ksb = k_s_bwd, best_c = 1e18;
    for (double ksf = 0.1; ksf <= 2.0; ksf += 0.1)
        for (double ksb = 0.1; ksb <= 3.0; ksb += 0.1) {
            double c = costo_forma(ksf, ksb);
            if (c < best_c) {best_c = c; best_ksf = ksf; best_ksb = ksb;} 
        }

    std::cout << "k_s_fwd ottimale: " << best_ksf << ", k_s_bwd ottimale: " <<
        best_ksb << std::endl;

 
    auto build_kernel_terma_gen2 = [&](double ks2f, double ks2b) {
        double w_el_s = mu_en_const_sc / mu_const_sc;
        double w_sc_s = 1.0 - w_el_s;
        double lambda_s2 = 1.0 / mu_const_sc;
        double a2f = ks2f / lambda_s2, a2b = ks2b / lambda_s2;
        double A2 = 1.0 / (1.0/a2f + 1.0/a2b);
        return [=](double d) {
            double el = w_el_s * kernel_el_scatter(d);
            double sc = (d>=0.0) ? (w_sc_s*A2)*std::exp(-a2f*d) 
                                 : (w_sc_s*A2)*std::exp(a2b*d);
            return el + sc;
        };
    };

    double integrale_dose_kt = integra_trapezi(z_check, dose_kt_check);

    auto costo_forma_tps = [&](double ks2f, double ks2b) {
        auto ker2 = build_kernel_terma_gen2(ks2f, ks2b);
        std::vector<double> dose_t_ps(z_check.size());
        for (size_t i = 0; i < z_check.size(); i++) {
            double z = z_check[i];
            auto integ = [&](double zp) {
                double kerma_p_prime = mu_en_const * E_beam * phi_p_continuo(zp);
                double terma_s_prime = mu_const_sc * E_beam_sc * phi_s_continuo(zp);
                return kerma_p_prime * kernel_el_primari(z - zp) 
                    + terma_s_prime * ker2(z - zp);
            };
            double fwd = (z>0.0)   ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ,0.0,z,1e-4,1e-6) : 0.0;
            double bwd = (z<z_max) ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ,z,z_max,1e-4,1e-6) : 0.0;
        dose_t_ps[i] = fwd + bwd;
        }
        double int_t_ps = integra_trapezi(z_check, dose_t_ps);
        double s2 = 0.0;
        for (size_t i = 0; i < z_check.size(); i++) {
            double forma_t_ps = dose_t_ps[i] / int_t_ps;
            double forma_kt = dose_kt_check[i] / integrale_dose_kt;
            double scarto = forma_t_ps - forma_kt;
            s2 += scarto * scarto;
        }
        return s2;
    };

    double best_ks2f = 1.0, best_ks2b = 1.0, best_c3 = 1e18;
    for (double ks2f = 0.2; ks2f <= 3.0; ks2f += 0.2)
        for (double ks2b = 0.2; ks2b <= 5.0; ks2b += 0.2) {
            double c = costo_forma_tps(ks2f, ks2b);
            if (c < best_c3) { best_c3 = c; best_ks2f = ks2f; best_ks2b = ks2b; }
    }
    std::cout << "k_s2_fwd (TERMA flusso Primario Scatter) ottimale: " << best_ks2f
        << ", k_s2_bwd (TERMA flusso Primario Scatter) ottimale: " << best_ks2b << std::endl;

    // fit su costante kernel propagazione forward dose elettronica usando come riferimento
    // distanza deposizione dose elettrone in acqua (valori)
    double me_c2 = 0.511; // MeV
    double alpha_edge = E_beam / me_c2;
    double E_max_electron = E_beam * (2.0*alpha_edge) / (1.0 + 2.0*alpha_edge); 
    // edge Compton, θ=π

    double R_max_electron = Rcsda_raw_loglog(E_max_electron);

    // criterio pratico: la coda in avanti del kernel deve annullarsi (≈95-99%
    // dell'energia gia' depositata) proprio a questa profondita'
    // elettroni si assorbono dopo range, integrale del flusso (exp neg) tra 0 e
    // il range deve essere circa 0 (\sim 1%) -> svolgendo il calcolo si ottiene exp 
    // che decade dopo dist = range con parametro calibrato
    double soglia = 0.96;  
    double k_fwd_el_calibrato = -std::log(1.0 - soglia) / R_max_electron; 

    // Derivazione empirico-fisica del parametro bwd
    // Gli elettroni backscatterati hanno molta meno energia,
    // quindi si fermano in uno spazio circa 3-4 volte più breve.
    double rapporto_bwd_fwd = 3.5;
    double k_bwd_el_calibrato = k_fwd_el_calibrato * rapporto_bwd_fwd;

    std::cout << "Calibrazione Kernel Elettronico (CSDA):" << std::endl;
    std::cout << " - k_fwd: " << k_fwd_el_calibrato << " cm^-1" << std::endl;
    std::cout << " - k_bwd: " << k_bwd_el_calibrato << " cm^-1" << std::endl;

    return 0;
}
