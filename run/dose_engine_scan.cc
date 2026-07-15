#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "Interpolation/interpolation.hh"
#include "runge_kutta.hh"

////////////////////////////////////////////////////////////////////////////////////
/// STRUTTURA RISULTATI: tutto quello che vogliamo confrontare fra energie diverse
////////////////////////////////////////////////////////////////////////////////////

struct RisultatiEnergia {
    double E_beam = 0.0, E_beam_sc = 0.0, f_fwd = 0.0;
    double mu = 0.0, mu_en = 0.0, mu_sc = 0.0, mu_en_sc = 0.0;
    double R_csda_p = 0.0, R_csda_s = 0.0;
    double k_fwd_el = 0.0, k_bwd_el = 0.0;
    double k_s2_fwd = 0.0, k_s2_bwd = 0.0;   // eventualmente dal fit, o riusati
    double d_max = 0.0;                       // profondita' del massimo di dose K_tot
    double integrale_riferimento = 0.0;       // K_col totale (p.+s.), nessun kernel
    double integrale_kcol = 0.0;              // K_col solo primari + kernel
    double integrale_terma = 0.0;             // TERMA con phi_p/phi_s separati
    double integrale_ktot = 0.0;              // K_tot da ODE, kernel separati
};

////////////////////////////////////////////////////////////////////////////////////
/// PIPELINE PARAMETRIZZATA: tutto quello che dipende da E_beam
///
/// mu_continuo/mu_en_continuo/Rcsda_raw_loglog sono passate come parametri
/// (dipendono solo dai file NIST/ESTAR, non da E_beam)
////////////////////////////////////////////////////////////////////////////////////

template <typename MuFunc, typename MuEnFunc, typename RcsdaFunc>
RisultatiEnergia esegui_pipeline(double E_beam, bool fit_completo,
                                  double k_s2_fwd_iniziale, double k_s2_bwd_iniziale,
                                  MuFunc&& mu_continuo, MuEnFunc&& mu_en_continuo,
                                  RcsdaFunc&& Rcsda_raw_loglog)
{
    RisultatiEnergia r{};
    r.E_beam = E_beam;

    ////////////////////////////////////////////////////////////////////////
    // 1. Klein-Nishina: E_beam_sc e f_fwd, derivati
    ////////////////////////////////////////////////////////////////////////
    double m_e_c2 = 0.511; // MeV
    double alpha = E_beam / m_e_c2;

    auto E_prime = [&](double theta) {
        return E_beam / (1.0 + alpha * (1.0 - std::cos(theta)));
    };
    auto dsigma_dOmega = [&](double theta) {
        double Ep = E_prime(theta);
        double rr = Ep / E_beam;
        double sin2 = std::sin(theta) * std::sin(theta);
        return rr*rr * (rr + 1.0/rr - sin2);
    };
    auto integranda_energia = [&](double theta) { return E_prime(theta) * dsigma_dOmega(theta) * std::sin(theta); };
    auto integranda_norm    = [&](double theta) { return dsigma_dOmega(theta) * std::sin(theta); };

    double tol = 1e-8;
    double num_E   = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_energia, 0.0, M_PI, tol, tol);
    double den_E   = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_norm,    0.0, M_PI, tol, tol);
    double num_fwd = Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_energia, 0.0, M_PI/2.0, tol, tol);

    r.E_beam_sc = num_E / den_E;
    r.f_fwd     = num_fwd / num_E;

    r.mu       = mu_continuo(E_beam);
    r.mu_en    = mu_en_continuo(E_beam);
    r.mu_sc    = mu_continuo(r.E_beam_sc);
    r.mu_en_sc = mu_en_continuo(r.E_beam_sc);

    double mu_scatt_create = r.mu - r.mu_en;

    ////////////////////////////////////////////////////////////////////////
    // 2. ODE a due stadi: phi_p, phi_s lungo z (stessa struttura del file originale)
    ////////////////////////////////////////////////////////////////////////
    rk::vd<2> phi_init;
    phi_init[0] = 1.0;
    phi_init[1] = 0.0;

    double E_beam_sc = r.E_beam_sc; // per catturare per valore nella rhs
    double f_fwd = r.f_fwd;
    double mu_const = r.mu, mu_const_sc = r.mu_sc;

    rk::rk_rhs_t<rk::vd<2>> rhs = [mu_const, mu_const_sc, mu_scatt_create, f_fwd,
        E_beam, E_beam_sc](double /*z*/, rk::vd<2>& state) {
        double phi_p = state(0);
        double phi_s = state(1);
        double dp_dz = -mu_const * phi_p;
        double ds_dz = (mu_scatt_create * (E_beam / E_beam_sc) * f_fwd) * phi_p - mu_const_sc * phi_s;
        state[0] = dp_dz;
        state[1] = ds_dz;
    };

    std::vector<double> z_points;
    for (int i = 0; i <= 50; i++) z_points.push_back(i * 0.5); // 0-25 cm, passo 0.5 cm
    rk::TimeInfo spatial_info(z_points);

    rk::RungeKutta<rk::vd<2>, 4> solver(rk::PreImplementedTableau::RKOriginal, phi_init, spatial_info, rhs);

    std::vector<double> z_history, phi_p_history, phi_s_history;

    rk::rk_callback_t<rk::vd<2>> callback = [&](double z, const rk::TimeInfo&, const rk::vd<2>& phi) {
        z_history.push_back(z);
        phi_p_history.push_back(phi(0));
        phi_s_history.push_back(phi(1));
    };
    solver.AddCallback(callback);
    solver.CallbackOnlyOnTimeStamp();
    solver();

    double z_max = z_history.back();

    ////////////////////////////////////////////////////////////////////////
    // 3. Interpolazione continua di phi_p, phi_s (raw lineare + grid Chebyshev)
    ////////////////////////////////////////////////////////////////////////
    auto phi_p_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return phi_p_history.front();
        if (it == z_history.end())   return phi_p_history.back();
        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double p0 = phi_p_history[i-1], p1 = phi_p_history[i];
        return p0 + (p1 - p0) * (z_req - z0) / (z1 - z0);
    };
    auto phi_s_raw_linear = [&](double z_req) {
        auto it = std::lower_bound(z_history.begin(), z_history.end(), z_req);
        if (it == z_history.begin()) return phi_s_history.front();
        if (it == z_history.end())   return phi_s_history.back();
        int i = std::distance(z_history.begin(), it);
        double z0 = z_history[i-1], z1 = z_history[i];
        double p0 = phi_s_history[i-1], p1 = phi_s_history[i];
        return p0 + (p1 - p0) * (z_req - z0) / (z1 - z0);
    };

    std::vector<double> tasselli_z = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0};
    std::vector<size_t> gradi_z = {15, 15, 15, 15, 15};
    Interpolation::Grid1D grid_terma(
        Interpolation::make_discretization_info<Interpolation::details::identity_maps>(tasselli_z, gradi_z));

    std::vector<double> fj_phi_p = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma, phi_p_raw_linear, [](size_t n) {return std::vector<double>(n, 0.);});
    std::vector<double> fj_phi_s = Interpolation::Discretize<std::vector<double>, double>(
        grid_terma, phi_s_raw_linear, [](size_t n) {return std::vector<double>(n, 0.);});

    auto phi_p_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(z_req, fj_phi_p, []() -> double {return 0.;});
    };
    auto phi_s_continuo = [&](double z_req) {
        return grid_terma.interpolate<double, std::vector<double>>(z_req, fj_phi_s, []() -> double {return 0.;});
    };

    ////////////////////////////////////////////////////////////////////////
    // 4. Kernel elettronici, calibrati da edge Compton + R_csda 
    ////////////////////////////////////////////////////////////////////////
    r.R_csda_p = Rcsda_raw_loglog(E_beam);
    r.R_csda_s = Rcsda_raw_loglog(r.E_beam_sc);

    double alpha_edge = E_beam / m_e_c2;
    double E_max_electron = E_beam * (2.0*alpha_edge) / (1.0 + 2.0*alpha_edge);
    double R_max_electron = Rcsda_raw_loglog(E_max_electron);

    double soglia = 0.96;
    r.k_fwd_el = -std::log(1.0 - soglia);  // costante di forma ADIMENSIONALE: a_fwd = k_fwd_el / R(E)
    double rapporto_bwd_fwd = 3.5;
    r.k_bwd_el = r.k_fwd_el * rapporto_bwd_fwd;

    struct KernelExp {
        double A, a_fwd, a_bwd;
        double operator()(double distanza) const {
            if (distanza >= 0.0) return A * std::exp(-a_fwd * distanza);
            else                 return A * std::exp( a_bwd * distanza);
        }
    };
    auto make_kernel_elettroni = [&](double R) {
        double a_fwd = r.k_fwd_el / R;
        double a_bwd = r.k_bwd_el / R;
        double A = 1.0 / (1.0/a_fwd + 1.0/a_bwd);
        return KernelExp{A, a_fwd, a_bwd};
    };

    // uso R_max_electron (range dell'elettrone all'edge Compton) 
    // per il kernel primario, non R_csda_p = Rcsda(E_beam),
    // e' la scelta coerente col criterio di calibrazione
    // (soglia sulla coda del kernel proprio a quella profondita').
    KernelExp kernel_el_primari = make_kernel_elettroni(R_max_electron);
    KernelExp kernel_el_scatter = make_kernel_elettroni(r.R_csda_s);

    std::cout << "E=" << E_beam << "  R_max_electron=" << R_max_electron
        << " cm  a_fwd=" << (r.k_fwd_el / R_max_electron) << " cm^-1"
        << "  1/a_fwd=" << (R_max_electron / r.k_fwd_el) << " cm (larghezza kernel)" << std::endl;

    auto kerma_p_continuo = [&](double z_req) { return r.mu_en * E_beam * phi_p_continuo(z_req); };
    auto kerma_s_continuo = [&](double z_req) { return r.mu_sc * r.E_beam_sc * phi_s_continuo(z_req); };

    ////////////////////////////////////////////////////////////////////////
    // 5. K_tot (ODE, kernel separati) -- serve come bersaglio per il fit e per il d_max
    ////////////////////////////////////////////////////////////////////////
    double toll_rel = 1.0e-5, toll_abs = 1.0e-8;

    std::vector<double> z_check, dose_kt_check, dose_kcol_check;
    for (double z = 0.0; z <= 15.0; z += 0.1) {
        auto integranda_kt = [&](double zp) {
            return kerma_p_continuo(zp) * kernel_el_primari(z - zp)
                 + kerma_s_continuo(zp) * kernel_el_scatter(z - zp);
        };
        double fwd = (z > 0.0)   ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_kt, 0.0, z, toll_rel, toll_abs) : 0.0;
        double bwd = (z < z_max) ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_kt, z, z_max, toll_rel, toll_abs) : 0.0;
        double dose_z_kt = fwd + bwd;

        auto integranda_k = [&](double zp) {
            return (r.mu_en * E_beam * phi_p_continuo(zp)) * kernel_el_primari(z - zp);
        };
        double fwd_k = (z > 0.0)   ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_k, 0.0, z, toll_rel, toll_abs) : 0.0;
        double bwd_k = (z < z_max) ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integranda_k, z, z_max, toll_rel, toll_abs) : 0.0;

        z_check.push_back(z);
        dose_kt_check.push_back(dose_z_kt);
        dose_kcol_check.push_back(fwd_k + bwd_k);
    }

    auto integra_trapezi = [](const std::vector<double>& z, const std::vector<double>& f) {
        double s = 0.0;
        for (size_t i = 0; i + 1 < z.size(); i++) s += 0.5 * (f[i] + f[i+1]) * (z[i+1] - z[i]);
        return s;
    };

    double integrale_dose_kt = integra_trapezi(z_check, dose_kt_check);
    r.integrale_ktot = integrale_dose_kt;
    r.integrale_kcol = integra_trapezi(z_check, dose_kcol_check);

    // d_max: punto di massimo del profilo K_tot (il piu' affidabile spazialmente)
    size_t idx_max = std::distance(dose_kt_check.begin(),
        std::max_element(dose_kt_check.begin(), dose_kt_check.end()));
    r.d_max = z_check[idx_max];

    // riferimento: K_col totale senza alcun kernel (energia grezza disponibile)
    std::vector<double> kcol_tot_check;
    for (double z : z_check) kcol_tot_check.push_back(kerma_p_continuo(z) + kerma_s_continuo(z));
    r.integrale_riferimento = integra_trapezi(z_check, kcol_tot_check);

    ////////////////////////////////////////////////////////////////////////
    // 6. Kernel gen2 (componente scatter del TERMA): fit oppure riuso
    ////////////////////////////////////////////////////////////////////////
    auto build_kernel_terma_gen2 = [&](double ks2f, double ks2b) {
        double w_el_s = r.mu_en_sc / r.mu_sc;
        double w_sc_s = 1.0 - w_el_s;
        double lambda_s2 = 1.0 / r.mu_sc;
        double a2f = ks2f / lambda_s2, a2b = ks2b / lambda_s2;
        double A2 = 1.0 / (1.0/a2f + 1.0/a2b);
        return [=](double d) {
            double el = w_el_s * kernel_el_scatter(d);
            double sc = (d >= 0.0) ? (w_sc_s*A2)*std::exp(-a2f*d) : (w_sc_s*A2)*std::exp(a2b*d);
            return el + sc;
        };
    };

    auto dose_terma_con_kernel = [&](double ks2f, double ks2b, std::vector<double>* out_profilo) {
        auto ker2 = build_kernel_terma_gen2(ks2f, ks2b);
        std::vector<double> dose_t(z_check.size());
        for (size_t i = 0; i < z_check.size(); i++) {
            double z = z_check[i];
            auto integ = [&](double zp) {
                double kerma_p_prime = r.mu_en * E_beam * phi_p_continuo(zp);
                double terma_s_prime = r.mu_sc * r.E_beam_sc * phi_s_continuo(zp);
                return kerma_p_prime * kernel_el_primari(z - zp) + terma_s_prime * ker2(z - zp);
            };
            double fwd = (z > 0.0)   ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ, 0.0, z, 1e-4, 1e-6) : 0.0;
            double bwd = (z < z_max) ? Interpolation::GaussKronrod<Interpolation::GK_61>::integrate(integ, z, z_max, 1e-4, 1e-6) : 0.0;
            dose_t[i] = fwd + bwd;
        }
        if (out_profilo) *out_profilo = dose_t;
        return integra_trapezi(z_check, dose_t);
    };

    if (fit_completo) {
        auto costo_forma = [&](double ks2f, double ks2b) {
            std::vector<double> dose_t;
            double int_t = dose_terma_con_kernel(ks2f, ks2b, &dose_t);
            double s2 = 0.0;
            for (size_t i = 0; i < z_check.size(); i++) {
                double forma_t = dose_t[i] / int_t;
                double forma_kt = dose_kt_check[i] / integrale_dose_kt;
                double scarto = forma_t - forma_kt;
                s2 += scarto * scarto;
            }
            return s2;
        };
        double best_ks2f = k_s2_fwd_iniziale, best_ks2b = k_s2_bwd_iniziale, best_c = 1e18;
        for (double ks2f = 0.2; ks2f <= 3.0; ks2f += 0.2)
            for (double ks2b = 0.2; ks2b <= 5.0; ks2b += 0.2) {
                double c = costo_forma(ks2f, ks2b);
                if (c < best_c) { best_c = c; best_ks2f = ks2f; best_ks2b = ks2b; }
            }
        r.k_s2_fwd = best_ks2f;
        r.k_s2_bwd = best_ks2b;
    } else {
        r.k_s2_fwd = k_s2_fwd_iniziale;
        r.k_s2_bwd = k_s2_bwd_iniziale;
    }

    r.integrale_terma = dose_terma_con_kernel(r.k_s2_fwd, r.k_s2_bwd, nullptr);

    return r;
}

////////////////////////////////////////////////////////////////////////////////////
/// MAIN: lettura NIST una sola volta, poi scan su piu' energie
////////////////////////////////////////////////////////////////////////////////////

int main()
{
    // --- lettura dati NIST/ESTAR, una sola volta, indipendente da E_beam ---
    std::vector<double> E_data, mu_data, mu_en_data;
    std::ifstream file("nist_water.txt");
    if (!file.is_open()) { std::cerr << "Errore: nist_water.txt non trovato" << std::endl; return 1; }
    double E_val, mu_val, mu_en_val;
    while (file >> E_val >> mu_val >> mu_en_val) {
        E_data.push_back(E_val); mu_data.push_back(mu_val); mu_en_data.push_back(mu_en_val);
    }
    file.close();

    std::vector<double> E_csda_data, Rcsda_data;
    std::ifstream file_csda("estar_water.txt");
    if (!file_csda.is_open()) { std::cerr << "Errore: estar_water.txt non trovato" << std::endl; return 1; }
    double Ec, Rc;
    while (file_csda >> Ec >> Rc) { E_csda_data.push_back(Ec); Rcsda_data.push_back(Rc); }
    file_csda.close();

    std::cout << "Letti " << E_data.size() << " punti NIST, " << E_csda_data.size() << " punti ESTAR." << std::endl;

    auto mu_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_data.begin(), E_data.end(), E_req);
        if (it == E_data.begin()) return mu_data.front();
        if (it == E_data.end())   return mu_data.back();
        int i = std::distance(E_data.begin(), it);
        double E1 = E_data[i-1], E2 = E_data[i], m1 = mu_data[i-1], m2 = mu_data[i];
        double lm = std::log(m1) + (std::log(m2)-std::log(m1)) * (std::log(E_req)-std::log(E1)) / (std::log(E2)-std::log(E1));
        return std::exp(lm);
    };
    auto mu_en_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_data.begin(), E_data.end(), E_req);
        if (it == E_data.begin()) return mu_en_data.front();
        if (it == E_data.end())   return mu_en_data.back();
        int i = std::distance(E_data.begin(), it);
        double E1 = E_data[i-1], E2 = E_data[i], m1 = mu_en_data[i-1], m2 = mu_en_data[i];
        double lm = std::log(m1) + (std::log(m2)-std::log(m1)) * (std::log(E_req)-std::log(E1)) / (std::log(E2)-std::log(E1));
        return std::exp(lm);
    };
    auto Rcsda_raw_loglog = [&](double E_req) {
        auto it = std::lower_bound(E_csda_data.begin(), E_csda_data.end(), E_req);
        if (it == E_csda_data.begin()) return Rcsda_data.front();
        if (it == E_csda_data.end())   return Rcsda_data.back();
        int i = std::distance(E_csda_data.begin(), it);
        double E1 = E_csda_data[i-1], E2 = E_csda_data[i], R1 = Rcsda_data[i-1], R2 = Rcsda_data[i];
        double lR = std::log(R1) + (std::log(R2)-std::log(R1)) * (std::log(E_req)-std::log(E1)) / (std::log(E2)-std::log(E1));
        return std::exp(lR);
    };

    // grid Chebyshev per mu/mu_en
    std::vector<double> tasselli = {E_data.front(), 0.05, 0.1, 1, 5, E_data.back()};
    std::vector<size_t> gradi = {20, 20, 20, 20, 20};
    Interpolation::Grid1D grid(
        Interpolation::make_discretization_info<Interpolation::details::log_0_maps>(tasselli, gradi));

    std::vector<double> fj = Discretize<std::vector<double>, double>(
        grid, [&](double E) { return std::log(mu_raw_loglog(E)); }, [](size_t n) {return std::vector<double>(n, 0.);});
    std::vector<double> fj_en = Discretize<std::vector<double>, double>(
        grid, [&](double E) { return std::log(mu_en_raw_loglog(E)); }, [](size_t n) {return std::vector<double>(n, 0.);});

    auto mu_continuo = [&](double E_req) {
        return std::exp(grid.interpolate<double, std::vector<double>>(E_req, fj, []() -> double {return 0.;}));
    };
    auto mu_en_continuo = [&](double E_req) {
        return std::exp(grid.interpolate<double, std::vector<double>>(E_req, fj_en, []() -> double {return 0.;}));
    };

    // --- passo 1: scan rapido su molte energie, fit disattivato ---
    std::vector<double> energie_scan = {0.3, 0.5, 0.849864, 1.0, 1.75, 3.0, 6.0, 10.0};
    std::vector<RisultatiEnergia> scan;
    double ks2f_riuso = 2.8, ks2b_riuso = 4.8; // ultimi valori noti, aggiornati dopo il fit sotto

    std::cout << "\n--- Scan rapido (fit disattivato, riuso ultimo k_s2 noto) ---\n";
    for (double E : energie_scan) {
        RisultatiEnergia r = esegui_pipeline(E, /*fit_completo=*/false, ks2f_riuso, ks2b_riuso,
                                              mu_continuo, mu_en_continuo, Rcsda_raw_loglog);
        scan.push_back(r);
        std::cout << "E=" << E << " MeV  d_max=" << r.d_max << " cm  Rif=" << r.integrale_riferimento
            << "  Kcol=" << r.integrale_kcol << "  TERMA=" << r.integrale_terma
            << "  Ktot=" << r.integrale_ktot << std::endl;
        std::cout << "\n" << std::endl;
    }

    // --- passo 2: fit completo solo su alcune energie di controllo (lento) ---
    std::cout << "\n--- Fit completo su energie di controllo ---\n";
    for (double E : {0.5, 1.75, 6.0}) {
        RisultatiEnergia r = esegui_pipeline(E, /*fit_completo=*/true, ks2f_riuso, ks2b_riuso,
                                              mu_continuo, mu_en_continuo, Rcsda_raw_loglog);
        std::cout << "E=" << E << " MeV -> k_s2_fwd=" << r.k_s2_fwd
            << "  k_s2_bwd=" << r.k_s2_bwd
            << "  (TERMA=" << r.integrale_terma << "  Ktot=" << r.integrale_ktot << ")" << std::endl;
        std::cout << "\n" << std::endl;
    }

    return 0;
}
