#include "Interpolation/interpolation.hh"
#include <iostream>
#include "Interpolation/ran2.hh"

struct DYGrid {
    DYGrid(size_t nOgata, size_t nXi) : nO(nOgata), nXi(nXi), data(nOgata*nXi, 0.)
    {
    }

    double &operator()(size_t i, size_t j)
    {
        return data[i * nXi + j];
    }

    const double &operator()(size_t i, size_t j) const
    {
        return data[i * nXi + j];
    } 

    size_t nO, nXi;
    std::vector<double> data;
};

struct Table {
    std::vector<double> qT;
    std::vector<double> xi_grid;
    std::vector<double> ogata_grid;

    std::vector<DYGrid> WDY;
};

Table compute(const std::vector<double> &qT, std::pair<double, double> bin, size_t nOgata)
{
    using namespace Interpolation;

    if (bin.first >= bin.second) {
       throw std::logic_error("bin.first cannot be >= bin.second");
    }

    std::vector<double> bounds;
    std::vector<size_t> sizes;

    if (bin.first < 1 && bin.second > 1) {
       bounds = {bin.first, 1., bin.second};
       sizes  = {10, 10};
    } else {
       bounds = {bin.first, bin.second};
       sizes  = {20};
    }

    Grid1D grid(make_discretization_info<details::log_0_maps>(bounds, sizes));
    Table table;
    table.qT = qT;
    table.xi_grid = grid._coord;

    Ogata::Integrator ogata_integrator;
    table.ogata_grid.resize(nOgata, 0.);

    for (size_t i=0; i<nOgata; i++) {
        table.ogata_grid[i] = ogata_integrator.nodes[i];
    }
    using integrator = GaussKonrod<GK_21>;

    for (size_t iqT = 0; iqT < qT.size(); iqT++) {
        table.WDY.emplace_back(DYGrid(nOgata, table.xi_grid.size()));

        for (size_t iO = 0; iO<nOgata; iO++) {
            const double b = ogata_integrator.nodes[iO] / qT[iqT];
            const double pref = exp(-b * b) * ogata_integrator.nodes[iO] * 
                ogata_integrator.weights[iO] / qT[iqT] / (0.1 + qT[iqT] * qT[iqT]);

            size_t j_xi = 0;
            auto integrand = [&j_xi, &grid](double u) {
                const double xi = grid._d_info.to_phys_space(u);
                const double jac = grid._d_info._to_phys_space_der(u);
                const double w = grid._weights[j_xi](u, grid.get_std_grid(j_xi));

                return pow(xi, -0.25 /*-1+0.75*/) * w * jac;
            };
            for (size_t iXi = 0; iXi < table.xi_grid.size(); iXi++) {
                j_xi = iXi;
                auto [umin, umax] = grid.get_support_weight_aj(jXi);
                dygrid(iO, grid._from_iw_to_ic[iXi]) += pref * 
                    integrator::integrate(integrand, umin, umax);

            }
        }
    }

    return table;
}

std::vector<double> convolve(const Table &tab, const std::function<double(double,double)> &fnc) 
{
    std::vector<double> result(tab.qT.size(), 0.);

    for (size_t iqT = 0; iqT < tab.qT.size(); iqT++) {
        double &r = result[iqT];

        const DYGrid &dygrid = tab.WDY[iqT];
        const double qT = tab.qT[iqT];

        for (size_t iO = 0; iO < tab.ogata_grid.size(); iO++) {
            const double b = tab.ogata_grid[iO] / qT;

        double acc = 0.;

        for (size_t iXi = 0; iXi < tab.xi_grid.size(); iXi++) {
            acc += dygrid(iO, iXi) * fnc(tab.xi_grid[iXi], b) * fnc(1. / tab    .xi_grid[iXi], b);
        }
 
        r += acc;
        if ((std::abs(acc / r) < 1.0e-5 && iO > 10) || (std::abs(r) < 1e-16     && iO > 10)) break;
        }
    }
    return result;
}

std::pair<std::vector<double>, std::vector<std::vector<double>>> convolve_der(const Table &tab, 
        const std::function<double(double, double)> &fnc, const std::function<void(double, double, 
        std::vector<double> &)> &fnc_der, size_t n_par)
{
    std::vector<std::vector<double>> result_der(tab.qT.size(), std::vector<double>(n_par, 0.));
    std::vector<double> result(tab.qT.size(), 0.);
 
    std::vector<double> cache_1(n_par, 0.);
    std::vector<double> cache_2(n_par, 0.);
 
    for (size_t iqT = 0; iqT < tab.qT.size(); iqT++) {
        const double qT = tab.qT[iqT];
                                                                             
        auto &r_der          = result_der[iqT];
        double &r            = result[iqT];
        const DYGrid &dygrid = tab.WDY[iqT];
        std::vector<int> should_stop(n_par + 1, false);
        for (size_t iO = 0; iO < tab.ogata_grid.size(); iO++) {
        const double b = tab.ogata_grid[iO] / qT;
 
        double acc = 0.;                                                   
        std::vector<double> acc_der(n_par, 0.);
 
        for (size_t iXi = 0; iXi < tab.xi_grid.size(); iXi++) {
            fnc_der(tab.xi_grid[iXi], b, cache_1);
            fnc_der(1. / tab.xi_grid[iXi], b, cache_2);
 
            double f1 = fnc(tab.xi_grid[iXi], b);
            double f2 = fnc(1. / tab.xi_grid[iXi], b);

            acc += dygrid(iO, iXi) * f1 * f2;
            for (size_t iPar = 0; iPar < n_par; iPar++) {
                acc_der[iPar] += dygrid(iO, iXi) * (cache_1[iPar] * f2 + f1 *     cache_2[iPar]);
            }
        }
 
        int stop_count = 0;
        for (size_t iPar = 0; iPar < n_par; iPar++) {
            if (!should_stop[iPar]) {
                r_der[iPar]       += acc_der[iPar];
                should_stop[iPar]  = ((std::abs(acc_der[iPar] / r_der[iPar])     < 1.0e-5 && iO > 10) 
                        || (std::abs(r_der[iPar]) < 1e-16 && iO     > 10));
            } else stop_count++;
        }
        if (!should_stop[n_par]) {
            r += acc;
            should_stop[n_par] = (std::abs(acc / r) < 1.0e-5 && iO > 10) || (std::abs(r) 
                    < 1e-16 && iO > 10);
        } else stop_count++;
        if ((size_t)stop_count == n_par + 1) break;
        }
    }
    return {result, result_der};
}     
