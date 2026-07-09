#include "Interpolation/gauss_kronrod.hh"
#include <Interpolation/grid_1d.hh>
#include <algorithm>

namespace Interpolation
{

namespace 
{
    // u = a <=> t=+1
    // u = b <=> t=-1
    // m1p1 means Minus 1 Plus 1
    // noexcept means that this function won't throw any exception, no error
    // will be in the running of the function. In order to be faster.

    double from_ab_to_m1p1(double u, double a, double b) noexcept
    {
        return -2 * (u - a) / (b - a) +1;
    };
    double from_m1p1_to_ab(double t, double a, double b) noexcept
    {
        return (b - a) * (1 - t) * 0.5 + a;
    };
    double from_ab_to_m1p1_der(double a, double b) noexcept
    {
        return -2 / (b - a);
    };

} // namespace

SingleDiscretizationInfo::SingleDiscretizationInfo(std::vector<double> inter, 
                            std::vector<size_t> g_size,
                            std::function<double(double)> to_i_space, 
                            std::function<double(double)> to_i_space_der,
                            std::function<double(double)> to_p_space,
                            std::function<double(double)> to_p_space_der)
    : intervals(inter.size()-1, {0, 0}),
      intervals_phys(inter.size()-1, {0, 0}),
      grid_sizes(g_size),
      to_inter_space(to_i_space),
      to_inter_space_der(to_i_space_der),
      to_phys_space(to_p_space),
      to_phys_space_der(to_p_space_der)
{
    if (inter.size() < 2) {
        throw std::invalid_argument("[SingleDiscretizationInfo::SingleDiscretizationInfo]: Need at least 2 points, found " + std::to_string(inter.size()) );
    }

    if (g_size.size() != (inter.size() - 1)) {
        throw std::invalid_argument("[SingleDiscretizationInfo::SingleDiscretizationInfo]: Degree (" + std::to_string(g_size.size()) + ") must match intervals (" + std::to_string(inter.size()-1) + ")");
    }

    if (std::any_of(g_size.begin(), g_size.end(), [](size_t p){ return p < 1; })) {
        throw std::invalid_argument("[SingleDiscretizationInfo::SingleDiscretizationInfo]: All polynomial degrees must be higher than 1");
    }

    for (size_t i=0; i<inter.size()-1; i++)
    {
        intervals[i] = {to_i_space(inter[i]), to_i_space(inter[i+1])};
        intervals_phys[i] = {inter[i], inter[i+1]};       
    }
} // SingleDiscretizationInfo::SingleDiscretizationInfo
  
double Grid1D::get_der_matrix(size_t a, size_t j, size_t b, size_t k) const
{
    if (a != b) return 0.0;
    return _stored_grids.at(_d_info.grid_sizes[a])._Dij[j][k];
    // at(key): access to the element

    // the matrix is a block-matrix, each block is the Dij of a single 
    // sub-interval, of a single sub-grid. So to access the element I have
    // to define both the sub interval and the point in the block matrix.
    // a and b define which block, j and k the position in the block.
    // a represent which row (block). b which column (block). a represent which 
    // interval. b represent which interval belong f you are calculating the 
    // derivative of. if they're not the same, the result is zero.
    // remember that ~ f'i = Dij fj. Then, if a and b are compatible, you choose
    // the correct element with i and j, which belong to the single block.
    //
    // stored_grids stores different chebyshev grid inside, one for each differnet 
    // degree. so you don't have to create a new one, you can take the one already
    // built from this map. the key is the degree. this is the reason behind we check
    // for grid_sizes of the a-th interval.
 
} // Grid1D::get_der_matrix

Grid1D::Grid1D(const SingleDiscretizationInfo &d_info) : _d_info(d_info)
{
    size = 0; // will represent N, total number of point also repeated
    if (d_info.intervals.size() == 0) {
        size_li = 0;
        c_size = 0;
        c_size_li = 0;
        return;
    }

    for (size_t i = 0; i < d_info.intervals.size(); i++) {
        // NOTE: this works ok because grid_sizes stores N, 
        // but the points are N+1 always
        // NOTE: this is implementation of non-overlapping grids, 
        // leading to simpler weights
        // TODO: interpolation was good in tests, but i need to check that 
        // the kernels are not screwed by this (related to quadrature & diff. eq.)
        size += d_info.grid_sizes[i]+1;
        if (_stored_grids.find(d_info.grid_sizes[i]) == _stored_grids.end()) 
        // find(key): search for an element, if not found return end()
        {
            _stored_grids.emplace(d_info.grid_sizes[i], Chebyshev::StandardGrid(d_info.grid_sizes[i]));
            // emplace(key, value): build the element in the map
        }
    }

    size_li = static_cast<index_t>(size); // static_cast to force size_t to be an integer with sign
                                          // index_t is integer with sign
    c_size = size - d_info.intervals.size() + 1;
    // coordinate size, number of unique points
    c_size_li = static_cast<index_t>(c_size);

    // resize create an object of dimension size, with default value.
    // i'm going to fill them later, but memory already allocated
    _weights.resize(size); // allocate dimension size bc for point on 
                           // border two contribute, from both intervals.
                           // so two weights, similar to point counted twice
    _weights_der.resize(size);
    _weights_sub.resize(size);
    _from_iw_to_ic.resize(size); // function to pass from index n to l
    _from_idx_to_inter.resize(size); // function from n to interval a

    _coord.resize(c_size); // i have to use c_size bc in real space points are not duplicated
    _coord_inter.resize(c_size);
    _delim_indexes.resize(d_info.intervals.size() + 1);

    _der_matrix.resize(size); // create the matrix, allocate space for rows of matrix
    for (size_t i=0; i<size; i++) {
        _der_matrix[i].resize(size); // for each element, allocate the correct quantity of memory
                                     // for the column 
    }

    // No check on order of arrival
    index_t index = 0;  // index n, all point with duplicates
    index_t index_coord = 0; // index l, all point without duplicates
    for (size_t a=0; a<d_info.intervals.size(); a++) { // for all intervals
        const Chebyshev::StandardGrid &sg = _stored_grids.at(d_info.grid_sizes[a]);
        _delim_indexes[a] = index;

        for (size_t j=0; j<=d_info.grid_sizes[a]; j++) {
            // _coord stores only the unique coordinates
            if (j!=d_info.grid_sizes[a] || (a==d_info.intervals.size()-1)) {
                // Compute coordinate in intermediate u space by linear map 
                _coord_inter[index_coord] = from_m1p1_to_ab(sg.t(j),
                        d_info.intervals[a].first, d_info.intervals[a].second);
               // Compute coordinate in physical x space
               _coord[index_coord] = d_info.to_phys_space(_coord_inter[index_coord]); 
            }

            // Fill map for n->l, from index_weight to index_count
            _from_iw_to_ic[index] = index_coord;
            // Fill map for n->i(n), i.e. from index to interval index
            _from_idx_to_inter[index] = a; // means: for a certain n I'm referring to
                                           // interval a

            if (j != d_info.grid_sizes[a]) index_coord++;

            // Note: in order to have proper assignement operator I cannot use
            // [this] in the lambdas. Hence, I need to pass the needed variables
            // by value to the lambdas.
            // lambda function: anonymous function that are characterized by a 
            // particular syntax -> [variable captured](parameter){proper function}
            // If a pass to a lambda as [var capt] something from another object
            // like d_info, e.g d_info.intervals.size, I'm passing a puntator to 
            // d_info. If I remove d_info, teh lambdas will crash. So I have to pass
            // by value the object. Hence I make copy of the value here.
            size_t cached_int_size = d_info.intervals.size(); 
            std::pair<double, double> cached_inter = d_info.intervals[a]; 

            // -----------------------
            // std::vector<std::function<double(double, const Chebyshev::StandardGrid &)>> _weights;
            // _weights are function, so I have to implement function. I'm using lambda function
            // [captures](parameters) -> returned_type { body }
            _weights[index] = [a, j, cached_int_size, cached_inter](double u, 
                    const Chebyshev::StandardGrid &sg) -> double {
                double res = 0;
                bool condition_x
                    = a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second;
                if (u >= cached_inter.first && condition_x) {
                    res += sg.poli_weight(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j);
                    // double poli_weight(double t, size_t j); 
                    // t the point in which to evaluate the weight & j the index of the weight to evaluate
                }
                return res;
            };
            // -----------------------
            
            _weights_der[index] = [a, j, cached_int_size, cached_inter](double u,
                    const Chebyshev::StandardGrid &sg) -> double {
                double res = 0;
                bool condition_x 
                    = a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second;
                if (u >= cached_inter.first && condition_x) {
                    const double dl_dx = from_ab_to_m1p1_der(cached_inter.first, cached_inter.second);
                    const double dw_dl = sg.poli_weight_der(from_ab_to_m1p1(u, cached_inter.first,
                                cached_inter.second), j);
                    // double poli_weight_der(double t, size_t j);
                    // From Notes, here we are trying to derivate bj(eta^{-1}(x)). In order to derivate 
                    // respect to x, we have to calculate d/d_eta^{-1} (bj) * d/dx (eta^{-1}). Here 
                    // eta^{-1} is the function from x (ab interval) to -1,+1 (m1p1).
                    // dw_dx = dw_dl * dl_dx. x is the variable in phys space, here represented by u.
                    res += dw_dl * dl_dx;
                }
                return res;
            };

            _weights_sub[index] = [a, j, cached_int_size, cached_inter](double u,
                    const Chebyshev::StandardGrid &sg) -> double {
                double res = 0;
                bool condition_x
                    = a == cached_int_size - 1 ? u <= cached_inter.second : u < cached_inter.second;
                if (u >= cached_inter.first && condition_x) {
                    res += sg.poli_weight(from_ab_to_m1p1(u, cached_inter.first, cached_inter.second), j)
                        - 1;
                }
                return res;
            }; // this will be used with integration of the function, if presence of singolarity
            
            index_t inner_index = 0;
            for (size_t b=0; b<d_info.intervals.size(); b++) {
                for (size_t k=0; k<=d_info.grid_sizes[b]; k++) {
                    _der_matrix[index][inner_index++] = get_der_matrix(a, j, b, k);
                }
            }
            index++;
        }
    }
    _delim_indexes[d_info.intervals.size()] = index;
    // Before I was associating for each a (from 0 to intervals.size-1) the value of the index n
    // associated to each new sub-intervals. I still have to assign the last limit, associated to 
    // the end of the last sub-intervals. The right limit of last sub-intervals is the last value 
    // of index (which represent the n index of the notes).

    _from_ic_to_iw.resize(c_size);
    // Resize the l->n function to the number of l index (c_size, index not repeated).

    for (size_t i=0; i<_from_iw_to_ic.size(); i++) {
        _from_ic_to_iw[_from_iw_to_ic[i]].push_back(i);
        // Get l from _from_iw_to_ic[n], then in _from_ic_to_iw[l] append the value n.
        // _form_ic_to_iw was defined as a vector of vector. In this way, I can append 
        // two different values at the same index l if I'm at the boundary and one if not.
    }

    {
      using integrator = GaussKronrod<GK_41>;
      // Decides to use Gauss Konrod as integrator

      _integral_weights.resize(c_size_li, 0.);
      // Prepare vector for the integral -> calc as a sum of f(node) * weights, where the
      // weight is the integral on the interval of the weight for the interpolation.
      // It's necessary to store a vector with number of values equal to the number
      // of physical nodes, without duplicates 

      size_t i_w;

      std::function<double(double)> full_integrand = [&](double v) -> double {
         return _d_info.to_phys_space_der(v) * _weights[i_w](v, get_std_grid(i_w));
      };
      // [&] means the function is able to see and use all the variables, by references,
      // so if something changes the function see it immediately
      // it takes as input variable v, we are integrating on interpolating space, not 
      // physical one. it returns to_phys_space_der(v), which is the jacobian due to the 
      // change of variable (to interpolating space, function is dx/dv), multiplied by
      // the value of the weight calc in v.

      for (size_t j = 0; j < size; j++) {
      // size number total of weights, also duplicates

         size_t j_c = _from_iw_to_ic[j];

         i_w = j;
         // assign to this variable the value of j, changing the value given to the 
         // full_integrand since given by references. In this way is possible to calc
         // integral of the correct weight

         auto [vmin, vmax] = get_support_weight_aj(j);
         // Get the support of the j-th weight in interpolation space, out from it 
         // integral goes to zero.

         if (std::fabs(vmax - vmin) < 1.0e-15) {
            continue;
         }
         _integral_weights[j_c]
             += integrator::integrate(full_integrand, vmin, vmax, 1.0e-10, 1.0e-10);
         // it calls gauss konrod library to integrate everything, take the result and 
         // sum in _integral_weights[j_c], using += and not only = bc if j_c is the 
         // junction between two intervals, the for will pass twice summing the "half"
         // areas

      } // j-loop
   }
}

vector_d Grid1D::integrate(const std::function<double(double)> &integrand, double eps) const
{
    using integrator = GaussKronrod<GK_41>;
    
    // Inizializza il vettore dei risultati a zero
    vector_d V(c_size, 0.0);

    for (size_t j = 0; j < size; j++) {
        size_t j_c = _from_iw_to_ic[j];
        auto [vmin, vmax] = get_support_weight_aj(j);

        if (std::fabs(vmax - vmin) < 1.0e-15) continue;

        // Nuova funzione da integrare: include integrand(y)
        std::function<double(double)> full_integrand = [&](double v) -> double {
            // 1. Trovo la coordinata fisica y
            double y = _d_info.to_phys_space(v);
            // 2. Moltiplico: K(y) * Jacobiano * w_j(v)
            return integrand(y) * _d_info.to_phys_space_der(v) * _weights[j](v, get_std_grid(j));
        };

        // Calcolo e sommo nel "cassetto" del nodo fisico
        V[j_c] += integrator::integrate(full_integrand, vmin, vmax, eps, eps);
    }
    return V;
}

matrix_d Grid1D::integrate(
    const std::function<double(double, double)> &integrand,
    const std::function<double(double)> &weight_fnc) const
{
    using integrator = GaussKronrod<GK_41>;
    
    // Assumendo che matrix_d abbia un costruttore (righe, colonne)
    matrix_d M(c_size, c_size); 
    // NOTA: Se matrix_d usa std::vector, potrebbe servire un resize manuale

    // Ciclo esterno: fisso il nodo fisico x_i
    for (size_t i = 0; i < c_size; i++) {
        double xi = _coord[i];

        // Ciclo interno: integro sulle funzioni di base w_j (uguale a prima)
        for (size_t j = 0; j < size; j++) {
            size_t j_c = _from_iw_to_ic[j];
            auto [vmin, vmax] = get_support_weight_aj(j);

            if (std::fabs(vmax - vmin) < 1.0e-15) continue;

            std::function<double(double)> full_integrand = [&](double v) -> double {
                double y = _d_info.to_phys_space(v);
                // Moltiplico: K(xi, y) * Q(y) * Jacobiano * w_j(v)
                return integrand(xi, y) * weight_fnc(y) * _d_info.to_phys_space_der(v) * _weights[j](v, get_std_grid(j));
            };

            // M(riga i, colonna j_c)
            M(i, j_c) += integrator::integrate(full_integrand, vmin, vmax, 1.0e-10, 1.0e-10); 
        }
    }
    return M;
}

vector_d Grid1D::integrate_subtr(const std::function<double(double)> &integrand) const
{
    using integrator = GaussKronrod<GK_41>;
    vector_d V(c_size, 0.0);

    for (size_t j = 0; j < size; j++) {
        size_t j_c = _from_iw_to_ic[j];
        auto [vmin, vmax] = get_support_weight_aj(j);

        if (std::fabs(vmax - vmin) < 1.0e-15) continue;

        std::function<double(double)> full_integrand = [&](double v) -> double {
            double y = _d_info.to_phys_space(v);
            
            // Il trucco è qui: usiamo _weights_sub invece di _weights!
            // _weights_sub contiene analiticamente (w - 1)
            double peso = _weights_sub[j](v, get_std_grid(j));
            
            return integrand(y) * _d_info.to_phys_space_der(v) * peso;
        };

        V[j_c] += integrator::integrate(full_integrand, vmin, vmax, 1.0e-10, 1.0e-10);
    }
    
    // Al risultato del nodo 0 (la singolarità), la matematica richiede 
    // di ri-sommare l'integrale puro di K(y) senza pesi
    // (vedi le tue dispense o la formula nel commento .hh)
    
    // Opzionale: calcolo integrale_puro = \int K(y) dy
    // V[0] += integrale_puro;

    return V;
}

} // namespace Interpolation

