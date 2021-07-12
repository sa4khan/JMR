// Generated by rstantools.  Do not edit by hand.

/*
    JMR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JMR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with JMR.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_eweibullsurv_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_eweibullsurv");
    reader.add_event(47, 45, "end", "model_eweibullsurv");
    return reader;
}
#include <stan_meta_header.hpp>
class model_eweibullsurv
  : public stan::model::model_base_crtp<model_eweibullsurv> {
private:
        int N;
        int ni;
        int pp1;
        int qpoints;
        double st;
        vector_d y;
        vector_d sigma;
        std::vector<vector_d> xalpha;
        matrix_d z_rand;
        vector_d pred0;
        vector_d phi;
        vector_d rho;
        vector_d kappa;
        vector_d gam;
        row_vector_d weights;
        matrix_d nodes;
        std::vector<matrix_d> chol_cov_matrix;
public:
    model_eweibullsurv(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_eweibullsurv(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_eweibullsurv_namespace::model_eweibullsurv";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 1);
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "ni", "int", context__.to_vec());
            ni = int(0);
            vals_i__ = context__.vals_i("ni");
            pos__ = 0;
            ni = vals_i__[pos__++];
            check_greater_or_equal(function__, "ni", ni, 1);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "pp1", "int", context__.to_vec());
            pp1 = int(0);
            vals_i__ = context__.vals_i("pp1");
            pos__ = 0;
            pp1 = vals_i__[pos__++];
            check_greater_or_equal(function__, "pp1", pp1, 2);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "qpoints", "int", context__.to_vec());
            qpoints = int(0);
            vals_i__ = context__.vals_i("qpoints");
            pos__ = 0;
            qpoints = vals_i__[pos__++];
            check_greater_or_equal(function__, "qpoints", qpoints, 3);
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "st", "double", context__.to_vec());
            st = double(0);
            vals_r__ = context__.vals_r("st");
            pos__ = 0;
            st = vals_r__[pos__++];
            check_greater_or_equal(function__, "st", st, 0);
            current_statement_begin__ = 7;
            validate_non_negative_index("y", "(ni * N)", (ni * N));
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec((ni * N)));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>((ni * N));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = (ni * N);
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("sigma", "(ni * N)", (ni * N));
            context__.validate_dims("data initialization", "sigma", "vector_d", context__.to_vec((ni * N)));
            sigma = Eigen::Matrix<double, Eigen::Dynamic, 1>((ni * N));
            vals_r__ = context__.vals_r("sigma");
            pos__ = 0;
            size_t sigma_j_1_max__ = (ni * N);
            for (size_t j_1__ = 0; j_1__ < sigma_j_1_max__; ++j_1__) {
                sigma(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("xalpha", "ni", ni);
            validate_non_negative_index("xalpha", "N", N);
            context__.validate_dims("data initialization", "xalpha", "vector_d", context__.to_vec(N,ni));
            xalpha = std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >(N, Eigen::Matrix<double, Eigen::Dynamic, 1>(ni));
            vals_r__ = context__.vals_r("xalpha");
            pos__ = 0;
            size_t xalpha_j_1_max__ = ni;
            size_t xalpha_k_0_max__ = N;
            for (size_t j_1__ = 0; j_1__ < xalpha_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < xalpha_k_0_max__; ++k_0__) {
                    xalpha[k_0__](j_1__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("z_rand", "ni", ni);
            validate_non_negative_index("z_rand", "pp1", pp1);
            context__.validate_dims("data initialization", "z_rand", "matrix_d", context__.to_vec(ni,pp1));
            z_rand = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(ni, pp1);
            vals_r__ = context__.vals_r("z_rand");
            pos__ = 0;
            size_t z_rand_j_2_max__ = pp1;
            size_t z_rand_j_1_max__ = ni;
            for (size_t j_2__ = 0; j_2__ < z_rand_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < z_rand_j_1_max__; ++j_1__) {
                    z_rand(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("pred0", "N", N);
            context__.validate_dims("data initialization", "pred0", "vector_d", context__.to_vec(N));
            pred0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("pred0");
            pos__ = 0;
            size_t pred0_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < pred0_j_1_max__; ++j_1__) {
                pred0(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("phi", "N", N);
            context__.validate_dims("data initialization", "phi", "vector_d", context__.to_vec(N));
            phi = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("phi");
            pos__ = 0;
            size_t phi_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < phi_j_1_max__; ++j_1__) {
                phi(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 13;
            validate_non_negative_index("rho", "N", N);
            context__.validate_dims("data initialization", "rho", "vector_d", context__.to_vec(N));
            rho = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("rho");
            pos__ = 0;
            size_t rho_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < rho_j_1_max__; ++j_1__) {
                rho(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 14;
            validate_non_negative_index("kappa", "N", N);
            context__.validate_dims("data initialization", "kappa", "vector_d", context__.to_vec(N));
            kappa = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("kappa");
            pos__ = 0;
            size_t kappa_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < kappa_j_1_max__; ++j_1__) {
                kappa(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 15;
            validate_non_negative_index("gam", "N", N);
            context__.validate_dims("data initialization", "gam", "vector_d", context__.to_vec(N));
            gam = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("gam");
            pos__ = 0;
            size_t gam_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < gam_j_1_max__; ++j_1__) {
                gam(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 16;
            validate_non_negative_index("weights", "qpoints", qpoints);
            context__.validate_dims("data initialization", "weights", "row_vector_d", context__.to_vec(qpoints));
            weights = Eigen::Matrix<double, 1, Eigen::Dynamic>(qpoints);
            vals_r__ = context__.vals_r("weights");
            pos__ = 0;
            size_t weights_j_1_max__ = qpoints;
            for (size_t j_1__ = 0; j_1__ < weights_j_1_max__; ++j_1__) {
                weights(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "weights", weights, 0);
            check_less_or_equal(function__, "weights", weights, 1);
            current_statement_begin__ = 17;
            validate_non_negative_index("nodes", "qpoints", qpoints);
            validate_non_negative_index("nodes", "pp1", pp1);
            context__.validate_dims("data initialization", "nodes", "matrix_d", context__.to_vec(qpoints,pp1));
            nodes = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(qpoints, pp1);
            vals_r__ = context__.vals_r("nodes");
            pos__ = 0;
            size_t nodes_j_2_max__ = pp1;
            size_t nodes_j_1_max__ = qpoints;
            for (size_t j_2__ = 0; j_2__ < nodes_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < nodes_j_1_max__; ++j_1__) {
                    nodes(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 18;
            validate_non_negative_index("chol_cov_matrix", "pp1", pp1);
            validate_non_negative_index("chol_cov_matrix", "pp1", pp1);
            validate_non_negative_index("chol_cov_matrix", "N", N);
            context__.validate_dims("data initialization", "chol_cov_matrix", "matrix_d", context__.to_vec(N,pp1,pp1));
            chol_cov_matrix = std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >(N, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(pp1, pp1));
            vals_r__ = context__.vals_r("chol_cov_matrix");
            pos__ = 0;
            size_t chol_cov_matrix_j_2_max__ = pp1;
            size_t chol_cov_matrix_j_1_max__ = pp1;
            size_t chol_cov_matrix_k_0_max__ = N;
            for (size_t j_2__ = 0; j_2__ < chol_cov_matrix_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < chol_cov_matrix_j_1_max__; ++j_1__) {
                    for (size_t k_0__ = 0; k_0__ < chol_cov_matrix_k_0_max__; ++k_0__) {
                        chol_cov_matrix[k_0__](j_1__, j_2__) = vals_r__[pos__++];
                    }
                }
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 21;
            validate_non_negative_index("b_unscaled", "pp1", pp1);
            validate_non_negative_index("b_unscaled", "N", N);
            num_params_r__ += (pp1 * N);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_eweibullsurv() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 21;
        if (!(context__.contains_r("b_unscaled")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable b_unscaled missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("b_unscaled");
        pos__ = 0U;
        validate_non_negative_index("b_unscaled", "pp1", pp1);
        validate_non_negative_index("b_unscaled", "N", N);
        context__.validate_dims("parameter initialization", "b_unscaled", "matrix_d", context__.to_vec(pp1,N));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_unscaled(pp1, N);
        size_t b_unscaled_j_2_max__ = N;
        size_t b_unscaled_j_1_max__ = pp1;
        for (size_t j_2__ = 0; j_2__ < b_unscaled_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < b_unscaled_j_1_max__; ++j_1__) {
                b_unscaled(j_1__, j_2__) = vals_r__[pos__++];
            }
        }
        try {
            writer__.matrix_unconstrain(b_unscaled);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable b_unscaled: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 21;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> b_unscaled;
            (void) b_unscaled;  // dummy to suppress unused var warning
            if (jacobian__)
                b_unscaled = in__.matrix_constrain(pp1, N, lp__);
            else
                b_unscaled = in__.matrix_constrain(pp1, N);
            // transformed parameters
            current_statement_begin__ = 24;
            validate_non_negative_index("b", "pp1", pp1);
            validate_non_negative_index("b", "N", N);
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> > b(N, Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1>(pp1));
            stan::math::initialize(b, DUMMY_VAR__);
            stan::math::fill(b, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 25;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 26;
                stan::model::assign(b, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            multiply(get_base1(chol_cov_matrix, i, "chol_cov_matrix", 1), col(b_unscaled, i)), 
                            "assigning variable b");
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 24;
            size_t b_k_0_max__ = N;
            size_t b_j_1_max__ = pp1;
            for (size_t k_0__ = 0; k_0__ < b_k_0_max__; ++k_0__) {
                for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
                    if (stan::math::is_uninitialized(b[k_0__](j_1__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: b" << "[" << k_0__ << "]" << "(" << j_1__ << ")";
                        stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable b: ") + msg__.str()), current_statement_begin__, prog_reader__());
                    }
                }
            }
            // model body
            {
            current_statement_begin__ = 30;
            int k(0);
            (void) k;  // dummy to suppress unused var warning
            stan::math::fill(k, std::numeric_limits<int>::min());
            stan::math::assign(k,1);
            current_statement_begin__ = 31;
            validate_non_negative_index("mu", "(ni * N)", (ni * N));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu((ni * N));
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 32;
            validate_non_negative_index("llik", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> llik(N);
            stan::math::initialize(llik, DUMMY_VAR__);
            stan::math::fill(llik, DUMMY_VAR__);
            current_statement_begin__ = 33;
            local_scalar_t__ psi(DUMMY_VAR__);
            (void) psi;  // dummy to suppress unused var warning
            stan::math::initialize(psi, DUMMY_VAR__);
            stan::math::fill(psi, DUMMY_VAR__);
            current_statement_begin__ = 34;
            local_scalar_t__ logs0(DUMMY_VAR__);
            (void) logs0;  // dummy to suppress unused var warning
            stan::math::initialize(logs0, DUMMY_VAR__);
            stan::math::fill(logs0, DUMMY_VAR__);
            current_statement_begin__ = 35;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 36;
                stan::model::assign(mu, 
                            stan::model::cons_list(stan::model::index_min_max(k, ((k + ni) - 1)), stan::model::nil_index_list()), 
                            add(get_base1(xalpha, i, "xalpha", 1), multiply(z_rand, get_base1(b, i, "b", 1))), 
                            "assigning variable mu");
                current_statement_begin__ = 37;
                stan::math::assign(k, (k + ni));
                current_statement_begin__ = 38;
                stan::math::assign(psi, (((multiply(weights, stan::math::exp(multiply(-(get_base1(phi, i, "phi", 1)), multiply(nodes, get_base1(b, i, "b", 1))))) * st) * stan::math::exp(-(get_base1(pred0, i, "pred0", 1)))) * 0.5));
                current_statement_begin__ = 39;
                stan::math::assign(logs0, log1m(stan::math::exp((get_base1(gam, i, "gam", 1) * log1m_exp(-(pow((get_base1(rho, i, "rho", 1) * psi), get_base1(kappa, i, "kappa", 1))))))));
                current_statement_begin__ = 40;
                stan::model::assign(llik, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (is_inf(logs0) ? stan::math::promote_scalar<local_scalar_t__>((stan::math::log(get_base1(gam, i, "gam", 1)) - pow((get_base1(rho, i, "rho", 1) * psi), get_base1(kappa, i, "kappa", 1)))) : stan::math::promote_scalar<local_scalar_t__>(logs0) ), 
                            "assigning variable llik");
            }
            current_statement_begin__ = 42;
            lp_accum__.add(normal_log(y, mu, sigma));
            current_statement_begin__ = 43;
            lp_accum__.add(llik);
            current_statement_begin__ = 44;
            lp_accum__.add(std_normal_log(to_vector(b_unscaled)));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("b_unscaled");
        names__.push_back("b");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(pp1);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dims__.push_back(pp1);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_eweibullsurv_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_unscaled = in__.matrix_constrain(pp1, N);
        size_t b_unscaled_j_2_max__ = N;
        size_t b_unscaled_j_1_max__ = pp1;
        for (size_t j_2__ = 0; j_2__ < b_unscaled_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < b_unscaled_j_1_max__; ++j_1__) {
                vars__.push_back(b_unscaled(j_1__, j_2__));
            }
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 24;
            validate_non_negative_index("b", "pp1", pp1);
            validate_non_negative_index("b", "N", N);
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > b(N, Eigen::Matrix<double, Eigen::Dynamic, 1>(pp1));
            stan::math::initialize(b, DUMMY_VAR__);
            stan::math::fill(b, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 25;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 26;
                stan::model::assign(b, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            multiply(get_base1(chol_cov_matrix, i, "chol_cov_matrix", 1), col(b_unscaled, i)), 
                            "assigning variable b");
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t b_j_1_max__ = pp1;
                size_t b_k_0_max__ = N;
                for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
                    for (size_t k_0__ = 0; k_0__ < b_k_0_max__; ++k_0__) {
                        vars__.push_back(b[k_0__](j_1__));
                    }
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_eweibullsurv";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t b_unscaled_j_2_max__ = N;
        size_t b_unscaled_j_1_max__ = pp1;
        for (size_t j_2__ = 0; j_2__ < b_unscaled_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < b_unscaled_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "b_unscaled" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t b_j_1_max__ = pp1;
            size_t b_k_0_max__ = N;
            for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < b_k_0_max__; ++k_0__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "b" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t b_unscaled_j_2_max__ = N;
        size_t b_unscaled_j_1_max__ = pp1;
        for (size_t j_2__ = 0; j_2__ < b_unscaled_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < b_unscaled_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "b_unscaled" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t b_j_1_max__ = pp1;
            size_t b_k_0_max__ = N;
            for (size_t j_1__ = 0; j_1__ < b_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < b_k_0_max__; ++k_0__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "b" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_eweibullsurv_namespace::model_eweibullsurv stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
