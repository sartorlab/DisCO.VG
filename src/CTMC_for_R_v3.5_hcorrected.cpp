#include <RcppArmadillo.h>
#include <set>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
using namespace Rcpp;

// =============================================================================
  // CTMC: Coupled Tensor-Matrix Completion (v3h — legacy engine, h corrections only)
//
  // Implements Equation 8 from Bagherian et al. (2021)
//
  // Copy of CTMC_for_R.cpp (v3.1) carrying ONLY the corrections that concern the
// metric transforms h / g = h^{-1}. Everything else is the legacy engine
// verbatim: input preprocessing (ensureSimilarity + scale + epsilon), the
// relative eigenvalue floor, the Cholesky ridge fallback, the sweep order, and
// the dual convergence rule. The full strict engine (no floors, descent
                                                      // rejection, require_converged, state-change stopping, KKT diagnostics) is the
// separate CTMC_for_R_v2corrected.cpp.
//
  // Notation: Q = g^T g = h^{-2} (precision metric), S = Q^{-1} = h h^T = h^2.
// For A = C(w) + M Q_other M^T with eigendecomposition V diag(λ) V^T and
// c = det(A)^{1/n}:  S = A / c,  Q = c A^{-1}, det(S) = det(Q) = 1.
//
  // h corrections:
  //   [1] Interpolation kernel. With M fixed-metric optimal, stationarity of
//       ||g_X M g_Y||² subject to M[Ω] = v gives M = S_X E(u) S_Y and
//       K[k,l] = S_X(x_k,x_l) S_Y(y_k,y_l). The legacy engine used h in place of
//       S = h² in both K and the reconstruction h_X E(u) h_Y, which is feasible
//       but not the minimiser (docs/v4_ctmc_theory_audit_2026-09-10.md).
//   [2] Metric updates see the corrected M. The update rule is unchanged,
//       A_X = C_X + M Q_Y M^T (legacy wrote this as h_X E E^T h_X, which equals
                                //       M Q_Y M^T only for M = h_X E h_Y). M is now kept explicitly (n0 x n1),
//       and the Y update uses A_Y = C_Y + M^T Q_X M with the NEW Q_X and the same
//       M, so each metric step is an exact block minimiser. (Under the legacy
                                                              //       M = h E h the old and new Q_X gave the same A_Y; under M = S E S they do
                                                              //       not, and the fixed-M version is the one that descends.) This also drops
//       the dense m x m sharing matrices J_x, J_y and the selection matrices.
//   [3] Slice weights balance Tr(C_i Q) = accu(C_i % g²), the quantity in the
//       objective, not h² (v3.1 fix, 2026-06-17).
//   Objective: F = Tr(M^T Q_X M Q_Y) + Σ w_x Tr(C_i Q_X) + Σ w_y Tr(C_j Q_Y).
//   For M = S_X E(u) S_Y the data term equals u^T K u = dot(u, M[Ω]); the
//   legacy sum(u²) held only for M = h_X E(u) h_Y.
//
  // Linear solve for K u = v (solver = "cholesky" default, or "cg"):
  //   cholesky  dense K (m x m), Cholesky; on failure legacy 1e-10 ridge + solve.
//   cg        matrix-free Jacobi-preconditioned CG, warm-started from the
//             previous u; (K z)_k = [S_X E(z) S_Y](x_k, y_k) in O(m (n0 + n1)).
//             Stops at ||v - K u|| <= cg_tol ||v|| (checked on the true
                                                     //             residual). cg_max = 0 means 10 m iterations.
//
  // Flagging (numerical workarounds are recorded, never silent):
  //   Legacy numerics are kept, so a fit can still pass through the eigenvalue
//   floor, the ridge fallback, or an approximate solve. Armadillo's solve()
//   would otherwise substitute an SVD least-squares answer when rcond < eps
//   with only a console message; here every solve runs with no_approx, and on
//   failure the same SVD answer is computed explicitly with force_approx, so
//   results are unchanged but labelled. Per sweep:
//     solve_status  "ok" | "approx" | "ridge" | "ridge+approx" |
//                   "cg_unconverged" | "cg_breakdown"
//     eig_floored   eigenvalues raised to the 1e-10 relative floor (X + Y)
//     relres        ||v - M[Ω]|| / ||v||
//   Per fit:
//     fit_status    "ok"      no workaround anywhere, every relres <= relres_tol
//                   "flagged" completed and reproduces the observations, but
//                             some sweep used a workaround (see flags)
//                   "failed"  final relres > relres_tol, or a metric update /
//                             solve produced nonfinite values (loop stops; the
//                             last completed sweep is returned; see failure)
//     flags         distinct reasons: approx, ridge, cg_unconverged,
//                   cg_breakdown, eig_floor, relres, nonfinite
//   Invalid arguments still raise R errors.
//
// Convergence (unchanged from legacy v3): after a 2-sweep burn-in, stop when
//   |F - F_prev| / |F_prev| < tol  OR  ||u - u_prev|| / ||u|| < tol.
//   Reaching iter without either is returned with converged = FALSE.
// =============================================================================

// ---------------------------------------------------------------------------
// generateOmegaAndV
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec generateOmegaAndV_ctmc(arma::mat& xy, arma::mat& Omega) {
  arma::sp_mat mxy = arma::sp_mat(xy);
  arma::vec ans(mxy.n_nonzero);
  Omega = arma::mat(mxy.n_nonzero, 2);

  typedef arma::sp_mat::const_row_iterator iter;
  int k = 0;
  for (unsigned int i = 0; i < mxy.n_rows; i++) {
    for (iter j = mxy.begin_row(i); j != mxy.end_row(i); ++j) {
      ans.at(k) = *j;
      Omega(k, 0) = (double)i;
      Omega(k, 1) = (double)j.col();
      k++;
    }
  }
  return ans;
}

// ---------------------------------------------------------------------------
// metricUpdate: S = A / c and Q = c A^{-1}, c = det(A)^{1/n}
//
// One eigendecomposition A = V Λ V^t, with the legacy relative floor on Λ, so
// S and Q are exact inverses with det 1. S = h² and Q = g² in legacy terms.
// Returns false, leaving S and Q untouched, if A or the result is nonfinite.
// ---------------------------------------------------------------------------
bool metricUpdate(const arma::mat& A, arma::mat& S, arma::mat& Q, int& n_floored) {
  n_floored = 0;
  arma::vec eigval;
  arma::mat eigvec;
  if (!A.is_finite() || !arma::eig_sym(eigval, eigvec, A)) return false;

  double tol = 1e-10 * max(eigval);
  if (!(tol > 0)) return false;
  for (arma::uword i = 0; i < eigval.n_elem; ++i) {
    if (eigval(i) < tol) {
      eigval(i) = tol;
      ++n_floored;
    }
  }

  arma::vec log_eig = arma::log(eigval);
  double log_c = arma::mean(log_eig);

  arma::mat S_new = eigvec * arma::diagmat(arma::exp(log_eig - log_c)) * eigvec.t();
  arma::mat Q_new = eigvec * arma::diagmat(arma::exp(log_c - log_eig)) * eigvec.t();
  if (!S_new.is_finite() || !Q_new.is_finite()) return false;

  S = S_new;
  Q = Q_new;
  return true;
}

// ---------------------------------------------------------------------------
// ensureSimilarity: enforce PSD + unit diagonal (valid similarity matrix)
//
// Similarity matrices by convention have diag = 1 (maximal self-similarity)
// and off-diagonal in [0, 1]. Eigenvalue clamping can break unit diagonal,
// so we re-normalize after clamping. This is the standard contract assumed
// by the Bagherian et al. formulation.
//
// If your input is a covariance matrix (diag != 1), normalize it to a
// correlation matrix upstream, or use a different preprocessing function.
// ---------------------------------------------------------------------------
arma::mat ensureSimilarity_ctmc(const arma::mat& A) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A);
  eigval = arma::clamp(eigval, 0.0, datum::inf);

  arma::mat result = eigvec * arma::diagmat(eigval) * eigvec.t();

  // Restore unit diagonal (may have drifted above 1 due to clamping)
  arma::vec d = arma::sqrt(result.diag());
  d.replace(0.0, 1.0);
  mat D_inv = arma::diagmat(1.0 / d);
  return D_inv * result * D_inv;
}

// ---------------------------------------------------------------------------
// leftProduct: W = S_X E(z), n0 x n1, E(z) holding z at the observed cells
// ---------------------------------------------------------------------------
arma::mat leftProduct(const arma::mat& SX, int n1, const arma::uvec& idx_x, const arma::uvec& idx_y,
                const arma::vec& z) {
  arma::mat W(SX.n_rows, n1, arma::fill::zeros);
  for (arma::uword k = 0; k < z.n_elem; ++k) {
    W.col(idx_y(k)) += z(k) * SX.col(idx_x(k));
  }
  return W;
}

// ---------------------------------------------------------------------------
// applyK: (K z)_k = Σ_l S_X(x_k,x_l) S_Y(y_k,y_l) z_l = [S_X E(z) S_Y](x_k,y_k)
// ---------------------------------------------------------------------------
arma::vec applyK(const arma::mat& SX, const arma::mat& SY, const arma::uvec& idx_x, const arma::uvec& idx_y,
           const arma::vec& z) {
  // Transposed so each observation's dot product reads contiguous columns
arma::mat Wt = leftProduct(SX, SY.n_rows, idx_x, idx_y, z).t();
arma::vec out(z.n_elem);
for (arma::uword k = 0; k < z.n_elem; ++k) {
  out(k) = arma::dot(Wt.col(idx_x(k)), SY.col(idx_y(k)));
}
return out;
}

// ---------------------------------------------------------------------------
  // solveCG: Jacobi-preconditioned CG for K u = v, warm-started from u
  //
    // Returns "ok" when the true residual meets cg_tol, "cg_unconverged" when the
  // budget runs out first, and "cg_breakdown" on nonpositive curvature. u holds
  // the last iterate in every case.
  // ---------------------------------------------------------------------------
    std::string solveCG(const arma::mat& SX, const arma::mat& SY, const arma::uvec& idx_x,
                        const arma::uvec& idx_y, const arma::vec& v, arma::vec& u,
                        double cg_tol, int cg_max, int& cg_iters) {
      int m = v.n_elem;
      int budget = (cg_max > 0) ? cg_max : 10 * m;
      double tol_abs = cg_tol * arma::norm(v, 2);
      
      arma::vec diagK(m);
      for (int k = 0; k < m; ++k) {
        diagK(k) = SX(idx_x(k), idx_x(k)) * SY(idx_y(k), idx_y(k));
      }
      
      arma::vec r = v - applyK(SX, SY, idx_x, idx_y, u);
      arma::vec z = r / diagK;
      arma::vec p = z;
      double rz = arma::dot(r, z);
      
      cg_iters = 0;
      while (arma::norm(r, 2) > tol_abs && cg_iters < budget) {
        vec Kp = applyK(SX, SY, idx_x, idx_y, p);
        double pKp = arma::dot(p, Kp);
        if (!std::isfinite(pKp) || pKp <= 0) return "cg_breakdown";
        double alpha = rz / pKp;
        u += alpha * p;
        r -= alpha * Kp;
        ++cg_iters;
        
        if (arma::norm(r, 2) <= tol_abs) {
          // Accept only on the true residual; restart if the recursion drifted
          r = v - applyK(SX, SY, idx_x, idx_y, u);
          if (arma::norm(r, 2) <= tol_abs) break;
          z = r / diagK;
          p = z;
          rz = arma::dot(r, z);
          continue;
        }
        
        z = r / diagK;
        double rz_new = arma::dot(r, z);
        p = z + (rz_new / rz) * p;
        rz = rz_new;
      }
      
      r = v - applyK(SX, SY, idx_x, idx_y, u);
      return (arma::norm(r, 2) <= tol_abs) ? "ok" : "cg_unconverged";
    }
  
  // ---------------------------------------------------------------------------
    // CoupledTensorMatrixCompletion (v3h: h-corrected, legacy numerics, flagged)
  // ---------------------------------------------------------------------------
    arma::mat CoupledTensorMatrixCompletion(arma::mat& Omega, const arma::vec& v,
                                            arma::field<arma::mat>& C_xx,
                                            arma::field<arma::mat>& C_yy,
                                            int iter,
                                            double tol_converge,
                                            const std::string& solver,
                                            double cg_tol,
                                            int cg_max,
                                            double relres_tol,
                                            arma::vec& w_x_out,
                                            arma::vec& w_y_out,
                                            int& iters_used,
                                            double& F_final,
                                            arma::vec& F_trace_out,
                                            bool& converged,
                                            std::string& stop_reason,
                                            std::vector<std::string>& solve_status_out,
                                            arma::ivec& cg_iters_out,
                                            arma::ivec& eig_floored_out,
                                            arma::vec& relres_out,
                                            std::string& fit_status,
                                            std::vector<std::string>& flags_out,
                                            std::string& failure,
                                            arma::mat& SX, arma::mat& SY, arma::mat& QX, arma::mat& QY,
                                            arma::vec& u) {
      int m = Omega.n_rows;
      int nU = C_xx.n_rows;
      int nZ = C_yy.n_rows;
      int n0 = C_xx(0, 0).n_rows;
      int n1 = C_yy(0, 0).n_rows;
      
      // Precompute index vectors
      arma::uvec idx_x = arma::conv_to<arma::uvec>::from(Omega.col(0));
      arma::uvec idx_y = arma::conv_to<arma::uvec>::from(Omega.col(1));
      
      // Initialize: S = Q = I, u = v, so M = S_X E(u) S_Y = E(v)
      SX = arma::eye<arma::mat>(n0, n0);
      SY = arma::eye<arma::mat>(n1, n1);
      QX = arma::eye<arma::mat>(n0, n0);
      QY = arma::eye<arma::mat>(n1, n1);
      arma::vec w_x(nU, arma::fill::ones);
      arma::vec w_y(nZ, arma::fill::ones);
      u = v;
      arma::vec u_prev = v;
      double F_prev = arma::datum::inf;
      
      arma::mat M(n0, n1, arma::fill::zeros);
      for (int k = 0; k < m; ++k) {
        M(idx_x(k), idx_y(k)) = v(k);
      }
      
      // Preallocate temporaries
      arma::mat C_X(n0, n0);
      arma::mat C_Y(n1, n1);
      arma::mat A_x(n0, n0);
      arma::mat A_y(n1, n1);
      arma::mat K, L;
      const double v_norm = arma::norm(v, 2);
      
      // Per-sweep diagnostics
      arma::vec F_trace(iter, arma::fill::zeros);
      arma::vec relres(iter, arma::fill::zeros);
      arma::ivec cg_iters(iter, arma::fill::zeros);
      arma::ivec eig_floored(iter, arma::fill::zeros);
      std::vector<std::string> solve_status(iter);
      std::set<std::string> flags;
      
      iters_used = 0;
      F_final = NA_REAL;
      converged = false;
      stop_reason = "max_iter";
      failure = "";
      
      // ---- Main loop ----
        for (int l = 0; l < iter; ++l) {
          int floored_x = 0, floored_y = 0;
          
          // --- Update metric X: A_X = C_X + M Q_Y M^T ---
            {
              C_X = w_x(0) * C_xx(0, 0);
              for (int i = 1; i < nU; ++i) {
                C_X += w_x(i) * C_xx(i, 0);
              }
              
              A_x = C_X + M * QY * M.t();
              A_x = 0.5 * (A_x + A_x.t());
              if (!metricUpdate(A_x, SX, QX, floored_x)) {
                failure = "sweep " + std::to_string(l + 1) + ": X metric update nonfinite";
                break;
              }
            }
          
          // --- Update metric Y: A_Y = C_Y + M^T Q_X M (same M, new Q_X) ---
            {
              C_Y = w_y(0) * C_yy(0, 0);
              for (int j = 1; j < nZ; ++j) {
                C_Y += w_y(j) * C_yy(j, 0);
              }
              
              A_y = C_Y + M.t() * QX * M;
              A_y = 0.5 * (A_y + A_y.t());
              if (!metricUpdate(A_y, SY, QY, floored_y)) {
                failure = "sweep " + std::to_string(l + 1) + ": Y metric update nonfinite";
                break;
              }
            }
          
          // --- Update g_U weights ---
            // The regularizer is Tr(g_X C_i g_X) = accu(C_i % Q_X), so the weights
            // balance Q_X = g_X² (v3.1; was h² before).
            if (nU > 1) {
              arma::vec traces_x(nU);
              for (int i = 0; i < nU; ++i) {
                traces_x(i) = arma::accu(C_xx(i, 0) % QX);
                if (traces_x(i) < 1e-12) traces_x(i) = 1e-12;
              }
              double geomean = std::exp(arma::mean(arma::log(traces_x)));
              for (int i = 0; i < nU; ++i) {
                w_x(i) = geomean / traces_x(i);
              }
            }
            
            // --- Update g_Z weights ---
              if (nZ > 1) {
                arma::vec traces_y(nZ);
                for (int j = 0; j < nZ; ++j) {
                  traces_y(j) = arma::accu(C_yy(j, 0) % QY);
                  if (traces_y(j) < 1e-12) traces_y(j) = 1e-12;
                }
                double geomean = std::exp(arma::mean(arma::log(traces_y)));
                for (int j = 0; j < nZ; ++j) {
                  w_y(j) = geomean / traces_y(j);
                }
              }
            
            // --- Update u: K u = v with K[k,l] = S_X(x_k,x_l) S_Y(y_k,y_l) ---
              // no_approx makes a rank-deficient solve report failure instead of quietly
            // returning Armadillo's SVD answer; force_approx then computes that same
    // answer explicitly, so the legacy result is kept and the sweep is labelled.
    arma::vec u_new = u;
    std::string status = "ok";
    bool solved = true;
    if (solver == "cholesky") {
      K = SX(idx_x, idx_x) % SY(idx_y, idx_y);

      if (arma::chol(L, K, "lower")) {
        arma::vec z;
        bool exact = arma::solve(z, arma::trimatl(L), v, arma::solve_opts::no_approx);
        if (!exact) solved = arma::solve(z, arma::trimatl(L), v, arma::solve_opts::force_approx);
        if (solved && !arma::solve(u_new, arma::trimatu(L.t()), z, arma::solve_opts::no_approx)) {
          exact = false;
          solved = arma::solve(u_new, arma::trimatu(L.t()), z, arma::solve_opts::force_approx);
        }
        if (!exact) status = "approx";
      } else {
        K.diag() += 1e-10;
        status = "ridge";
        if (!arma::solve(u_new, K, v, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) {
          status = "ridge+approx";
          solved = arma::solve(u_new, K, v, arma::solve_opts::force_approx);
        }
      }
    } else {
      int steps = 0;
      status = solveCG(SX, SY, idx_x, idx_y, v, u_new, cg_tol, cg_max, steps);
      cg_iters(l) = steps;
    }

    if (!solved || !u_new.is_finite()) {
      failure = "sweep " + std::to_string(l + 1) + ": linear solve (" + status +
                ") failed or returned nonfinite coefficients";
      break;
    }

    // --- Reconstruct: M = S_X E(u) S_Y ---
    mat M_new = leftProduct(SX, n1, idx_x, idx_y, u_new) * SY;
    if (!M_new.is_finite()) {
      failure = "sweep " + std::to_string(l + 1) + ": reconstruction nonfinite";
      break;
    }
    u = u_new;
    M = M_new;

    // Sweep l completed: record its diagnostics
    iters_used = l + 1;
    solve_status[l] = status;
    eig_floored(l) = floored_x + floored_y;
    if (status != "ok") {
      if (status == "ridge+approx") {
        flags.insert("ridge");
        flags.insert("approx");
      } else {
        flags.insert(status);
      }
    }
    if (eig_floored(l) > 0) flags.insert("eig_floor");

    // ==============================================================
    // Objective F and dual convergence criteria
    //
    // F = Tr(M^T Q_X M Q_Y) + Σ_i w_x(i) Tr(C_i Q_X) + Σ_j w_y(j) Tr(C_j Q_Y)
    //
    // With M = S_X E(u) S_Y the data term is u^T K u, and M[Ω] = K u, so it
    // is dot(u, M[Ω]) exactly, whatever the accuracy of the linear solve.
    //
    // Dual criteria:
    //   1. F converged (relative change < tol)
    //   2. u stabilizes (relative change < tol)
    // Stop when EITHER triggers (after burn-in of 2 iterations).
    // ==============================================================
    {
      arma::vec M_obs(m);
      for (int k = 0; k < m; ++k) {
        M_obs(k) = M(idx_x(k), idx_y(k));
      }
      relres(l) = arma::norm(v - M_obs, 2) / (v_norm + 1e-30);
      if (relres(l) > relres_tol) flags.insert("relres");

      double F_data = arma::dot(u, M_obs);

      double F_reg_x = 0.0;
      for (int i = 0; i < nU; ++i) {
        F_reg_x += w_x(i) * arma::accu(C_xx(i, 0) % QX);
      }

      double F_reg_y = 0.0;
      for (int j = 0; j < nZ; ++j) {
        F_reg_y += w_y(j) * arma::accu(C_yy(j, 0) % QY);
      }

      double F_total = F_data + F_reg_x + F_reg_y;
      F_trace(l) = F_total;
      F_final = F_total;

      // Convergence check (dual criteria, skip first 2 iterations)
      if (l > 1 && tol_converge > 0) {
        // Criterion 1: objective effectively stopped decreasing
        double f_rel_change = fabs(F_total - F_prev) / (fabs(F_prev) + 1e-30);
        bool f_converged = (f_rel_change < tol_converge);

        // Criterion 2: u vector stabilized
        double u_change = arma::norm(u - u_prev, 2) / (arma::norm(u, 2) + 1e-30);
        bool u_stable = (u_change < tol_converge);

        if (f_converged || u_stable) {
          converged = true;
          stop_reason = (f_converged && u_stable) ? "objective+u"
                        : (f_converged ? "objective" : "u");
          break;
        }
      }

      F_prev = F_total;
      u_prev = u;
    }
  }

  // Fit status
  if (!failure.empty()) {
    flags.insert("nonfinite");
    stop_reason = "failed";
  } else if (relres(iters_used - 1) > relres_tol) {
    failure = tfm::format("final fit misses the observations: relative residual %.3e > relres_tol %.1e",
                          relres(iters_used - 1), relres_tol);
  }
  fit_status = !failure.empty() ? "failed" : (flags.empty() ? "ok" : "flagged");
  flags_out.assign(flags.begin(), flags.end());

  // Output (completed sweeps only)
  w_x_out = w_x;
  w_y_out = w_y;
  solve_status_out.assign(solve_status.begin(), solve_status.begin() + iters_used);
  if (iters_used > 0) {
    F_trace_out = F_trace(arma::span(0, iters_used - 1));
    relres_out = relres(arma::span(0, iters_used - 1));
    cg_iters_out = cg_iters(arma::span(0, iters_used - 1));
    eig_floored_out = eig_floored(arma::span(0, iters_used - 1));
  }

  return M;
}

// ---------------------------------------------------------------------------
// CTMC: Main R entry point
//
// Returns: xy_hat, w_x, w_y, nU, nZ, iters_used, F_final, F_trace (as legacy)
//   plus fit_status, flags, failure, converged, stop_reason, solver, and per
//   sweep solve_status, relres, eig_floored, cg_iters; then u, SX, SY (= h²),
//   QX, QY (= g²). See the header for the flagging rules.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List CTMC(arma::mat& xy,
          List xx_list,
          List yy_list,
          NumericVector xscales,
          NumericVector yscales,
          double xepsilon,
          double yepsilon,
          int iter,
          double tol = 1e-4,
          double cg_tol = 1e-10,
          int cg_max = 0,
          std::string solver = "cholesky",
          double relres_tol = 1e-8) {

  int nU = xx_list.size();
  int nZ = yy_list.size();

  if (nU == 0) stop("xx_list must contain at least one similarity matrix");
  if (nZ == 0) stop("yy_list must contain at least one similarity matrix");
  if (iter < 1) stop("iter must be at least 1");
  if (solver != "cholesky" && solver != "cg")
    stop("solver must be \"cholesky\" or \"cg\"");
  if (!(cg_tol > 0) || cg_max < 0) stop("cg_tol must be positive and cg_max nonnegative");
  if (!(relres_tol > 0)) stop("relres_tol must be positive");
  if (solver == "cg" && cg_tol > relres_tol)
    stop("cg_tol (%g) must not exceed relres_tol (%g)", cg_tol, relres_tol);

  arma::mat xx_first = as<arma::mat>(xx_list[0]);
  arma::mat yy_first = as<arma::mat>(yy_list[0]);
  int a = xx_first.n_rows;
  int b = yy_first.n_rows;

  // Handle scalar vs vector scales (with bounds check)
  arma::vec xscale_vec(nU);
  arma::vec yscale_vec(nZ);
  if (xscales.size() == 1) {
    xscale_vec.fill(xscales[0]);
  } else if (xscales.size() == nU) {
    for (int i = 0; i < nU; ++i) xscale_vec(i) = xscales[i];
  } else {
    stop("xscales must be length 1 (scalar) or length %d (one per VV slice), got %d",
         nU, xscales.size());
  }
  if (yscales.size() == 1) {
    yscale_vec.fill(yscales[0]);
  } else if (yscales.size() == nZ) {
    for (int j = 0; j < nZ; ++j) yscale_vec(j) = yscales[j];
  } else {
    stop("yscales must be length 1 (scalar) or length %d (one per GG slice), got %d",
         nZ, yscales.size());
  }

  // Extract Omega and v
  arma::mat Omega;
  arma::vec v = generateOmegaAndV_ctmc(xy, Omega);
  if (v.n_elem == 0) stop("Mxy has no nonzero (observed) entries");

  // Preprocess VV slices
  arma::field<arma::mat> C_xx(nU, 1);
  for (int i = 0; i < nU; ++i) {
    arma::mat xx_i = as<arma::mat>(xx_list[i]);
    xx_i = ensureSimilarity_ctmc(xx_i);
    C_xx(i, 0) = xscale_vec(i) * xx_i + xepsilon * arma::eye<arma::mat>(a, a);
  }

  // Preprocess GG slices
  arma::field<arma::mat> C_yy(nZ, 1);
  for (int j = 0; j < nZ; ++j) {
    arma::mat yy_j = as<arma::mat>(yy_list[j]);
    yy_j = ensureSimilarity_ctmc(yy_j);
    C_yy(j, 0) = yscale_vec(j) * yy_j + yepsilon * arma::eye<arma::mat>(b, b);
  }

  // Run
  arma::vec w_x_final(nU);
  arma::vec w_y_final(nZ);
  int iters_used = 0;
  double F_final = 0.0;
  arma::vec F_trace;
  bool converged = false;
  std::string stop_reason;
  std::vector<std::string> solve_status;
  arma::ivec cg_iters;
  arma::ivec eig_floored;
  arma::vec relres;
  std::string fit_status;
  std::vector<std::string> flags;
  std::string failure;
  arma::mat SX, SY, QX, QY;
  arma::vec u;
  arma::mat xy_hat = CoupledTensorMatrixCompletion(Omega, v, C_xx, C_yy, iter,
                                                   tol, solver, cg_tol, cg_max,
                                                   relres_tol,
                                                   w_x_final, w_y_final,
                                                   iters_used, F_final, F_trace,
                                                   converged, stop_reason,
                                                   solve_status, cg_iters,
                                                   eig_floored, relres,
                                                   fit_status, flags, failure,
                                                   SX, SY, QX, QY, u);

  List out = List::create(
    Named("xy_hat") = xy_hat,
    Named("w_x") = w_x_final,
    Named("w_y") = w_y_final,
    Named("nU") = nU,
    Named("nZ") = nZ,
    Named("iters_used") = iters_used,
    Named("F_final") = F_final,
    Named("F_trace") = F_trace,
    Named("fit_status") = fit_status,
    Named("flags") = CharacterVector(flags.begin(), flags.end()),
    Named("failure") = failure,
    Named("converged") = converged,
    Named("stop_reason") = stop_reason,
    Named("solver") = solver,
    Named("solve_status") = CharacterVector(solve_status.begin(), solve_status.end()),
    Named("relres") = relres,
    Named("eig_floored") = eig_floored,
    Named("cg_iters") = cg_iters
  );
  out["u"] = u;
  out["SX"] = SX;
  out["SY"] = SY;
  out["QX"] = QX;
  out["QY"] = QY;
  return out;
}

// ---------------------------------------------------------------------------
// CMMC: Backward-compatible wrapper
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List CMMC(arma::mat& xy, arma::mat& xx, arma::mat& yy,
          double xscale, double yscale,
          double xepsilon, double yepsilon,
          int iter,
          std::string solver = "cholesky") {

  List xx_list = List::create(xx);
  List yy_list = List::create(yy);
  NumericVector xscales = NumericVector::create(xscale);
  NumericVector yscales = NumericVector::create(yscale);

  return CTMC(xy, xx_list, yy_list, xscales, yscales,
              xepsilon, yepsilon, iter, 1e-4, 1e-10, 0, solver, 1e-8);
}

