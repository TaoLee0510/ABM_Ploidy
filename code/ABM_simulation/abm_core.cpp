#include <Rcpp.h>
#include <limits>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//
// --------------- Helper: fetch CHR_W from the global environment ----------------
//

// Get NumericVector CHR_W from R_GlobalEnv
NumericVector get_chr_w() {
  SEXP chrw = Rf_findVarInFrame(R_GlobalEnv, Rf_install("CHR_W"));
  if (chrw == R_UnboundValue) {
    stop("CHR_W not found in global environment.");
  }
  return as<NumericVector>(chrw);
}

//
// --------------- Weighted ploidy helpers ----------------
//

// [[Rcpp::export]]
double mean_copy_number_weighted_cpp(IntegerVector kt, NumericVector w){
  int n = kt.size();
  if (w.size() != n) stop("mean_copy_number_weighted_cpp: weight length mismatch");
  long double num = 0.0L, den = 0.0L;
  for(int i=0;i<n;++i){
    num += (long double)kt[i] * (long double)w[i];
    den += (long double)w[i];
  }
  if (den <= 0.0L) return NA_REAL;
  return (double)(num / den);
}

// [[Rcpp::export]]
NumericVector mean_copy_number_rows_weighted_cpp(IntegerMatrix kt_mat, NumericVector w){
  int n = kt_mat.nrow();
  int m = kt_mat.ncol();
  if (w.size() != m) stop("mean_copy_number_rows_weighted_cpp: weight length mismatch");
  NumericVector out(n);
  long double den = 0.0L;
  for(int j=0;j<m;++j) den += (long double)w[j];
  if (den <= 0.0L) {
    for(int i=0;i<n;++i) out[i] = NA_REAL;
    return out;
  }
  for(int i=0;i<n;++i){
    long double num = 0.0L;
    for(int j=0;j<m;++j){
      num += (long double)kt_mat(i,j) * (long double)w[j];
    }
    out[i] = (double)(num / den);
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector karyotype_from_P_weighted_cpp(double P, NumericVector w,
                                            double tol=1e-3, int max_iter=20000, int max_copy=-1){
  const int chr_n = w.size();
  if (chr_n != 22) stop("karyotype_from_P_weighted_cpp expects 22 weights (autosomes 1..22)");
  IntegerVector kt(chr_n);
  for(int i=0;i<chr_n;++i) kt[i]=2;

  auto wsum = [&]()->long double{
    long double s=0.0L;
    for(int t=0;t<chr_n;++t) s += (long double)w[t];
    return s;
  };
  auto meanW = [&](const IntegerVector& v)->double{
    long double num=0.0L, den=wsum();
    for(int t=0;t<chr_n;++t) num += (long double)v[t] * (long double)w[t];
    if (den<=0.0L) return NA_REAL;
    return (double)(num/den);
  };

  long double SW = wsum();
  if (SW <= 0.0L) return kt;

  double currP = meanW(kt);
  double err = std::fabs(currP - P);
  int iter=0;

  while (iter < max_iter && err > tol){
    int dir = (currP < P) ? +1 : -1;
    int best_j = -1;
    double best_err = err;
    for(int j=0;j<chr_n;++j){
      int nv = kt[j] + dir;
      if (nv < 0) continue;
      if (max_copy > 0 && nv > max_copy) continue;
      double newP = currP + dir * ((double)w[j] / (double)SW);
      double e2 = std::fabs(newP - P);
      if (e2 + 1e-12 < best_err){
        best_err = e2;
        best_j = j;
      }
    }
    if (best_j < 0) break;
    kt[best_j] = kt[best_j] + dir;
    currP = currP + dir * ((double)w[best_j] / (double)SW);
    err = std::fabs(currP - P);
    ++iter;
  }
  return kt;
}

// [[Rcpp::export]]
IntegerMatrix karyotype_from_P_batch_weighted_cpp(NumericVector Pvec, NumericVector w,
                                                  double tol=1e-3, int max_iter=20000, int max_copy=-1){
  const int chr_n = w.size();
  if (chr_n != 22) stop("karyotype_from_P_batch_weighted_cpp expects 22 weights");
  int n = Pvec.size();
  IntegerMatrix out(n, chr_n);
  long double SW = 0.0L;
  for(int t=0;t<chr_n;++t) SW += (long double)w[t];
  if (SW <= 0.0L) return out;

  for (int i=0;i<n;++i){
    double P = Pvec[i];
    IntegerVector kt(chr_n);
    for(int c=0;c<chr_n;++c) kt[c]=2;
    double currP = 0.0;
    {
      long double num=0.0L;
      for(int t=0;t<chr_n;++t) num += (long double)kt[t]*(long double)w[t];
      currP = (double)(num/SW);
    }
    double err = std::fabs(currP - P);
    int iter=0;
    while (iter < max_iter && err > tol){
      int dir = (currP < P) ? +1 : -1;
      int best_j = -1;
      double best_err = err;
      for(int j=0;j<chr_n;++j){
        int nv = kt[j] + dir;
        if (nv < 0) continue;
        if (max_copy > 0 && nv > max_copy) continue;
        double newP = currP + dir * ((double)w[j] / (double)SW);
        double e2 = std::fabs(newP - P);
        if (e2 + 1e-12 < best_err){
          best_err = e2;
          best_j = j;
        }
      }
      if (best_j < 0) break;
      kt[best_j] = kt[best_j] + dir;
      currP = currP + dir * ((double)w[best_j] / (double)SW);
      err = std::fabs(currP - P);
      ++iter;
    }
    for (int c=0;c<chr_n;++c) out(i,c) = kt[c];
  }
  return out;
}

//
// --------------- Oxygen field PDE update ----------------
//

// [[Rcpp::export]]
NumericMatrix update_oxygen_field_cpp(NumericMatrix O2, IntegerMatrix grid,
                                      NumericVector Pv,
                                      double Coxy, double diff_rate,
                                      double supply_rate, double consume,
                                      int iters) {
  int N = O2.nrow();
  int M = O2.ncol();
  int nP = Pv.size();
  for (int it = 0; it < iters; ++it) {
    NumericMatrix nextO2(N, M);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        int iu = (i-1<0)?0:i-1;
        int id = (i+1>=N)?N-1:i+1;
        int jl = (j-1<0)?0:j-1;
        int jr = (j+1>=M)?M-1:j+1;
        double up    = O2(iu,j);
        double down  = O2(id,j);
        double left  = O2(i,jl);
        double right = O2(i,jr);
        double center = O2(i,j);
        double lap = (up + down + left + right - 4.0*center);

        double cons = 0.0;
        int idx = grid(i,j);
        if (idx != NA_INTEGER) {
          int k = idx - 1;
          if (k >= 0 && k < nP) {
            double P = Pv[k];
            if (!(R_finite(P)) || P <= 0.0) P = 1e-6;
            cons = consume * (P / 2.0);
          } else {
            cons = consume;
          }
        }
        double val = center + diff_rate * lap + supply_rate * (Coxy - center) - cons;
        if (val < 0.0) val = 0.0;
        if (val > 1.0) val = 1.0;
        nextO2(i,j) = val;
      }
    }
    O2 = nextO2;
  }
  return O2;
}

// [[Rcpp::export]]
NumericMatrix update_oxygen_field_vessels_cpp(NumericMatrix O2, IntegerMatrix grid,
                                              NumericVector Pv, IntegerMatrix vessel_mask,
                                              double Coxy, double diff_rate,
                                              double supply_rate, double consume,
                                              int iters, int use_soft_supply) {
  int N = O2.nrow();
  int M = O2.ncol();
  int nP = Pv.size();
  for (int it = 0; it < iters; ++it) {
    NumericMatrix nextO2(N, M);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        int iu = (i-1<0)?0:i-1;
        int id = (i+1>=N)?N-1:i+1;
        int jl = (j-1<0)?0:j-1;
        int jr = (j+1>=M)?M-1:j+1;
        double up    = O2(iu,j);
        double down  = O2(id,j);
        double left  = O2(i,jl);
        double right = O2(i,jr);
        double center = O2(i,j);
        double lap = (up + down + left + right - 4.0*center);

        double cons = 0.0;
        int idx = grid(i,j);
        if (idx != NA_INTEGER) {
          int k = idx - 1;
          if (k >= 0 && k < nP) {
            double P = Pv[k];
            if (!(R_finite(P)) || P <= 0.0) P = 1e-6;
            cons = consume * (P / 2.0);
          } else {
            cons = consume;
          }
        }

        if (vessel_mask(i,j) == 1) {
          nextO2(i,j) = Coxy;  // Vessel cell value fixed to Coxy
        } else {
          double val = center + diff_rate * lap - cons;
          if (use_soft_supply == 1) {
            val += supply_rate * (Coxy - center);
          }
          if (val < 0.0) val = 0.0;
          if (val > 1.0) val = 1.0;
          nextO2(i,j) = val;
        }
      }
    }
    O2 = nextO2;
  }
  return O2;
}

// [[Rcpp::export]]
NumericMatrix update_oxygen_field_bc_cpp(NumericMatrix O2, IntegerMatrix grid,
                                         NumericVector Pv,
                                         double Coxy, double diff_rate,
                                         double supply_rate, double consume,
                                         int iters, int bc_mode, double bc_value) {
  // bc_mode: 0 = Neumann-like (reflect), 1 = Dirichlet edges (fixed edge value = bc_value)
  int N = O2.nrow();
  int M = O2.ncol();
  int nP = Pv.size();
  for (int it = 0; it < iters; ++it) {
    NumericMatrix nextO2(N, M);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        double up    = (i-1 >= 0) ? O2(i-1,j) : (bc_mode==1 ? bc_value : O2(i,j));
        double down  = (i+1 <  N) ? O2(i+1,j) : (bc_mode==1 ? bc_value : O2(i,j));
        double left  = (j-1 >= 0) ? O2(i,j-1) : (bc_mode==1 ? bc_value : O2(i,j));
        double right = (j+1 <  M) ? O2(i,j+1) : (bc_mode==1 ? bc_value : O2(i,j));
        double center = O2(i,j);
        double lap = (up + down + left + right - 4.0*center);

        double cons = 0.0;
        int idx = grid(i,j);
        if (idx != NA_INTEGER) {
          int k = idx - 1;
          if (k >= 0 && k < nP) {
            double P = Pv[k];
            if (!(R_finite(P)) || P <= 0.0) P = 1e-6;
            cons = consume * (P / 2.0);
          } else {
            cons = consume;
          }
        }
        double val = center + diff_rate * lap +
          supply_rate * (Coxy - center) - cons;
        if (val < 0.0) val = 0.0;
        if (val > 1.0) val = 1.0;

        if (bc_mode == 1 && (i==0 || j==0 || i==N-1 || j==M-1)) {
          nextO2(i,j) = bc_value;
        } else {
          nextO2(i,j) = val;
        }
      }
    }
    O2 = nextO2;
  }
  return O2;
}

//
// --------------- Batched computation of G and DivisionTime ----------------
//

// [[Rcpp::export]]
List compute_G_and_Div_cpp(IntegerMatrix grid, NumericMatrix O2,
                           IntegerVector X, IntegerVector Y,
                           NumericVector Rv, NumericVector Pv,
                           double beta, int Kcap, int radius) {
  int N = grid.nrow();
  int M = grid.ncol();
  int n = X.size();
  NumericVector G(n);
  NumericVector Div(n);
  for (int idx = 0; idx < n; ++idx) {
    int xi = X[idx] - 1;
    int yj = Y[idx] - 1;
    if (xi < 0 || xi >= N || yj < 0 || yj >= M) {
      G[idx] = 0.0;
      Div[idx] = R_PosInf;
      continue;
    }
    int rmin = std::max(0, xi - radius);
    int rmax = std::min(N - 1, xi + radius);
    int cmin = std::max(0, yj - radius);
    int cmax = std::min(M - 1, yj + radius);
    int n_occ = 0;
    for (int ii = rmin; ii <= rmax; ++ii) {
      for (int jj = cmin; jj <= cmax; ++jj) {
        if (grid(ii,jj) != NA_INTEGER) ++n_occ;
      }
    }
    double Rcell = Rv[idx];
    double th = 1.0 - O2(xi, yj);
    double Pval = Pv[idx];
    if (Pval < 0.0) Pval = 0.0;
    double base = Rcell * (1.0 - th * (beta * Pval));
    if (base < 0.0) base = 0.0;
    double dens = 1.0 - ((double)n_occ / (double)Kcap);
    if (dens < 0.0) dens = 0.0;
    double g = base * dens;
    G[idx] = g;
    if (g > 0.0) Div[idx] = (std::log(2.0) / g) * 24.0;
    else Div[idx] = R_PosInf;
  }
  return List::create(_["G"] = G, _["DivisionTime"] = Div);
}

// [[Rcpp::export]]
List compute_G_Div_quiescence_cpp(int N, NumericMatrix O2,
                                  IntegerVector X, IntegerVector Y,
                                  NumericVector Rv, NumericVector Pv,
                                  double beta, int Kcap, int radius,
                                  IntegerVector Quiescent, NumericVector QuiescentTime,
                                  double q_death_thresh_hours,
                                  double softness) {
  int n = X.size();
  if (Y.size()!=n || Rv.size()!=n || Pv.size()!=n ||
      Quiescent.size()!=n || QuiescentTime.size()!=n) {
    stop("compute_G_Div_quiescence_cpp: Input vector lengths must match");
  }

  // Occupancy matrix: 1 if occupied, 0 otherwise
  IntegerMatrix occ(N, N);
  for (int i = 0; i < n; ++i) {
    int xi = X[i] - 1;
    int yj = Y[i] - 1;
    if (xi >= 0 && xi < N && yj >= 0 && yj < N) {
      occ(xi, yj) = 1;
    }
  }

  NumericVector G(n);
  NumericVector Div(n);
  IntegerVector Qout = clone(Quiescent);
  NumericVector QTout = clone(QuiescentTime);
  IntegerVector Alive(n, 1);

  for (int idx=0; idx<n; ++idx){
    int xi = X[idx];
    int yj = Y[idx];
    int r2 = radius * radius;
    int n_occ = 0;
    int local_sites = 0;
    for (int di = -radius; di <= radius; ++di) {
      int ii = xi + di;
      if (ii < 1 || ii > N) continue;
      for (int dj = -radius; dj <= radius; ++dj) {
        int jj = yj + dj;
        if (jj < 1 || jj > N) continue;
        if (di * di + dj * dj > r2) continue; // outside circle
        ++local_sites;
        if (occ(ii - 1, jj - 1) == 1) {
          ++n_occ;
        }
      }
    }

    double Rcell = Rv[idx];
    double th = 1.0 - O2(xi - 1, yj - 1);
    double Pval = Pv[idx];
    if (!(R_finite(Pval)) || Pval < 0.0) Pval = 0.0;

    double base = Rcell * (1.0 - th * (beta * Pval));
    if (base < 0.0) base = 0.0;

    // Smooth density factor based on local circular neighborhood occupancy.
    // occ_frac = 1.0 corresponds to a fully packed local neighborhood.
    double occ_frac = 0.0;
    if (local_sites > 0) {
      occ_frac = static_cast<double>(n_occ) / static_cast<double>(local_sites);
    }

    double dens;
    if (occ_frac >= 1.0) {
      // Geometric saturation: no remaining growth space.
      dens = 0.0;
    } else {
      double x = (occ_frac - 1.0) / softness;
      dens = 1.0 / (1.0 + std::exp(x));
    }

    double g = base * dens;

    // Unified rule: any non‑positive growth (g <= 0) is treated as quiescent.
    // Such cells get DivisionTime = Inf and are processed by the quiescence timer
    // rather than remaining alive forever without dividing.
    if (g <= 0.0) {
      g = 0.0;
      Qout[idx] = 1;
    }

    G[idx] = g;
    Div[idx] = (g > 0.0) ? ((std::log(2.0) / g) * 24.0) : R_PosInf;

    if (Qout[idx] == 1) {
      if (g > 0.0) {
        // Cell recovers from quiescence when growth becomes positive again.
        Qout[idx] = 0;
        QTout[idx] = 0.0;
      } else {
        // Persisting quiescence with G == 0: accumulate quiescent time
        // and trigger death if the threshold is exceeded.
        QTout[idx] = QTout[idx] + 1.0;
        if (QTout[idx] > q_death_thresh_hours) {
          Alive[idx] = 0;
        }
      }
    }
  }

  return List::create(_["G"]=G,
                      _["DivisionTime"]=Div,
                      _["Quiescent"]=Qout,
                      _["QuiescentTime"]=QTout,
                      _["AliveStatus"]=Alive);
}

//
// --------------- Various mask / neighborhood / MS / WGD helpers ----------------
//

// [[Rcpp::export]]
LogicalVector lethal_zero_copy_mask_cpp(IntegerMatrix kt_mat){
  int n = kt_mat.nrow();
  int m = kt_mat.ncol();
  LogicalVector out(n, false);
  for (int i=0; i<n; ++i){
    for (int j=0; j<m; ++j){
      int v = kt_mat(i,j);
      if (v == NA_INTEGER || v < 1){
        out[i] = true;
        break;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector stochastic_death_mask_cpp(NumericVector D){
  int n = D.size();
  LogicalVector out(n);
  for (int i=0; i<n; ++i){
    double p = D[i];
    if (!R_finite(p) || p <= 0.0) {
      out[i] = false;
      continue;
    }
    if (p >= 1.0) {
      out[i] = true;
      continue;
    }
    out[i] = (R::runif(0.0, 1.0) < p);
  }
  return out;
}

// [[Rcpp::export]]
LogicalVector ploidy_over_cap_mask_cpp(NumericVector P, double thr){
  int n = P.size();
  LogicalVector out(n, false);
  for (int i=0; i<n; ++i){
    double v = P[i];
    if (R_finite(v) && v > thr) out[i] = true;
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix empty_neighbors_cpp(IntegerMatrix grid, int i, int j) {
  int N = grid.nrow();
  int M = grid.ncol();
  static const int di[8] = {-1,-1,-1, 0,0, 1,1,1};
  static const int dj[8] = {-1, 0, 1,-1,1,-1,0,1};
  std::vector<int> ris, cjs;
  for (int k=0;k<8;++k){
    int ii = i + di[k];
    int jj = j + dj[k];
    if (ii>=1 && ii<=N && jj>=1 && jj<=M){
      if (IntegerVector::is_na(grid(ii-1,jj-1))) {
        ris.push_back(ii);
        cjs.push_back(jj);
      }
    }
  }
  int n = ris.size();
  IntegerMatrix out(n,2);
  for (int t=0;t<n;++t){
    out(t,0)=ris[t];
    out(t,1)=cjs[t];
  }
  return out;
}

// [[Rcpp::export]]
List apply_MS_to_two_cpp(IntegerVector parent_vec, int k_events, int max_copy){
  int n = parent_vec.size();
  IntegerVector a = clone(parent_vec);
  IntegerVector b = clone(parent_vec);
  if (k_events <= 0) return List::create(_["a"]=a, _["b"]=b);

  int upper = (max_copy > 0) ? max_copy : std::numeric_limits<int>::max();
  std::vector<int> valid;
  valid.reserve(n);
  for(int j=0;j<n;++j){
    int v = parent_vec[j];
    if (v >= 1 && v <= upper - 1) valid.push_back(j);
  }
  if (valid.empty()) return List::create(_["a"]=a, _["b"]=b);

  int k = std::min((int)valid.size(), k_events);
  for (int t = (int)valid.size()-1; t>0; --t){
    int u = (int)std::floor(R::runif(0.0, t+1));
    std::swap(valid[t], valid[u]);
  }
  for (int t=0; t<k; ++t){
    int j = valid[t];
    if (R::runif(0.0,1.0) < 0.5){
      a[j] = a[j] + 1;
      if (a[j] > upper) a[j] = upper;
      b[j] = std::max(0,  b[j]-1);
    } else {
      a[j] = std::max(0,  a[j]-1);
      b[j] = b[j] + 1;
      if (b[j] > upper) b[j] = upper;
    }
  }
  return List::create(_["a"]=a, _["b"]=b);
}

// [[Rcpp::export]]
double mean_copy_number_cpp(IntegerVector kt){
  int n = kt.size();
  double s = 0.0;
  for(int i=0;i<n;++i) s += kt[i];
  return s / (double)n;
}

// [[Rcpp::export]]
NumericVector mean_copy_number_rows_cpp(IntegerMatrix kt_mat){
  int n = kt_mat.nrow();
  int m = kt_mat.ncol();
  NumericVector out(n);
  for(int i=0;i<n;++i){
    double s=0.0;
    for(int j=0;j<m;++j) s += kt_mat(i,j);
    out[i] = s / (double)m;
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector apply_WGD_to_row_vec_cpp(IntegerVector kt, int max_copy){
  int n = kt.size();
  IntegerVector out = clone(kt);
  int upper = (max_copy > 0) ? max_copy : std::numeric_limits<int>::max();
  for(int i=0;i<n;++i){
    long v = (long)out[i] * 2L;
    if (v > upper) v = upper;
    out[i] = (int)v;
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector karyotype_from_P_cpp(double P, double tol=1e-3,
                                   int max_iter=20000, int max_copy=-1){
  const int chr_n = 22;
  IntegerVector kt(chr_n);
  for(int i=0;i<chr_n;++i) kt[i]=2;
  auto mean22 = [&](const IntegerVector& v)->double{
    double s=0.0;
    for(int i=0;i<chr_n;++i) s += v[i];
    return s/22.0;
  };
  double currP = mean22(kt);
  double err = std::fabs(currP - P);
  int iter=0;
  while (iter < max_iter && err > tol){
    int dir = (currP < P) ? +1 : -1;
    int best_j = -1;
    double best_err = err;
    for(int j=0;j<chr_n;++j){
      int nv = kt[j] + dir;
      if (nv < 0) continue;
      if (max_copy > 0 && nv > max_copy) continue;
      double newP = currP + (nv - kt[j]) / 22.0;
      double e2 = std::fabs(newP - P);
      if (e2 + 1e-12 < best_err){
        best_err = e2;
        best_j = j;
      }
    }
    if (best_j < 0) break;
    kt[best_j] = kt[best_j] + dir;
    currP = mean22(kt);
    err = best_err;
    ++iter;
  }
  return kt;
}

// [[Rcpp::export]]
IntegerMatrix karyotype_from_P_batch_cpp(NumericVector Pvec, double tol=1e-3,
                                         int max_iter=20000, int max_copy=-1){
  const int chr_n = 22;
  int n = Pvec.size();
  IntegerMatrix out(n,22);
  for (int i=0;i<n;++i){
    double P = Pvec[i];
    IntegerVector kt(chr_n);
    for(int c=0;c<chr_n;++c) kt[c]=2;
    auto mean22 = [&](const IntegerVector& v)->double{
      double s=0.0;
      for(int t=0;t<chr_n;++t) s += v[t];
      return s/22.0;
    };
    double currP = mean22(kt);
    double err = std::fabs(currP - P);
    int iter=0;
    while (iter < max_iter && err > tol){
      int dir = (currP < P) ? +1 : -1;
      int best_j = -1;
      double best_err = err;
      for(int j=0;j<chr_n;++j){
        int nv = kt[j] + dir;
        if (nv < 0) continue;
        if (max_copy > 0 && nv > max_copy) continue;
        double newP = currP + (nv - kt[j]) / 22.0;
        double e2 = std::fabs(newP - P);
        if (e2 + 1e-12 < best_err){
          best_err = e2;
          best_j = j;
        }
      }
      if (best_j < 0) break;
      kt[best_j] = kt[best_j] + dir;
      currP = mean22(kt);
      err = best_err;
      ++iter;
    }
    for (int c=0;c<chr_n;++c) out(i,c) = kt[c];
  }
  return out;
}

// [[Rcpp::export]]
IntegerMatrix rebuild_grid_from_cells_cpp(IntegerMatrix grid,
                                          IntegerVector X, IntegerVector Y){
  int N = grid.nrow();
  int M = grid.ncol();
  for (int i=0;i<N;++i)
    for (int j=0;j<M;++j)
      grid(i,j) = NA_INTEGER;
  int n = X.size();
  for (int r=0;r<n;++r){
    int xi = X[r]-1;
    int yj = Y[r]-1;
    if (xi>=0 && xi<N && yj>=0 && yj<M){
      grid(xi,yj) = r+1;
    }
  }
  return grid;
}

//
// --------------- Main core: simulate_step_cpp ----------------
//

// [[Rcpp::export]]
List simulate_step_cpp(List state) {
  // Unpack state
  IntegerMatrix grid = as<IntegerMatrix>(state["grid"]);
  NumericMatrix O2   = as<NumericMatrix>(state["O2"]);
  DataFrame cells_df = as<DataFrame>(state["cells"]);
  List cfg           = state["cfg"];
  int step           = as<int>(state["step"]);

  SEXP vessel_mask_sexp =
    state.containsElementNamed("vessel_mask") ? state["vessel_mask"] : R_NilValue;
  IntegerMatrix vessel_mask;
  bool has_vessel = false;
  if (!Rf_isNull(vessel_mask_sexp)) {
    vessel_mask = as<IntegerMatrix>(vessel_mask_sexp);
    has_vessel = true;
  }

  int N = grid.nrow();
  int n_cells = cells_df.nrows();

  // Cell-level columns
  IntegerVector X = cells_df["X"];
  IntegerVector Y = cells_df["Y"];
  NumericVector R = cells_df["R"];
  NumericVector D = cells_df["D"];
  NumericVector P = cells_df["P"];
  NumericVector G = cells_df["G"];
  NumericVector DivisionTime = cells_df["DivisionTime"];
  NumericVector Time = cells_df["Time"];
  IntegerVector Status = cells_df["Status"];

  IntegerVector Quiescent, WGDLabel, WGDCount;
  NumericVector QuiescentTime;
  if (cells_df.containsElementNamed("Quiescent"))
    Quiescent = cells_df["Quiescent"];
  else
    Quiescent = IntegerVector(n_cells, 0);

  if (cells_df.containsElementNamed("QuiescentTime"))
    QuiescentTime = cells_df["QuiescentTime"];
  else
    QuiescentTime = NumericVector(n_cells, 0.0);

  if (cells_df.containsElementNamed("WGDLabel"))
    WGDLabel = cells_df["WGDLabel"];
  else
    WGDLabel = IntegerVector(n_cells, 0);

  if (cells_df.containsElementNamed("WGDCount"))
    WGDCount = cells_df["WGDCount"];
  else
    WGDCount = IntegerVector(n_cells, 0);

  CharacterVector DeathReason =
    cells_df.containsElementNamed("DeathReason")
      ? cells_df["DeathReason"]
      : CharacterVector(n_cells, NA_STRING);

  // Karyotype chromosome columns (C1..C22)
  std::vector<std::string> c_names;
  for (int c = 1; c <= 22; ++c)
    c_names.push_back("C" + std::to_string(c));
  std::vector<IntegerVector> C(22);
  for (int c = 0; c < 22; ++c) {
    C[c] = as<IntegerVector>(cells_df[c_names[c]]);
  }

  // -------- 0) Oxygen update --------
  double Coxy_run = 1.0;  // default
  if (cfg.containsElementNamed("Coxy_scalar")) {
    SEXP coxy_s = cfg["Coxy_scalar"];
    if (!Rf_isNull(coxy_s)) {
      Coxy_run = as<double>(coxy_s);
    }
  } else if (cfg.containsElementNamed("Coxy")) {
    SEXP coxy_s = cfg["Coxy"];
    if (!Rf_isNull(coxy_s) && Rf_length(coxy_s) >= 1) {
      Coxy_run = as<double>(coxy_s);
    }
  }

  double dr = cfg.containsElementNamed("O2DiffRate")
                ? as<double>(cfg["O2DiffRate"]) : 0.2;
  double sr = cfg.containsElementNamed("O2SupplyRate")
                ? as<double>(cfg["O2SupplyRate"]) : 0.05;
  double cs = cfg.containsElementNamed("O2Consume")
                ? as<double>(cfg["O2Consume"]) : 0.02;
  int iters = cfg.containsElementNamed("O2JacobiIters")
                ? as<int>(cfg["O2JacobiIters"]) : 1;
  if (iters < 1) iters = 1;

  NumericVector Pv = P;

  if (has_vessel) {
    int use_soft = (cfg.containsElementNamed("O2VesselSoftSupply") &&
                    as<bool>(cfg["O2VesselSoftSupply"])) ? 1 : 0;
    O2 = update_oxygen_field_vessels_cpp(
      O2, grid, Pv, vessel_mask,
      Coxy_run, dr, sr, cs, iters, use_soft
    );
  } else {
    bool use_boundary = false;
    if (cfg.containsElementNamed("O2UseBoundary")) {
      SEXP ub = cfg["O2UseBoundary"];
      if (!Rf_isNull(ub)) {
        use_boundary = as<bool>(ub);
      }
    }

    if (use_boundary) {
      int bc_mode = 0;
      if (cfg.containsElementNamed("O2BoundaryMode")) {
        SEXP bm = cfg["O2BoundaryMode"];
        if (!Rf_isNull(bm)) {
          std::string mode = as<std::string>(bm);
          if (mode == "dirichlet" || mode == "dirichlet_edges") {
            bc_mode = 1;
          }
        }
      }
      double bc_val = Coxy_run;
      O2 = update_oxygen_field_bc_cpp(
        O2, grid, Pv,
        Coxy_run, dr, sr, cs,
        1, bc_mode, bc_val
      );
    } else {
      O2 = update_oxygen_field_cpp(
        O2, grid, Pv,
        Coxy_run, dr, sr, cs,
        1
      );
    }
  }

  // -------- 1) Time++ for alive --------
  for (int i = 0; i < n_cells; ++i) {
    if (Status[i] == 1) {
      Time[i] += 1.0;
    }
  }

  // -------- 2) Division + WGD --------
  IntegerVector alive_rows;
  for (int i = 0; i < n_cells; ++i) {
    if (Status[i] == 1) alive_rows.push_back(i + 1); // 1-based
  }
  int alive_n = alive_rows.size();

  std::vector<int> snapshot(alive_n);
  for (int k = 0; k < alive_n; ++k) {
    snapshot[k] = alive_rows[k] - 1;
  }

  NumericVector CHR_W = get_chr_w();

  double WGDr   = cfg.containsElementNamed("WGDr")   ? as<double>(cfg["WGDr"])   : 0.0;
  double MSR    = cfg.containsElementNamed("MSR")    ? as<double>(cfg["MSR"])    : 0.0;
  double alpha  = cfg.containsElementNamed("alpha")  ? as<double>(cfg["alpha"])  : 0.0;
  double Dv     = cfg.containsElementNamed("Dv")     ? as<double>(cfg["Dv"])     : 0.0;
  double beta   = cfg.containsElementNamed("beta")   ? as<double>(cfg["beta"])   : 1.0;
  double PloidyMax =
      cfg.containsElementNamed("PloidyMax") ? as<double>(cfg["PloidyMax"]) : R_PosInf;

  double MigProb =
      cfg.containsElementNamed("MigProb") ? as<double>(cfg["MigProb"]) : 0.2;
  if (!R_finite(MigProb)) MigProb = 0.2;
  if (MigProb < 0.0) MigProb = 0.0;
  if (MigProb > 1.0) MigProb = 1.0;

  bool MigBiasToVessel =
      (cfg.containsElementNamed("MigBiasToVessel") &&
       as<bool>(cfg["MigBiasToVessel"]));

  int max_copy_cpp =
    (R_finite(PloidyMax)) ? (int)PloidyMax : std::numeric_limits<int>::max();
  int ploidy_cap_death =
    (R_finite(PloidyMax)) ? (int)PloidyMax : 10;

  if (WGDCount.size() != n_cells)
    WGDCount = IntegerVector(n_cells, 0);

  for (size_t snap_idx = 0; snap_idx < snapshot.size(); ++snap_idx) {
    int irow = snapshot[snap_idx];
    if (irow >= n_cells) continue;
    if (Status[irow] != 1) continue;

    // Current coordinates
    int i = X[irow];
    int j = Y[irow];

    // Random migration: with probability MigProb, move one step into a random empty neighbor cell.
    if (R::runif(0.0, 1.0) < MigProb) {
      IntegerMatrix empty_m = empty_neighbors_cpp(grid, i, j);
      std::vector<int> mig_idx;
      if (empty_m.nrow() > 0) {
        for (int k = 0; k < empty_m.nrow(); ++k) {
          int ni = empty_m(k, 0);
          int nj = empty_m(k, 1);
          if (has_vessel && vessel_mask(ni - 1, nj - 1) == 1) continue;
          mig_idx.push_back(k);
        }
      }
      if (!mig_idx.empty()) {
        int choice_m;
        // If no bias is requested or no vessel mask is available, keep uniform choice
        if (!MigBiasToVessel || !has_vessel) {
          if (mig_idx.size() == 1) {
            choice_m = mig_idx[0];
          } else {
            choice_m = mig_idx[(int)std::floor(R::runif(0.0, 1.0) * mig_idx.size())];
          }
        } else {
          // Bias migration toward neighbors with higher oxygen (proximal to vessels).
          int m = (int)mig_idx.size();
          std::vector<double> w(m);
          double sumw = 0.0;
          for (int t = 0; t < m; ++t) {
            int row_idx = mig_idx[t];
            int ni = empty_m(row_idx, 0);
            int nj = empty_m(row_idx, 1);
            double o2_val = O2(ni - 1, nj - 1);
            if (!R_finite(o2_val) || o2_val < 0.0) o2_val = 0.0;
            if (o2_val > 1.0) o2_val = 1.0;
            // Simple linear weight: neighbors with higher O2 more likely.
            double wt = 1e-6 + o2_val;
            w[t] = wt;
            sumw += wt;
          }
          if (sumw <= 0.0) {
            // Fallback to uniform if for some reason no positive weight
            if (m == 1) {
              choice_m = mig_idx[0];
            } else {
              choice_m = mig_idx[(int)std::floor(R::runif(0.0, 1.0) * mig_idx.size())];
            }
          } else {
            double r = R::runif(0.0, 1.0) * sumw;
            double acc = 0.0;
            int chosen_idx = m - 1;
            for (int t = 0; t < m; ++t) {
              acc += w[t];
              if (r <= acc) {
                chosen_idx = t;
                break;
              }
            }
            choice_m = mig_idx[chosen_idx];
          }
        }

        int xi_new = empty_m(choice_m, 0);
        int yj_new = empty_m(choice_m, 1);

        // Clear old position if it still references this cell
        if (i >= 1 && i <= N && j >= 1 && j <= N) {
          if (grid(i - 1, j - 1) == (irow + 1)) {
            grid(i - 1, j - 1) = NA_INTEGER;
          }
        }
        // Update coordinates and grid with new position
        X[irow] = xi_new;
        Y[irow] = yj_new;
        grid(xi_new - 1, yj_new - 1) = irow + 1;

        // Refresh local variables
        i = xi_new;
        j = yj_new;
      }
    }

    // Only attempt division if the cell has reached its division time
    if (!R_finite(DivisionTime[irow]) || DivisionTime[irow] <= 0.0) {
      // Non-proliferative or undefined division schedule: skip division this step
      continue;
    }
    if (Time[irow] < DivisionTime[irow]) {
      // Not yet ready to divide
      continue;
    }

    // Quiescence: full 11x11 window
    int rmin = std::max(1, i - 5);
    int rmax = std::min(N, i + 5);
    int cmin = std::max(1, j - 5);
    int cmax = std::min(N, j + 5);
    int n_occ = 0;
    for (int ii = rmin; ii <= rmax; ++ii)
      for (int jj = cmin; jj <= cmax; ++jj)
        if (!IntegerVector::is_na(grid(ii - 1, jj - 1))) ++n_occ;

    if (n_occ >= 121) {
      G[irow] = 0.0;
      Quiescent[irow] = 1;
      continue;
    }

    // 8-neighborhood empty sites, excluding vessel cells
    IntegerMatrix empty = empty_neighbors_cpp(grid, i, j);
    std::vector<int> keep_idx;
    if (empty.nrow() > 0) {
      for (int k = 0; k < empty.nrow(); ++k) {
        int ni = empty(k, 0);
        int nj = empty(k, 1);
        if (has_vessel && vessel_mask(ni - 1, nj - 1) == 1) continue;
        keep_idx.push_back(k);
      }
    }

    if (keep_idx.empty()) {
      if (R::runif(0.0, 1.0) < Dv) {
        Status[irow] = 0;
        if (DeathReason[irow] == NA_STRING ||
            DeathReason[irow] == "" ||
            DeathReason[irow] == std::string("")) {
          DeathReason[irow] = "division_death_Dv";
        }
      } else {
        Quiescent[irow] = 1;
      }
      continue;
    }

    int choice =
      (keep_idx.size() == 1)
        ? keep_idx[0]
        : keep_idx[(int)std::floor(R::runif(0.0, 1.0) * keep_idx.size())];

    int xi = empty(choice, 0);
    int yj = empty(choice, 1);

    // WGD branch
    if (R::runif(0.0, 1.0) < WGDr) {
      IntegerVector v(22);
      for (int c = 0; c < 22; ++c) v[c] = C[c][irow];
      v = apply_WGD_to_row_vec_cpp(v, max_copy_cpp);
      for (int c = 0; c < 22; ++c) C[c][irow] = v[c];
      WGDLabel[irow] = 1;
      WGDCount[irow] = WGDCount[irow] + 1;
      P[irow] = mean_copy_number_weighted_cpp(v, CHR_W);
      Time[irow] = 0.0;
      G[irow] = NA_REAL;
      DivisionTime[irow] = R_PosInf;
      continue;
    }

    // mitosis branch
    double Pval = std::max(P[irow], 0.0);
    double val_o2 = O2(i - 1, j - 1);
    double theta_local = 1.0 - val_o2;
    double MSR_eff = MSR * (1.0 + theta_local * alpha);
    if (!R_finite(MSR_eff) || MSR_eff < 0) MSR_eff = 0;
    double lambda_ms = MSR_eff * 22.0 * Pval;
    if (!R_finite(lambda_ms) || lambda_ms < 0) lambda_ms = 0;
    int k_ms = R::rpois(lambda_ms);

    IntegerVector parent_kt(22);
    for (int c = 0; c < 22; ++c) parent_kt[c] = C[c][irow];

    List ab = apply_MS_to_two_cpp(parent_kt, k_ms, max_copy_cpp);
    IntegerVector a_vec = ab["a"];
    IntegerVector b_vec = ab["b"];

    bool lethal_a = false, lethal_b = false;
    for (int c = 0; c < 22; ++c) {
      if (!R_finite(a_vec[c]) || a_vec[c] < 1) lethal_a = true;
      if (!R_finite(b_vec[c]) || b_vec[c] < 1) lethal_b = true;
    }

    if (lethal_a && lethal_b) {
      Status[irow] = 0;
      if (DeathReason[irow] == NA_STRING ||
          DeathReason[irow] == "" ||
          DeathReason[irow] == std::string("")) {
        DeathReason[irow] = "division_death_Dv";
      }
      continue;
    }

    if (lethal_a && !lethal_b) {
      // Only generate daughter cell
      X.push_back(xi);
      Y.push_back(yj);
      for (int c = 0; c < 22; ++c) {
        C[c].push_back(b_vec[c]);
      }
      R.push_back(R[irow]);
      D.push_back(D[irow]);
      P.push_back(mean_copy_number_weighted_cpp(b_vec, CHR_W));
      G.push_back(NA_REAL);
      WGDLabel.push_back(0);
      WGDCount.push_back(WGDCount[irow]);
      DivisionTime.push_back(R_PosInf);
      Time.push_back(0.0);
      Status.push_back(1);
      Quiescent.push_back(0);
      QuiescentTime.push_back(0.0);
      DeathReason.push_back(NA_STRING);
      grid(xi - 1, yj - 1) = X.size();
      Status[irow] = 0;
      if (DeathReason[irow] == NA_STRING ||
          DeathReason[irow] == "" ||
          DeathReason[irow] == std::string("")) {
        DeathReason[irow] = "division_death_Dv";
      }
      continue;
    }

    if (!lethal_a && lethal_b) {
      for (int c = 0; c < 22; ++c) C[c][irow] = a_vec[c];
      P[irow] = mean_copy_number_weighted_cpp(a_vec, CHR_W);
      Time[irow] = 0.0;
      G[irow] = NA_REAL;
      DivisionTime[irow] = R_PosInf;
      continue;
    }

    // both viable
    for (int c = 0; c < 22; ++c) C[c][irow] = a_vec[c];

    X.push_back(xi);
    Y.push_back(yj);
    for (int c = 0; c < 22; ++c) {
      C[c].push_back(b_vec[c]);
    }
    R.push_back(R[irow]);
    D.push_back(D[irow]);

    P[irow] = mean_copy_number_weighted_cpp(a_vec, CHR_W);
    P.push_back(mean_copy_number_weighted_cpp(b_vec, CHR_W));

    G[irow] = NA_REAL;
    G.push_back(NA_REAL);

    WGDLabel.push_back(0);
    WGDCount.push_back(WGDCount[irow]);

    DivisionTime[irow] = R_PosInf;
    DivisionTime.push_back(R_PosInf);

    Time[irow] = 0.0;
    Time.push_back(0.0);

    Status.push_back(1);

    Quiescent.push_back(0);
    QuiescentTime.push_back(0.0);
    DeathReason.push_back(NA_STRING);

    grid(xi - 1, yj - 1) = X.size();
  }

  // -------- 3) Death rules --------
  int nrow2 = X.size();

  IntegerMatrix cur_kt_mat(nrow2, 22);
  for (int c = 0; c < 22; ++c) {
    IntegerVector col = C[c];
    for (int i = 0; i < nrow2; ++i) {
      cur_kt_mat(i, c) = col[i];
    }
  }

  // A) lethal 0-copy
  LogicalVector lethal_mask = lethal_zero_copy_mask_cpp(cur_kt_mat);
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 1 && lethal_mask[i]) {
      Status[i] = 0;
      if (DeathReason[i] == NA_STRING ||
          DeathReason[i] == "" ||
          DeathReason[i] == std::string("")) {
        DeathReason[i] = "lethal_zero_copy";
      }
    }
  }

  // B) ploidy over cap
  NumericVector P2 = P;
  LogicalVector over_mask = ploidy_over_cap_mask_cpp(P2, (double)ploidy_cap_death);
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 1 && over_mask[i]) {
      Status[i] = 0;
      if (DeathReason[i] == NA_STRING ||
          DeathReason[i] == "" ||
          DeathReason[i] == std::string("")) {
        DeathReason[i] = "ploidy_over_cap";
      }
    }
  }

  // C) Time > Dtime
  double Dtime = cfg.containsElementNamed("Dtime")
                   ? as<double>(cfg["Dtime"]) : R_PosInf;
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 1 && R_finite(Time[i]) && Time[i] > Dtime) {
      Status[i] = 0;
      if (DeathReason[i] == NA_STRING ||
          DeathReason[i] == "" ||
          DeathReason[i] == std::string("")) {
        DeathReason[i] = "over_Dtime";
      }
    }
  }

  // D) baseline stochastic death
  LogicalVector kill_mask = stochastic_death_mask_cpp(D);
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 1 && kill_mask[i]) {
      Status[i] = 0;
      if (DeathReason[i] == NA_STRING ||
          DeathReason[i] == "" ||
          DeathReason[i] == std::string("")) {
        DeathReason[i] = "baseline_death_Db";
      }
    }
  }

  // -------- 4) Recompute G / DivisionTime via compute_G_Div_quiescence_cpp --------
  IntegerVector alive_idx;
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 1) alive_idx.push_back(i + 1);
  }

  if (alive_idx.size() > 0) {
    IntegerVector X_a(alive_idx.size()), Y_a(alive_idx.size()), Q_a(alive_idx.size());
    NumericVector R_a(alive_idx.size()), P_a(alive_idx.size()), QT_a(alive_idx.size());

    for (int k = 0; k < alive_idx.size(); ++k) {
      int i = alive_idx[k] - 1;
      X_a[k] = X[i];
      Y_a[k] = Y[i];
      R_a[k] = R[i];
      P_a[k] = P[i];
      Q_a[k] = Quiescent[i];
      QT_a[k] = QuiescentTime[i];
    }

    int Kcap = 121;
    int radius = 5;
    double qthresh =
      cfg.containsElementNamed("QuiescentDeathHours")
        ? as<double>(cfg["QuiescentDeathHours"]) : 72.0;

    double softness =
      cfg.containsElementNamed("CrowdingSoftness")
        ? as<double>(cfg["CrowdingSoftness"]) : 0.15;
    if (!R_finite(softness) || softness <= 0.0) {
      softness = 0.15;
    }

    List res = compute_G_Div_quiescence_cpp(
      N, O2,
      X_a, Y_a,
      R_a, P_a,
      beta, Kcap, radius,
      Q_a, QT_a, qthresh,
      softness
    );

    NumericVector G2  = res["G"];
    NumericVector Div2 = res["DivisionTime"];
    IntegerVector Q2  = res["Quiescent"];
    IntegerVector AliveFlag = res["AliveStatus"];
    NumericVector QT2 = res["QuiescentTime"];

    for (int k = 0; k < alive_idx.size(); ++k) {
      int i = alive_idx[k] - 1;
      G[i] = G2[k];
      DivisionTime[i] = Div2[k];
      Quiescent[i] = Q2[k];
      QuiescentTime[i] = QT2[k];
      if (AliveFlag[k] == 0) {
        Status[i] = 0;
        if (DeathReason[i] == NA_STRING ||
            DeathReason[i] == "" ||
            DeathReason[i] == std::string("")) {
          DeathReason[i] = "quiescence_timeout";
        }
      }
    }
  }

  // -------- 5) Build per-step dead_log, remove dead cells, and rebuild grid --------
  IntegerVector dead_idx;
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 0) dead_idx.push_back(i + 1);
  }

  DataFrame dead_log;

  if (dead_idx.size() > 0) {
    for (int k = 0; k < dead_idx.size(); ++k) {
      int i = dead_idx[k] - 1;
      if (DeathReason[i] == NA_STRING ||
          DeathReason[i] == "" ||
          DeathReason[i] == std::string("")) {
        DeathReason[i] = "unspecified";
      }
    }

    IntegerVector Step(dead_idx.size(), step + 1);

    List dead_cells;
    dead_cells["X"] = X[dead_idx - 1];
    dead_cells["Y"] = Y[dead_idx - 1];

    for (int c = 0; c < 22; ++c) {
      IntegerVector col(dead_idx.size());
      IntegerVector cc = C[c];
      for (int k = 0; k < dead_idx.size(); ++k) {
        col[k] = cc[dead_idx[k] - 1];
      }
      dead_cells[c_names[c]] = col;
    }

    dead_cells["R"] = R[dead_idx - 1];
    dead_cells["D"] = D[dead_idx - 1];
    dead_cells["P"] = P[dead_idx - 1];
    dead_cells["G"] = G[dead_idx - 1];
    dead_cells["WGDLabel"] = WGDLabel[dead_idx - 1];
    dead_cells["WGDCount"] = WGDCount[dead_idx - 1];
    dead_cells["DivisionTime"] = DivisionTime[dead_idx - 1];
    dead_cells["Time"] = Time[dead_idx - 1];
    dead_cells["Status"] = Status[dead_idx - 1];
    dead_cells["Quiescent"] = Quiescent[dead_idx - 1];
    dead_cells["QuiescentTime"] = QuiescentTime[dead_idx - 1];
    dead_cells["DeathReason"] = DeathReason[dead_idx - 1];
    dead_cells["Step"] = Step;

    dead_log = DataFrame(dead_cells);
  } else {
    // Construct an empty dead_log (same columns as cells plus Step)
    List empty_cells;
    empty_cells["X"] = IntegerVector(0);
    empty_cells["Y"] = IntegerVector(0);
    for (int c = 0; c < 22; ++c) {
      empty_cells[c_names[c]] = IntegerVector(0);
    }
    empty_cells["R"] = NumericVector(0);
    empty_cells["D"] = NumericVector(0);
    empty_cells["P"] = NumericVector(0);
    empty_cells["G"] = NumericVector(0);
    empty_cells["WGDLabel"] = IntegerVector(0);
    empty_cells["WGDCount"] = IntegerVector(0);
    empty_cells["DivisionTime"] = NumericVector(0);
    empty_cells["Time"] = NumericVector(0);
    empty_cells["Status"] = IntegerVector(0);
    empty_cells["Quiescent"] = IntegerVector(0);
    empty_cells["QuiescentTime"] = NumericVector(0);
    empty_cells["DeathReason"] = CharacterVector(0);
    empty_cells["Step"] = IntegerVector(0);
    dead_log = DataFrame(empty_cells);
  }

  // Keep only alive rows
  IntegerVector keep_idx;
  for (int i = 0; i < nrow2; ++i) {
    if (Status[i] == 1) keep_idx.push_back(i + 1);
  }
  int n_alive = keep_idx.size();

  List keep_cells;
  keep_cells["X"] = X[keep_idx - 1];
  keep_cells["Y"] = Y[keep_idx - 1];
  for (int c = 0; c < 22; ++c) {
    IntegerVector col(n_alive);
    IntegerVector cc = C[c];
    for (int k = 0; k < n_alive; ++k) {
      col[k] = cc[keep_idx[k] - 1];
    }
    keep_cells[c_names[c]] = col;
  }
  keep_cells["R"] = R[keep_idx - 1];
  keep_cells["D"] = D[keep_idx - 1];
  keep_cells["P"] = P[keep_idx - 1];
  keep_cells["G"] = G[keep_idx - 1];
  keep_cells["WGDLabel"] = WGDLabel[keep_idx - 1];
  keep_cells["WGDCount"] = WGDCount[keep_idx - 1];
  keep_cells["DivisionTime"] = DivisionTime[keep_idx - 1];
  keep_cells["Time"] = Time[keep_idx - 1];
  keep_cells["Status"] = Status[keep_idx - 1];
  keep_cells["Quiescent"] = Quiescent[keep_idx - 1];
  keep_cells["QuiescentTime"] = QuiescentTime[keep_idx - 1];
  keep_cells["DeathReason"] = DeathReason[keep_idx - 1];

  DataFrame alive_df(keep_cells);

  IntegerVector X_final = as<IntegerVector>(alive_df["X"]);
  IntegerVector Y_final = as<IntegerVector>(alive_df["Y"]);

  IntegerMatrix new_grid = clone(grid);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      new_grid(i, j) = NA_INTEGER;
  for (int r = 0; r < X_final.size(); ++r) {
    int xi = X_final[r] - 1;
    int yj = Y_final[r] - 1;
    if (xi >= 0 && xi < N && yj >= 0 && yj < N) {
      new_grid(xi, yj) = r + 1;
    }
  }

  // Write updates back into state
  state["cells"]    = alive_df;
  state["grid"]     = new_grid;
  state["O2"]       = O2;
  state["dead_log"] = dead_log;
  state["step"]     = step + 1;

  return state;
}
