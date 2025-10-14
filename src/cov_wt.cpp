#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Helper: weighted or unweighted crossproduct
arma::mat crossprod_mat(const arma::mat& X, const arma::vec& wt = arma::vec()) {
  arma::mat Y = X;
  if (wt.n_elem > 0) {
    Y.each_col() %= arma::sqrt(wt);
  }
  return arma::trans(Y) * Y;
}

// [[Rcpp::export]]
arma::mat crossprod_cpp(const arma::mat& X) {
  return arma::trans(X) * X;
}

// [[Rcpp::export]]
List cov_wt_cpp(const arma::mat& x,
                Rcpp::Nullable<Rcpp::NumericVector> wt_ = R_NilValue,
                bool cor = false,
                Rcpp::Nullable<Rcpp::NumericVector> center_ = R_NilValue,
                std::string method = "unbiased") {

  int n = x.n_rows;

  // 1. weights
  arma::vec wt;
  bool with_wt = !wt_.isNull();
  if (with_wt) {
    wt = Rcpp::as<arma::vec>(wt_);
    wt /= arma::sum(wt);
  } else {
    wt = arma::vec(n, arma::fill::ones);
    wt /= static_cast<double>(n);
  }

  // 2. center
  arma::rowvec center;
  if (center_.isNull()) {
    center = arma::sum(x.each_col() % wt, 0);
  } else {
    center = Rcpp::as<arma::rowvec>(center_);
  }

  // 3. subtract center
  arma::mat X = x;
  X.each_row() -= center;

  // 4. compute covariance using helper
  arma::mat cov;
  if (method == "ML") {
    cov = crossprod_mat(X, wt);
  } else { // unbiased
    double denom = 1.0 - arma::dot(wt, wt);
    cov = crossprod_mat(X, wt) / denom;
  }

  // 5. build output
  List out = List::create(
    Named("cov") = cov,
    Named("center") = center,
    Named("n.obs") = n
  );

  if (with_wt) {
    out["wt"] = Rcpp::NumericVector(wt.begin(), wt.end());
  }

  // 6. correlation if requested
  if (cor) {
    arma::vec invsd = 1.0 / arma::sqrt(cov.diag());
    arma::mat R = cov;
    R.each_col() %= invsd;
    R.each_row() %= invsd.t();
    out["cor"] = R;
  }

  return out;
}

