#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;


// Symphony utils----------------------------------

// [[Rcpp::depends(RcppArmadillo)]]

// Computes the soft cluster assignments of query cells across reference clusters.
// 
// Y: Reference cluster centroid locations
// Z: Query cells projected into PC space (d x m)
// sigma: Soft k-means "fuzziness" parameter, sigma = 0 is hard clustering
// [[Rcpp::export]]
arma::mat soft_cluster(const arma::mat& Y, const arma::mat& Z, float sigma) {
    arma::mat Y_cos = arma::normalise(Y, 2, 0); // L2 normalize the columns
    arma::mat Z_cos = arma::normalise(Z, 2, 0); // L2 normalize the columns
    arma::mat R = -2 * (1 - Y_cos.t() * Z_cos) / sigma; // dist_mat 

    R.each_row() -= arma::max(R, 0);  
    R = exp(R);
    R.each_row() /= arma::sum(R, 0);    
    return R;
}


// Computes the Symphony reference compression terms, Nr and C.
// 
// Rr: Soft cluster assignments of reference cells (cols) across clusters (rows).
// Zr: Corrected embedding for reference cells (cols) in harmonized PCs (rows).
// [[Rcpp::export]]
List compute_ref_cache(
    const arma::mat& Rr, 
    const arma::mat& Zr
) {
    List result(2);

    result[0] = arma::sum(Rr, 1); // Nr (k x 1)
    result[1] = Rr * Zr.t();      // C (k x d)
    return result;
}

// Computes the corrected query cell embedding.
// 
// Zq: Query cells projected into PC space (d x m)
// Xq: Query design matrix ((c + 1) x m)
// Rq: Query soft cluster assignments across reference clusters (k x m)
// Nr: Reference cluster sizes (first compression term) (length k)
// RrZtr: Second reference compression term (C) (k x d)
// [[Rcpp::export]]
arma::mat moe_correct_ref(
    const arma::mat& Zq, // query cells projected into PC space
    const arma::mat& Xq, // query design matrix
    const arma::mat& Rq, // query soft cluster assignments
    const arma::vec& Nr, // ref cluster sizes
    const arma::mat& RrZtr // ref matrix cached
) {
    unsigned K = Rq.n_rows;
    arma::mat Zq_corr = Zq;
    arma::mat Xq_Rk, beta;
    arma::mat mat1, mat2, lambda_I;
    for (unsigned k = 0; k < K; k++) { 
        Xq_Rk = Xq * arma::diagmat(Rq.row(k));

        // (B+1) x (B+1)
        mat1 = Xq_Rk * Xq.t();
        mat1(0,0) += Nr[k];
        
        // ridge
        unsigned nrows = Xq.n_rows;
        lambda_I = arma::eye(nrows, nrows); 
        lambda_I(0,0) -= 1; //do not penalize the intercept
        mat1 += lambda_I;
        
        // (B+1) x d
        mat2 = Xq_Rk * Zq.t();
        mat2.row(0) += RrZtr.row(k);
        
        beta = arma::inv(mat1) * (mat2);
        beta.row(0).zeros(); //do not correct the intercept terms
        Zq_corr -= beta.t() * Xq_Rk;
    } 
    return Zq_corr; //(d x m)
}

// Returns the batch coefficients of the linear mixture model as a 3D tensor.
// [[Rcpp::export]]
arma::cube get_betas(const arma::mat& R, const arma::mat& Z, const arma::mat& lambda, const arma::mat& design) {
  unsigned K = R.n_rows;
  unsigned B = design.n_rows;
  unsigned D = Z.n_rows;
  unsigned N = Z.n_cols;
  arma::cube W_cube(B, D, K); // rows, cols, slices
  arma::mat Phi_Rk(B, N);
  for (unsigned k = 0; k < K; k++) { 
    Phi_Rk = design * arma::diagmat(R.row(k));
    W_cube.slice(k) = arma::inv(Phi_Rk * design.t() + lambda) * Phi_Rk * Z.t();
  }
  return W_cube;
}
