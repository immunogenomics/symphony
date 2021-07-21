#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//typedef arma::mat MATTYPE;
//typedef arma::vec VECTYPE;
//typedef arma::fmat MATTYPE;
//typedef arma::fvec VECTYPE;


// [[Rcpp::export]]
arma::mat exp_mean(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, int nrow, const arma::uvec& groups, const arma::uvec& group_sizes) {
    int ngroups = group_sizes.n_elem;
    arma::mat res = arma::zeros<arma::mat>(nrow, ngroups);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            res(i[j], groups[c]) += std::expm1(x[j]);
        }
    }
    
    for (int c = 0; c < ngroups; c++) {
        for (int r = 0; r < nrow; r++) {
            res(r, c) /= group_sizes[c];
        }
    }
        
    return(res);
}



// [[Rcpp::export]]
arma::mat log_vmr(const arma::vec& x, const arma::vec& p, const arma::vec& i, 
                  int ncol, int nrow, const arma::mat& means,
                  const arma::uvec& groups, const arma::uvec& group_sizes) {
    
    int ngroups = group_sizes.n_elem;
    arma::mat res = arma::zeros<arma::mat>(nrow, ngroups);
    arma::mat nnzero = arma::zeros<arma::mat>(nrow, ngroups);
    double tmp;
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            tmp = std::expm1(x[j]) - means(i[j], groups(c));
            res(i[j], groups[c]) += tmp * tmp;
            nnzero(i[j], groups(c))++;
        }
    }
    
    for (int c = 0; c < ngroups; c++) {
        for (int r = 0; r < nrow; r++) {
            res(r, c) += (group_sizes[c] - nnzero(r, c)) * means(r, c) * means(r, c);
            res(r, c) /= (group_sizes[c] - 1);
        }
    }
    
    res = log(res / means);
    res.replace(arma::datum::nan, 0);
    
    return(res);
}

// [[Rcpp::export]]
arma::vec normalizeCLR_dgc(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, int nrow, int margin) {    
    arma::vec res = x;
    if (margin == 1) {
        // first compute scaling factors for each row
        arma::vec geo_mean = arma::zeros<arma::vec>(nrow);
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                // i[j] gives the row num
                geo_mean(i[j]) += std::log1p(x[j]);
            }
        }
        for (int i = 0; i < nrow; i++) {
//            geo_mean(i) = (geo_mean(i) / (1 + ncol));    
            geo_mean(i) = std::exp(geo_mean(i) / ncol);    
        }
        // then  scale data
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                res(j) = std::log1p(res(j) / geo_mean(i[j]));
            }
        }        
    } else {
        // first compute scaling factors for each column
        arma::vec geo_mean = arma::zeros<arma::vec>(ncol);
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                geo_mean(c) += std::log1p(x[j]);
            }
            geo_mean(c) = std::exp(geo_mean(c) / nrow);
        }
        
        // then  scale data
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                res(j) = std::log1p(res(j) / geo_mean(c));
            }
        }        
        
    }
    
    return res;
}



// [[Rcpp::export]]
arma::mat scaleRowsWithStats_dgc(const arma::vec& x, const arma::vec& p, 
                                 const arma::vec& i, const arma::vec& mean_vec,
                                 const arma::vec& sd_vec, int ncol, int nrow, 
                                 float thresh) {
    // fill in non-zero elements
    arma::mat res = arma::zeros<arma::mat>(nrow, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res(i[j], c) = x(j);
        }
    }
    // scale rows with given means and SDs
    res.each_col() -= mean_vec;
    res.each_col() /= sd_vec;
    res.elem(find(res > thresh)).fill(thresh);
    res.elem(find(res < -thresh)).fill(-thresh);
    return res;
}


// [[Rcpp::export]]
arma::mat scaleRows_dgc(const arma::vec& x, const arma::vec& p, const arma::vec& i, 
                        int ncol, int nrow, float thresh) {
    // (0) fill in non-zero elements
    arma::mat res = arma::zeros<arma::mat>(nrow, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res(i[j], c) = x(j);
        }
    }

    // (1) compute means
    arma::vec mean_vec = arma::zeros<arma::vec>(nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            mean_vec(i[j]) += x[j];
        }
    }
    mean_vec /= ncol;
    
    // (2) compute SDs
    arma::vec sd_vec = arma::zeros<arma::vec>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
            nz(i[j])--;
        }
    }
        
    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }
    
    sd_vec = arma::sqrt(sd_vec / (ncol - 1));
    
    // (3) scale values
    res.each_col() -= mean_vec;
    res.each_col() /= sd_vec;
    res.elem(find(res > thresh)).fill(thresh);
    res.elem(find(res < -thresh)).fill(-thresh);
    return res;
}


// [[Rcpp::export]]
arma::vec rowMeansWeighted_dgc(const arma::vec& x, const arma::vec& p, 
                     const arma::vec& i, const arma::vec& weights,
                     int ncol, int nrow) {

    arma::vec res = arma::zeros<arma::vec>(nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res[i[j]] += x[j] * weights[c];
        }
    }
    
    res /= arma::accu(weights);    
    return res;
}

// [[Rcpp::export]]
arma::vec rowSDs_dgc(const arma::vec& x, const arma::vec& p, 
                     const arma::vec& i, const arma::vec& mean_vec, 
                     int ncol, int nrow, bool do_sqrt) {

    arma::vec sd_vec = arma::zeros<arma::vec>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
            nz(i[j])--;
        }
    }
    
    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }
    
    sd_vec = sd_vec / (ncol - 1);
    if (do_sqrt) {
        sd_vec = arma::sqrt(sd_vec);
    }
    
    return sd_vec;    
}


// [[Rcpp::export]]
arma::vec rowVarSDs_dgc(
    const arma::vec& x, const arma::vec& p, 
    const arma::vec& i, const arma::vec& mean_vec, const arma::vec& sd_vec, 
    double vmax, int ncol, int nrow, bool do_sqrt) {

    arma::vec res = arma::zeros<arma::vec>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    double val;
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            val = std::min(vmax, (x[j] - mean_vec(i[j])) / sd_vec(i[j]));
            res(i[j]) += val * val; // [(x - mu)/sig]^2
            nz(i[j])--;
        }
    }
    
    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        res(r) += nz(r) * mean_vec(r) * mean_vec(r) / (sd_vec(r) * sd_vec(r));
    }
    
    res = res / (ncol - 1);
    if (do_sqrt) {
        res = arma::sqrt(res);
    }
    
    return res;    
}


// [[Rcpp::export]]
arma::vec rowSDsWeighted_dgc(const arma::vec& x, const arma::vec& p, 
                     const arma::vec& i, const arma::vec& mean_vec, 
                     const arma::vec& weights, 
                     int ncol, int nrow, bool do_sqrt) {

    arma::vec sd_vec = arma::zeros<arma::vec>(nrow);
    double sum_weights = arma::accu(weights);
    arma::vec nz = arma::zeros<arma::vec>(nrow);
    nz.fill(sum_weights);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            sd_vec(i[j]) += weights[c] * (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // w * (x - mu)^2
            nz(i[j]) -= weights[c];
        }
    }
    
    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }

    sd_vec *= sum_weights / (sum_weights * sum_weights - arma::accu(weights % weights));
    if (do_sqrt) {
        sd_vec = arma::sqrt(sd_vec);
    }
    return sd_vec;
}


// [[Rcpp::export]]
arma::mat cosine_normalize_cpp(arma::mat & V, int dim) {
  // norm rows: dim=1
  // norm cols: dim=0 or dim=2
  if (dim == 2) dim = 0;
  return arma::normalise(V, 2, dim);
}

// [[Rcpp::export]]
List soft_kmeans_cpp(arma::mat Y, arma::mat Z, unsigned max_iter, float sigma) {
    Y = arma::normalise(Y, 2, 0); // L2 normalize the columns
    Z = arma::normalise(Z, 2, 0); // L2 normalize the columns
    arma::mat R;// = -2 * (1 - Y.t() * Z) / sigma; // dist_mat 
    for (unsigned i = 0; i < max_iter; i++) {
        R = -2 * (1 - Y.t() * Z) / sigma; // dist_mat 
        R.each_row() -= arma::max(R, 0);  
        R = exp(R);
        R.each_row() /= arma::sum(R, 0);
        Y = arma::normalise(Z * R.t(), 2, 0); 
    }

    List result = List::create(Named("R") = R , _["Y"] = Y);
    return result;
    
}
