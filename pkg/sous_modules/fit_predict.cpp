#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

double matern_kernel_normal(const Eigen::VectorXd& ligne_quanti1,
                            const Eigen::VectorXd& ligne_post1,
                            const Eigen::VectorXd& ligne_quanti2,
                            const Eigen::VectorXd& ligne_post2,
                            const Eigen::VectorXd& w_vec) {
  double coef1 = sqrt(5);
  double coef2 = 3.0;
  int p = ligne_quanti1.size();
  int q = ligne_post1.size();
  double prod = 1.0;
  for (int i = 0; i < p; ++i) {
    double b = (coef1 * abs(ligne_quanti1[i] - ligne_quanti2[i])) / w_vec[i];
    prod *= (1 + b + (1/coef2) * pow(b, 2)) * exp(-b);
  }
  for (int j = 0; j < q; ++j) {
    double b = coef1 * abs(ligne_post1[j] - ligne_post2[j]);
    prod *= (1 + b + (1/coef2) * pow(b, 2)) * exp(-b);
  }
  
  return prod;
}

double matern_kernel_quanti(const Eigen::VectorXd& ligne_quanti1,
                            const Eigen::VectorXd& ligne_quanti2,
                            const Eigen::VectorXd& w_vec) {
  double coef1 = sqrt(5);
  double coef2 = 3.0;
  int p = ligne_quanti1.size();
  double prod = 1.0;
  for (int i = 0; i < p; ++i) {
    double b = (coef1 * abs(ligne_quanti1[i] - ligne_quanti2[i])) / w_vec[i];
    prod *= (1 + b + (1/coef2) * pow(b, 2)) * exp(-b);
  }
  return prod;
}

double matern_kernel_quali(const Eigen::VectorXd& ligne_post1,
                            const Eigen::VectorXd& ligne_post2) {
  double coef1 = sqrt(5);
  double coef2 = 3.0;
  int q = ligne_post1.size();
  double prod = 1.0;
  for (int j = 0; j < q; ++j) {
    double b = coef1 * abs(ligne_post1[j] - ligne_post2[j]);
    prod *= (1 + b + (1/coef2) * pow(b, 2)) * exp(-b);
  }
  
  return prod;
}

double exp_kernel_normal(const Eigen::VectorXd& ligne_quanti1, const Eigen::VectorXd& ligne_post1, const Eigen::VectorXd& ligne_quanti2, const Eigen::VectorXd& ligne_post2, const Eigen::VectorXd& w_vec) {
  int p = ligne_quanti1.size();
  double sum_squared_diff = (ligne_post1 - ligne_post2).squaredNorm();
  
  double sum_weighted_diff = 0.0;
  for (int k = 0; k < p; ++k) {
    sum_weighted_diff += w_vec[k] * pow(ligne_quanti1[k] - ligne_quanti2[k], 2);
  }
  
  return exp(-sum_squared_diff - sum_weighted_diff);
}

double exp_kernel_quanti(const Eigen::VectorXd& ligne_quanti1, const Eigen::VectorXd& ligne_quanti2, const Eigen::VectorXd& w_vec) {
  int p = ligne_quanti1.size();
  double sum_weighted_diff = 0.0;
  
  for (int k = 0; k < p; ++k) {
    sum_weighted_diff += w_vec[k] * pow(ligne_quanti1[k] - ligne_quanti2[k], 2);
  }
  
  return exp(-sum_weighted_diff);
}

double exp_kernel_quali(const Eigen::VectorXd& ligne_post1, const Eigen::VectorXd& ligne_post2) {
  int p = ligne_post1.size();
  double sum_squared_diff = (ligne_post1 - ligne_post2).squaredNorm();
  
  return exp(-sum_squared_diff);
}

Eigen::MatrixXd compute_R_normal_matern(const Eigen::MatrixXd& quanti_matrix,
                                     const Eigen::MatrixXd& post_matrix,
                                     const Eigen::VectorXd& w_vec,
                                     const double& reg) {
  int n = quanti_matrix.rows();
  Eigen::MatrixXd R(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R(i, j) = matern_kernel_normal(quanti_matrix.row(i), post_matrix.row(i),
        quanti_matrix.row(j), post_matrix.row(j), w_vec);
    }
    R(i,i) += reg; 
  }
  
  return R;
}

Eigen::MatrixXd compute_R_quanti_matern(const Eigen::MatrixXd& quanti_matrix, const Eigen::VectorXd& w_vec, const double& reg) {
  int n = quanti_matrix.rows();
  Eigen::MatrixXd R(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R(i, j) = matern_kernel_quanti(quanti_matrix.row(i), quanti_matrix.row(j), w_vec);
    }
    R(i,i) += reg; 
  }
  
  return R;
}

Eigen::MatrixXd compute_R_quali_matern(const Eigen::MatrixXd& post_matrix, const double& reg) {
  int n = post_matrix.rows();
  Eigen::MatrixXd R(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R(i, j) = matern_kernel_quali(post_matrix.row(i), post_matrix.row(j));
    }
    R(i,i) += reg; 
  }
  return R;
}

Eigen::MatrixXd compute_R_normal_exp(const Eigen::MatrixXd& quanti_matrix,
                                     const Eigen::MatrixXd& post_matrix,
                                     const Eigen::VectorXd& w_vec,
                                     const double& reg) {
  int n = quanti_matrix.rows();
  Eigen::MatrixXd R(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R(i, j) = exp_kernel_normal(quanti_matrix.row(i), post_matrix.row(i),
        quanti_matrix.row(j), post_matrix.row(j), w_vec);
    }
    R(i,i) += reg; 
  }
  
  return R;
}

Eigen::MatrixXd compute_R_quanti_exp(const Eigen::MatrixXd& quanti_matrix, const Eigen::VectorXd& w_vec, const double& reg) {
  int n = quanti_matrix.rows();
  Eigen::MatrixXd R(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R(i, j) = exp_kernel_quanti(quanti_matrix.row(i), quanti_matrix.row(j), w_vec);
    }
    R(i,i) += reg; 
  }
  
  return R;
}

Eigen::MatrixXd compute_R_quali_exp(const Eigen::MatrixXd& post_matrix, const double& reg) {
  int n = post_matrix.rows();
  Eigen::MatrixXd R(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      R(i, j) = exp_kernel_quali(post_matrix.row(i), post_matrix.row(j));
    }
    R(i,i) += reg; 
  }
  return R;
}

Eigen::MatrixXd compute_g_normal_matern(const Eigen::MatrixXd& new_quanti_matrix,
                                     const Eigen::MatrixXd& new_post_matrix,
                                     const Eigen::MatrixXd& quanti_matrix,
                                     const Eigen::MatrixXd& post_matrix,
                                     const Eigen::VectorXd& w_vec,
                                     const double& sigma2) {
  int n = quanti_matrix.rows();
  int n2 = new_quanti_matrix.rows();
  Eigen::MatrixXd g(n, n2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n2; ++j) {
      g(i, j) = sigma2*matern_kernel_normal(quanti_matrix.row(i), post_matrix.row(i),
        new_quanti_matrix.row(j), new_post_matrix.row(j), w_vec);
    }
  }
  return g;
} 

Eigen::MatrixXd compute_g_quanti_matern(const Eigen::MatrixXd& new_quanti_matrix, 
                                     const Eigen::MatrixXd& quanti_matrix,
                                     const Eigen::VectorXd& w_vec,
                                     const double& sigma2) {
  int n = quanti_matrix.rows();
  int n2 = new_quanti_matrix.rows();
  Eigen::MatrixXd g(n, n2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n2; ++j) {
      g(i, j) = sigma2*matern_kernel_quanti(quanti_matrix.row(i), new_quanti_matrix.row(j), w_vec);
    }
  }
  return g;
}

Eigen::MatrixXd compute_g_quali_matern(const Eigen::MatrixXd& new_post_matrix,
                                    const Eigen::MatrixXd& post_matrix,
                                    const double& sigma2) {
  int n = post_matrix.rows();
  int n2 = new_post_matrix.rows();
  Eigen::MatrixXd g(n, n2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n2; ++j) {
      g(i, j) = sigma2*matern_kernel_quali(post_matrix.row(i), new_post_matrix.row(j));
    }
  }
  return g;
}

Eigen::MatrixXd compute_g_normal_exp(const Eigen::MatrixXd& new_quanti_matrix,
                                     const Eigen::MatrixXd& new_post_matrix,
                                     const Eigen::MatrixXd& quanti_matrix,
                                     const Eigen::MatrixXd& post_matrix,
                                     const Eigen::VectorXd& w_vec,
                                     const double& sigma2) {
  int n = quanti_matrix.rows();
  int n2 = new_quanti_matrix.rows();
  Eigen::MatrixXd g(n, n2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n2; ++j) {
      g(i, j) = sigma2*exp_kernel_normal(quanti_matrix.row(i), post_matrix.row(i),
        new_quanti_matrix.row(j), new_post_matrix.row(j), w_vec);
    }
  }
  return g;
} 

Eigen::MatrixXd compute_g_quanti_exp(const Eigen::MatrixXd& new_quanti_matrix, 
                                     const Eigen::MatrixXd& quanti_matrix,
                                     const Eigen::VectorXd& w_vec,
                                     const double& sigma2) {
  int n = quanti_matrix.rows();
  int n2 = new_quanti_matrix.rows();
  Eigen::MatrixXd g(n, n2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n2; ++j) {
      g(i, j) = sigma2*exp_kernel_quanti(quanti_matrix.row(i), new_quanti_matrix.row(j), w_vec);
    }
  }
  return g;
}

Eigen::MatrixXd compute_g_quali_exp(const Eigen::MatrixXd& new_post_matrix,
                                    const Eigen::MatrixXd& post_matrix,
                                    const double& sigma2) {
  int n = post_matrix.rows();
  int n2 = new_post_matrix.rows();
  Eigen::MatrixXd g(n, n2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n2; ++j) {
      g(i, j) = sigma2*exp_kernel_quali(post_matrix.row(i), new_post_matrix.row(j));
    }
  }
  return g;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd invertMatrixChol(const Eigen::MatrixXd& mat) {
  Eigen::LLT<Eigen::MatrixXd> lltOfA(mat);
  Eigen::MatrixXd invMat = lltOfA.solve(Eigen::MatrixXd::Identity(mat.rows(), mat.cols()));
  return invMat;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd cpp_beta_hat(const Eigen::MatrixXd& F_matrix,
                             const Eigen::MatrixXd& inv_R_matrix,
                             const Eigen::VectorXd& Y_vec) {
  Eigen::MatrixXd t_F_inv_R = F_matrix.transpose() * inv_R_matrix;
  Eigen::MatrixXd temp = invertMatrixChol(t_F_inv_R * F_matrix);
  Eigen::VectorXd beta_hat = temp * (t_F_inv_R * Y_vec);
  return beta_hat;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double cpp_sigma_hat(const Eigen::MatrixXd& F_matrix,
                     const Eigen::MatrixXd& inv_R_matrix,
                     const Eigen::VectorXd& beta_vec,
                     const Eigen::VectorXd& Y_vec) {
  int n = Y_vec.size();
  Eigen::VectorXd M = Y_vec - F_matrix * beta_vec;
  double sigma_hat = (1.0 / n) * M.transpose() * inv_R_matrix * M;
  return sigma_hat;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpp_compute_R_normal(const Eigen::MatrixXd& quanti_matrix, const Eigen::MatrixXd& post_matrix, const Eigen::VectorXd& w_vec, const std::string& kernel, const double& reg) {
  Eigen::MatrixXd res;
  if (kernel == "exp") {
    res = compute_R_normal_exp(quanti_matrix, post_matrix, w_vec, reg);
  } else if (kernel == "mat") {
    res = compute_R_normal_matern(quanti_matrix, post_matrix, w_vec, reg);
  } else {
    Rcpp::stop("Invalid kernel type. Supported types are 'exp' and 'mat'.");
  }
  return res;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpp_compute_R_quanti(const Eigen::MatrixXd& quanti_matrix, const Eigen::VectorXd& w_vec, const std::string& kernel, const double& reg) {
  Eigen::MatrixXd res;
  if (kernel == "exp") {
    res = compute_R_quanti_exp(quanti_matrix, w_vec, reg);
  } else if (kernel == "mat") {
    res = compute_R_quanti_matern(quanti_matrix, w_vec, reg);
  } else {
    Rcpp::stop("Invalid kernel type. Supported types are 'normal', 'quanti', and 'quali'.");
  }
  
  return res;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpp_compute_R_quali(const Eigen::MatrixXd& post_matrix, const std::string& kernel, const double& reg) {
  Eigen::MatrixXd res;
  if (kernel == "exp") {
    res = compute_R_quali_exp(post_matrix, reg);
  } else if (kernel == "mat") {
    res = compute_R_quali_matern(post_matrix, reg);
  } else {
    Rcpp::stop("Invalid kernel type. Supported types are 'normal', 'quanti', and 'quali'.");
  }
  return res;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double L_normal_cpp(const Eigen::MatrixXd& quanti_matrix,
                    const Eigen::MatrixXd& prior_matrix,
                    const Eigen::MatrixXd& F_matrix,
                    const std::string& kernel_str,
                    const Eigen::MatrixXd& A_matrix,
                    const Eigen::VectorXd& w_vec,
                    const Eigen::VectorXd& Y_vec,
                    const double& reg) {
  int n = Y_vec.size();
  Eigen::MatrixXd post_matrix = prior_matrix * A_matrix;
  Eigen::MatrixXd R_matrix = cpp_compute_R_normal(quanti_matrix, post_matrix, w_vec, kernel_str, reg);
  Eigen::LLT<Eigen::MatrixXd> lltOfA(R_matrix);
  Eigen::MatrixXd inv_R_matrix = lltOfA.solve(Eigen::MatrixXd::Identity(R_matrix.rows(), R_matrix.cols()));
  double detR = lltOfA.matrixLLT().diagonal().array().square().prod();
  Eigen::VectorXd beta_vec = cpp_beta_hat(F_matrix, inv_R_matrix, Y_vec);
  double sig = cpp_sigma_hat(F_matrix, inv_R_matrix, beta_vec, Y_vec);
  double lik = n * log(sig) + log(detR);
  return lik;
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double L_quanti_cpp(const Eigen::MatrixXd& quanti_matrix,
                    const Eigen::MatrixXd& F_matrix,
                    const std::string& kernel_str,
                    const Eigen::VectorXd& w_vec,
                    const Eigen::VectorXd& Y_vec,
                    const double& reg) {
  int n = Y_vec.size();
  Eigen::MatrixXd R_matrix = cpp_compute_R_quanti(quanti_matrix, w_vec, kernel_str, reg);
  Eigen::LLT<Eigen::MatrixXd> lltOfA(R_matrix);
  Eigen::MatrixXd inv_R_matrix = lltOfA.solve(Eigen::MatrixXd::Identity(R_matrix.rows(), R_matrix.cols()));
  double detR = lltOfA.matrixLLT().diagonal().array().square().prod();
  Eigen::VectorXd beta_vec = cpp_beta_hat(F_matrix, inv_R_matrix, Y_vec);
  double sig = cpp_sigma_hat(F_matrix, inv_R_matrix, beta_vec, Y_vec);
  double lik = n * log(sig) + log(detR);
  return lik;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double L_quali_cpp(const Eigen::MatrixXd& prior_matrix,
                   const Eigen::MatrixXd& F_matrix,
                   const std::string& kernel_str,
                   const Eigen::MatrixXd& A_matrix,
                   const Eigen::VectorXd& Y_vec,
                   const double& reg) {
  int n = Y_vec.size();
  Eigen::MatrixXd post_matrix = prior_matrix * A_matrix;
  Eigen::MatrixXd R_matrix = cpp_compute_R_quali(post_matrix, kernel_str, reg);
  Eigen::LLT<Eigen::MatrixXd> lltOfA(R_matrix);
  Eigen::MatrixXd inv_R_matrix = lltOfA.solve(Eigen::MatrixXd::Identity(R_matrix.rows(), R_matrix.cols()));
  double detR = lltOfA.matrixLLT().diagonal().array().square().prod();
  Eigen::VectorXd beta_vec = cpp_beta_hat(F_matrix, inv_R_matrix, Y_vec);
  double sig = cpp_sigma_hat(F_matrix, inv_R_matrix, beta_vec, Y_vec);
  double lik = n * log(sig) + log(detR);
  return lik;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpp_compute_g_normal(const Eigen::MatrixXd& new_quanti_matrix, 
                                     const Eigen::MatrixXd& new_post_matrix, 
                                     const Eigen::MatrixXd& quanti_matrix, 
                                     const Eigen::MatrixXd& post_matrix,
                                     const Eigen::VectorXd& w_vec, 
                                     const std::string& kernel,
                                     const double& sigma2) {
  Eigen::MatrixXd res;
  if (kernel == "exp") {
    res = compute_g_normal_exp(new_quanti_matrix, 
                               new_post_matrix, 
                               quanti_matrix, 
                               post_matrix,
                               w_vec,
                               sigma2);
  } else if (kernel == "mat") {
    res = compute_g_normal_matern(new_quanti_matrix, 
                               new_post_matrix, 
                               quanti_matrix, 
                               post_matrix,
                               w_vec,
                               sigma2);
  } else { 
    Rcpp::stop("Invalid kernel type. Supported types are 'exp' and 'mat'.");
  } 
  return res;
} 

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpp_compute_g_quanti(const Eigen::MatrixXd& new_quanti_matrix, 
                                     const Eigen::MatrixXd& quanti_matrix, 
                                     const Eigen::VectorXd& w_vec, 
                                     const std::string& kernel,
                                     const double& sigma2) {
  Eigen::MatrixXd res;
  if (kernel == "exp") {
    res = compute_g_quanti_exp(new_quanti_matrix,quanti_matrix, w_vec, sigma2);
  } else if (kernel == "mat") { 
    res = compute_g_quanti_matern(new_quanti_matrix,quanti_matrix, w_vec, sigma2);
  } else { 
    Rcpp::stop("Invalid kernel type. Supported types are 'normal', 'quanti', and 'quali'.");
  } 
  
  return res;
} 

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd cpp_compute_g_quali(const Eigen::MatrixXd& new_post_matrix,
                                    const Eigen::MatrixXd& post_matrix,
                                    const std::string& kernel,
                                    const double& sigma2) {
  Eigen::MatrixXd res;
  if (kernel == "exp") {
    res = compute_g_quali_exp(new_post_matrix,post_matrix, sigma2);
  } else if (kernel == "mat") { 
    res = compute_g_quali_matern(new_post_matrix,post_matrix, sigma2);
  } else { 
    Rcpp::stop("Invalid kernel type. Supported types are 'normal', 'quanti', and 'quali'.");
  } 
  return res;
} 


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericVector cpp_var_gp(const Eigen::MatrixXd& G,
                         const Eigen::MatrixXd& inv_V,
                         const Eigen::MatrixXd& Fnew,
                         const Eigen::MatrixXd& f_v,
                         const Eigen::MatrixXd& inv_f_v_f,
                         const double& sigma2) {
  int p = G.cols();
  NumericVector var_vec(p);
  
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd G_col = G.col(j);
    Eigen::VectorXd Fnew_row = Fnew.row(j);
    Eigen::VectorXd h = Fnew_row - f_v * G_col;
    
    double term1 = G_col.transpose() * inv_V * G_col;
    double term2 = h.transpose() * inv_f_v_f * h;
    
    var_vec[j] = sigma2 - term1 + term2;
  }
  
  return var_vec;
}