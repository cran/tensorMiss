#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <math.h>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
Eigen::MatrixXd K3_D_nu(int k, Eigen::MatrixXd D, Eigen::MatrixXd Q,
                        NumericVector dim, NumericVector C, NumericVector Y,
                        int jj, int nu = 0) {
  // this function is adapted to data with K=3 for each time point
  // kk: mode to work on, starting from 2 (using the way that R counts) (as mode 1 is time mode)
  // D: r_k by r_k eigenvalue matrix of sample covariance matrix
  // Q: estimated mode-k loading matrix
  // dim: vector specifying the dimension of C (and hence Y)
  // C: vectorisation of estimated common component
  // Y: vectorisation of data
  // j: row index for D_{k,nu,j}, start from 0
  // nu: nu in D_{k,nu,j}

  int r_k = D.rows();
  
  int dim1 = dim[0];
  int dim2 = dim[1];
  int dim3 = dim[2];
  int dim4 = dim[3];
  int dimkm1 = dim[k-1];
  
  Eigen::Tensor<double, 4> tensor(dim1, dim2, dim3, dim4);
  
  for (int i=0; i < dim1; i++){
    for (int j=0; j < dim2; j++){
      for (int s=0; s < dim3; s++){
        for (int l=0; l < dim4; l++){
          tensor(i,j,s,l) = Y(i + j*dim1 + s*dim1*dim2 + l*dim1*dim2*dim3);
        }
      }
    }
  }
  
  Eigen::Tensor<double, 4> C_tensor(dim1, dim2, dim3, dim4);
  
  for (int i=0; i < dim1; i++){
    for (int j=0; j < dim2; j++){
      for (int s=0; s < dim3; s++){
        for (int l=0; l < dim4; l++){
          C_tensor(i,j,s,l) = C(i + j*dim1 + s*dim1*dim2 + l*dim1*dim2*dim3);
        }
      }
    }
  }
  // finish constructing the data tensor by now -----------------------------------------------
  // finish constructing the estimated common component tensor by now -------------------------
  
  
  // create unfolding first, for each t -----------------------------------------------
  int rows = dimkm1;
  long cols = tensor.size() / rows;
  cols /= dim1;
  
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> matrixList(dim1);

  for (int t = 0; t < dim1; t++) {
    // Create an Eigen Matrix to store the unfolded data
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> k_mat(rows, cols);

    // Unfold the tensor into a matrix
    for (int u1 = 0; u1 < rows; u1++) {
      for (long u2 = 0; u2 < cols; u2++) {
        long index = u1 * cols + u2;
        Eigen::array<Eigen::Index, 4> indices;
        indices[0] = t;
        indices[k-1] = u1;
        for (int w = 1; w < 4; w++) {
          if ( w != (k-1)) {
            indices[w] = int( fmod(index, dim[w]) );
            index /= dim[w];
          }
        }

        k_mat(u1, u2) = tensor(indices);
      }
    }
    //matrixList.push_back(k_mat);
    matrixList[t] = k_mat;
  }
  
  
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> C_matrixList(dim1);
  
  for (int t = 0; t < dim1; t++) {
    // Create an Eigen Matrix to store the unfolded data
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> k_mat(rows, cols);
    
    // Unfold the tensor into a matrix
    for (int u1 = 0; u1 < rows; u1++) {
      for (long u2 = 0; u2 < cols; u2++) {
        long index = u1 * cols + u2;
        Eigen::array<Eigen::Index, 4> indices;
        indices[0] = t;
        indices[k-1] = u1;
        for (int w = 1; w < 4; w++) {
          if ( w != (k-1)) {
            indices[w] = int( fmod(index, dim[w]) );
            index /= dim[w];
          }
        }
        
        k_mat(u1, u2) = C_tensor(indices);
      }
    }
    C_matrixList[t] = k_mat;
  }
  

  // create psi list -----------------------------------------------------------
  std::vector<std::vector<std::vector<bool> > > psi_list(dimkm1);

  for (int i=0; i < dimkm1; i++){
    std::vector<std::vector<bool>> i_list(cols);
    for (long h=0; h < cols; h++){
      std::vector<bool> psi_ih(dim1);
      for (int t=0; t < dim1; t++){
        bool t_ih_obs = !( NumericVector::is_na( matrixList[t](i,h) ) );
        bool t_jh_obs = !( NumericVector::is_na( matrixList[t](jj,h) ) );

        psi_ih[t] = (t_ih_obs and t_jh_obs);
      }
      i_list[h] = psi_ih;
    }
    psi_list[i] = i_list;
  }

  // create psi matrix for nu = 0, to obtain |psi_{k,ij,h}|
  //std::vector< Eigen::MatrixXd > psi_mat_list(dim1);
  Eigen::MatrixXd psi_mat = Eigen::MatrixXd::Constant(dimkm1, cols, 0.0);
  for (int i=0; i < dimkm1; i++){
    for (long h=0; h < cols; h++){
      int psi_count = 0;
      for (int t=0; t < dim1; t++){
        bool t_ih_obs = !( NumericVector::is_na( matrixList[t](i,h) ) );
        bool t_jh_obs = !( NumericVector::is_na( matrixList[t](jj,h) ) );

        if (t_ih_obs and t_jh_obs){
          psi_count += 1;
        }
      }
      psi_mat(i,h) = psi_count;
    }
  }

  
  
  Eigen::MatrixXd DQ;
  DQ = D.inverse() * Q.transpose();
  
  // construct H part in D_{k,nu,j} for each i----------------------------------
  std::vector< Eigen::VectorXd > H_i_list(dimkm1);
  for (int i=0; i < dimkm1; i++){
    Eigen::VectorXd C_start = Eigen::VectorXd::Constant(dimkm1, 0.0);
    for (int t = 0; t < dim1; t++){
      Eigen::Map<Eigen::MatrixXd> C_hat_t(C_matrixList[t].data(), dimkm1, cols);
      Eigen::VectorXd C_hat_t_row(cols);
      C_hat_t_row = C_hat_t.row(i);
      
      C_start += C_hat_t * C_hat_t_row;
    }
    C_start /= dim1;
    
    Eigen::VectorXd H_i = Eigen::VectorXd::Constant(r_k, 0.0);
    H_i = DQ * C_start;
    
    H_i_list[i] = H_i;
  }
  

  
  // construct final D_{k,nu,j}-------------------------------------------------
  Eigen::MatrixXd D_nu_j = Eigen::MatrixXd::Constant(r_k, r_k, 0.0);
  //Eigen::VectorXd D_nu_j = Eigen::VectorXd::Constant(r_k, 0.0);

  for (int t=0; t < (dim1-nu); t++){
    // initialise first part of D_{k,nu,j} for time t
    Eigen::VectorXd D_start1 = Eigen::VectorXd::Constant(r_k, 0.0);
    for (int i=0; i < dimkm1; i++){
      Eigen::VectorXd H_i;
      H_i = H_i_list[i];
      
      Eigen::VectorXd C_hat_t_row;
      C_hat_t_row = C_matrixList[t].row(i);

      Eigen::VectorXd Y_hat_t_row;
      Eigen::VectorXd YC_hat_t_row;
      Y_hat_t_row = matrixList[t].row(jj);
      YC_hat_t_row = C_matrixList[t].row(jj);

      Eigen::VectorXd E_hat_t_row;
      E_hat_t_row = Y_hat_t_row - YC_hat_t_row;
      
      double i_sum = 0.0;
      for (int h=0; h < cols; h++){
        double h_sum = 0.0;
        
        if ( !(psi_list[i][h][t]) ){
          i_sum += h_sum;
        } else {
          h_sum += C_hat_t_row[h] * E_hat_t_row[h];
          h_sum /= psi_mat(i,h);
          i_sum += h_sum;
        }
      }

      D_start1 += H_i * i_sum;
    }
    

    // initialise second part of D_{k,nu,j} for time t+nu
    int tnu;
    tnu = t+nu;
    Eigen::VectorXd D_start2 = Eigen::VectorXd::Constant(r_k, 0.0);
    for (int i=0; i < dimkm1; i++){
      Eigen::VectorXd H_i;
      H_i = H_i_list[i];
      
      Eigen::VectorXd C_hat_t_row;
      C_hat_t_row = C_matrixList[tnu].row(i);
      
      Eigen::VectorXd Y_hat_t_row;
      Eigen::VectorXd YC_hat_t_row;
      Y_hat_t_row = matrixList[tnu].row(jj);
      YC_hat_t_row = C_matrixList[tnu].row(jj);
      
      Eigen::VectorXd E_hat_t_row;
      E_hat_t_row = Y_hat_t_row - YC_hat_t_row;
      
      double i_sum = 0.0;
      for (int h=0; h < cols; h++){
        double h_sum = 0.0;
        
        if ( !(psi_list[i][h][tnu]) ){
          i_sum += h_sum;
        } else {
          h_sum += C_hat_t_row[h] * E_hat_t_row[h];
          h_sum /= psi_mat(i,h);
          i_sum += h_sum;
        }
      }
      
      D_start2 += H_i * i_sum;
    }

    D_nu_j += D_start1 * D_start2.transpose();
  }
  
  
  return D_nu_j;
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//



/*** R
#a <- K3_D_nu(2, diag((svd(est_loading[[1]])$d)[(1:r1)]), A1.est, dim(Y.imp), c(Y.imp), c(ten_t_re), 1, 1)

#a <- K3_D_nu(2, diag((svd(est_loading[[1]])$d)[(1:r1)]), A1.est, dim(Y.imp), c(Y.imp), c(ten_comp), 0, 0)

*/
