// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>

using namespace boost::math;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::mat mydigamma(const arma::mat x) {
  unsigned int N = x.n_rows;
  unsigned int M = x.n_cols;
  arma::mat A(N,M);
  for (unsigned int i=0;i<M;i++){
    for(unsigned int j=0;j<N;j++){
      A(j,i) = digamma(x(j,i));
    }
  }
  return A;
}

// [[Rcpp::export]]
arma::mat mytrigamma(const arma::mat x) {
  unsigned int N = x.n_rows;
  unsigned int M = x.n_cols;
  arma::mat A(N,M);
  for (unsigned int i=0;i<M;i++){
    for(unsigned int j=0;j<N;j++){
      A(j,i) = Rf_trigamma(x(j,i));
    }
  }
  return A;
}

// [[Rcpp::export]]
arma::mat mylgamma(const arma::mat x) {
  unsigned int N = x.n_rows;
  unsigned int M = x.n_cols;
  arma::mat A(N,M);
  for (unsigned int i=0;i<M;i++){
    for(unsigned int j=0;j<N;j++){
      A(j,i) = lgamma(x(j,i));
    }
  }
  return A;
}


//[[Rcpp::export]]
arma::mat N_ab(const arma::mat x, const arma::rowvec ab){
    unsigned int N = x.n_rows;
    unsigned int Nrep = x.n_cols;
    unsigned int K = ab.n_elem;
    unsigned int M = Nrep*K;
    arma::mat A(N,K*Nrep);
    for (unsigned int i=0;i<M;i++){
        A.col(i) = x.col(i%Nrep) + ab(i/Nrep);
    }
    return A;
}


  //[[Rcpp::export]]
  arma::mat cmpLgm( arma::mat x, unsigned int K){
      unsigned int N = x.n_rows;
      unsigned int M = x.n_cols;
      unsigned int Nrep = M/K;
      x = mylgamma(x);
      arma::mat A(N,K);
      A.zeros();
     for (unsigned int i=0;i<M;i++){
         A.col(i/Nrep) = A.col(i/Nrep) + x.col(i);
     }
      return A;
  }

//[[Rcpp::export]]
arma::mat cmpDgm( arma::mat x, unsigned int K){
    unsigned int N = x.n_rows;
    unsigned int M = x.n_cols;
    unsigned int Nrep = M/K;
    x = mydigamma(x);
    arma::mat A(N,K);
    A.zeros();
    for (unsigned int i=0;i<M;i++){
        A.col(i/Nrep) = A.col(i/Nrep) + x.col(i);
    }
    return A;
}

//[[Rcpp::export]]
arma::mat cmpTgm( arma::mat x, unsigned int K){
    unsigned int N = x.n_rows;
    unsigned int M = x.n_cols;
    unsigned int Nrep = M/K;
    x = mytrigamma(x);
    arma::mat A(N,K);
    A.zeros();
    for (unsigned int i=0;i<M;i++){
        A.col(i/Nrep) = A.col(i/Nrep) + x.col(i);
    }
    return A;
}


//[[Rcpp::export]]
arma::mat ab_sv(arma::mat ab, double s, arma::vec v, unsigned int K){
  arma::mat out = ab;
  for(unsigned int i=0;i<2*K;i++){
    out(i%2, i/2) += v(i)*s;
  }
  return out;
}

//[[Rcpp::export]]
arma::vec logvec(arma::vec v, unsigned int N){
  arma::vec out = v;
  for (unsigned int i=0;i<N;i++){
    out(i) = log(v(i));
  }
  return out;
}


struct CVXOBJ {
   arma::vec g;
   arma::mat H;
   double f;
};

CVXOBJ cvxObj(double t, arma::mat Postprob, arma::rowvec wght, arma::mat x, arma::mat y, arma::mat n, arma::mat ab, 
                 arma::mat F, arma::vec h, unsigned int N, unsigned int Nrep ){
   unsigned int m = F.n_rows;
   unsigned int K = Postprob.n_cols;
   arma::vec d(m);
   arma::mat D(m,m);
   arma::mat J(2,K);
   arma::mat H_d(2,K);
   arma::mat H_o(2,K);
   arma::mat n_ab(N,Nrep*K);
   arma::mat x_a(N,Nrep*K);
   arma::mat y_b(N,Nrep*K);
   arma::mat dn_ab(K,K);
   arma::mat dx_a(K,K);
   arma::mat dy_b(K,K);
   arma::mat px(N,K);
   // List out(3);
   double f;
   arma::vec g(2*K);
   arma::mat H(2*K,2*K);
   J.zeros();
   H_d.zeros();
   H_o.zeros();

   d = F * vectorise(ab) - h;
   D = diagmat(1 / d);
   n_ab = N_ab(n,sum(ab,0));
   x_a = N_ab(x,ab.row(0));
   y_b = N_ab(y,ab.row(1));

   J =mydigamma( ones(2,1) * sum(ab,0) ) % (ones(2,1)*wght) * Nrep;
   cmpDgm(N_ab(n,sum(ab,0)),K).t() * Postprob;
  mydigamma(ab) % (ones(2,1)*wght) *Nrep;
   dn_ab = cmpDgm(n_ab,K).t() * Postprob;
   dx_a = cmpDgm(x_a,K) .t() * Postprob;
   dy_b = cmpDgm(y_b,K) .t() * Postprob;

   J =mydigamma( ones(2,1) * sum(ab,0) ) % (ones(2,1)*wght) * Nrep -  
      mydigamma(ab) % (ones(2,1)*wght) *Nrep;
   for (unsigned int i=0;i<K;i++){
     J(0,i) += dx_a(i,i) - dn_ab(i,i);
     J(1,i) += dy_b(i,i) - dn_ab(i,i);
   }
   // compute the trigamma for hessian
   // repeatedly use dn_ab, dx_a, dy_b; H_d represents the diagnal; H_o off diagnal of H
   dn_ab = cmpTgm(n_ab,K).t() * Postprob;
   dx_a = cmpTgm(x_a,K) .t() * Postprob;
   dy_b = cmpTgm(y_b,K) .t() * Postprob;

   H_d = mytrigamma( ones(2,1) * sum(ab,0) ) % (ones(2,1)*wght) * Nrep -  
       mytrigamma(ab) % (ones(2,1)*wght) *Nrep;
   for (unsigned int i=0;i<K;i++){
     H_d(0,i) += dx_a(i,i) - dn_ab(i,i); //
     H_d(1,i) += dy_b(i,i) - dn_ab(i,i);
   }
   H_o = ( mytrigamma( ones(2,1) * sum(ab,0) ) % (ones(2,1)*wght) * Nrep ) - ones(2,1) * dn_ab.diag().t();
   H = diagmat(vectorise(H_d));
   for(unsigned int i=0;i<K;i++){
       H(2*i,1+2*i) = H_o(1,i);
       H(1+2*i,2*i) = H_o(1,i);
   }

   // likelihood computation 
   px = Nrep*mylgamma(ones(N,1)*sum(ab,0)) - Nrep*mylgamma(ones(N,1)*ab.row(0)) - Nrep*mylgamma(ones(N,1)*ab.row(1))
         - cmpLgm(n_ab,K) + cmpLgm(x_a,K) + cmpLgm(y_b,K);

   g = -t * vectorise(J) - F.t() * D * ones(m,1); // g first dierivatives of the function
   H = -t * H;        // H Hessian
   for (unsigned int i=0;i<m;i++){
      d(i) =  log(-d(i));
   }
   f = -t * accu(px % Postprob) - sum(d); // f value for output function

   return {g,H,f};

}

//[[Rcpp::export]]
List cmpHmm(const arma::mat x, arma::mat y, arma::mat trans, arma::mat ab, arma::mat F, arma::vec h){
 unsigned int Nrep = x.n_cols;
 unsigned int K = ab.n_cols;
 unsigned int N = x.n_rows;
 unsigned int Nit = 30;
 double t_step;
 double mu = 1e6;
 double PRECISION = 1e8;
 double t_INITIAL = 1;
 double ALPHA = 0.1;
 double BETA = 0.7;
 double NTTOL = 1e-6;
 double s;
 double lambda;
 List out(4);

 CVXOBJ cvx;
 CVXOBJ cvx_new;

 arma::mat n = x + y;
 arma::mat Dens(N,K);
 arma::mat Forward(N,K);
 arma::mat Backward(N,K);
 arma::mat Postprob(N,K);
 arma::rowvec wght(K);
 arma::vec Scale(N);
 Backward.row(N-1).fill(0);
 Scale(0) = 1.0;
 arma::vec logl;
 logl.zeros(Nit+1);
 arma::vec v(2*K); // inverse(H)*J

 for(unsigned int nit=0;nit<Nit;nit++){
      
      Dens = exp(  Nrep*mylgamma(ones(N,1)*sum(ab,0)) - Nrep*mylgamma(ones(N,1)*ab.row(0)) - Nrep*mylgamma(ones(N,1)*ab.row(1))
         - cmpLgm(N_ab(n,sum(ab,0)),K) + cmpLgm(N_ab(x,ab.row(0)),K) + cmpLgm(N_ab(y,ab.row(1)),K)
         + mylgamma(n+1.0) * ones(Nrep,K) - mylgamma(x+1.0) * ones(Nrep,K) - mylgamma(y+1.0) * ones(Nrep,K) 
         );

      Forward.row(0) =  Dens.row(0)/K;
      for (unsigned int n=1;n<N;n++){
         Forward.row(n) = ( Forward.row(n-1) * trans ) % Dens.row(n);
         Scale(n) = sum(Forward.row(n));
         Forward.row(n) = Forward.row(n) / Scale(n);
      }
      logl(nit+1) = log(sum(Forward.row(N-1))) + sum(logvec(Scale,N)); // caculate the likelhood remember end at N-1

      Backward.row(N-1) += 1;
      for (int n=N-2;n>=0;n--){
         Backward.row(n) = ( Backward.row(n+1) % Dens.row(n+1) ) * trans.t();
         Backward.row(n) = Backward.row(n) / Scale(n);
      }
      
      trans = trans % ( Forward.rows(0,N-2).t() * ( Dens.rows(1,N-1) % Backward.rows(1,N-1) ) );
      trans = trans / ( sum(trans,1) * ones(1,K) );
      
      Postprob = Forward % Backward;
      Postprob = Postprob / ( sum(Postprob,1) * ones(1,K) );
      wght = sum(Postprob,0);
      
      t_step = t_INITIAL;
      while (t_step < PRECISION ){
      // for (unsigned int i=0;i<2;i++){

        cvx = cvxObj(t_step,Postprob,wght,x,y,n,ab,F,h,N,Nrep);
        v = -solve( cvx.H,cvx.g );
        lambda = sum(v % cvx.g);
        
        s = 1;
        while(min(h - F * (vectorise(ab) + s*v)) < 0) {
          s = BETA*s; 
        }
        
        cvx_new = cvxObj(t_step,Postprob,wght,x,y,n,ab_sv(ab,s,v,K),F,h,N,Nrep);
        while (cvx_new.f > cvx.f + ALPHA*s*lambda) {
          s = BETA * s;
          cvx_new = cvxObj(t_step,Postprob,wght,x,y,n,ab_sv(ab,s,v,K),F,h,N,Nrep);
        }

        if (std::abs(lambda/2) < NTTOL) {break;}
        ab = ab_sv(ab,s,v,K);
       t_step = mu * t_step; 
      }
      
      if (std::abs( (logl(nit+1) - logl(nit)) / logl(nit+1) ) < NTTOL) {break;}
 }

out[0] = ab;
out[1] = trans;
out[2] = Postprob;
out[3] = logl;

 return out;
}



// //[[Rcpp::export]]
// arma::mat dens1(const arma::mat x, const arma::mat y, const::mat n, const arma::mat ab){
//  unsigned int N = x.n_rows;
//  unsigned int Nrep = x.n_cols;
//  unsigned int K = ab.n_cols;
//  arma::mat A(N,K); 
//  A = exp(  
//        Nrep*mylgamma(ones(N,1)*sum(ab,0)) - Nrep*mylgamma(ones(N,1)*ab.row(0)) - Nrep*mylgamma(ones(N,1)*ab.row(1))
//          - cmpLgm(N_ab(n,sum(ab,0)),K) 
//           + cmpLgm(N_ab(x,ab.row(0)),K) 
//           + cmpLgm(N_ab(y,ab.row(1)),K)
//          + mylgamma(n+1.0) * ones(Nrep,K) - mylgamma(x+1.0) * ones(Nrep,K) - mylgamma(y+1.0) * ones(Nrep,K) 
//          );
//  //sum(ab,0): colSums; sum(ab,1): rowSums;
//  return A;
// }


// //[[Rcpp::export]]
// arma::mat forward1(const arma::mat dens, arma::mat trans){
//   unsigned int N = dens.n_rows;
//   int K = trans.n_cols;
//   arma::mat Forward(N,K);
//   arma::mat Postprob(N,K);
//   arma::vec Scale(N);
//   Scale(0) = 1.0;
//   Forward.row(0) =  dens.row(0)/K;

//   for (unsigned int n=1;n<N;n++){
//     Forward.row(n) = ( Forward.row(n-1) * trans ) % dens.row(n);
//     Scale(n) = sum(Forward.row(n));
//     Forward.row(n) = Forward.row(n) / Scale(n);
//   }
//   return Scale;
// }


// //[[Rcpp::export]]
// arma::mat frbk1(const arma::mat dens, arma::mat trans){
//   unsigned int N = dens.n_rows;
//   int K = trans.n_cols;
//   arma::mat Forward(N,K);
//   arma::mat Backward(N,K);
//   arma::mat Postprob(N,K);
//   arma::vec Scale(N);
//   Scale(0) = 1.0;

//   Forward.row(0) =  dens.row(0)/K;
//   for (unsigned int n=1;n<N;n++){
//     Forward.row(n) = ( Forward.row(n-1) * trans ) % dens.row(n);
//     Scale(n) = sum(Forward.row(n));
//     Forward.row(n) = Forward.row(n) / Scale(n);
//   }

//   Backward.row(N-1) += 1;
//   for (int n=N-2;n>=0;n--){
//     Backward.row(n) = ( Backward.row(n+1) % dens.row(n+1) ) * trans.t();
//     Backward.row(n) = Backward.row(n) / Scale(n);
//   }

//   trans = trans % ( Forward.rows(0,N-2).t() * ( dens.rows(1,N-1) % Backward.rows(1,N-1) ) );
//   trans = trans / ( sum(trans,1) * ones(1,K) );

//   Postprob = Forward % Backward;
//   Postprob = Postprob / ( sum(Postprob,1) * ones(1,K) );

//   return Postprob;
// }

