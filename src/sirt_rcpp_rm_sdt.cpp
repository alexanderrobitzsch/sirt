//// File Name: sirt_rcpp_rm_sdt.cpp
//// File Version: 0.423


// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
// #include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// user includes


///********************************************************************
///** array3_index
// [[RcppNOexport]]
int array3_index( int ii1, int ii2, int ii3, int I1, int I2 )
{
    int index = ii1+ii2*I1+ii3*I1*I2;
    return index;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_rm_sdt_posterior
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_rm_sdt_posterior( Rcpp::NumericVector prob_item,
        Rcpp::NumericVector prob_rater, int I0, int K1, int TP,
        Rcpp::IntegerVector item_index, Rcpp::IntegerMatrix dat2,
        Rcpp::LogicalMatrix dat2_resp, int I, Rcpp::NumericVector pi_k0)
{
    int N = dat2.nrow();
    int NU = I0*K1*TP;
    Rcpp::NumericVector hjikl(NU);
    Rcpp::NumericVector xjikl(NU);
    Rcpp::NumericMatrix gjil(I0,TP);
    Rcpp::NumericVector wjl(TP);
    Rcpp::NumericVector nik_item(NU);
    Rcpp::NumericVector nik_rater(I*K1*K1);
    Rcpp::NumericMatrix fqkyi(N,TP);
    Rcpp::NumericMatrix fyiqk(N,TP);
    Rcpp::NumericVector pi_k(TP);
    Rcpp::NumericVector like(N);

    pi_k.fill(0);
    like.fill(0);
    double eps=1e-30;
    double temp1=0;
    double temp2=0;
    int vv_temp=0;
    int hh_temp=0;
    double p1=0;
    double temp=0;
    double ll=0;
    double eps_like=1e-300;

    for (int nn=0; nn<N; nn++){
        //*** hjikl
            //    xjikl <- hjikl <- array(0, dim=c(I0,K1,TP) )
            //    for (ii in 1:I0){
            //        hjikl[ii,,] <- prob.item[ii,,]
            //    }
        for (int uu=0; uu<NU; uu++){
            hjikl[uu] = prob_item[uu];
        }
            // # prob.rater[ item, observed_category, true_category]
            // for (jj in 1:I0){
            //    if ( dat2_resp[nn,jj] ){
            //        for (kk in 1:K1){
            //            p1 <- prob.rater[ jj, dat2[nn,jj] + 1, kk]
            //            hjikl[ item.index[jj], kk,] <- hjikl[ item.index[jj], kk,] * p1
            //    }
            // }
        for( int jj=0; jj<I; jj++){
            if ( dat2_resp(nn,jj) ){
                for (int kk=0; kk<K1; kk++){
                    hh_temp = array3_index( jj, dat2(nn,jj), kk, I, K1);
                    p1 = prob_rater[hh_temp];
                    for (int tt=0; tt<TP; tt++){
                        vv_temp = array3_index( item_index[jj], kk, tt, I0, K1);
                        hjikl[vv_temp] = hjikl[vv_temp] * p1;
                    }
                }
            }
        }
            // gjil <- array(hjikl[,1,], dim=c(I0,TP) )
            // for (kk in 2:K1){
            //    gjil <- gjil + hjikl[,kk,]
            // }
            // hjikl <- array(0, dim=c(I0,K1,TP) )
        gjil.fill(0);
        for (int ii=0; ii<I0; ii++){
            for (int tt=0; tt<TP; tt++){
                for (int kk=0; kk<K1; kk++){
                    gjil(ii,tt) += hjikl[ ii+kk*I0+tt*I0*K1];
                }
            }
        }
            // wjl <- rep(0,TP)
            // for (tt in 1:TP){    wjl[tt] <- pi.k[tt] * prod(gjil[,tt])}
            // wjl <- wjl / sum(wjl)
        wjl.fill(0);
        temp2=0;
        for (int tt=0; tt<TP; tt++){
            temp1=1;
            for (int ii=0; ii<I0; ii++){
                temp1 = temp1*gjil(ii,tt);
            }
            wjl[tt] = pi_k0[tt] * temp1;
            temp2 += wjl[tt];
        }
        like[nn] = temp2;
        ll += std::log( like[nn]+eps_like );
        for (int tt=0; tt<TP; tt++){
            fyiqk(nn,tt) = wjl[tt];
            wjl[tt] = wjl[tt] / ( temp2 + eps );
            fqkyi(nn,tt) = wjl[tt];
            pi_k[tt] += wjl[tt];
        }
            //    for (kk in 1:K1){
            //         xjikl[,kk,] <- hjikl[,kk,] / gjil
            //    }
            //    nik_item <- array(0, dim=c(I0, K1, TP) )
            // for (tt in 1:TP){
            //    nik_item[,,tt] <- nik_item[,,tt] + wjl[tt] * xjikl[,,tt]
            // }
        for (int ii=0; ii<I0; ii++){
            for (int tt=0; tt<TP; tt++){
                for (int kk=0; kk<K1; kk++){
                    vv_temp = array3_index( ii, kk, tt, I0, K1);
                    xjikl[ vv_temp ] = hjikl[ vv_temp ] / ( gjil(ii,tt) + eps );
                    nik_item[ vv_temp ] += wjl[tt] * xjikl[ vv_temp ];
                }
            }
        }
            //for (jj in 1:I){
            //    for (cc in 1:K1){
            //        temp <- 0
            //        for (tt in 1:TP){
            //            temp <- temp + wjl[tt] * xjikl[ item.index[jj],cc,tt]
            //        }
            //        for (oo in 1:K1){
            //            nik_rater[jj,oo,cc] <- nik_rater[jj,oo,cc] +
            //                                dat2.ind.resp[nn,jj,oo] * temp
            //        }
            //     }
            // }
        for (int jj=0; jj<I; jj++){
            if ( dat2_resp(nn,jj) ){
                for (int cc=0; cc<K1; cc++){
                    temp=0;
                    for (int tt=0; tt<TP; tt++){
                        temp += wjl[tt]*xjikl[ item_index[jj]+cc*I0+tt*I0*K1 ];
                    }
                    nik_rater[ jj+dat2(nn,jj)*I+cc*I*K1 ] += temp;
                }
            }
        }

    } // end nn
    for (int tt=0; tt<TP; tt++){
        pi_k[tt] = pi_k[tt] / N;
    }
    //--- output
    return Rcpp::List::create(
            Rcpp::Named("nik_item") = nik_item,
            Rcpp::Named("nik_rater") = nik_rater,
            Rcpp::Named("fqkyi") = fqkyi,
            Rcpp::Named("fyiqk") = fyiqk,
            Rcpp::Named("pi_k") = pi_k,
            Rcpp::Named("like") = like,
            Rcpp::Named("ll") = ll
        );
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_rm_sdt_calc_gradient_item_deriv_a
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_rm_sdt_calc_gradient_item_deriv_a(
        Rcpp::NumericVector prob0, Rcpp::IntegerVector prob_D1_dim,
        Rcpp::NumericVector theta_k )
{
    int VV = prob_D1_dim[0];
    int K1 = prob_D1_dim[1];
    int TP = prob_D1_dim[2];
    Rcpp::NumericVector prob_D1(VV*K1*TP);
    Rcpp::NumericMatrix exp_val(VV,TP);
    //    exp_val <- prob0[,1,]
    for (int ii=0; ii<VV; ii++){
        //    for (kk in 2:(K+1) ){
        //        exp_val <- exp_val + prob0[,kk,]*kk
        //    }
        for (int kk=1; kk<K1; kk++){
            for (int tt=0; tt<TP; tt++){
                exp_val(ii,tt) += kk*prob0[ ii + kk*VV + tt*VV*K1 ];
            }
        }
        //    for (kk in 1:(K+1)){
        //        prob_D1[,kk,] <- kk - exp_val
        //    }
        for (int tt=0; tt<TP; tt++){
            for (int kk=0; kk<K1; kk++){
                prob_D1[ ii + kk*VV + tt*VV*K1 ] = theta_k[tt] * ( kk - exp_val(ii,tt) );
            }
        }
        //    for (tt in 1:TP){
        //        prob_D1[,,tt] <- prob_D1[,,tt]*theta.k[tt]
        //    }
    }
    return prob_D1;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_rm_sdt_calc_gradient_item_deriv_tau
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_rm_sdt_calc_gradient_item_deriv_tau(
            Rcpp::NumericVector prob0, Rcpp::IntegerVector prob_D1_dim,
            int cat_pp)
{
    int VV = prob_D1_dim[0];
    int K1 = prob_D1_dim[1];
    int TP = prob_D1_dim[2];
    Rcpp::NumericVector prob_D1(VV*K1*TP);
        //     for (kk in 1:K1){
        //        prob_D1[,kk,] <- prob0[,cat_pp,]
        //        if ( kk == cat_pp ){
        //            prob_D1[,kk,] <- -1 + prob_D1[,kk,]
        //        }
        //    }
    int index=0;
    for (int vv=0; vv<VV; vv++){
        for (int tt=0; tt<TP; tt++){
            index = vv + cat_pp*VV + tt*VV*K1;
            for (int kk=0; kk<K1; kk++){
                prob_D1[ vv + kk*VV + tt*VV*K1 ] = prob0[ index ];
            }
            prob_D1[ index ] += -1;
        }
    }
    return prob_D1;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_rm_sdt_calc_gradient_likelihood_item_llgrad
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_rm_sdt_calc_gradient_likelihood_item_llgrad(
        Rcpp::NumericVector logprob_D1, Rcpp::IntegerVector prob_D1_dim,
        Rcpp::NumericVector nik_item )
{
    int VV = prob_D1_dim[0];
    int K1 = prob_D1_dim[1];
    int TP = prob_D1_dim[2];
    Rcpp::NumericVector ll_grad(VV);
        // ll <- rowSums( logprob_D1[,1,] * nik.item[,1,] )
        // for (kk in 2:(K+1) ){
        //    ll <- ll + rowSums( logprob_D1[,kk,] * nik.item[,kk,] )
        // }
    int ind=0;
    for (int vv=0; vv<VV; vv++){
        for (int kk=0; kk<K1; kk++){
            for (int tt=0; tt<TP; tt++){
                ind = vv + kk*VV + tt*VV*K1;
                ll_grad[vv] += logprob_D1[ind] * nik_item[ind];
            }
        }
    }
    return ll_grad;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_squeeze_eps
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_squeeze_eps( Rcpp::NumericVector x, double eps)
{
    int NH = x.size();
    Rcpp::NumericVector y = Rcpp::clone(x);
    for (int hh=0; hh<NH; hh++){
        if ( x[hh] < eps ){
            y[hh] = eps;
        }
    }
    //---- output
    return y;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_log
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_log( Rcpp::NumericVector x)
{
    int NH = x.size();
    Rcpp::NumericVector y(NH);
    for (int hh=0; hh<NH; hh++){
        y[hh] = std::log( x[hh] );
    }
    //---- output
    return y;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_plogis
// [[Rcpp::export]]
double sirt_rcpp_plogis( double x)
{
    double y=0;
    double z=0;
    if (x >= 0){
        z = std::exp(-x);
        y = 1 / ( 1 + z );
    } else {
        z = std::exp(x);
        y = z / ( 1 + z );
    }
    //---- output
    return y;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_rm_sdt_calc_probs_gpcm
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_rm_sdt_calc_probs_gpcm( Rcpp::NumericVector a,
        Rcpp::NumericMatrix tau, Rcpp::NumericVector theta_k, int VV, int K1,
        int TP, double eps, bool use_log )
{
    int NH = VV*K1*TP;
    Rcpp::NumericVector probs(NH);
    Rcpp::NumericVector p1(K1);
    double max_temp=0;
    double val=0;
    double sumval=0;
    for (int vv=0; vv<VV; vv++){
        for (int tt=0; tt<TP; tt++){
            max_temp=0;
            p1[0]=0;
            for (int kk=1; kk<K1; kk++){
                val = a[vv] * kk * theta_k[tt] - tau(vv,kk-1);
                p1[kk] = val;
                if (val > max_temp){
                    max_temp = val;
                }
            }
            sumval=0;
            for (int kk=0; kk<K1; kk++){
                p1[kk] = std::exp(p1[kk] - max_temp);
                sumval += p1[kk];
            }
            for (int kk=0; kk<K1; kk++){
                probs[ vv + kk*VV + tt*VV*K1 ] = p1[kk] / sumval;
            }
        }
    }
    if ( eps > 0 ){
        probs = sirt_rcpp_squeeze_eps( probs, eps);
    }
    if ( use_log ){
        probs = sirt_rcpp_log(probs);
    }
    return probs;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_rm_sdt_calc_probs_grm_rater
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_rm_sdt_calc_probs_grm_rater( Rcpp::NumericMatrix c_rater,
        Rcpp::NumericVector d_rater, int I, int K, double eps, bool use_log )
{
    int K1 = K + 1;
    int NH = I*K1*K1;
    Rcpp::NumericVector prob_cum(NH);
    Rcpp::NumericVector prob_rater(NH);
    double temp=0;
    //    for (eta in 0:K){
    //        for (k in 1:K){
    //            prob_cum2[ k, eta+1 ] <- plogis(c.rater[ii,k] - d.rater[ii]*eta)
    //        }
    //        prob_cum2[K1,] <- 1
    //    }
    for (int ii=0; ii<I; ii++){
        for (int cc=0; cc<K1; cc++){
            for (int kk=0; kk<K; kk++){
                temp=c_rater(ii,kk) - d_rater[ii]*cc;
                prob_cum[ ii + kk*I + cc*I*K1 ] = sirt_rcpp_plogis(temp);
            }
            prob_cum[ ii + K*I + cc*I*K1 ] = 1.0;
        }
    }
    // category probabilities
    int index = 0;
    int index0 = 0;
    for (int ii=0; ii<I; ii++){
        for (int cc=0; cc<K1; cc++){
            index = ii + cc*I*K1;
            prob_rater[ index ] = prob_cum[ index ];
            index0 = index;
            for (int kk=1; kk<K1; kk++){
                index = ii+kk*I + cc*I*K1;
                prob_rater[ index ] = prob_cum[ index ] - prob_cum[ index0 ];
                index0 = index;
            }
        }
    }
    if ( eps > 0 ){
        prob_rater = sirt_rcpp_squeeze_eps( prob_rater, eps);
    }
    if ( use_log ){
        prob_rater = sirt_rcpp_log(prob_rater);
    }
    //---- output
    return prob_rater;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_rm_sdt_calc_probs_grm_item
// [[Rcpp::export]]
Rcpp::NumericVector sirt_rcpp_rm_sdt_calc_probs_grm_item( Rcpp::NumericMatrix tau_item,
        Rcpp::NumericVector a_item, Rcpp::NumericVector theta_k, int VV, int K, int TP,
        double eps, bool use_log )
{
    int K1 = K+1;
    int NH = VV*K1*TP;
    Rcpp::NumericVector prob_cum(NH);
    Rcpp::NumericVector prob_item(NH);
    double val_ii_tt = 0;
    double temp=0;
    for (int ii=0; ii<VV; ii++){
        for (int tt=0; tt<TP; tt++){
            val_ii_tt = a_item[ii]*theta_k[tt];
            for (int kk=0; kk<K; kk++){
                temp=tau_item(ii,kk) - val_ii_tt;
                prob_cum[ ii + kk*VV + tt*VV*K1 ] = sirt_rcpp_plogis(temp);
            }
            prob_cum[ ii + K*VV + tt*VV*K1 ] = 1.0;
        }
    }
    int index=0;
    int index2=0;
    for (int ii=0; ii<VV; ii++){
        for (int tt=0; tt<TP; tt++){
            index = ii + tt*VV*K1;
            prob_item[ index ] = prob_cum[ index ];
            for (int kk=1; kk<K1; kk++){
                index2 = ii + kk*VV + tt*VV*K1;
                prob_item[ index2 ] = prob_cum[ index2 ] - prob_cum[ index ];
                index = index2;
            }
        }
    }
    if ( eps > 0 ){
        prob_item = sirt_rcpp_squeeze_eps( prob_item, eps);
    }
    if ( use_log ){
        prob_item = sirt_rcpp_log(prob_item);
    }

    //---- output
    return prob_item;
}
///********************************************************************
