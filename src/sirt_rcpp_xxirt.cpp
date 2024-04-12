//// File Name: sirt_rcpp_xxirt.cpp
//// File Version: 0.485



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


///********************************************************************
///** sirt_rcpp_xxirt_compute_posterior_expected_counts
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_xxirt_compute_posterior_expected_counts(
        Rcpp::LogicalMatrix dat1_resp_gg, Rcpp::NumericMatrix p_aj_xi_gg,
        Rcpp::NumericVector weights_gg)
{
    int N = dat1_resp_gg.nrow();
    int I = dat1_resp_gg.ncol();
    int TP = p_aj_xi_gg.ncol();

    Rcpp::NumericMatrix nij(I, TP);
    double val=0;

    //*** loop over items and categories
    for (int ii=0;ii<I;ii++){
        for (int tt=0;tt<TP;tt++){
            val=0;
            for (int nn=0;nn<N;nn++){
                if ( dat1_resp_gg(nn,ii) ){
                    val +=  weights_gg(nn)*p_aj_xi_gg(nn,tt);
                }
            }  // end nn
            nij(ii,tt) = val;
        }   // end tt
    }   // end ii

    //--- OUTPUT
    return nij;
}
///********************************************************************

///********************************************************************
///** sirt_rcpp_xxirt_compute_likelihood
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_xxirt_compute_likelihood(
        Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
        Rcpp::NumericMatrix probs, int TP, int maxK )
{
    int N = dat.nrow();
    int I = dat.ncol();
    Rcpp::NumericMatrix p_xi_aj(N, TP);

    for (int nn=0;nn<N;nn++){
        for (int tt=0;tt<TP;tt++){
            p_xi_aj(nn,tt) = 1;
        }
        for (int ii=0;ii<I;ii++){
            if ( dat_resp_bool(nn,ii) ){
                for (int tt=0;tt<TP;tt++){
                    p_xi_aj(nn,tt) = p_xi_aj(nn,tt) * probs(ii, dat(nn,ii) + tt*maxK );
                }
            }
        }
    }

    //---- OUTPUT
    return p_xi_aj;
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_xxirt_hessian_reduced_probs
// [[Rcpp::export]]
Rcpp::NumericMatrix sirt_rcpp_xxirt_hessian_reduced_probs(
        Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
        Rcpp::NumericMatrix probs_ratio, int TP, int maxK,
        int itemnr, int itemnr2, bool use_itemnr2, Rcpp::NumericMatrix p_xi_aj )
{
    int N = dat.nrow();
    Rcpp::NumericMatrix t1(N, TP);

    int ii=itemnr;

    for (int nn=0;nn<N;nn++){
        if (dat_resp_bool(nn,itemnr)){
            for (int tt=0; tt<TP; tt++){
                t1(nn,tt) = probs_ratio(ii, dat(nn,ii) + tt*maxK );
            }
        } else {
            for (int tt=0; tt<TP; tt++){
                t1(nn,tt) = 1;
            }
        }
    }
    if (use_itemnr2){
        ii = itemnr2;
        for (int nn=0;nn<N;nn++){
            if (dat_resp_bool(nn,itemnr2)){
                for (int tt=0; tt<TP; tt++){
                    t1(nn,tt) = t1(nn,tt) * probs_ratio(ii, dat(nn,ii) + tt*maxK );
                }
            }
        }
    }

    for (int nn=0;nn<N;nn++){
        for (int tt=0; tt<TP; tt++){
            t1(nn,tt) = t1(nn,tt)*p_xi_aj(nn,tt);
        }
    }

    //---- OUTPUT
    return t1;
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_xxirt_newton_raphson_reduced_probs
// [[Rcpp::export]]
double sirt_rcpp_xxirt_newton_raphson_derivative_par(
        Rcpp::IntegerMatrix dat, Rcpp::LogicalMatrix dat_resp_bool,
        Rcpp::NumericMatrix ratio, Rcpp::NumericMatrix p_xi_aj,
        int item, Rcpp::NumericMatrix prior_Theta,
        Rcpp::IntegerVector group0, Rcpp::NumericVector weights,
        Rcpp::NumericVector ll_case0, double eps)
{
    int N = dat.nrow();
    int TP = p_xi_aj.ncol();
    double temp = 0;
    int ii = item-1;
    double ll_case_der_nn=0;
    double grad_pp=0;

    for (int nn=0; nn<N; nn++){
        if (dat_resp_bool(nn,ii)){
            ll_case_der_nn=0;
            for (int tt=0; tt<TP; tt++){
                temp=p_xi_aj(nn,tt)*ratio( dat(nn,ii), tt);
                ll_case_der_nn += prior_Theta(tt,group0[nn])*temp;
            }
            grad_pp -= weights[nn] * ll_case_der_nn / ( ll_case0[nn] + eps );
        }
    }

    //---- OUTPUT
    return grad_pp;
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_xxirt_nr_pml_opt_fun
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_xxirt_nr_pml_opt_fun(
        Rcpp::NumericVector prior_Theta, Rcpp::NumericVector probs_items,
        Rcpp::NumericVector freq1, Rcpp::NumericVector freq2,
        Rcpp::NumericVector W1, Rcpp::NumericMatrix W2_long,
        int G, int K, int I, int TP, int NI2, double eps )
{
    double val=0;
    int ind=0;
    Rcpp::NumericVector val1(I*K*G);
    Rcpp::NumericVector val2(NI2*K*K*G);

    // first-order frequencies
    for (int gg=0; gg<G; gg++){
        for (int ii=0; ii<I; ii++){
            for (int hh=0; hh<K; hh++){
                double p1=0;
                for (int tt=0; tt<TP; tt++){
                    // p1 <- sum( prior_Theta[,gg] * probs_items[ii,hh,] )
                    p1 += prior_Theta[tt+gg*TP]*probs_items[ii+hh*I+tt*I*K];
                } // end tt
                // t1 <- t1 + freq1[ii,hh,gg]*log(p1+eps)
                ind=ii+hh*I+gg*I*K;
                val1[ind]=p1;
                val += W1[ii]*freq1[ind]*std::log(p1+eps);
            } // end hh
        } // end ii
    }  // end gg

    // second-order frequencies
    double t1=0;
    double w2ii=0;
    for (int gg=0; gg<G; gg++){
        for (int nn=0; nn<NI2; nn++){
            // ii1 <- W2_long[nn,"item1"]
            //    ii2 <- W2_long[nn,"item2"]
            //    w2ii <- W2_long[nn,"w2"]
            int ii1=W2_long(nn,0)-1;
            int ii2=W2_long(nn,1)-1;
            w2ii=W2_long(nn,2);
            for (int hh=0; hh<K; hh++){
                for (int kk=0; kk<K; kk++){
                    t1=0;
                    for (int tt=0; tt<TP; tt++){
                        //            p1 <- sum( prior_Theta[,gg] * probs_items[ii1,hh,]*
                        //                                        probs_items[ii2,kk,])
                        t1+=prior_Theta[tt+gg*TP]*probs_items[ii1+hh*I+tt*I*K]*
                                                    probs_items[ii2+kk*I+tt*I*K];
                    } // end tt
                    //            t1 <- t1 + freq2[nn,hh,kk,gg]*log(p1+eps)
                    ind=nn+hh*NI2+kk*NI2*K+gg*NI2*K*K;
                    val2[ind]=t1;
                    val += w2ii*freq2[ind]*std::log(t1+eps);
                } // end kk
            } // end hh
        } // end nn
    }  // end gg


    //---- OUTPUT
    return Rcpp::List::create(
                Rcpp::Named("val")=val,
                Rcpp::Named("val1")=val1,
                Rcpp::Named("val2")=val2
            );
}
///********************************************************************



///********************************************************************
///** sirt_rcpp_xxirt_nr_pml_grad_fun_eval
// [[Rcpp::export]]
double sirt_rcpp_xxirt_nr_pml_grad_fun_eval(
        Rcpp::NumericVector prior_Theta, Rcpp::NumericVector probs_items,
        Rcpp::NumericVector freq1, Rcpp::NumericVector freq2,
        Rcpp::NumericVector W1, Rcpp::NumericMatrix W2_long,
        int G, int K, int I, int TP, int NI2, double eps, int NP,
        Rcpp::NumericVector der_prior_Theta, Rcpp::NumericVector val1,
        Rcpp::NumericVector val2, bool pp_Theta,
        Rcpp::NumericVector der_probs_items, Rcpp::IntegerVector index_freq1,
        Rcpp::IntegerVector index_freq2)
{
    double val=0;
    double pt=0;
    double pi=0;
    int ind=0;
    int ind1=0;
    int ind2=0;
    int jnd=0;
    double grad=0;
    int I1=0;
    int ii=0;
    double p1=0;

    // first-order frequencies
    for (int gg=0; gg<G; gg++){
        I1 = index_freq1.size();
        for (int iii=0; iii<I1; iii++){
            ii = index_freq1[iii];
            for (int hh=0; hh<K; hh++){
                p1=0;
                for (int tt=0; tt<TP; tt++){
                    ind1=tt+gg*TP;
                    ind2=ii+hh*I+tt*I*K;
                    if (pp_Theta){
                        pt=der_prior_Theta[ind1];
                        pi=probs_items[ind2];
                    } else {
                        pt=prior_Theta[ind1];
                        pi=der_probs_items[ind2];
                    }
                    p1 += pt*pi;
                } // end tt
                ind=ii+hh*I+gg*I*K;
                val += W1[ii]*freq1[ind]*p1/(val1[ind]+eps);
            } // end hh
        } // end ii
    }  // end gg

    // second-order frequencies
    double t1=0;
    double w2ii=0;
    double pi1=0;
    double pi2=0;
    int jnd1=0;
    int jnd2=0;
    int N2=0;
    int nn=0;

    for (int gg=0; gg<G; gg++){
        N2 = index_freq2.size();
        for (int nnn=0; nnn<N2; nnn++){
            nn = index_freq2[nnn];
            int ii1=W2_long(nn,0)-1;
            int ii2=W2_long(nn,1)-1;
            w2ii=W2_long(nn,2);
            for (int hh=0; hh<K; hh++){
                for (int kk=0; kk<K; kk++){
                    t1=0;
                    for (int tt=0; tt<TP; tt++){
                        jnd=tt+gg*TP;
                        jnd1=ii1+hh*I+tt*I*K;
                        jnd2=ii2+kk*I+tt*I*K;
                        if (pp_Theta){
                            pt=der_prior_Theta[jnd];
                            pi=probs_items[jnd1]*probs_items[jnd2];
                        } else {
                            pt=prior_Theta[jnd];
                            pi1=der_probs_items[jnd1]*probs_items[jnd2];
                            pi2=probs_items[jnd1]*der_probs_items[jnd2];
                            pi=pi1+pi2;
                        }
                        t1+=pt*pi;
                    } // end tt
                    ind=nn+hh*NI2+kk*NI2*K+gg*NI2*K*K;
                    val += w2ii*freq2[ind]*t1/(val2[ind]+eps);
                } // end kk
            } // end hh
        } // end nn
    }  // end gg

    grad=-val;

    return grad;
}
///********************************************************************


///********************************************************************
///** sirt_rcpp_xxirt_nr_pml_casewise_opt_fun
// [[Rcpp::export]]
Rcpp::List sirt_rcpp_xxirt_nr_pml_casewise_opt_fun(
        Rcpp::NumericVector prior_Theta, Rcpp::NumericVector probs_items,
        Rcpp::NumericVector W1, Rcpp::NumericMatrix W2_long,
        int G, int K, int I, int TP, int NI2, double eps,
        Rcpp::IntegerVector group0, Rcpp::NumericVector weights,
        Rcpp::IntegerMatrix dat1, Rcpp::LogicalMatrix dat_resp)
{
    double val=0;
    int NC=dat1.nrow();
    int gg=0;
    double weight_cc=1.0;
    double p1=0;
    double temp=0;
    double input=0;
    int ind=0;
    Rcpp::NumericVector case_val(NC);
    Rcpp::NumericVector model_probs1(I*G*K);

    //-- first-order probabilities
    for (int gg=0; gg<G; gg++){
        for (int ii=0; ii<I; ii++){
            for (int hh=0; hh<K; hh++){
                p1=0;
                for (int tt=0; tt<TP; tt++){
                    p1 += prior_Theta[tt+gg*TP]*probs_items[ii+hh*I+tt*I*K];
                } // end tt
                ind=ii+hh*I+gg*I*K;
                model_probs1[ind]=std::log(p1+eps);
            } // end hh
        } // end ii
    }  // end gg

    //-- opt fun first-order frequencies
    for (int cc=0; cc<NC; cc++){
        gg=group0[cc];
        weight_cc=weights[cc];
        for (int ii=0; ii<I; ii++){
            if (dat_resp(cc,ii)){
                for (int hh=0; hh<K; hh++){
                    if (dat1(cc,ii)==hh){
                        ind=ii+hh*I+gg*I*K;
                        input=model_probs1[ind];
                        temp = W1[ii]*weight_cc*input;
                        val += temp;
                        case_val[cc] += temp;
                    } // end hh
                } // end if    dat1(cc,ii)==hh
            }  // end if dat_resp(cc,ii) = TRUE
        } // end ii
    }  // end cc

    double t1=0;
    double w2ii=0;
    Rcpp::NumericVector model_probs2(NI2*K*K*G);

    //--- log probs second-order frequencies
    for (int gg=0; gg<G; gg++){
        for (int nn=0; nn<NI2; nn++){
            int ii1=W2_long(nn,0)-1;
            int ii2=W2_long(nn,1)-1;
            w2ii=W2_long(nn,2);
            for (int hh=0; hh<K; hh++){
                for (int kk=0; kk<K; kk++){
                    t1=0;
                    for (int tt=0; tt<TP; tt++){
                        t1+=prior_Theta[tt+gg*TP]*probs_items[ii1+hh*I+tt*I*K]*
                                                    probs_items[ii2+kk*I+tt*I*K];
                    } // end tt
                    ind=nn+hh*NI2+kk*NI2*K+gg*NI2*K*K;
                    model_probs2[ind]=std::log(t1+eps);
                } // end kk
            } // end hh
        } // end nn
    }  // end gg

    //-- opt fun second-order frequencies
    for (int cc=0; cc<NC; cc++){
        weight_cc=weights[cc];
        gg=group0[cc];
        for (int nn=0; nn<NI2; nn++){
            int ii1=W2_long(nn,0)-1;
            int ii2=W2_long(nn,1)-1;
            w2ii=W2_long(nn,2);
            if (dat_resp(cc,ii1)&dat_resp(cc,ii2)){
                for (int hh=0; hh<K; hh++){
                    for (int kk=0; kk<K; kk++){
                        if ( (dat1(cc,ii1)==hh) & (dat1(cc,ii2)==kk) ){
                            ind=nn+hh*NI2+kk*NI2*K+gg*NI2*K*K;
                            input=model_probs2[ind];
                            // input=std::log(t1+eps);
                            temp = w2ii*weight_cc*input;
                            val += temp;
                            case_val[cc] += temp;
                        }  // end if
                    } // end kk
                } // end hh
            } // if
        } // end nn
    }  // end cc


    //---- OUTPUT
    return Rcpp::List::create(
                Rcpp::Named("val")=val,
                Rcpp::Named("case_val")=case_val
            );
}
///********************************************************************



// if ( ! R_IsNA( resp(nn,ii) ) ){
