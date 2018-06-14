//// File Name: noharm_sirt_rcpp.cpp
//// File Version: 3.32

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

////************************************
///// noharm auxiliary functions

// trim increment
double noharm_trim_increment( double x, double maxincrement)
{
    double increment = x;
    if (increment < -maxincrement){
        increment = - maxincrement;
    }
    if (increment > maxincrement){
        increment = maxincrement;
    }
    return increment;
}

// absolute value function
double noharm_abs(double x)
{
    double x0=0;
    if (x < 0){
        x0 = - x;
    }
    return x0;
}


// avoid division by zero
double noharm_avoid_zero(double x, double eps)
{
    double x0=x;
    if ( ( x < 0) & ( x > - eps ) ) {
        x0 = - eps;
    }
    if ( ( x > 0) & ( x < eps ) ) {
        x0 = eps;
    }
    return x0;
}

//---- calculate dj
///** noharm_compute_dj
// [[Rcpp::export]]
Rcpp::NumericVector noharm_compute_dj(Rcpp::NumericMatrix Fval, Rcpp::NumericMatrix Pval,
    int I, int D)
{
    Rcpp::NumericVector dj(I);
    //*************** calculate dj
    // dj <- sqrt( diag( Fval_ %*% Pval %*% t(Fval_) ) )
    for ( int jj = 0; jj <I;jj++){
        for( int dd1 = 0; dd1<D;dd1++){
            for( int dd2 = 0; dd2<D;dd2++){
                dj[jj] += Fval(jj,dd1) * Pval(dd1,dd2) * Fval(jj,dd2) ;
            }
        }
        dj[jj] = sqrt( dj[jj] ) ;
    }
    return dj;
}


//---- calculate ej
///** noharm_compute_ej
// [[Rcpp::export]]
Rcpp::List noharm_compute_ej(Rcpp::NumericVector dj, int I)
{
    Rcpp::NumericVector ej(I) ;
    Rcpp::NumericMatrix ej_ek(I,I) ;
    for( int jj = 0 ; jj<I;jj++){
        ej[jj] = sqrt( 1 + pow( dj[jj] , 2.0) ) ;
    }
    for ( int jj = 0 ; jj<I;jj++){
        for (int kk=0;kk<I;kk++){
            if ( jj != kk){
                ej_ek(jj,kk) = 1/( ej[jj] * ej[kk] ) ;
            }
        }
    }
    //--- output
    return Rcpp::List::create(
            Rcpp::Named("ej") = ej ,
            Rcpp::Named("ej_ek") = ej_ek
        ) ;
}



//---- calculate vjk
///** noharm_compute_vjk
// [[Rcpp::export]]
Rcpp::List noharm_compute_vjk(Rcpp::NumericMatrix b0jk, Rcpp::NumericMatrix b1jk,
        Rcpp::NumericMatrix b2jk, Rcpp::NumericMatrix b3jk, Rcpp::NumericMatrix ej_ek,
        int I, int modtype )
{
    Rcpp::NumericMatrix v0jk(I,I);
    Rcpp::NumericMatrix v1jk(I,I);
    Rcpp::NumericMatrix v2jk(I,I);
    Rcpp::NumericMatrix v3jk(I,I);
    for (int ii=0; ii<I; ii++){
        v0jk(_,ii) = b0jk(_,ii);
    }

    for( int jj=0; jj < I ; jj ++ ){
        for( int kk=0; kk < I ; kk ++ ){
            if (modtype != 2){
                v1jk(jj,kk) = b1jk(jj,kk) * ej_ek(jj,kk) ;
                v2jk(jj,kk) = b2jk(jj,kk) * pow( ej_ek(jj,kk) ,2.0);
                v3jk(jj,kk) = b3jk(jj,kk) * pow( ej_ek(jj,kk) ,3.0);
            } else {
                v0jk(jj,kk)=0;
                v1jk(jj,kk)=1;
                v2jk(jj,kk)=0;
                v3jk(jj,kk)=0;
            }
        }
    }

    //--- output
    return Rcpp::List::create(
            Rcpp::Named("v0jk") = v0jk ,
            Rcpp::Named("v1jk") = v1jk ,
            Rcpp::Named("v2jk") = v2jk ,
            Rcpp::Named("v3jk") = v3jk
            ) ;
}


//---- calculate gammajk
///** noharm_compute_gammajk
// [[Rcpp::export]]
Rcpp::List noharm_compute_gammajk(Rcpp::NumericMatrix Fval, Rcpp::NumericMatrix Pval,
        Rcpp::NumericMatrix Psival, int I, int D )
{

    Rcpp::NumericMatrix gammajk(I,I);
    Rcpp::NumericMatrix gammajk2(I,I);
    Rcpp::NumericMatrix gammajk3(I,I);
    for (int jj=0;jj<I;jj++){
        for (int kk=0;kk<I;kk++){
            for( int dd1 = 0 ; dd1<D;dd1++){
                for( int dd2 = 0 ; dd2<D;dd2++){
                    gammajk(jj,kk) += Fval(jj,dd1) * Pval(dd1,dd2) * Fval(kk,dd2) + Psival(jj,kk) ;
                }
            }
            gammajk2(jj,kk) = gammajk(jj,kk)*gammajk(jj,kk) ;
            gammajk3(jj,kk) = gammajk2(jj,kk)*gammajk(jj,kk) ;
        }
    }
    //--- output
    return Rcpp::List::create(
            Rcpp::Named("gammajk") = gammajk ,
            Rcpp::Named("gammajk2") = gammajk2 ,
            Rcpp::Named("gammajk3") = gammajk3
            ) ;
}


//---- calculate pdfk
///** noharm_compute_pdfk
// [[Rcpp::export]]
Rcpp::List noharm_compute_pdfk(Rcpp::NumericMatrix Fval, Rcpp::NumericMatrix Pval,
        int I, int D )
{
    Rcpp::NumericMatrix pdfk(I,D);
    Rcpp::NumericMatrix pdfk2(I,D);
    for (int jj=0 ;jj<I;jj++){
        for (int dd=0 ;dd<D;dd++){
            for (int ee=0;ee<D;ee++){
                pdfk(jj,dd) += Fval(jj,ee)*Pval(ee,dd) ;
            }
            pdfk2(jj,dd) = pdfk(jj,dd)*pdfk(jj,dd) ;
        }
    }
    //--- output
    return Rcpp::List::create(
            Rcpp::Named("pdfk") = pdfk ,
            Rcpp::Named("pdfk2") = pdfk2
            ) ;
}



/////////////////***********************************////////////////
/// estimation subfunction for F

///********************************************************************
///** noharm_estFcpp
// [[Rcpp::export]]
Rcpp::List noharm_estFcpp( Rcpp::NumericMatrix Fval, Rcpp::NumericMatrix Pval, Rcpp::NumericMatrix Fpatt,
    Rcpp::NumericMatrix Ppatt, int I, int D, Rcpp::NumericMatrix b0jk, Rcpp::NumericMatrix b1jk,
    Rcpp::NumericMatrix b2jk, Rcpp::NumericMatrix b3jk, Rcpp::NumericMatrix wgtm, Rcpp::NumericMatrix pm,
    Rcpp::NumericMatrix Psival, Rcpp::NumericMatrix Psipatt, double maxincrement , int modtype )
{

    int FR = Fval.nrow();
    int FC = Fval.ncol();
    Rcpp::NumericMatrix Fval_(FR,FC);
    for (int cc=0; cc<FC; cc++){
        Fval_(_,cc) = Fval(_,cc);
    }

    //*************** calculate dj
    // dj <- sqrt( diag( Fval_ %*% Pval %*% t(Fval_) ) )
    Rcpp::NumericVector dj = noharm_compute_dj( Fval_, Pval , I, D );

    //*************** calculate ej
    //     ej <- sqrt( 1 + dj^2 )
    //     ej.ek <- 1 / outer( ej , ej )
    //     diag( ej.ek ) <- 0
    Rcpp::List res1 = noharm_compute_ej( dj, I);
    Rcpp::NumericVector ej = res1["ej"];
    Rcpp::NumericMatrix ej_ek = res1["ej_ek"];

    //*************** calculate vjk
    Rcpp::List res2 = noharm_compute_vjk(b0jk, b1jk, b2jk, b3jk, ej_ek, I, modtype );
    Rcpp::NumericMatrix v0jk = res2["v0jk"];
    Rcpp::NumericMatrix v1jk = res2["v1jk"];
    Rcpp::NumericMatrix v2jk = res2["v2jk"];
    Rcpp::NumericMatrix v3jk = res2["v3jk"];

    //***********  compute gamma.jk = f_j' P f_k
    //    gamma.jk <- Fval_ %*% Pval %*% t(Fval_ )
    Rcpp::List res3 = noharm_compute_gammajk(Fval, Pval, Psival, I, D );
    Rcpp::NumericMatrix gammajk = res3["gammajk"];
    Rcpp::NumericMatrix gammajk2 = res3["gammajk2"];
    Rcpp::NumericMatrix gammajk3 = res3["gammajk3"];

    //******* compute  pd.fk <- Fval_ %*% Pval
    Rcpp::List res4 = noharm_compute_pdfk( Fval, Pval, I, D );
    Rcpp::NumericMatrix pdfk = res4["pdfk"];
    Rcpp::NumericMatrix pdfk2 = res4["pdfk2"];

    //*!!!!!!!!!!!!!!!!!!!!!!!
    //*!!!!! one iteration
    Rcpp::NumericVector eps0jj(I);
    Rcpp::NumericVector eps1jj(I);
    Rcpp::NumericVector eps2jj(I);
    double f1jj=0;
    double f2jj=0;
    double increment=0;
    double parchange=0;
    double aincr=0;
    double eps1 = 1e-7;

    for (int jj=0;jj<I;jj++){
        for (int dd=0 ;dd<D;dd++){
            if ( Fpatt(jj,dd) > 0 ){
                f1jj=0;
                f2jj=0;
                for (int kk=0; kk<I;kk++){
                    // eps0.jj <- ( wgtm[jj,] * ( pm[jj,] - v0.jk[jj,] -
                    //  v1.jk[jj,]*gamma.jk[jj,] - v2.jk[jj,]*gamma.jk[jj,]^2 -
                    //        v3.jk[jj,]*gamma.jk[jj,]^3 ) )
                    eps0jj[kk] =  pm(jj,kk) - v0jk(jj,kk) -
                            v1jk(jj,kk) * gammajk(jj,kk) - v2jk(jj,kk) * gammajk2(jj,kk) -
                            v3jk(jj,kk) * gammajk3(jj,kk) ;
                    eps0jj[kk] = wgtm(jj,kk)* eps0jj[kk] ;
                    //  eps1.jj <- - ( v1.jk[jj,]  * pd.fk[,dd] ) -
                    //    2 * ( v2.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd] ) -
                    //    3 * ( v3.jk[jj,]  * gamma.jk[jj,]^2  * pd.fk[ , dd] )
                    eps1jj[kk] = - v1jk(jj,kk) * pdfk(kk,dd) -
                            2*v2jk(jj,kk) * gammajk(jj,kk)*pdfk(kk,dd) -
                            3*v3jk(jj,kk) * gammajk2(jj,kk)*pdfk(kk,dd) ;
                    eps1jj[kk] = wgtm(jj,kk)* eps1jj[kk] ;
                    //  eps2.jj <- - 2 * ( v2.jk[jj,]  * pd.fk[ , dd]^2 ) -
                    //    6 * ( v3.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd]^2 )
                    eps2jj[kk] = - 2*v2jk(jj,kk) * pdfk2(kk,dd) -
                            6*v3jk(jj,kk) * gammajk(jj,kk) * pdfk2(kk,dd) ;
                    eps2jj[kk] = wgtm(jj,kk)* eps2jj[kk] ;
                    // f1.jj <- 2* eps0.jj * eps1.jj
                    f1jj += 2 * eps0jj[kk] * eps1jj[kk] ;
                    // f2.jj <- 2* eps1.jj^2 + 2*eps0.jj * eps2.jj
                    f2jj += 2 * eps1jj[kk] * eps1jj[kk] + 2*eps0jj[kk]* eps2jj[kk] ;
                }
                increment = - f1jj / noharm_avoid_zero( f2jj, eps1) ;
                // trim increment
                increment = noharm_trim_increment( increment, maxincrement);
                // absolute value
                aincr = noharm_abs( increment );
                if ( aincr > parchange ){ parchange = aincr ; }
                Fval_(jj,dd) = Fval_(jj,dd) + increment ;
            }
        }
    }
    //*************************************************
    // OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("Fval_") = Fval_ ,
            Rcpp::Named("change") = parchange
        ) ;
}


/////////////////***********************************////////////////
/// estimation subfunction for P

///********************************************************************
///** noharm_estPcpp
// [[Rcpp::export]]
Rcpp::List noharm_estPcpp( Rcpp::NumericMatrix Fval, Rcpp::NumericMatrix Pval, Rcpp::NumericMatrix Fpatt,
    Rcpp::NumericMatrix Ppatt, int I, int D, Rcpp::NumericMatrix b0jk, Rcpp::NumericMatrix b1jk,
    Rcpp::NumericMatrix b2jk, Rcpp::NumericMatrix b3jk, Rcpp::NumericMatrix wgtm, Rcpp::NumericMatrix pm,
    Rcpp::NumericMatrix Psival, Rcpp::NumericMatrix Psipatt, double maxincrement , int modtype )
{
    int FR = Pval.nrow();
    int FC = Pval.ncol();
    Rcpp::NumericMatrix Pval_(FR,FC);
    for (int cc=0; cc<FC; cc++){
        Pval_(_,cc) = Pval(_,cc);
    }

     //*************** calculate dj
     //     dj <- sqrt( diag( Fval %*% Pval_ %*% t(Fval) ) )
    Rcpp::NumericVector dj = noharm_compute_dj( Fval, Pval , I, D );

     //*************** calculate ej
     //     ej <- sqrt( 1 + dj^2 )
     //     ej.ek <- 1 / outer( ej , ej )
     //     diag( ej.ek ) <- 0
    Rcpp::List res1 = noharm_compute_ej( dj, I);
    Rcpp::NumericVector ej = res1["ej"];
    Rcpp::NumericMatrix ej_ek = res1["ej_ek"];

     //*************** calculate vjk
    Rcpp::List res2 = noharm_compute_vjk(b0jk, b1jk, b2jk, b3jk, ej_ek, I, modtype );
    Rcpp::NumericMatrix v0jk = res2["v0jk"];
    Rcpp::NumericMatrix v1jk = res2["v1jk"];
    Rcpp::NumericMatrix v2jk = res2["v2jk"];
    Rcpp::NumericMatrix v3jk = res2["v3jk"];

     //***********  compute gamma.jk = f_j' P f_k
     //    gamma.jk <- Fval %*% Pval_ %*% t(Fval )
    Rcpp::List res3 = noharm_compute_gammajk(Fval, Pval, Psival, I, D );
    Rcpp::NumericMatrix gammajk = res3["gammajk"];
    Rcpp::NumericMatrix gammajk2 = res3["gammajk2"];
    Rcpp::NumericMatrix gammajk3 = res3["gammajk3"];

     //******* compute  pd.fk <- Fval %*% Pval_
    Rcpp::List res4 = noharm_compute_pdfk( Fval, Pval, I, D );
    Rcpp::NumericMatrix pdfk = res4["pdfk"];
    Rcpp::NumericMatrix pdfk2 = res4["pdfk2"];

     //*!!!!!!!!!!!!!!!!!!!!!!!
     //*!!!!! one iteration
     Rcpp::NumericMatrix eps0jj(I,I);
     Rcpp::NumericMatrix eps1jj(I,I);
     Rcpp::NumericMatrix eps2jj(I,I);
     double f1jj=0;
     double f2jj=0;
     double increment=0;
     double parchange=0;
     double gammajk1 =0;
     double gammajk12=0;
     double aincr=0;
     double eps1 = 1e-7;

     for ( int dd=0;dd<D;dd++){
     for (int ee=0 ;ee<D;ee++){

     if (dd >= ee ){
      if ( Ppatt(dd,ee) > 0 ){

      f1jj=0;
      f2jj=0;
     for (int jj=0;jj<I;jj++){
     for (int kk=0;kk<I;kk++){

     // gammajk1 <- outer( Fval[ ,dd] , Fval[ ,ee] )
     gammajk1 = Fval(jj,dd) * Fval( kk,ee) ;
     if (dd==ee){ gammajk1 = 2*gammajk1 ; }
     gammajk12 = gammajk1 * gammajk1 ;
     // eps0.jj <- ( wgtm * ( pm - v0.jk - v1.jk*gamma.jk -
     //     v2.jk*gamma.jk^2 - v3.jk*gamma.jk^3 ) )
         eps0jj(jj,kk) =  pm(jj,kk) - v0jk(jj,kk) -
             v1jk(jj,kk) * gammajk(jj,kk) - v2jk(jj,kk) * gammajk2(jj,kk) -
             v3jk(jj,kk) * gammajk3(jj,kk) ;
         eps0jj(jj,kk) = wgtm(jj,kk)* eps0jj(jj,kk) ;
     //    eps1.jj <- - ( v1.jk  * gammajk1 ) -
     //        2 * ( v2.jk  * gamma.jk  * gammajk1 ) -
     //            3 * ( v3.jk  * gamma.jk^2  * gammajk1 )
         eps1jj(jj,kk) = - v1jk(jj,kk) * gammajk1 -
             2*v2jk(jj,kk) * gammajk(jj,kk)*gammajk1 -
             3*v3jk(jj,kk) * gammajk2(jj,kk)*gammajk1 ;
         eps1jj(jj,kk) = wgtm(jj,kk)* eps1jj(jj,kk) ;
     //    eps2.jj <- - 2 * ( v2.jk  * gammajk1^2 ) -
     //            6 * ( v3.jk  * gamma.jk  * gammajk1^2 )
         eps2jj(jj,kk) = - 2*v2jk(jj,kk) * gammajk12 -
             6*v3jk(jj,kk) * gammajk(jj,kk) * gammajk12 ;
         eps2jj(jj,kk) = wgtm(jj,kk)* eps2jj(jj,kk) ;
         // f1.jj <- 2* eps0.jj * eps1.jj
         f1jj += 2 * eps0jj(jj,kk) * eps1jj(jj,kk) ;
         // f2.jj <- 2* eps1.jj^2 + 2*eps0.jj * eps2.jj
         f2jj += 2 * eps1jj(jj,kk) * eps1jj(jj,kk) +
             2*eps0jj(jj,kk)* eps2jj(jj,kk) ;
     }
     }

                increment = - f1jj / noharm_avoid_zero( f2jj, eps1) ;
                // trim increment
                increment = noharm_trim_increment( increment, maxincrement);
                // absolute value
                aincr = noharm_abs( increment );
     if ( aincr > parchange ){ parchange = aincr ; }

     //   Pval_[dd,ee] <- Pval_[ dd,ee] + increment
     //   if ( dd > ee ){ Pval_[ee,dd] <- Pval_[dd,ee] }
     Pval_(dd,ee) = Pval_(dd,ee) + increment ;
     if (dd>ee){ Pval_(ee,dd)=Pval_(dd,ee) ; }

      }
     }
     }
     }

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["Pval_"] = Pval_ ,
         Rcpp::_["change"] = parchange ,
         Rcpp::_["residuals"] = eps0jj
         ) ;
}




/////////////////***********************************////////////////
/// estimation subfunction for Psi

///********************************************************************
///** noharm_estPcpp
// [[Rcpp::export]]
Rcpp::List noharm_estPsicpp( Rcpp::NumericMatrix Fval, Rcpp::NumericMatrix Pval, Rcpp::NumericMatrix Fpatt,
    Rcpp::NumericMatrix Ppatt, int I, int D, Rcpp::NumericMatrix b0jk, Rcpp::NumericMatrix b1jk,
    Rcpp::NumericMatrix b2jk, Rcpp::NumericMatrix b3jk, Rcpp::NumericMatrix wgtm, Rcpp::NumericMatrix pm,
    Rcpp::NumericMatrix Psival, Rcpp::NumericMatrix Psipatt, double maxincrement , int modtype )
{

    int FR = Psival.nrow();
    int FC = Psival.ncol();
    Rcpp::NumericMatrix Psival_(FR,FC);
    for (int cc=0; cc<FC; cc++){
        Psival_(_,cc) = Psival(_,cc);
    }

     //*************** calculate dj
     //     dj <- sqrt( diag( Fval %*% Pval %*% t(Fval) ) )
    Rcpp::NumericVector dj = noharm_compute_dj( Fval, Pval , I, D );

     //*************** calculate ej
     //     ej <- sqrt( 1 + dj^2 )
     //     ej.ek <- 1 / outer( ej , ej )
     //     diag( ej.ek ) <- 0
    Rcpp::List res1 = noharm_compute_ej( dj, I);
    Rcpp::NumericVector ej = res1["ej"];
    Rcpp::NumericMatrix ej_ek = res1["ej_ek"];

     //*************** calculate vjk
    Rcpp::List res2 = noharm_compute_vjk(b0jk, b1jk, b2jk, b3jk, ej_ek, I, modtype );
    Rcpp::NumericMatrix v0jk = res2["v0jk"];
    Rcpp::NumericMatrix v1jk = res2["v1jk"];
    Rcpp::NumericMatrix v2jk = res2["v2jk"];
    Rcpp::NumericMatrix v3jk = res2["v3jk"];

     //***********  compute gamma.jk = f_j' P f_k
     //    gamma.jk <- Fval %*% Pval %*% t(Fval )
    Rcpp::List res3 = noharm_compute_gammajk(Fval, Pval, Psival, I, D );
    Rcpp::NumericMatrix gammajk = res3["gammajk"];
    Rcpp::NumericMatrix gammajk2 = res3["gammajk2"];
    Rcpp::NumericMatrix gammajk3 = res3["gammajk3"];

     //******* compute  pd.fk <- Fval %*% Pval
    Rcpp::List res4 = noharm_compute_pdfk( Fval, Pval, I, D );
    Rcpp::NumericMatrix pdfk = res4["pdfk"];
    Rcpp::NumericMatrix pdfk2 = res4["pdfk2"];

     //*!!!!!!!!!!!!!!!!!!!!!!!
     //*!!!!! one iteration
     Rcpp::NumericMatrix eps0jj(I,I);
     Rcpp::NumericMatrix eps1jj(I,I);
     Rcpp::NumericMatrix eps2jj(I,I);
     double f1jj=0;
     double f2jj=0;
     double increment=0;
     double parchange=0;
     double aincr=0;
     double eps1 = 1e-7;

     for (int jj=0;jj<I;jj++){
     for (int kk=0 ;kk<I;kk++){

     if ( Psipatt(jj,kk) > 0 ){
     if (jj>kk){

         // eps0.jj <- ( wgtm[jj,] * ( pm[jj,] - v0.jk[jj,] -
         //  v1.jk[jj,]*gamma.jk[jj,] - v2.jk[jj,]*gamma.jk[jj,]^2 -
         //        v3.jk[jj,]*gamma.jk[jj,]^3 ) )
         eps0jj(jj,kk) =  pm(jj,kk) - v0jk(jj,kk) -
             v1jk(jj,kk) * gammajk(jj,kk) - v2jk(jj,kk) * gammajk2(jj,kk) -
             v3jk(jj,kk) * gammajk3(jj,kk)  ;
         eps0jj(jj,kk) = wgtm(jj,kk)* eps0jj(jj,kk) ;
         //  eps1.jj <- - ( v1.jk[jj,]  * pd.fk[,dd] ) -
         //    2 * ( v2.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd] ) -
         //    3 * ( v3.jk[jj,]  * gamma.jk[jj,]^2  * pd.fk[ , dd] )
         eps1jj(jj,kk) = - v1jk(jj,kk)  -
             2*v2jk(jj,kk) * gammajk(jj,kk) -
             3*v3jk(jj,kk) * gammajk2(jj,kk) ;
         eps1jj(jj,kk) = wgtm(jj,kk)* eps1jj(jj,kk) ;
         //  eps2.jj <- - 2 * ( v2.jk[jj,]  * pd.fk[ , dd]^2 ) -
         //    6 * ( v3.jk[jj,]  * gamma.jk[jj,]  * pd.fk[ , dd]^2 )
         eps2jj(jj,kk) = 0 ;
         // f1.jj <- 2* eps0.jj * eps1.jj
         f1jj = 2 * eps0jj(jj,kk) * eps1jj(jj,kk) ;
         // f2.jj <- 2* eps1.jj^2 + 2*eps0.jj * eps2.jj
         f2jj = 2 * eps1jj(jj,kk) * eps1jj(jj,kk) + 2*eps0jj(jj,kk)* eps2jj(jj,kk) ;
        increment = - f1jj / noharm_avoid_zero( f2jj, eps1) ;
        // trim increment
        increment = noharm_trim_increment( increment, maxincrement);
        // absolute value
        aincr = noharm_abs( increment );
        if ( aincr > parchange ){ parchange = aincr ; }

        Psival_(jj,kk) = Psival_(jj,kk) + increment ;
        Psival_(kk,jj) = Psival_(jj,kk) ;
     }
     }
     }
     }

     //*************************************************
     // OUTPUT
     return Rcpp::List::create(
         Rcpp::_["Psival_"] = Psival_ ,
         Rcpp::_["change"] = parchange
         ) ;
}


