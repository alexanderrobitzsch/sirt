//// File Name: gooijer_isop_rcpp.cpp
//// File Version: 5.33


// [[Rcpp::depends(RcppArmadillo)]]


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;



//**********************************************************************
// compute means and covariances of estimators obtained by Jackknife
Rcpp::NumericVector gta( Rcpp::NumericMatrix dat )
{
    int N=dat.nrow();
    int I=dat.ncol(); // number of items
    int VV=I*(I-1)/ 2;
    int K=I+1;

    Rcpp::NumericVector score(N);
    Rcpp::NumericVector nk(I+1);
    Rcpp::NumericMatrix xik(I+1,I);
    Rcpp::NumericMatrix rijk(I+1,VV);
    Rcpp::NumericVector rijkmax(I+1);
    Rcpp::NumericVector stat(1);

    //****************
    // compute score distribution
    for( int nn=0;nn<N;nn++){
        for (int ii=0;ii<I;ii++){
            score[nn] += dat(nn,ii);
        }
        nk[ score[nn] ] ++;
    }
    nk = nk + 1E-20;

    //********
    // compute scorewise item means and counts
    for (int ii=0;ii<I;ii++){
        for (int nn=0;nn<N;nn++){
            xik( score[nn], ii ) += dat(nn,ii);
        }  // end nn
        for (int kk=0;kk<I+1;kk++){
            xik(kk,ii ) = xik(kk,ii)/nk[kk];
        } // end kk
    }   // end ii

    //****************
    // compute conditional covariance
    int zz=0;
    for (int ii =0;ii<I+1;ii++){
        rijkmax[ii] = rijkmax[ii] - 100;
    }

    for (int ii=0;ii<I;ii++){
        for (int jj=ii+1;jj<I;jj++){
            for (int nn=0;nn<N;nn++){
                rijk( score[nn], zz ) += dat(nn,ii)*dat(nn,jj);
            } // end nn
            for (int kk=0;kk<K;kk++){
                rijk( kk, zz ) = rijk( kk, zz ) - nk[kk]*xik(kk,ii)*xik(kk,jj);
                rijk( kk, zz ) = rijk( kk, zz ) / nk[kk];
                if ( rijk(kk,zz) > rijkmax[kk] ){
                    rijkmax[kk] = rijk(kk,zz);
                }
            }  // end kk
        zz ++;
        } // end jj
    }  // end ii

    for (int kk=0;kk<K;kk++){
        stat[0] += nk[kk] / N * rijkmax[kk];
    }
    return stat;
}

//***************************************************
// permutation of elements in dataset
Rcpp::NumericMatrix gooijer_permutation(Rcpp::NumericMatrix sampleM, int NS, int N,
    Rcpp::NumericMatrix score_index )
{
    Rcpp::NumericVector i1(2);
    int I = sampleM.ncol();
    int NS2 = sampleM.nrow();

    Rcpp::NumericMatrix sampleM_(NS2,I);
    for (int ii=0;ii<I; ii++){
        sampleM_(_,ii) = sampleM(_,ii);
    }
    int SS=score_index.nrow();
    int t1=0;
    int NS1=0;
    int n0=0;
    int n1=0;
    for (int ss=0; ss<SS; ss++){
        NS1 = score_index(ss,3);
        n0 = score_index(ss,0) - 1;
        n1 = score_index(ss,1) - 1;
        for (int ii =0; ii<I;ii++){ // variable ii
            for ( int zz=0;zz<NS1;zz++){
                i1 = Rcpp::floor( runif(2)*(n1-n0)+n0 );
                t1 = sampleM_( i1[0],ii );
                sampleM_( i1[0], ii)  = sampleM_( i1[1],ii);
                sampleM_( i1[1], ii) = t1;
            }  // end zz
        }  // end ii
    } // end ss

    return sampleM_;
}



//********************************************************************
//********************************************************************
//********************************************************************
// unidimensionality test of Gooijer
//********************************************************************
//********************************************************************
//********************************************************************


///********************************************************************
///** gooijer_csn_table
// [[Rcpp::export]]
Rcpp::List gooijer_csn_table( Rcpp::NumericMatrix dat, Rcpp::NumericMatrix dat_perm, int RR, int NS,
    int progress, Rcpp::NumericVector progress_vec, Rcpp::NumericMatrix score_index )
{
    int N=dat.nrow();
    Rcpp::NumericVector stat_perm(RR);
    int RP = dat_perm.nrow();
    int CP = dat_perm.ncol();
    Rcpp::NumericMatrix sampleM(RP, CP);
    for (int cc=0; cc<CP; cc++){
        sampleM(_, cc) = dat_perm(_,cc);
    }
    // compute Gooijer statistic for original data
    Rcpp::NumericVector stat = gta(dat);
    Rcpp::NumericVector s1(1);
    int zz=0;
    if ( progress==1){
        Rcpp::Rcout << "|" <<   std::flush;
    }
    for (int rr=0;rr<RR;rr++){
        // permutation of original data;
        sampleM = gooijer_permutation(sampleM, NS, N, score_index );
        // compute statistic for permuted data;
        s1 = gta(sampleM );
        stat_perm[rr] = s1[0];
        if ( (progress==1) & ( rr==progress_vec[zz] ) ){
            zz = zz+1;
            if ( zz==10){ zz = 9; }
            Rcpp::Rcout << "-" <<   std::flush;
        }
    }
    if ( progress==1){
        Rcpp::Rcout << "|" <<   std::flush << std::endl;
    }
    //*************************************************
    // OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("stat") = stat,
            Rcpp::Named("stat_perm") = stat_perm
        );
}
//********************************************************************
//********************************************************************
//********************************************************************




//********************************************************************
//********************************************************************
//********************************************************************
// ISOP test
//********************************************************************
//********************************************************************
//********************************************************************



///********************************************************************
///** isop_tests_C
// [[Rcpp::export]]
Rcpp::List isop_tests_C( Rcpp::NumericMatrix dat,
    Rcpp::NumericMatrix dat_resp, Rcpp::NumericVector weights,
    Rcpp::NumericVector jackunits, int JJ )
{

    int N=dat.nrow();
    int I=dat.ncol();
    Rcpp::NumericMatrix Esi(I,JJ+1);
    Rcpp::NumericMatrix Edi(I,JJ+1);
    Rcpp::NumericMatrix W1i(I,JJ+1);
    Rcpp::NumericVector W1test(JJ+1);

    for (int ii=0;ii<I;ii++){
        for (int nn=0;nn<N;nn++){
            for (int mm=0;mm<N;mm++){
                if (mm!=nn){
                    for (int jj=0;jj<I;jj++){
                        if ( jj != ii ){
                            if ( ( dat(nn,ii) > dat(mm,ii) ) & ( dat(nn,jj) > dat(mm,jj) ) ){
                                Esi.row(ii) = Esi.row(ii) + weights[nn];
                                Esi(ii, jackunits[nn]+1 ) +=  - weights[nn];
                                if ( jackunits[nn] != jackunits[mm] ){
                                    Esi( ii, jackunits[mm] + 1 ) += - weights[nn];
                                }
                            }
                            if ( ( dat(nn,ii) < dat(mm,ii) ) & ( dat(nn,jj) < dat(mm,jj) ) ){
                                Esi.row(ii) = Esi.row(ii) + weights[nn];
                                Esi(ii, jackunits[nn]+1 ) +=  - weights[nn];
                                if ( jackunits[nn] != jackunits[mm] ){
                                    Esi( ii, jackunits[mm] + 1 ) += - weights[nn];
                                }
                            }
                            if ( ( dat(nn,ii) < dat(mm,ii) ) & ( dat(nn,jj) > dat(mm,jj) ) ){
                                Edi.row(ii) = Edi.row(ii) + weights[nn];
                                Edi(ii, jackunits[nn]+1 ) +=  - weights[nn];
                                if ( jackunits[nn] != jackunits[mm] ){
                                    Edi( ii, jackunits[mm] + 1 ) += - weights[nn];
                                }
                            }
                            if ( ( dat(nn,ii) > dat(mm,ii) ) & ( dat(nn,jj) < dat(mm,jj) ) ){
                                Edi.row(ii) = Edi.row(ii) + weights[nn];
                                Edi(ii, jackunits[nn]+1 ) +=  - weights[nn];
                                if ( jackunits[nn] != jackunits[mm] ){
                                    Edi( ii, jackunits[mm] + 1 ) += - weights[nn];
                                }
                            }
                        }  // end if jj != ii
                    }  // end for jj
                }  // end if mm!= nn
            }   // end for mm
        } // end for nn
    } // end for ii

    // compute W1i
    for (int ii=0;ii<I;ii++){
        W1i.row(ii) = ( Esi.row(ii) - Edi.row(ii)  ) /  ( Esi.row(ii) + Edi.row(ii)  );
    }

    // compute statistic for the whole test
    double tmp1, tmp2;
    for (int jj=0;jj<JJ+1;jj++){
        tmp2 = 0;
        for (int ii=0;ii<I;ii++){
            tmp1 = ( Esi(ii,jj) + Edi(ii,jj) );
            tmp2 += tmp1;
            W1test[jj] += tmp1 * W1i(ii,jj);
        }
        W1test[jj] = W1test[jj] / tmp2;
    }

    ///////////////////////////////////////
    /// OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("W1test") = W1test,
            Rcpp::Named("W1i") = W1i,
            Rcpp::Named("Esi") = Esi,
            Rcpp::Named("Edi") = Edi
        );
}
//********************************************************************
//********************************************************************
//********************************************************************


