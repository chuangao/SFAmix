#include <stdlib.h>

#include <math.h>
#include <Eigen/Dense>
#include <limits>
#include <stdio.h>
#include <ostream>
#include <string>
#include <iostream>
#include <float.h>
//#include <boost/math/special_functions/gamma.hpp>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

const double PI  =3.141592653589793238462;

void range_colwise(MatrixXd& range,MatrixXd& x,int n,int m){
	int i,j;
	for(i=0;i<m;i++){
		//cout << "i" << i << endl;
		int s;
		double min=1E300,max=0;
		// for(s=0;s<n;s++){
		// 	if(x(s,i)!=0){
		// 		min=abs(x(s,i));
		// 		max=abs(x(s,i));
		// 		continue;		
		// 	}
		// 	//cout << "s" << s << endl;
		// }
		for(j=s;j<n;j++){
			if(x(j,i)!=0){
				if(min>abs(x(j,i))){
					min=abs(x(j,i));
				}
				if(max<abs(x(j,i))){
					max=abs(x(j,i));
				}
			}
		}
		range(0,i)=min;
		range(1,i)=max;
	}
	//cout << "range_coleise_finished" << endl;
}

void range_rowwise(MatrixXd& range,MatrixXd& x,int n,int m){
	int i,j;
	for(j=0;j<n;j++){
		int s;
		double min=1E300,max=0;
		// for(s=0;s<m;s++){
		// 	if(x(j,s)!=0){
		// 		min=abs(x(j,s));
		// 		max=abs(x(j,s));
		// 		continue;
		// 	}
			
		// }
		for(i=s;i<m;i++){
			if(x(j,i)!=0){
				if(min>abs(x(j,i))){
					min=abs(x(j,i));
				}
				if(max<abs(x(j,i))){
					max=abs(x(j,i));
				}
			}
		}
		range(0,j)=min;
		range(1,j)=max;
	}
}

void cpy_row_matrix(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int row){
    for(int i=0;i<n;i++){
        M1(0,i)=M2(row,index(i));
    }
}
void cpy_row_matrix_bak(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int row){
    for(int i=0;i<n;i++){
        M1(row,index(i))=M2(0,i);
    }
}

void inv_psi(Eigen::MatrixXd& psi,Eigen::MatrixXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i,i)=double(1)/psi(i,i);
    }
}
void inv_psi_vec(Eigen::VectorXd& psi,Eigen::VectorXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i)=double(1)/psi(i);
    }
}
double log_norm(double x,double e, double v){
  //const double PI  =3.141592653589793238462;
    double a = std::numeric_limits<double>::infinity();
    //cout << "a " << a << endl;
	
    if(v==0){
        return a;
    }else{
        return -0.5*(x-e)*(x-e)/v-0.5*log(v)-0.5*log(2*PI);
		//return(log(gsl_ran_gaussian_pdf(x-e,sqrt(v))));
    } 
	
	
}
double log_gamma(double x, double a, double beta){
    //return a*log(beta)-log(gsl_sf_gamma(a))+(a-1)*log(x)-beta*x;
	if(a==1&&x==0){
		return(beta);
	}else{
		return a*log(beta)-lgamma(a)+(a-1)*log(x)-beta*x;
	}
	//return(log(gsl_ran_gamma_pdf(x,a,double(1)/beta)));
    //boost::lambda
}
void cumsum(Eigen::VectorXd& S, Eigen::VectorXd& D, int n){
  for(int i=0;i<n;i++){
    if(i==0){
      D(i)=S(i);
    }
    D(i)=D(i-1)+S(i);
  }
}

double fx(double x,double c){
    return -1*log(x)+0.5/x+c;
}

double dfx(double x){
    return -1*(double(1)/x+0.5/x/x);
}

double NR(double c){
    double x=1e-10;
    for(int i=0;i<500;i++){
        x=x-fx(x,c)/dfx(x);
    }
    return x;
}


void cal_lam(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV, int s_n,int nf){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
		for(int i=0;i<nf;i++){              
			if(THETA(j,i)!=0 && PHI(i)!=0){
				indexALL(count_indexALL)=i;
				count_indexALL++;
			}
		}
            
		if(count_indexALL==0){
			LAM.row(j).setZero();
			continue;
		}
            
		partV.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}
							
		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			
		MatrixXd vl = partVI.lu().solve(IDNF);
		LAMI=partLI*vl;

		for(int i=0;i<count_indexALL;i++){
			if(THETA(j,indexALL(i))==0){
				LAMI(i)=0;
			}
		}
            
		cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
			
		LLI=LAMI.transpose()*LAMI;


		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				//vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
			
	}

}


double like(Eigen::MatrixXd& PSI,Eigen::MatrixXd& EXX,Eigen::MatrixXd& LAM,Eigen::MatrixXd& THETA,Eigen::MatrixXd& DELTA,Eigen::VectorXd& PHI,Eigen::VectorXd& TAU,Eigen::MatrixXd& Z,Eigen::MatrixXd& V,double ETA, double GAMMA,double alpha, double beta, int n,int p, int nf,double a, double b, double c, double d, double e, double f, double nu){
  double det_psi=0;
  
  for(int i=0;i<n;i++){
    det_psi = det_psi+log(PSI(i,i));
  }
  
  double like=(-1)*0.5*n*p*log(2*PI)-0.5*p*det_psi;
  
  double sum_x=0;
  for(int i=0;i<nf;i++){
    sum_x = sum_x + (-1)*0.5*EXX(i,i);
  }
  
  like = like - 0.5*nf*p*log(2*PI) + sum_x;
  
  for(int i=0;i<nf;i++){
    like=like + (Z(0,i)+alpha)*V(0,i) + (Z(1,i)+beta)*V(1,i);
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<nf;j++){
      if(THETA(i,j)!=0){
        like=like+Z(0,j)*log_norm(LAM(i,j),0,THETA(i,j));
        //if(DELTA(i,j)!=0){
          like=like+Z(0,j)*log_gamma(THETA(i,j),a,DELTA(i,j));
          //}
      }
      like=like+Z(0,j)*log_gamma(DELTA(i,j),b,PHI(j));
      like=like+Z(1,j)*log_norm(LAM(i,j),0,PHI(j));
    }
  }
 
 for(int i=0;i<nf;i++){
   like=like+log_gamma(PHI(i),c,TAU(i));
   like=like+log_gamma(TAU(i),d,ETA);
 }
 
 like=like+log_gamma(ETA,e,GAMMA);
 like=like+log_gamma(GAMMA,f,nu);
  
  return like;
 
}



/*
#include <iostream>     // cout
#include <math.h>       // acos
#include <float.h>      // DBL_MAX
#include <limits>       // numeric_limits
*/
template<typename T>
bool is_infinite( const T &value )
{
    // Since we're a template, it's wise to use std::numeric_limits<T>
    //
    // Note: std::numeric_limits<T>::min() behaves like DBL_MIN, and is the smallest absolute value possible.
    //
 
    T max_value = std::numeric_limits<T>::max();
    T min_value = - max_value;
 
    return ! ( min_value <= value && value <= max_value );
}
 
template<typename T>
bool is_nan( const T &value )
{
    // True if NAN
    return value != value;
}
 
template<typename T>
bool is_valid( const T &value )
{
    return ! is_infinite(value) && ! is_nan(value);
}




/*************************************
 An ANSI-C implementation of the digamma-function for real arguments based
 on the Chebyshev expansion proposed in appendix E of
 http://arXiv.org/abs/math.CA/0403344 . This is identical to the implementation
 by Jet Wimp, Math. Comp. vol 15 no 74 (1961) pp 174 (see Table 1).
 For other implementations see
 the GSL implementation for Psi(Digamma) in
 http://www.gnu.org/software/gsl/manual/html_node/Psi-_0028Digamma_0029-Function.html
 
 Richard J. Mathar, 2005-11-24
 **************************************/
#include <math.h>

#ifndef M_PIl
/** The constant Pi in high precision */
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

/** The digamma function in long double precision.
 * @param x the real value of the argument
 * @return the value of the digamma (psi) function at that point
 * @author Richard J. Mathar
 * @since 2005-11-24
 */
long double digammal(long double x)
{
    /* force into the interval 1..3 */
    if( x < 0.0L )
        return digammal(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;    /* reflection formula */
    else if( x < 1.0L )
        return digammal(1.0L+x)-1.0L/x ;
    else if ( x == 1.0L)
        return -M_GAMMAl ;
    else if ( x == 2.0L)
        return 1.0L-M_GAMMAl ;
    else if ( x == 3.0L)
        return 1.5L-M_GAMMAl ;
    else if ( x > 3.0L)
    /* duplication formula */
        return 0.5L*(digammal(x/2.0L)+digammal((x+1.0L)/2.0L))+M_LN2l ;
    else
    {
        /* Just for your information, the following lines contain
         * the Maple source code to re-generate the table that is
         * eventually becoming the Kncoe[] array below
         * interface(prettyprint=0) :
         * Digits := 63 :
         * r := 0 :
         *
         * for l from 1 to 60 do
         *     d := binomial(-1/2,l) :
         *     r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
         *     evalf(r) ;
         *     print(%,evalf(1+Psi(1)-r)) ;
         *o d :
         *
         * for N from 1 to 28 do
         *     r := 0 :
         *     n := N-1 :
         *
         *    for l from iquo(n+3,2) to 70 do
         *        d := 0 :
         *        for s from 0 to n+1 do
         *         d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
         *        od :
         *        if 2*l-n > 1 then
         *        r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
         *        fi :
         *    od :
         *    print(evalf((-1)^n*2*r)) ;
         *od :
         *quit :
         */
        static long double Kncoe[] = { .30459198558715155634315638246624251L,
            .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
            .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
            .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
            .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
            .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
            .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
            .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
            .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
            .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
            .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
            .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
            .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
            .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
            .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;
        
        register long double Tn_1 = 1.0L ;    /* T_{n-1}(x), started at n=1 */
        register long double Tn = x-2.0L ;    /* T_{n}(x) , started at n=1 */
        register long double resul = Kncoe[0] + Kncoe[1]*Tn ;
        
        x -= 2.0L ;
        
        for(int n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++)
        {
            const long double Tn1 = 2.0L * x * Tn - Tn_1 ;    /* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
            resul += Kncoe[n]*Tn1 ;
            Tn_1 = Tn ;
            Tn = Tn1 ;
        }
        return resul ;
    }
}

/** The optional interface to CREASO's IDL is added if someone has defined
 * the cpp macro export_IDL_REF, which typically has been done by including the
 * files stdio.h and idl_export.h before this one here.
 */
#ifdef export_IDL_REF
/** CALL_EXTERNAL interface.
 * A template of calling this C function from IDL  is
 * @verbatim
 * dg = CALL_EXTERNAL('digamma.so',X)
 * @endverbatim
 * @param argc the number of arguments. This is supposed to be 1 and not
 *    checked further because that might have negative impact on performance.
 * @param argv the parameter list. The first element is the parameter x
 *    supposed to be of type DOUBLE in IDL
 * @return the return value, again of IDL-type DOUBLE
 * @since 2007-01-16
 * @author Richard J. Mathar
 */
double digamma_idl(int argc, void *argv[])
{
    long double x = *(double*)argv[0] ;
    return (double)digammal(x) ;
}
#endif /* export_IDL_REF */

#ifdef TEST

/* an alternate implementation for test purposes, using formula 6.3.16 of Abramowitz/Stegun with the
 first n terms */
#include <stdio.h>

long double digammalAlt(long double x, int n)
{
    /* force into the interval 1..3 */
    if( x < 0.0L )
        return digammalAlt(1.0L-x,n)+M_PIl/tanl(M_PIl*(1.0L-x)) ;    /* reflection formula */
    else if( x < 1.0L )
        return digammalAlt(1.0L+x,n)-1.0L/x ;
    else if ( x == 1.0L)
        return -M_GAMMAl ;
    else if ( x == 2.0L)
        return 1.0L-M_GAMMAl ;
    else if ( x == 3.0L)
        return 1.5L-M_GAMMAl ;
    else if ( x > 3.0L)
        return digammalAlt(x-1.0L,n)+1.0L/(x-1.0L) ;
    else
    {
        x -= 1.0L ;
        register long double resul = -M_GAMMAl ;
        
        for( ; n >= 1 ;n--)
            resul += x/(n*(n+x)) ;
        return resul ;
    }
}
/*
int main(int argc, char *argv[])
{
    for( long double x=0.01 ; x < 5. ; x += 0.02)
        printf("%.2Lf %.30Lf %.30Lf %.30Lf\n",x, digammal(x), digammalAlt(x,100), digammalAlt(x,200) ) ;
}
*/
#endif /* TEST */









