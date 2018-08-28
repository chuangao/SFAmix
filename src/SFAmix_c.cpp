#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>


#include <iostream>
#include <Eigen/Dense>
#include "myHeader.cpp"


#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
//#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

//const double PI  =3.141592653589793238462;
/*
  Chuan Gao C++
*/

/* 
   usage
   ./SFAmix_conditional_psi_vec --nf 50 --y sim_data/Y_sparse_noise1_1.txt --out result --sep tab --a 0.1 --b 0.1 --c 0.5 --d 0.5 --e 1 --f 1
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 0.5 --b 0.5 --c 0.5 --d 0.5 --e 0.5 --f 0.5
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 1 --b 0.5 --c 1 --d 0.5 --e 1 --f 0.5
   ./SFAmix_conditional --nf 50 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 1 --b 1 --c 1 --d 1 --e 1 --f 1
   ./SFAmix --nf 20 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 1 --b 1 --c 1 --d 1 --e 1 --f 1
   ./SFAmix --nf 20 --y ../data/sim_data_SFA/Y_1.txt --out result --sep tab --a 0.5 --b 0.5 --c 0.5 --d 0.5 --e 0.5 --f 0.5
   ./SFAmix --y /gpfs/fs0/data/engelhardtlab/cg148/data/CAP/S480RjQ2N_t_noCrossHyb.dat --nf 50 --sep space --out result
*/







extern "C" void SFAmix(double *Y_TMP_param ,int *nrow_param, int *ncol_param, double *a_param,double *b_param, int *nf_param, int *itr_param, double *LAM_out, double *EX_out, double *Z_out, double *EXX_out, int *nf_out, int *out_itr, char **output_dir,int *itr_final){
    
    double a = *a_param;
    double b = *b_param;
    int nf = *nf_param;
	int s_n = *nrow_param;
	int d_y = *ncol_param;
	int n_itr = *itr_param;
    
    int itr_conv = 0;
    
    int interval = 500;
    int write_itr = *out_itr;
    string out_dir = *output_dir;
    std::replace( out_dir.begin(), out_dir.end(), '%', '/');
    
    stringstream ss;
    

	//cout << "a " << a << endl;

   
	//cout << "nf " <<  nf << endl;
	//cout << "s_n " << s_n << endl;
   //cout << "a " << a << endl;
   //cout << "a " << a << endl;
   

	double c=0.5,d=0.5,g=0.5,h=0.5,alpha=1,beta=1;

    MatrixXd Y=MatrixXd::Constant(s_n,d_y,0);
    
    for(int i = 0; i < d_y; i++){
        for(int j = 0; j < s_n; j++){
            Y(j,i) = Y_TMP_param[i*s_n+j];
			//cout << "y_i_j" << Y(j,i) << endl;
        }
    }
    //cout << 'Y' << endl;
	//cout << Y.block(0,0,3,3) << endl;
    
    VectorXd psi_v = VectorXd::Constant(s_n,1);
    VectorXd PSI=VectorXd::Constant(s_n,1);
    VectorXd PSI_INV=VectorXd::Constant(s_n,1);
    MatrixXd LX=MatrixXd::Constant(s_n,d_y,0);
    MatrixXd YLX=MatrixXd::Constant(s_n,d_y,0);
    MatrixXd YLXi=VectorXd::Constant(s_n,0);
    MatrixXd YLX2=MatrixXd::Constant(s_n,d_y,0);
    VectorXd TOP=VectorXd::Constant(s_n,0);
    
    // Declare variables that are dependent of the factor number
    int nf2 = nf;
    int nt=nf;
    
    MatrixXd EX=MatrixXd::Constant(nt,d_y,0);
    MatrixXd TEX=MatrixXd::Constant(d_y,nt,0);
    MatrixXd VX=MatrixXd::Constant(nt,nt,0);
    MatrixXd EXX=MatrixXd::Constant(nt,nt,0);
    MatrixXd LP=MatrixXd::Constant(nt,s_n,0);
    MatrixXd ID=MatrixXd::Constant(nt,nt,0);
    VectorXd id_v = VectorXd::Constant(nt,1);
    
    MatrixXd LAM=MatrixXd::Constant(s_n,nt,0);
    //MatrixXd LAM_T=MatrixXd::Constant(nt,s_n,0);
    MatrixXd LAM_BAK=MatrixXd::Constant(s_n,nt,0);
    MatrixXd THETA=MatrixXd::Constant(s_n,nf,0);
    MatrixXd DELTA=MatrixXd::Constant(s_n,nf,0);
    VectorXd PHI = VectorXd::Constant(nf,0);
    VectorXd TAU = VectorXd::Constant(nf,0);
    double nu = 1;
    double ETA = 1;
    double GAMMA = 1;
    
    VectorXd count_lam = VectorXd::Constant(nf,0);
    VectorXd index = VectorXd::Constant(nf,0);
    
    
    MatrixXd LAM_TOP=MatrixXd::Constant(s_n,nf,0);
    MatrixXd LAM_BOT=MatrixXd::Constant(s_n,nf,0);
    
    double nmix=2;
    double zi = double(1)/nmix;
    MatrixXd Z = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGV = MatrixXd::Constant(nmix,1,log(0.5));
    MatrixXd Zm = MatrixXd::Constant(nmix,nf,log(zi));
    
    MatrixXd LPL = MatrixXd::Constant(nf,nf,0);

	
	long seed;
	
	gsl_rng *r;  // random number generator
	r=gsl_rng_alloc(gsl_rng_mt19937);
	
	seed = time (NULL) * getpid();
	
	//seed=20877551893284;
	gsl_rng_set (r, seed);                  // set seed

	//gsl_rng_set (r, i_seed);
	
    // Initialize parameters 
    ID.diagonal()=id_v;
   
    //PSI.diagonal() = psi_v;
    //inv_psi(PSI,PSI_INV,s_n);

    // fill in the lambda matrix  
    for (int i=0; i<s_n; i++) {
        for(int j=0;j<nt;j++){
            LAM(i,j)=gsl_ran_gaussian(r,1);
        }
    }
	
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            THETA(i,j)=1;
            DELTA(i,j)=1;
        }
    }
    
    // continue initializing 
    for(int i=0;i<nf;i++){
        PHI(i)=1;
        TAU(i)=1;
    }
    VectorXd lam_count_v = VectorXd::Constant(n_itr,0);

	MatrixXd LAM_T=LAM.transpose();
    for(int i=0;i<s_n;i++){
        LP.col(i)=LAM_T.col(i)*PSI_INV(i);
    }
	VX=(LP*LAM+ID).inverse();
	EX=VX*LP*Y;
	EXX=EX*EX.transpose()+VX*d_y;
	
	TEX=EX.transpose();

	// int gene_start=0;
	// int gene_stop=0;
	
    for(int itr=0;itr<(n_itr-1);itr++){

		if(itr%10==0){
		  //cout << out_dir << endl;
            //cout << "after reduction" << endl;
            cout << "itr " << itr << endl;
            cout << "number of factors " << nf << endl;
            
            cout << "count_lam" << endl << count_lam.transpose() << endl;
			//cout << "out_dir" << endl << out_dir << endl;

			/*
			cout << "TAU " << TAU(1) << endl;
			cout << "PHI " << PHI(1) << endl;
			cout << "DELTA " << DELTA.block(0,0,2,2) << endl;
			cout << "THETA " << THETA.block(0,0,2,2) << endl;
			cout << "LAM " << LAM.block(0,0,2,2) << endl;
			cout << "EX " << EX.block(0,0,2,2) << endl;
			cout << "Z " << Z.block(0,0,2,2) << endl;
			cout << "EXX " << EXX.block(0,0,2,2) << endl;
			*/
         }
        
	
        // Expectation step

     	MatrixXd LAM_T=LAM.transpose();
		for(int i=0;i<s_n;i++){
			LP.col(i)=LAM_T.col(i)*PSI_INV(i);
		}
        //VX=(LPL+ID).lu().solve(ID);
		VX=(LP*LAM+ID).lu().solve(ID);
        EX=VX*LP*Y;
        EXX=EX*EX.transpose()+VX*d_y;		
        TEX=EX.transpose();
        
        // Mazimization step 
        double alpha_sum=TAU.sum();
        ETA=double((d*nf+g-1))/(GAMMA+alpha_sum);
        GAMMA=double(g+h)/(ETA+nu);
        for(int i=0;i<nf;i++){
            TAU(i)=double(c+d)/(PHI(i)+ETA);
			
		}
		for(int i=0;i<nf;i++){
			double sum_c=s_n*b*Z(0,i)+c-1-0.5*s_n*Z(1,i);
			double at = 2*(TAU(i)+Z(0,i)*(DELTA.col(i).sum()));
			double lam_sum=LAM.col(i).dot(LAM.col(i));
			double bt = Z(1,i)*lam_sum;
			PHI(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
			
			// if(PHI(i)<1e-300){
			// 	PHI(i)=1e-300;
			// }
			
		}

		
		double min_phi=1e100;
		for(int i=1;i<nf;i++){
			if(PHI(i)!=0&PHI(i)<min_phi){
				min_phi=PHI(i);
			}
		}
		for(int i=1;i<nf;i++){
			if(PHI(i)==0){
				PHI(i)=min_phi;
			}
		}
		     
		// for(int i=1;i<nf;i++){
		// 	if(PHI(i)<1e-100){
		// 		PHI(i)=1e-100;
		// 	}
		// }
	
		//MatrixXd range_LAM_cov = MatrixXd::Constant(2,nfc,0);
		//MatrixXd range_EX_cov = MatrixXd::Constant(2,nfc,0);	
		 
		// range_colwise(range_LAM,LAM,s_n,nf);
		// range_rowwise(range_EX,EX,nf,d_y);
		// range_colwise(range_THETA,THETA,s_n,nf);
		
		// count the number of loadings that are set to all zero
		count_lam.setZero();
		for(int i=0;i<nf;i++){
			for(int j=0;j<s_n;j++){
				if(LAM(j,i)!=0){
					count_lam(i) +=  1;
				}
			}
		}

		//exhaust the gene list, assign each gene to an empty list
                
		/*       
				 for(int i=0;i<nf;i++){
				 if(count_lam(i)==0){
				 if(gene_start<s_n){
				 //int gene_index = rand()%s_n;
				 LAM(gene_start,i)=gsl_ran_gaussian(r,1);
				 THETA(gene_start,i)=2;
				 DELTA(gene_start,i)=2;
				 count_lam(i)=1;
				 gene_start ++;
				 }
				 }
				 }
		*/
        
		// Count the number of loadings that are active, either all zero or PHI_k zero will kill 
		int count_nonzero = 0;
		for(int i=0;i<nf;i++){
			if(count_lam(i)>1&&PHI(i)!=0){
				index(count_nonzero)=i;
				count_nonzero++;
			}
		}
		// remove inactive factors, loadings etc. and assign to new matrix
		//cout << "number of factors before " << nf << endl;
		if(count_nonzero != nf){

			// MatrixXd cbind=MatrixXd::Constant(nf,4,0);
			// cbind.col(0)=PHI;
			// cbind.col(1)=count_lam;
			// cbind.col(2)=Z.row(0);
			// cbind.col(3)=index;

			//cout << "cbind " << count_nonzero <<endl << cbind <<endl;

			
			nf=count_nonzero;
			nt=nf;
			MatrixXd EX2=MatrixXd::Constant(nt,d_y,0);
			MatrixXd TEX2=MatrixXd::Constant(d_y,nt,0);
			MatrixXd VX2=MatrixXd::Constant(nt,nt,0);
			MatrixXd EXX2=MatrixXd::Constant(nt,nt,0);
			MatrixXd LP2=MatrixXd::Constant(nt,s_n,0);
			MatrixXd ID2=MatrixXd::Constant(nt,nt,0);
			VectorXd id_v2 = VectorXd::Constant(nt,1);
			MatrixXd LAM2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd LAM_BAK2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd THETA2=MatrixXd::Constant(s_n,nf,0);
			MatrixXd DELTA2=MatrixXd::Constant(s_n,nf,0);
			VectorXd PHI2 = VectorXd::Constant(nf,0);
			VectorXd TAU2 = VectorXd::Constant(nf,0);
			MatrixXd Z2 = MatrixXd::Constant(nmix,nf,zi);
			MatrixXd logZ2 = MatrixXd::Constant(nmix,nf,log(zi));
			VectorXd count_lam2 = VectorXd::Constant(nf,0);
			VectorXd index2 = VectorXd::Constant(nf,0);
			MatrixXd LAM_TOP2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd LAM_BOT2=MatrixXd::Constant(s_n,nt,0);
			MatrixXd LPL2 = MatrixXd::Constant(nf,nf,0);
			ID2.diagonal()=id_v2;
                     
			for(int i=0;i<nf;i++){
				EX2.row(i)=EX.row(index(i));
				//TEX2=EX2.transpose();
				for(int j=0;j<nf;j++){
					VX2(i,j)=VX(index(i),index(j));
					EXX2(i,j)=EXX(index(i),index(j));
					LPL2(i,j)=LPL(index(i),index(j));
				}
				LP2.row(i)=LP.row(index(i));
				LAM2.col(i)=LAM.col(index(i));
				LAM_BAK2.col(i)=LAM_BAK.col(index(i));
				THETA2.col(i)=THETA.col(index(i));
				DELTA2.col(i)=DELTA.col(index(i));
				PHI2(i)=PHI(index(i));
				TAU2(i)=TAU(index(i));
				Z2.col(i)=Z.col(index(i));
				logZ2.col(i)=logZ.col(index(i));
				LAM_TOP2.col(i)=LAM_TOP.col(index(i));
				LAM_BOT2.col(i)=LAM_BOT.col(index(i));
				count_lam2(i)=count_lam(index(i));
				index2(i)=index(i);
			}
            
			// Assign the new parameters back 
			EX=EX2;
			//TEX=TEX2;
			VX=VX2;
			EXX=EXX2;
			LP=LP2;
			ID=ID2;
			LAM=LAM2;
			LAM_BAK=LAM_BAK2;
			THETA=THETA2;
			DELTA=DELTA2;
			PHI=PHI2;
			TAU=TAU2;
			Z=Z2;
			logZ=logZ2;
			LAM_TOP=LAM_TOP2;
			LAM_BOT=LAM_BOT2;
			count_lam=count_lam2;
			index=index2;
			LPL=LPL2;
            
		}
		//cout << "number of factors after" << nf << endl;
        
		// continue updating parameters 
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				DELTA(i,j)=double((a+b))/(THETA(i,j)+PHI(j));
			}
		}
		
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				double a23=(2*a-3);
				THETA(i,j)=double(a23+sqrt(a23*a23+8*LAM(i,j)*LAM(i,j)*DELTA(i,j)))/4/DELTA(i,j);
			}
		}

        
		// Need to update separately because 0/0=NA
		/*
		  VectorXd indexALL = VectorXd::Constant(nf,0);
		  MatrixXd partV=MatrixXd::Constant(nf,nf,0);
		  MatrixXd partL = MatrixXd::Constant(s_n,nf,0);

		  cal_lam(LAM,indexALL, partL, Y, EX, PSI_INV, EXX, Z, LPL, THETA,PHI, partV, s_n,nf);
		*/

				
		YLX=Y-LAM*EX;
		for(int i=0;i<nf;i++){
			YLX=YLX+LAM.col(i)*(EX.row(i));
			//YLX=YLX+LAM.col(i)*(TEX.col(i).transpose());
			for(int j=0;j<s_n;j++){
				TOP(j) = YLX.row(j).dot(EX.row(i));
			}
                    
			for(int j=0;j<s_n;j++){
				if(Z(0,i)==0){
					LAM(j,i) = TOP(j)/(EXX(i,i)+PSI(j)*Z(1,i)/PHI(i));
				}
				else if(Z(1,i)==0){
					LAM(j,i) = TOP(j)/(EXX(i,i)+PSI(j)*Z(0,i)/THETA(j,i));
				}
				else{
					LAM(j,i) = TOP(j)/(EXX(i,i)+PSI(j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
				}
			}
			YLX=YLX-LAM.col(i)*EX.row(i);                   
		}



		/*
		  YLX=Y*TEX-LAM*EXX;
		  for(int i=0;i<nf;i++){
		  YLX.col(i)=YLX.col(i)+LAM.col(i)*(EXX(i,i));
		  for(int j=0;j<s_n;j++){
		  if(Z(0,i)==0){
		  LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*Z(1,i)/PHI(i));
		  }
		  else if(Z(1,i)==0){
		  LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*Z(0,i)/THETA(j,i));
		  }
		  else{
		  LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
		  }
		  }
		  YLX.col(i)=YLX.col(i)-LAM.col(i)*(EXX(i,i));
		  }
		*/		
		// If theta=0, then LAM=0

		/*
		  VectorXd maxEX = VectorXd::Constant(nf,0);
		
		  for(int k=0;k<nf;k++){
		  maxEX(k)=EX(k,0);
		  for(int i=1;i<d_y;i++){
		  if(abs(maxEX(k))<abs(EX(k,i))){
		  maxEX(k)=EX(k,i);
		  }
				
		  }
	
		  for(int i=1;i<s_n;i++){
			
		  if(abs(maxEX(k)*LAM(i,k))<1e-10){
		  LAM(i,k)=0;
		  }
		  }
	
		  }
		
		*/
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				if(THETA(i,j)==0){
					LAM(i,j)=0;
				}
				if(abs(LAM(i,j))<1e-10){
					LAM(i,j)=0;
				}
			}
		}
		
		// logZ 
		for(int i=0;i<nf;i++){
			logZ(0,i)=LOGV(0,0);
			logZ(1,i)=LOGV(1,0);
			for(int j=0;j<s_n;j++){
				logZ(0,i)=logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+log_gamma(THETA(j,i),a,DELTA(j,i))+log_gamma(DELTA(j,i),b,PHI(i));
				logZ(1,i)=logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
			}

		}	
        
		// Z 
		for(int i=0;i<nf;i++){
			Z(0,i)=double(1)/(1+exp(logZ(1,i)-logZ(0,i)));
			Z(1,i)=1-Z(0,i);
		}
        
		double ps1=alpha;
		double ps2=beta;
        
		// Probability 
		for(int i=0;i<nf;i++){
			ps1=ps1+Z(0,i);
			ps2=ps2+Z(1,i);
		}
		double dgama = gsl_sf_psi(ps1+ps2);
		LOGV(0,0)=gsl_sf_psi(ps1)-dgama;
		LOGV(1,0)=gsl_sf_psi(ps2)-dgama;
               
		// PSI 
		LX=LAM*EX;
		for(int i=0;i<s_n;i++){
			PSI(i)=Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i));
		}
		PSI=PSI/d_y;
		inv_psi_vec(PSI,PSI_INV,s_n);
        
		// LAM active 
		int lam_count=0;
		for(int i=0;i<s_n;i++){
			for(int j=0;j<nf;j++){
				if(LAM(i,j)!=0){
					lam_count++;
				}
			}
		}

		MatrixXd range_LAM = MatrixXd::Constant(2,nf,0);
		MatrixXd range_EX = MatrixXd::Constant(2,nf,0);
     
        lam_count_v(itr+1)=lam_count;
  
        // claim convergence if the number of values is stable for 200 iterations.
        if(itr>=interval){
            //if(lam_count_v(itr+1)!=(s_n*nf)&&(lam_count_v(itr-interval)-lam_count_v(itr+1))<(0.05*lam_count_v(itr-interval))){
            if(lam_count_v(itr)!=(s_n*nf)&&(lam_count_v(itr-interval)-lam_count_v(itr))==0){
                itr_conv = itr;
				break; 
            }   
        }
        if(out_dir.compare("NULL") != 0 & itr % write_itr == 0){
            
            ss.str("");
            ss.clear();
            ss << out_dir << "/LAM_" << itr;
            ofstream f_lam (ss.str().c_str());
            if (f_lam.is_open()){
                f_lam << LAM << endl;
            }
            f_lam.close();
            
            ss.str("");
            ss.clear();
            ss << out_dir << "/EX_" << itr;
            ofstream f_ex (ss.str().c_str());
            if (f_ex.is_open()){
                f_ex << EX << endl;
            }
            f_ex.close();
            
            ss.str("");
            ss.clear();
            ss << out_dir << "/Z_" << itr;
            ofstream f_z (ss.str().c_str());
            if (f_z.is_open()){
                f_z << Z.row(1) << endl;
            }
            f_z.close();
            
            ss.str("");
            ss.clear();
            ss << out_dir << "/EXX_" << itr;
            ofstream f_exx (ss.str().c_str());
            if (f_exx.is_open()){
                f_exx << EXX << endl;
            }
            f_exx.close();
            
            ss.str("");
            ss.clear();
            ss << out_dir << "/nf_" << itr;
            ofstream f_nf (ss.str().c_str());
            if (f_nf.is_open()){
                f_nf << nf << endl;
            }
            f_nf.close();
        }
    }

	cout << "nf " << nf << endl;
	//LAM_out[0] = 1.5;
	
	for(int i = 0; i < nf; i++){
		for(int j = 0; j < s_n; j++){
			LAM_out[i*s_n + j] = LAM(j,i);
		}
	}
	
	for(int i = 0; i < d_y; i++){
		for(int j = 0; j < nf; j++){
			EX_out[i*nf + j] = EX(j,i);
		}
	}
	
    for(int i = 0; i < nf; i++){
		for(int j = 0; j < nf; j++){
			EXX_out[i*nf + j] = EXX(j,i);
		}
	}
    
    for(int i = 0; i < nf; i++){
        Z_out[i] = Z(1,i); 
    }
    
	nf_out[0] = nf;
    
    itr_final[0] = itr_conv;
	
}



