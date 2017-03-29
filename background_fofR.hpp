//////////////////////////
// background_fofR.hpp
//////////////////////////
//
// code components related to background evolution for the f(R) model
//
// Author: David Daverio (Cambridge University)
//         Lorenzo Reverberi (Cape Town University)
//
// Last modified: February 2016
//
// Note: Some methods have been copy pasted from background.hpp
//
//////////////////////////

#ifndef BACKGROUND_FOFR_HEADER
#define BACKGROUND_FOFR_HEADER

#include <cstdio>
#include <sstream>
#include <string>
#include <iostream>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

inline double Tbar(const double a,const cosmology cosmo)
{
	return Omega_m(a,cosmo) + 4.0*Omega_Lambda(a,cosmo);
}

void loadBackground(gsl_spline *& a_spline,
										gsl_spline *& H_spline,
										gsl_spline *& Rbar_spline,
	  								const gsl_interp_type *t, const char * filename)
{
		ifstream file;
		string line;
		long lineCount;
		istringstream iss;


		if(parallel.grid_rank()[0]==0)
		{
			file.open(filename);
			if(!file)
			{
					COUT<<"Background file error: cannot open file"<<endl;
					if(parallel.rank()==0)parallel.abortForce();
					return;
			}

			COUT<<"loading background file: "<< filename <<endl;

			getline(file, line);
			getline(file, line);
			iss.str(line);
			iss>>lineCount;


		}
		parallel.broadcast_dim0(lineCount,0);

		double * tau = new double[lineCount];
		double * a = new double[lineCount];
		double * H = new double[lineCount];
		double * Rbar = new double[lineCount];



		if(parallel.grid_rank()[0]==0)
		{
			for(int i=0;i<lineCount;i++)
			{
				getline(file, line);
				iss.str(line);
				//COUT<<line<<endl;
				//COUT<< iss.str() <<endl;
				iss.seekg (0, iss.beg);
				iss >> tau[i] >> a[i] >> H[i] >> Rbar[i];
				//COUT << tau[i] <<" "<< a[i] <<" "<< H[i] <<" "<< phi[i] <<" "<< sigma[i] <<" "<< sigmaDot[i]<<endl;
			}
		}


		parallel.broadcast_dim0(tau,lineCount,0);
		parallel.broadcast_dim0(a,lineCount,0);
		parallel.broadcast_dim0(H,lineCount,0);
		parallel.broadcast_dim0(Rbar,lineCount,0);




		a_spline  = gsl_spline_alloc (t, lineCount);
		H_spline  = gsl_spline_alloc (t, lineCount);
		Rbar_spline  = gsl_spline_alloc (t, lineCount);



		gsl_spline_init (a_spline, tau, a, lineCount);
		gsl_spline_init (H_spline, tau, H, lineCount);
		gsl_spline_init (Rbar_spline, tau, Rbar, lineCount);


		if(parallel.grid_rank()[0]==0)file.close();

		delete[] tau;
		delete[] a;
		delete[] H;
		delete[] Rbar;
}

inline double Hconf(gsl_spline * spline, double tau, gsl_interp_accel * acc)
{
		return gsl_spline_eval(spline, tau, acc);
}

inline double getBackground(gsl_spline * spline, double tau, gsl_interp_accel * acc)
{
	  return gsl_spline_eval(spline, tau, acc);
}




#endif
