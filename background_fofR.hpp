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

#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

void loadBackground(gsl_spline *& a_spline,
										gsl_spline *& H_spline,
										gsl_spline *& phi_spline,
									 	gsl_spline *& sigma_spline,
										gsl_spline *& sigmaDot_spline,
	  								const gsl_interp_type *t, const char * filename)
{
		ifstream file;
		string line;
		long lineCount;
		istringstream iss;

		file.open(filename);
		if(!file)
		{
			COUT<<"Background file error: cannot open file"<<endl;
			if(parallel.rank()==0)parallel.abortForce();
			return;
		}

		COUT<<"loading background file: "<< filename <<endl;

		getline(file, line);
		COUT<<"line"<<endl;
		getline(file, line);
		iss.str(line);
		iss>>lineCount;



		double * tau = new double[lineCount];
		double * a = new double[lineCount];
		double * H = new double[lineCount];
		double * phi = new double[lineCount];
		double * sigma = new double[lineCount];
		double * sigmaDot = new double[lineCount];

		for(int i=0;i<lineCount;i++)
		{
			getline(file, line);
			iss.str(line);
			iss >> tau[i] >> a[i] >> H[i] >> phi[i] >> sigma[i] >> sigmaDot[i];
		}


		a_spline  = gsl_spline_alloc (t, lineCount);
		H_spline  = gsl_spline_alloc (t, lineCount);
		phi_spline  = gsl_spline_alloc (t, lineCount);
		sigma_spline  = gsl_spline_alloc (t, lineCount);
		sigmaDot_spline  = gsl_spline_alloc (t, lineCount);


		gsl_spline_init (a_spline, tau, a, lineCount);
		gsl_spline_init (H_spline, tau, H, lineCount);
		gsl_spline_init (phi_spline, tau, phi, lineCount);
		gsl_spline_init (sigma_spline, tau, sigma, lineCount);
		gsl_spline_init (sigmaDot_spline, tau, sigmaDot, lineCount);

		file.close();

		delete[] tau;
		delete[] a;
		delete[] H;
		delete[] phi;
		delete[] sigma;
		delete[] sigmaDot;
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
