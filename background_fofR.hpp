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

inline double Tbar(const double a, const cosmology cosmo)
{
	return - ( (cosmo.Omega_cdm + cosmo.Omega_b)/a/a/a + 4.*cosmo.Omega_Lambda );
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

			COUT << "loading background file: "<< filename <<endl;

			getline(file, line);
			getline(file, line);
			iss.str(line);
			iss >> lineCount;


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


inline double H_initial_fR(const double a, const double H, const double R, const double f, const double fr)
{
	return sqrt( (H*H + 0*a*a*(fr*R - f)/6.) / (1. + 0*fr) );
}


inline double R_initial_fR(const double a, const double eightpiG, const cosmology cosmo)
{
	return eightpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4. * (1 - cosmo.Omega_cdm - cosmo.Omega_b - bg_ncdm(a, cosmo)) ); // TODO: Maybe not accurate enough at this stage, but let's try this first
}


///////////////////////////////////////////////////////////
// Method when computing the background with RK4 solver
//////////////////////////////////////////////////////////

// System for f(R) background:
// a' = func_RK_adot(...)
// H' = func_RK_Hdot(...)
// R' = func_RK_Rdot(...)

inline double func_RK_adot(const double a, const double H)
{
	return a * H;
}

inline double func_RK_Hdot(const double a, const double H, const double R)
{
	return a*a*R/6. - H*H;
}

double func_RK_Rdot(const double a,
										const double H,
										const double R,
										const cosmology cosmo,
										const double fourpiG,
										double * params,
										const int fofR_type)
{
	double frr = FRR(R, params, fofR_type);
	double fr = FR(R, params, fofR_type);
	double rho = Hconf(a, fourpiG, cosmo);
	rho *= 3. * rho; // Actually computes 8 pi G * rho * a**2

	// return - 6 * fourpiG * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * H / a / a / a;

	if(3. * frr * H) return - (rho - 3.*H*H*(1. + fr) + 0.5*(fr * R - F(R, params, fofR_type))*a*a) / (3. * frr * H);
	else
	{
		COUT << " FRR evaluates to zero. Closing..." << endl
				 << " R = " << R << endl
				 << " H = " << H << endl
				 << " FRR = " << frr << endl;
		exit(3);
	}
}


// Runge-Kutta (4) solver for f(R) background
void rungekutta_fR(double &a,
									 double &H,
									 double &R,
									 const double fourpiG,
									 const cosmology cosmo,
									 const double dtau,
									 double * params,
									 const int fofR_type)
{
	double a1, a2, a3, a4;
	double H1, H2, H3, H4;
	double R1, R2, R3, R4;

	a1 = func_RK_adot(a, H);
	H1 = func_RK_Hdot(a, H, R);
	R1 = func_RK_Rdot(a, H, R, cosmo, fourpiG, params, fofR_type);

	// COUT << " a1 = " << a1 << "  H1 = " << H1 << "  R1 = " << R1 << endl;

	a2 = func_RK_adot(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = func_RK_Hdot(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = func_RK_Rdot(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, cosmo, fourpiG, params, fofR_type);

	// COUT << " a2 = " << a2 << "  H2 = " << H2 << "  R2 = " << R2 << endl;

	a3 = func_RK_adot(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = func_RK_Hdot(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = func_RK_Rdot(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, cosmo, fourpiG, params, fofR_type);

	// COUT << " a3 = " << a3 << "  H3 = " << H3 << "  R3 = " << R3 << endl;

	a4 = func_RK_adot(a + a3 * dtau, H + H3 * dtau);
	H4 = func_RK_Hdot(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = func_RK_Rdot(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, cosmo, fourpiG, params, fofR_type);

	// COUT << " a4 = " << a4 << "  H4 = " << H4 << "  R4 = " << R4 << endl;

	a += dtau * (a1 + 2.*a2 + 2.*a3 + a4) / 6.;
	H += dtau * (H1 + 2.*H2 + 2.*H3 + H4) / 6.;
	R += dtau * (R1 + 2.*R2 + 2.*R3 + R4) / 6.;

}




#endif
