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
	return -(cosmo.Omega_cdm + cosmo.Omega_b)/a/a/a - 4.*cosmo.Omega_Lambda;// TODO Put neutrinos in, if necessary
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
					if(parallel.rank()==0) parallel.abortForce();
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
				iss.seekg (0, iss.beg);
				iss >> tau[i] >> a[i] >> H[i] >> Rbar[i];
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

// TODO: At the moment, initial conditions are such that:
//	R_in = R_in(GR) = 8piG(rho_m + 4*rho_Lambda)
//	R'_in = R'_in(GR) = -24piG*H*rho_m
//	This gives H = (16piG rho + FR*R - F) * a^2 / 6 / (1 + FR - 24piG rho_m FRR)  -- see RK4 equation for R'
inline double R_initial_fR(const double a, const double eightpiG, const cosmology cosmo)
{
	return eightpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda ); // TODO: Maybe not accurate enough at this stage, but let's try this first.
}

inline double H_initial_fR(const double a, const double H, const double R, const double f, const double fr, const double frr_term)
{
	return H;
	// return sqrt( (H*H + a*a*(fr*R - f)/6.) / (1. + fr - frr_term) );
}




///////////////////////////////////////////////////////////
// Computing the background with RK4 solver
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

// For background_only option -- T00_background is computed from cosmological parameters Omega_m, Omega_Lambda etc
double func_RK_Rdot(const double a,
										const double H,
										const double R,
										const double rho, // Actually 8 pi G * rho_background * a**2 -- TODO: Should this be T00hom instead?
										const double fourpiG,
										double * params,
										const int fofR_type,
										const double t00hom = 1.)
{
	double frr = FRR(R, params, fofR_type, 25);
	double fr = FR(R, params, fofR_type);
	double denom = 3. * frr * H;

	if(denom)
	{
		double rdot = (rho - 3.*H*H*(1. + fr) + 0.5*(fr * R - F(R, params, fofR_type))*a*a) / denom;
		return rdot;
	}
	else
	{
		COUT << " FRR evaluates to zero. Closing..." << endl
				 << " R = " << R << endl
				 << " H = " << H << endl
				 << " FRR = " << frr << endl;
		exit(3);
	}
}




////////////////////////
// Runge-Kutta (4) solver for f(R) background
// TODO: more info here
////////////////////////
double rungekutta_fR(double & a,
									   double & H,
									   double & R,
									   const double fourpiG,
									   const cosmology cosmo,
									   const double dtau,
									   double * params,
									   const int fofR_type)
{
	double a1, a2, a3, a4;
	double H1, H2, H3, H4;
	double R1, R2, R3, R4;
	double rho0 = 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo);
	// the density appearing in the equation for R is the WHOLE background density: particles + radiation + Lambda

	a1 = func_RK_adot(a, H);
	H1 = func_RK_Hdot(a, H, R);
	R1 = func_RK_Rdot(a, H, R, rho0, fourpiG, params, fofR_type);

	// COUT << " a1 = " << a1 << "  H1 = " << H1 << "  R1 = " << R1 << endl;

	a2 = func_RK_adot(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = func_RK_Hdot(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = func_RK_Rdot(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, 3. * Hconf(a + 0.5 * a1 * dtau, fourpiG, cosmo) * Hconf(a + 0.5 * a1 * dtau, fourpiG, cosmo), fourpiG, params, fofR_type);

	// COUT << " a2 = " << a2 << "  H2 = " << H2 << "  R2 = " << R2 << endl;

	a3 = func_RK_adot(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = func_RK_Hdot(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = func_RK_Rdot(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, 3. * Hconf(a + 0.5 * a2 * dtau, fourpiG, cosmo) * Hconf(a + 0.5 * a2 * dtau, fourpiG, cosmo), fourpiG, params, fofR_type);

	// COUT << " a3 = " << a3 << "  H3 = " << H3 << "  R3 = " << R3 << endl;

	a4 = func_RK_adot(a + a3 * dtau, H + H3 * dtau);
	H4 = func_RK_Hdot(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = func_RK_Rdot(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, 3. * Hconf(a + a3 * dtau, fourpiG, cosmo) * Hconf(a + a3 * dtau, fourpiG, cosmo), fourpiG, params, fofR_type);

	// COUT << " a4 = " << a4 << "  H4 = " << H4 << "  R4 = " << R4 << endl;
	a += dtau * (a1 + 2.*a2 + 2.*a3 + a4) / 6.;
	H += dtau * (H1 + 2.*H2 + 2.*H3 + H4) / 6.;
	R += dtau * (R1 + 2.*R2 + 2.*R3 + R4) / 6.;

	return dtau;
}



////////////////////////
// Runge-Kutta-Fehlberg solver for f(R) background
// TODO: more info here
////////////////////////
double rungekutta_fR_45(double & a,
									 double & H,
									 double & R,
									 const double fourpiG,
									 const cosmology cosmo,
									 const double dtau,
									 double * params,
									 const int fofR_type)
{
	double ac, af, Hc, Hf, Rc, Rf;
	double a1, a2, a3, a4, a5, a6;
	double H1, H2, H3, H4, H5, H6;
	double R1, R2, R3, R4, R5, R6;
	double s, dtau_temp = dtau;
	double tol = 1.E-2;
	double rho0 = 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo);

	a1 = dtau * func_RK_adot(a, H);
	H1 = dtau * func_RK_Hdot(a, H, R);
	R1 = dtau * func_RK_Rdot(a, H, R, rho0, fourpiG, params, fofR_type);

	a2 = dtau * func_RK_adot(a + 0.25 * a1, H + 0.25 * H1);
	H2 = dtau * func_RK_Hdot(a + 0.25 * a1, H + 0.25 * H1, R + 0.25 * R1);
	R2 = dtau * func_RK_Rdot(a + 0.25 * a1, H + 0.25 * H1, R + 0.25 * R1, 3. * Hconf(a + 0.25 * a1, fourpiG, cosmo) * Hconf(a + 0.25 * a1, fourpiG, cosmo), fourpiG, params, fofR_type);

	a3 = dtau * func_RK_adot(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2);
	H3 = dtau * func_RK_Hdot(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2, R + 3./32. * R1 + 9./32. * R2);
	R3 = dtau * func_RK_Rdot(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2, R + 3./32. * R1 + 9./32. * R2, 3. * Hconf(a + 3./32. * a1 + 9./32. * a2, fourpiG, cosmo) * Hconf(a + 3./32. * a1 + 9./32. * a2, fourpiG, cosmo), fourpiG, params, fofR_type);

	a4 = dtau * func_RK_adot(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3);
	H4 = dtau * func_RK_Hdot(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3, R + 1932./2197. * R1 - 7200./2197. * R2 + 7296./2197. * R3);
	R4 = dtau * func_RK_Rdot(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3, R + 1932./2197. * R1 - 7200./2197. * R2 + 7296./2197. * R3, 3. * Hconf(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, fourpiG, cosmo) * Hconf(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, fourpiG, cosmo), fourpiG, params, fofR_type);

	a5 = dtau * func_RK_adot(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4);
	H5 = dtau * func_RK_Hdot(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4, R + 439./216. * R1 - 8. * R2 + 3680./513. * R3 - 845./4104. * R4);
	R5 = dtau * func_RK_Rdot(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4, R + 439./216. * R1 - 8. * R2 + 3680./513. * R3 - 845./4104. * R4, 3. * Hconf(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, fourpiG, cosmo) * Hconf(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, fourpiG, cosmo), fourpiG, params, fofR_type);

	a6 = dtau * func_RK_adot(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5);
	H6 = dtau * func_RK_Hdot(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5, R - 8./27. * R1 + 2. * R2 - 3544./2565. * R3 + 1859./4104. * R4 - 11./40. * R5);
	R6 = dtau * func_RK_Rdot(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5, R - 8./27. * R1 + 2. * R2 - 3544./2565. * R3 + 1859./4104. * R4 - 11./40. * R5, 3. * Hconf(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, fourpiG, cosmo) * Hconf(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, fourpiG, cosmo), fourpiG, params, fofR_type);

	ac = a + 25./216. * a1 + 1408./2565. * a3 + 2197./4101. * a4 - 1./5. * a5;
	Hc = H + 25./216. * H1 + 1408./2565. * H3 + 2197./4101. * H4 - 1./5. * H5;
	Rc = R + 25./216. * R1 + 1408./2565. * R3 + 2197./4101. * R4 - 1./5. * R5;

	a = af = a + 16./135. * a1 + 6656./12825. * a3 + 28561./56430. * a4 - 9./50. * a5 + 2./55. * a6;
	H = Hf = H + 16./135. * H1 + 6656./12825. * H3 + 28561./56430. * H4 - 9./50. * H5 + 2./55. * H6;
	R = Rf = R + 16./135. * R1 + 6656./12825. * R3 + 28561./56430. * R4 - 9./50. * R5 + 2./55. * R6;

	if(Rf - Rc)
	{
		s = 0.84 * sqrt(
			sqrt(
				tol * dtau / (fabs(Rf - Rc))
			)
		);
		dtau_temp = s*dtau;
	}

	if(Hf - Hc)
	{
		s = 0.84 * sqrt(
			sqrt(
				tol * dtau / (fabs(Hf - Hc))
			)
		);
		if(s*dtau < dtau_temp) dtau_temp = s*dtau;
	}

	if(af - ac)
	{
		s = 0.84 * sqrt(
			sqrt(
				tol * dtau / (fabs(af - ac))
			)
		);
		if(s*dtau < dtau_temp) dtau_temp = s*dtau;
	}

	return dtau;

}









#endif
