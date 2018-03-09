//////////////////////////
// background_fR.hpp
//////////////////////////
//
// code components related to background evolution for the f(R) model
//
// Author: David Daverio (Cambridge University)
//         Lorenzo Reverberi (Cape Town University, Czech Academy of Sciences & CEICO, Prague)
//
// Last modified: February 2016
//
// Note: Some methods have been copy pasted from background.hpp
//
//////////////////////////

#ifndef BACKGROUND_FR_HEADER
#define BACKGROUND_FR_HEADER

#include <cstdio>
#include <sstream>
#include <string>
#include <iostream>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

double R_dot_RungeKutta(const double, const double, const double,	const double,	const metadata &);

inline double Tbar(const double a, const cosmology & cosmo)
{
	return - cosmo.Omega_m/a/a/a - 4. * cosmo.Omega_Lambda;// TODO Put massive neutrinos in, if necessary
}

double R_GR(const double a, const double fourpiG, const cosmology & cosmo)
{
	return 2. * fourpiG * ( cosmo.Omega_m/a/a/a + 4. * (1. - cosmo.Omega_m - cosmo.Omega_rad) );
}

// TODO: At the moment, initial conditions are such that:
//	H_in = H_in(GR)
//	R_in = R_in(GR) = 8piG(rho_m + 4*rho_Lambda)
//	R'_in = R'_in(R_GR, H_GR) -- see differential equation for R'
inline double H_initial_fR(const double a, const double H, const double R, const double f, const double fr, const double frr_term)
{
	// return sqrt( (H*H + a*a*(fr*R - f)/6.) / (1. + fr - frr_term) );
	return H;
}

inline double R_initial_fR(const double a, const double fourpiG, const cosmology & cosmo)
{
	// TODO: Maybe not accurate enough at this stage, but let's try this first.
	return 2. * fourpiG * ( cosmo.Omega_m / a / a / a + 4. * (1. - cosmo.Omega_m - cosmo.Omega_rad) );
}

inline double Rbar_GR(const double a, const double fourpiG, const cosmology & cosmo)
{
	return 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4. * cosmo.Omega_Lambda);
}

inline double dot_R_initial_fR(const double a, const double H, const double fourpiG, const cosmology & cosmo, const metadata & sim)
{
	// Corresponds to background as computed using the modified Friedmann equation:
	// return R_dot_RungeKutta(a, H, R_initial_fR(a, fourpiG, cosmo), 2.*fourpiG*a*a*rho(a, cosmo), sim);

	// GR value:
	return -6. * fourpiG * H * cosmo.Omega_m / a / a / a; // TODO: add ncdm species
}

inline double dot_Rbar_GR(const double a, const double H, const double fourpiG, const cosmology & cosmo)
{
	return -6. * fourpiG * H * cosmo.Omega_m / a / a / a; // TODO: add ncdm species
}

///////////////////////////////////////////////////////////
// Computing the background with RK4 solver
//////////////////////////////////////////////////////////
// System for f(R) background:
// a' = a_dot_RungeKutta(...)
// H' = H_dot_RungeKutta(...)
// R' = R_dot_RungeKutta(...)

inline double a_dot_RungeKutta(const double a, const double H)
{
	return a * H;
}

inline double a_dot_RungeKutta_trace(const double a, const double H)
{
	return a * H;
}

inline double H_dot_RungeKutta(const double a, const double H, const double R)
{
	return a*a*R/6. - H*H;
}

inline double H_dot_RungeKutta_trace(const double a, const double H, const double R)
{
	return a*a*R/6. - H*H;
}

// For background_only option -- T00_background is computed from cosmological parameters Omega_m, Omega_Lambda etc
double R_dot_RungeKutta(const double a,
												const double H,
												const double R,
												const double eightpiG_rho_a2, // Actually 8 pi G * rho_background * a**2 -- TODO: Should this be T00_hom instead?
												const metadata & sim)
{
	double frr = fRR(R, sim, 210),
				 fr = fR(R, sim, 211),
				 denom = 3. * frr * H;

	if(denom)
	{
		double rdot = (eightpiG_rho_a2 - 3.*H*H*(1. + fr) + 0.5*(fr * R - f(R, sim, 212))*a*a) / denom;
		return rdot;
	}
	else
	{
		COUT << " fRR evaluates to zero. Closing..." << endl
				 << " R = " << R << endl
				 << " H = " << H << endl
				 << " fRR = " << frr << endl;
		exit(3);
	}
}

inline double R_dot_RungeKutta_trace(const double Y)
{
	return Y;
}

double Y_dot_RungeKutta_trace(const double a,
	 														const double H,
															const double R,
															const double Y,
															const double Trace_hom,
															const double fourpiG,
															const cosmology & cosmo,
															const metadata & sim)
{
	double res,
				 f0 = f(R, sim, 310),
				 fr = fR(R, sim, 311),
				 frr = fRR(R, sim, 312),
				 frrr = fRRR(R, sim, 313);

	if(frr)
	{
		res = fr * R - 2.*f0 - R - 2. * fourpiG * Trace_hom;
		res *= a * a / 3.;
		res -= frrr * Y * Y;
		res /= frr;
		res -= 2. * H * Y;
	}
	else
	{
		res = H / a;
		res *= 24. * res;
		res -= R;
		res *= fourpiG * cosmo.Omega_m / a;
	}
	return res;
}

////////////////////////
// Runge-Kutta (4) solver for f(R) background
// TODO: more info here
////////////////////////
double rungekutta_fR(double & a,
									   double & H,
									   double & R,
									   const double fourpiG,
									   const cosmology & cosmo,
									   const double dtau,
									   const metadata & sim)
{
	double a1, a2, a3, a4,
			   H1, H2, H3, H4,
				 R1, R2, R3, R4;

	a1 = a_dot_RungeKutta(a, H);
	H1 = H_dot_RungeKutta(a, H, R);
	R1 = R_dot_RungeKutta(a, H, R, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo), sim);

	a2 = a_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = H_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = R_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, 3. * Hconf(a + 0.5 * a1 * dtau, fourpiG, cosmo) * Hconf(a + 0.5 * a1 * dtau, fourpiG, cosmo), sim);

	a3 = a_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = H_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = R_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, 3. * Hconf(a + 0.5 * a2 * dtau, fourpiG, cosmo) * Hconf(a + 0.5 * a2 * dtau, fourpiG, cosmo), sim);

	a4 = a_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau);
	H4 = H_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = R_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, 3. * Hconf(a + a3 * dtau, fourpiG, cosmo) * Hconf(a + a3 * dtau, fourpiG, cosmo), sim);

	a += dtau * (a1 + 2.*a2 + 2.*a3 + a4) / 6.;
	H += dtau * (H1 + 2.*H2 + 2.*H3 + H4) / 6.;
	R += dtau * (R1 + 2.*R2 + 2.*R3 + R4) / 6.;

	return dtau;
}

////////////////////////////////////////////
// With explicit T00_hom instead of theoretical value
////////////////////////////////////////////
double rungekutta_fR(double & a,
									   double & H,
									   double & R,
									   const double fourpiG,
									   const cosmology & cosmo,
										 const double T00_hom,
									   const double dtau,
									   const metadata & sim)
{
	double a1, a2, a3, a4,
			   H1, H2, H3, H4,
				 R1, R2, R3, R4;

	a1 = a_dot_RungeKutta(a, H);
	H1 = H_dot_RungeKutta(a, H, R);
	R1 = R_dot_RungeKutta(a, H, R, 3. * Hconf(a, fourpiG, cosmo, T00_hom) * Hconf(a, fourpiG, cosmo, T00_hom), sim);

	a2 = a_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = H_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = R_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, 3. * Hconf(a + 0.5 * a1 * dtau, fourpiG, cosmo, T00_hom) * Hconf(a + 0.5 * a1 * dtau, fourpiG, cosmo, T00_hom), sim);

	a3 = a_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = H_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = R_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, 3. * Hconf(a + 0.5 * a2 * dtau, fourpiG, cosmo, T00_hom) * Hconf(a + 0.5 * a2 * dtau, fourpiG, cosmo, T00_hom), sim);

	a4 = a_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau);
	H4 = H_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = R_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, 3. * Hconf(a + a3 * dtau, fourpiG, cosmo, T00_hom) * Hconf(a + a3 * dtau, fourpiG, cosmo, T00_hom), sim);

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
									 const cosmology & cosmo,
									 const double dtau,
									 const metadata & sim)
{
	double ac, af, a1, a2, a3, a4, a5, a6,
				 Hc, Hf, H1, H2, H3, H4, H5, H6,
				 Rc, Rf, R1, R2, R3, R4, R5, R6,
				 s,
				 dtau_temp = dtau;

	a1 = dtau * a_dot_RungeKutta(a, H);
	H1 = dtau * H_dot_RungeKutta(a, H, R);
	R1 = dtau * R_dot_RungeKutta(a, H, R, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo), sim);

	a2 = dtau * a_dot_RungeKutta(a + 0.25 * a1, H + 0.25 * H1);
	H2 = dtau * H_dot_RungeKutta(a + 0.25 * a1, H + 0.25 * H1, R + 0.25 * R1);
	R2 = dtau * R_dot_RungeKutta(a + 0.25 * a1, H + 0.25 * H1, R + 0.25 * R1, 3. * Hconf(a + 0.25 * a1, fourpiG, cosmo) * Hconf(a + 0.25 * a1, fourpiG, cosmo), sim);

	a3 = dtau * a_dot_RungeKutta(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2);
	H3 = dtau * H_dot_RungeKutta(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2, R + 3./32. * R1 + 9./32. * R2);
	R3 = dtau * R_dot_RungeKutta(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2, R + 3./32. * R1 + 9./32. * R2, 3. * Hconf(a + 3./32. * a1 + 9./32. * a2, fourpiG, cosmo) * Hconf(a + 3./32. * a1 + 9./32. * a2, fourpiG, cosmo), sim);

	a4 = dtau * a_dot_RungeKutta(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3);
	H4 = dtau * H_dot_RungeKutta(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3, R + 1932./2197. * R1 - 7200./2197. * R2 + 7296./2197. * R3);
	R4 = dtau * R_dot_RungeKutta(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3, R + 1932./2197. * R1 - 7200./2197. * R2 + 7296./2197. * R3, 3. * Hconf(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, fourpiG, cosmo) * Hconf(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, fourpiG, cosmo), sim);

	a5 = dtau * a_dot_RungeKutta(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4);
	H5 = dtau * H_dot_RungeKutta(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4, R + 439./216. * R1 - 8. * R2 + 3680./513. * R3 - 845./4104. * R4);
	R5 = dtau * R_dot_RungeKutta(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4, R + 439./216. * R1 - 8. * R2 + 3680./513. * R3 - 845./4104. * R4, 3. * Hconf(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, fourpiG, cosmo) * Hconf(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, fourpiG, cosmo), sim);

	a6 = dtau * a_dot_RungeKutta(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5);
	H6 = dtau * H_dot_RungeKutta(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5, R - 8./27. * R1 + 2. * R2 - 3544./2565. * R3 + 1859./4104. * R4 - 11./40. * R5);
	R6 = dtau * R_dot_RungeKutta(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5, R - 8./27. * R1 + 2. * R2 - 3544./2565. * R3 + 1859./4104. * R4 - 11./40. * R5, 3. * Hconf(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, fourpiG, cosmo) * Hconf(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, fourpiG, cosmo), sim);

	ac = a + 25./216. * a1 + 1408./2565. * a3 + 2197./4101. * a4 - 1./5. * a5;
	Hc = H + 25./216. * H1 + 1408./2565. * H3 + 2197./4101. * H4 - 1./5. * H5;
	Rc = R + 25./216. * R1 + 1408./2565. * R3 + 2197./4101. * R4 - 1./5. * R5;

	a = af = a + 16./135. * a1 + 6656./12825. * a3 + 28561./56430. * a4 - 9./50. * a5 + 2./55. * a6;
	H = Hf = H + 16./135. * H1 + 6656./12825. * H3 + 28561./56430. * H4 - 9./50. * H5 + 2./55. * H6;
	R = Rf = R + 16./135. * R1 + 6656./12825. * R3 + 28561./56430. * R4 - 9./50. * R5 + 2./55. * R6;

	return dtau;
}

////////////////////////////////////////////
// With explicit T00_hom instead of theoretical value
// TODO Comments here
////////////////////////////////////////////
double rungekutta_fR_45(double & a,
									 double & H,
									 double & R,
									 const double fourpiG,
									 const cosmology & cosmo,
									 const double T00_hom_a3,
									 const double dtau,
									 const metadata & sim)
{
	double ac, af, a1, a2, a3, a4, a5, a6,
				 Hc, Hf, H1, H2, H3, H4, H5, H6,
				 Rc, Rf, R1, R2, R3, R4, R5, R6,
				 s,
				 dtau_temp = dtau;

	a1 = dtau * a_dot_RungeKutta(a, H);
	H1 = dtau * H_dot_RungeKutta(a, H, R);
	R1 = dtau * R_dot_RungeKutta(a, H, R, 3. * Hconf(a, fourpiG, cosmo, T00_hom_a3) * Hconf(a, fourpiG, cosmo, T00_hom_a3), sim);

	a2 = dtau * a_dot_RungeKutta(a + 0.25 * a1, H + 0.25 * H1);
	H2 = dtau * H_dot_RungeKutta(a + 0.25 * a1, H + 0.25 * H1, R + 0.25 * R1);
	R2 = dtau * R_dot_RungeKutta(a + 0.25 * a1, H + 0.25 * H1, R + 0.25 * R1, 3. * Hconf(a + 0.25 * a1, fourpiG, cosmo, T00_hom_a3) * Hconf(a + 0.25 * a1, fourpiG, cosmo, T00_hom_a3), sim);

	a3 = dtau * a_dot_RungeKutta(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2);
	H3 = dtau * H_dot_RungeKutta(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2, R + 3./32. * R1 + 9./32. * R2);
	R3 = dtau * R_dot_RungeKutta(a + 3./32. * a1 + 9./32. * a2, H + 3./32. * H1 + 9./32. * H2, R + 3./32. * R1 + 9./32. * R2, 3. * Hconf(a + 3./32. * a1 + 9./32. * a2, fourpiG, cosmo, T00_hom_a3) * Hconf(a + 3./32. * a1 + 9./32. * a2, fourpiG, cosmo, T00_hom_a3), sim);

	a4 = dtau * a_dot_RungeKutta(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3);
	H4 = dtau * H_dot_RungeKutta(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3, R + 1932./2197. * R1 - 7200./2197. * R2 + 7296./2197. * R3);
	R4 = dtau * R_dot_RungeKutta(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, H + 1932./2197. * H1 - 7200./2197. * H2 + 7296./2197. * H3, R + 1932./2197. * R1 - 7200./2197. * R2 + 7296./2197. * R3, 3. * Hconf(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, fourpiG, cosmo, T00_hom_a3) * Hconf(a + 1932./2197. * a1 - 7200./2197. * a2 + 7296./2197. * a3, fourpiG, cosmo, T00_hom_a3), sim);

	a5 = dtau * a_dot_RungeKutta(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4);
	H5 = dtau * H_dot_RungeKutta(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4, R + 439./216. * R1 - 8. * R2 + 3680./513. * R3 - 845./4104. * R4);
	R5 = dtau * R_dot_RungeKutta(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, H + 439./216. * H1 - 8. * H2 + 3680./513. * H3 - 845./4104. * H4, R + 439./216. * R1 - 8. * R2 + 3680./513. * R3 - 845./4104. * R4, 3. * Hconf(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, fourpiG, cosmo, T00_hom_a3) * Hconf(a + 439./216. * a1 - 8. * a2 + 3680./513. * a3 - 845./4104. * a4, fourpiG, cosmo, T00_hom_a3), sim);

	a6 = dtau * a_dot_RungeKutta(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5);
	H6 = dtau * H_dot_RungeKutta(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5, R - 8./27. * R1 + 2. * R2 - 3544./2565. * R3 + 1859./4104. * R4 - 11./40. * R5);
	R6 = dtau * R_dot_RungeKutta(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, H - 8./27. * H1 + 2. * H2 - 3544./2565. * H3 + 1859./4104. * H4 - 11./40. * H5, R - 8./27. * R1 + 2. * R2 - 3544./2565. * R3 + 1859./4104. * R4 - 11./40. * R5, 3. * Hconf(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, fourpiG, cosmo, T00_hom_a3) * Hconf(a - 8./27. * a1 + 2. * a2 - 3544./2565. * a3 + 1859./4104. * a4 - 11./40. * a5, fourpiG, cosmo, T00_hom_a3), sim);

	ac = a + 25./216. * a1 + 1408./2565. * a3 + 2197./4101. * a4 - 1./5. * a5;
	Hc = H + 25./216. * H1 + 1408./2565. * H3 + 2197./4101. * H4 - 1./5. * H5;
	Rc = R + 25./216. * R1 + 1408./2565. * R3 + 2197./4101. * R4 - 1./5. * R5;

	a = af = a + 16./135. * a1 + 6656./12825. * a3 + 28561./56430. * a4 - 9./50. * a5 + 2./55. * a6;
	H = Hf = H + 16./135. * H1 + 6656./12825. * H3 + 28561./56430. * H4 - 9./50. * H5 + 2./55. * H6;
	R = Rf = R + 16./135. * R1 + 6656./12825. * R3 + 28561./56430. * R4 - 9./50. * R5 + 2./55. * R6;

	return dtau;
}


////////////////////////////////////////////
// Runge-Kutta background evolution with trace equation (theoretical background)
// TODO: Details here
////////////////////////////////////////////
double rungekutta_fR_trace(double & a,
									         double & H,
									         double & R,
													 double & Y, // := dot_R
									         const double fourpiG,
									         const cosmology & cosmo,
									         const double dtau,
									         const metadata & sim)
{
	double a1, a2, a3, a4,
				 H1, H2, H3, H4,
				 R1, R2, R3, R4,
				 Y1, Y2, Y3, Y4;

	a1 = a_dot_RungeKutta_trace(a, H);
	H1 = H_dot_RungeKutta_trace(a, H, R);
	R1 = R_dot_RungeKutta_trace(Y);
	Y1 = Y_dot_RungeKutta_trace(a, H, R, Y, Tbar(a, cosmo), fourpiG, cosmo, sim);

	a2 = a_dot_RungeKutta_trace(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = H_dot_RungeKutta_trace(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = R_dot_RungeKutta_trace(Y + 0.5 * Y1 * dtau);
	Y2 = Y_dot_RungeKutta_trace(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, Y + 0.5 * Y1 * dtau, Tbar(a + 0.5 * a1 * dtau, cosmo), fourpiG, cosmo, sim);

	a3 = a_dot_RungeKutta_trace(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = H_dot_RungeKutta_trace(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = R_dot_RungeKutta_trace(Y + 0.5 * Y2 * dtau);
	Y3 = Y_dot_RungeKutta_trace(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, Y + 0.5 * Y2 * dtau, Tbar(a + 0.5 * a2 * dtau, cosmo), fourpiG, cosmo, sim);

	a4 = a_dot_RungeKutta_trace(a + a3 * dtau, H + H3 * dtau);
	H4 = H_dot_RungeKutta_trace(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = R_dot_RungeKutta_trace(Y + Y3 * dtau);
	Y4 = Y_dot_RungeKutta_trace(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, Y + Y3 * dtau, Tbar(a + a3 * dtau, cosmo), fourpiG, cosmo, sim);

	a += dtau * (a1 + 2.*a2 + 2.*a3 + a4) / 6.;
	H += dtau * (H1 + 2.*H2 + 2.*H3 + H4) / 6.;
	R += dtau * (R1 + 2.*R2 + 2.*R3 + R4) / 6.;
	Y += dtau * (Y1 + 2.*Y2 + 2.*Y3 + Y4) / 6.;

	return dtau;
}

////////////////////////////////////////////
// Runge-Kutta background evolution with trace equation (theoretical background)
// TODO: Details here
////////////////////////////////////////////
double rungekutta_fR_trace(double & a,
									         double & H,
									         double & R,
													 double & Y, // := dot_R
									         const double fourpiG,
									         const cosmology & cosmo,
													 const double T_hom,
									         const double dtau,
									         const metadata & sim)
{
	double a1, a2, a3, a4,
				 H1, H2, H3, H4,
				 R1, R2, R3, R4,
				 Y1, Y2, Y3, Y4,
				 Tbar_0 = Tbar(a, cosmo);

	a1 = a_dot_RungeKutta_trace(a, H);
	H1 = H_dot_RungeKutta_trace(a, H, R);
	R1 = R_dot_RungeKutta_trace(Y);
	Y1 = Y_dot_RungeKutta_trace(a, H, R, Y, T_hom, fourpiG, cosmo, sim);

	a2 = a_dot_RungeKutta_trace(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = H_dot_RungeKutta_trace(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = R_dot_RungeKutta_trace(Y + 0.5 * Y1 * dtau);
	Y2 = Y_dot_RungeKutta_trace(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, Y + 0.5 * Y1 * dtau, T_hom * Tbar(a + 0.5 * a1 * dtau, cosmo) / Tbar_0, fourpiG, cosmo, sim);

	a3 = a_dot_RungeKutta_trace(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = H_dot_RungeKutta_trace(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = R_dot_RungeKutta_trace(Y + 0.5 * Y2 * dtau);
	Y3 = Y_dot_RungeKutta_trace(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, Y + 0.5 * Y2 * dtau, T_hom * Tbar(a + 0.5 * a2 * dtau, cosmo) / Tbar_0, fourpiG, cosmo, sim);

	a4 = a_dot_RungeKutta_trace(a + a3 * dtau, H + H3 * dtau);
	H4 = H_dot_RungeKutta_trace(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = R_dot_RungeKutta_trace(Y + Y3 * dtau);
	Y4 = Y_dot_RungeKutta_trace(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, Y + Y3 * dtau, T_hom * Tbar(a + a3 * dtau, cosmo) / Tbar_0, fourpiG, cosmo, sim);

	a += dtau * (a1 + 2.*a2 + 2.*a3 + a4) / 6.;
	H += dtau * (H1 + 2.*H2 + 2.*H3 + H4) / 6.;
	R += dtau * (R1 + 2.*R2 + 2.*R3 + R4) / 6.;
	Y += dtau * (Y1 + 2.*Y2 + 2.*Y3 + Y4) / 6.;

	return dtau;
}



#endif
