//////////////////////////
// background_fR.hpp
//////////////////////////
//
// code components related to background evolution for the f(R) model
//
// Authors:
// David Daverio (Cambridge University)
// Lorenzo Reverberi (CEICO - Czech Academy of Sciences, Prague)
//
// Note: Some methods have been copied and pasted from background.hpp
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


inline double Tbar(const double a, const cosmology & cosmo)
{
	return - cosmo.Omega_m/a/a/a - 4. * cosmo.Omega_Lambda;
}

double R_GR(const double a, const double fourpiG, const cosmology & cosmo)
{
	return 2. * fourpiG * ( cosmo.Omega_m/a/a/a + 4. * (1. - cosmo.Omega_m - cosmo.Omega_rad) );
}

inline double H_initial_fR(const double a, const double H, const double R, const double f, const double fr, const double frr_term)
{
	// return sqrt( (H*H + a*a*(fr*R - f)/6.) / (1. + fr - frr_term) );
	return H;
}

inline double R_initial_fR(const double a, const double fourpiG, const cosmology & cosmo)
{
	return 2. * fourpiG * ( cosmo.Omega_m / a / a / a + 4. * (1. - cosmo.Omega_m - cosmo.Omega_rad) );
}

inline double Rbar_GR(const double a, const double fourpiG, const cosmology & cosmo)
{
	return 2. * fourpiG * ( cosmo.Omega_m / a / a / a + 4. * cosmo.Omega_Lambda );
}

inline double dot_R_initial_fR(const double a, const double H, const double fourpiG, const cosmology & cosmo, const metadata & sim)
{
	return -6. * fourpiG * H * cosmo.Omega_m / a / a / a;
}

inline double dot_Rbar_GR(const double a, const double H, const double fourpiG, const cosmology & cosmo)
{
	return -6. * fourpiG * H * cosmo.Omega_m / a / a / a;
}

///////////////////////////////////////////
// Computing the background with RK4 solver
///////////////////////////////////////////
// System for f(R) background:
// a' = a_dot_RungeKutta(...)
// H' = H_dot_RungeKutta(...)
// R' = R_dot_RungeKutta(...)

inline double a_dot_RungeKutta(const double a, const double H)
{
	return a * H;
}

inline double H_dot_RungeKutta(const double a, const double H, const double R)
{
	return a * a * R / 6. - H * H;
}

inline double R_dot_RungeKutta(const double Y)
{
	return Y;
}

double Y_dot_RungeKutta(
	const double a,
	const double H,
	const double R,
	const double Y,
	const double Trace_hom,
	const double fourpiG,
	const cosmology & cosmo,
	const metadata & sim
)
{
	double
	result,
	f0 = f(R, sim, 310),
	fr = fR(R, sim, 311),
	frr = fRR(R, sim, 312),
	frrr = fRRR(R, sim, 313);

	if(frr && false)
	{
		result = 2. * f0 + (1. - fr) * R + 2. * fourpiG * Trace_hom;
		result *= - a * a / 3.;
		result -= frrr * Y * Y;
		result /= frr;
		result -= 2. * H * Y;
	}
	else
	{
		result = H / a;
		result *= 24. * result;
		result -= R;
		result *= fourpiG * cosmo.Omega_m / a;
	}

	return result;
}

///////////////////////////////////////////////////////
// Runge-Kutta background evolution with trace equation
///////////////////////////////////////////////////////
double rungekutta_fR(
	double & a,
	double & H,
	double & R,
	double & Y, // := dot_R
	const double fourpiG,
	const double dtau,
	const double Trace_hom,
	const cosmology & cosmo,
	const metadata & sim
)
{
	double
	trace_new,
	a1, a2, a3, a4,
	H1, H2, H3, H4,
	R1, R2, R3, R4,
	Y1, Y2, Y3, Y4,
	Tbar_0 = Tbar(a, cosmo);

	a1 = a_dot_RungeKutta(a, H);
	H1 = H_dot_RungeKutta(a, H, R);
	R1 = R_dot_RungeKutta(Y);
	Y1 = Y_dot_RungeKutta(a, H, R, Y, Trace_hom, fourpiG, cosmo, sim);

	trace_new = Trace_hom * Tbar(a + 0.5 * a1 * dtau, cosmo) / Tbar_0;
	a2 = a_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau);
	H2 = H_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau);
	R2 = R_dot_RungeKutta(Y + 0.5 * Y1 * dtau);
	Y2 = Y_dot_RungeKutta(a + 0.5 * a1 * dtau, H + 0.5 * H1 * dtau, R + 0.5 * R1 * dtau, Y + 0.5 * Y1 * dtau, trace_new, fourpiG, cosmo, sim);

	trace_new =  Trace_hom * Tbar(a + 0.5 * a2 * dtau, cosmo) / Tbar_0;
	a3 = a_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau);
	H3 = H_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau);
	R3 = R_dot_RungeKutta(Y + 0.5 * Y2 * dtau);
	Y3 = Y_dot_RungeKutta(a + 0.5 * a2 * dtau, H + 0.5 * H2 * dtau, R + 0.5 * R2 * dtau, Y + 0.5 * Y2 * dtau, trace_new, fourpiG, cosmo, sim);

	trace_new =  Trace_hom * Tbar(a + a3 * dtau, cosmo) / Tbar_0;
	a4 = a_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau);
	H4 = H_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau);
	R4 = R_dot_RungeKutta(Y + Y3 * dtau);
	Y4 = Y_dot_RungeKutta(a + a3 * dtau, H + H3 * dtau, R + R3 * dtau, Y + Y3 * dtau, trace_new, fourpiG, cosmo, sim);

	a += dtau * (a1 + 2.*a2 + 2.*a3 + a4) / 6.;
	H += dtau * (H1 + 2.*H2 + 2.*H3 + H4) / 6.;
	R += dtau * (R1 + 2.*R2 + 2.*R3 + R4) / 6.;
	Y += dtau * (Y1 + 2.*Y2 + 2.*Y3 + Y4) / 6.;

	return dtau;
}


///////////////////////////
// Computing the background
///////////////////////////
void rungekutta_background(
	double & a,
	double & Hubble,
	double & Rbar,
	double & dot_Rbar,
	const double fourpiG,
	const double dtau,
	const double Trace_hom,
	const cosmology & cosmo,
	const metadata & sim,
	const double dtau_bg,
	const int numsteps_bg
)
{
	// GR or Newtonian evolution
	if(sim.modified_gravity_flag != MODIFIED_GRAVITY_FLAG_FR || sim.lcdm_background)
	{
		rungekutta4bg(a, fourpiG, cosmo, dtau);
		Hubble = Hconf(a, fourpiG, cosmo);
		Rbar = Rbar_GR(a, fourpiG, cosmo);
		dot_Rbar = dot_Rbar_GR(a, Hubble, fourpiG, cosmo);
	}
	else // f(R) gravity, non-LCDM background
	{
		if(numsteps_bg == 1)
		{
			rungekutta_fR(a, Hubble, Rbar, dot_Rbar, fourpiG, dtau, Trace_hom, cosmo, sim);
		}
		else
		{
			for(int g=0; g<numsteps_bg; g++) // TODO Rescale for numpsteps_ncdm[i]?
			{
				rungekutta_fR(a, Hubble, Rbar, dot_Rbar, fourpiG, dtau_bg, Trace_hom, cosmo, sim);
			}
		}
	}

	return;
}



#endif
