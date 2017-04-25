//////////////////////////
// fofR_tools.hpp
//////////////////////////
//
// code components related to f(R)
//
// Author: David Daverio (Cambridge University)
//         Lorenzo Reverberi (Cape Town University)
//
// Last modified: February 2016
//
//////////////////////////

#ifndef FOFR_TOOLS_HEADER
#define FOFR_TOOLS_HEADER

// Model: f(R) = R + F(R)
// Hu-Sawicki Model: arXiv:0705.1158 -- See paper for details
// m2 := m^2 (comparison between these values and the nomenclature of the original paper)
// F(R) = -m2 * c1 * (R/M)^n / (1 + c2*(R/M)^n)
// In settings.ini:
// params[0] = m^2 is given in units of (matter energy density today) * 8piG / 3.
// params[1] = c2 (normally >> 1)
// params[2] = n
// Constructed parameters:
// params[3] = c1, calibrated to give c1*m2/c2 = 2 * 8piG * rho_Lambda [see 0705.1158 equations (25,26)]

void rescale_params_Hu_Sawicki(const cosmology cosmo, const double fourpiG, metadata * sim)
{
	sim->fofR_params[3] *= sim->fofR_params[1] * 6. * (1 - cosmo.Omega_m - cosmo.Omega_rad) / cosmo.Omega_m; //Omega_m > 0, checked in parseMetadata()
	sim->fofR_params[0] *= fourpiG * cosmo.Omega_m / 1.5;
}


///////////////////// F(R) FUNCTIONS /////////////////////

// F(R) = f(R) - R
Real F(double R, double * params, const int fofR_type)
{
	double output;

	switch(fofR_type)
	{
		case(FOFR_TYPE_RN):
		{
			output = params[0]*pow(R,params[1]);
			break;
		}

		case(FOFR_TYPE_HU_SAWICKI):
		{
			double m2, c1, c2, n;
			m2 = params[0];
			c2 = params[1];
			n  = params[2];
			c1 = params[3];
			double Rpow = pow(R/m2, n);
			output = - m2 * c1 * Rpow / (1. + c2 * Rpow);
			break;
		}

		default:
			exit(2);
	}

	return output;
}



// F_R(R) -- First derivative with respect to R
Real FR(double R, double * params, const int fofR_type)
{
	double output;

	switch(fofR_type)
	{
		case(FOFR_TYPE_RN):
		{
			output = params[0]*params[1]*pow(R,params[1]-1.);
			break;
		}

		case(FOFR_TYPE_HU_SAWICKI):
		{
			double m2, c1, c2, n;
			m2 = params[0];
			c2 = params[1];
			n  = params[2];
			c1 = params[3];
			double Rpow = pow(R/m2, n);
			output = - c1 * n * Rpow * m2 / R / ( (1. + c2*Rpow) * (1. + c2*Rpow) );
			break;
		}

		default:
			exit(2);
	}

	return output;
}


// F_{RR}(R) -- Second derivative with respect to R
Real FRR(double R, double * params, const int fofR_type, int n = -1)
{
	double output;

	switch(fofR_type)
	{
		case(FOFR_TYPE_RN):
		{
			output = params[0] * params[1] * (params[1]-1.) * pow(R,params[1]-2.);
			break;
		}

		case(FOFR_TYPE_HU_SAWICKI):
		{
			double m2, c1, c2, n;
			m2 = params[0];
			c2 = params[1];
			n  = params[2];
			c1 = params[3];
			double Rpow = pow(R/m2, n);
			output = c1 * m2 * n * Rpow * (1. - n + c2*(1. + n)*Rpow ) / ( R * R * (1. + c2*Rpow) * (1. + c2*Rpow) * (1. + c2*Rpow) );
			break;
		}

		default:
			exit(2);
	}

	if(!output)
	{
		COUT << " Something went wrong when computing FRR, with code " << n << ". Closing...";
		exit(-11);
	}

	return output;
}






#endif
