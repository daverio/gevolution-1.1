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
// M := m^2, B := c1, D := c2 (comparison between these values and the nomenclature of the original paper)
// F(R) = -M * B * (R/M)**N / (1 + D*(R/M)**N)
// params[0] is given in settings.ini in units of (matter energy density) * 8piG / 3.
void rescale_params_Hu_Sawicki(const cosmology cosmo, const double fourpiG, double (&params)[MAX_FOFR_PARAMS])
{
	params[0] *= fourpiG * (cosmo.Omega_m + cosmo.Omega_rad) / 1.5;
	params[2] *= params[1] * 6. * (1 - cosmo.Omega_m - cosmo.Omega_rad) / (cosmo.Omega_m + cosmo.Omega_rad); //Omega_m + Omega_rad > 0, it is checked in parseMetadata()
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
			double M, B, D, N;
			M = params[0];
			D = params[1];
			B = params[2];
			N = params[3];
			output = - M * B * pow(1.*R/M, N) / (1. + D * pow(R/M, N));
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
			double M, B, D, N;
			M = params[0];
			D = params[1];
			B = params[2];
			N = params[3];
			output = - B * N * pow(R/M, N-1.) / ( (1. + D*pow(R/M, N)) * (1. + D*pow(R/M, N)) );
			break;
		}

		default:
			exit(2);
	}

	return output;
}


// F_{RR}(R) -- Second derivative with respect to R
Real FRR(double R, double * params, const int fofR_type)
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
			double M, B, D, N;
			M = params[0];
			D = params[1];
			B = params[2];
			N = params[3];
			output = B * M * N * pow(R/M, N) * (1. - N + D*(1. + N)*pow(R/M, N) ) / ( R * R * (1. + D*pow(R/M, N)) * (1. + D*pow(R/M, N)) * (1. + D*pow(R/M, N)) );
			break;
		}

		default:
			exit(2);
	}

	return output;
}




#endif
