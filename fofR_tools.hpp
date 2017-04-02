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
// Hu-Sawicki Model: arXiv:0705.1158
// M := m^2, B := c1, D := c2 (comparison between these values and the nomenclature of the original paper)
// F(R) = -M * B * (R/M)**N / (1 + D*(R/M)**N)


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
			B = D * params[2];
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
			B = D * params[2];
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
			B = D * params[2];
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
