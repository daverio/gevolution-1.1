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

Real fofr(double R, double * params, const double fofR_type)
{
	double output;
	if(fofR_type==FOFR_TYPE_RN)
	{
		output = params[0]*pow(R,params[0]);
	}
	else exit(2);

	return output;
}

Real fofr_p(double R, double * params, const double fofR_type)
{
	double output;
	if(fofR_type==FOFR_TYPE_RN)
	{
		if(params[1]==0) output;
		else output = params[0]*params[1]*pow(R,params[1]-1.0);
	}
	else exit(2);

	return output;
}

Real fofr_pp(double R, double * params, const double fofR_type)
{
	double output;
	if(fofR_type==FOFR_TYPE_RN)
	{
		if(params[1]==1) output;
		else output = params[0]*params[1]*params[1]*pow(R,params[1]-2.0);
	}
	else exit(2);

	return output;
}

Real fofr_ppp(double R, double * params, const double fofR_type)
{
	double output;
	if(fofR_type==FOFR_TYPE_RN)
	{
		if(params[1]==2) output;
		else output = params[0]*params[1]*params[1]*params[1]*pow(R,params[1]-3.0);
	}
	else exit(2);

	return output;
}

#endif
