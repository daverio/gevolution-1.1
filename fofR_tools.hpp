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

inline Real fofr_potential(double phi, double m)
{
	return m*m*phi*phi/2.0
}

#endif
