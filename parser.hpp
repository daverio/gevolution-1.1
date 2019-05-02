//////////////////////////
// parser.hpp
//////////////////////////
//
// Parser for settings file
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: November 2016
//
//////////////////////////

#ifndef PARSER_HEADER
#define PARSER_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "metadata.hpp"

using namespace std;

struct parameter
{
	char name[PARAM_MAX_LENGTH];
	char value[PARAM_MAX_LENGTH];
	bool used;
};

int sort_descending(const void * z1, const void * z2)
{
	return (* (double *) z1 > * (double *) z2) ? -1 : ((* (double *) z1 < * (double *) z2) ? 1 : 0);
}

//////////////////////////
// readline
//////////////////////////
// Description:
//   reads a line of characters and checks if it declares a parameter; if yes,
//   i.e. the line has the format "<parameter name> = <parameter value>" (with
//   an optional comment added, preceded by a hash-symbol, '#'), the parameter
//   name and value are copied to the corresponding arrays, and 'true' is returned.
//   If the format is not recognized (or the line is commented using the hash-symbol)
//   'false' is returned instead.
//
// Arguments:
//   line       string containing the line to be read
//   pname      will contain the name of the declared parameter (if found)
//   pvalue     will contain the value of the declared parameter (if found)
//
// Returns:
//   'true' if a parameter is declared in the line, 'false' otherwise.
//
//////////////////////////

bool readline(char * line, char * pname, char * pvalue)
{
	char * pequal;
	char * phash;
	char * l;
	char * r;

	pequal = strchr(line, '=');

	if(pequal == NULL || pequal == line) return false;

	phash = strchr(line, '#');

	if(phash != NULL && phash < pequal) return false;

	l = line;
	while(*l == ' ' || *l == '\t') l++;

	r = pequal-1;
	while((*r == ' ' || *r == '\t') && r > line) r--;

	if(r < l) return false;

	if(r-l+1 >= PARAM_MAX_LENGTH) return false;

	strncpy(pname, l, r-l+1);
	pname[r-l+1] = '\0';

	l = pequal+1;
	while(*l == ' ' || *l == '\t') l++;

	if(phash == NULL)
		r = line+strlen(line)-1;
	else
		r = phash-1;

	while(*r == ' ' || *r == '\t' || *r == '\n' || *r == '\r') r--;

	if(r < l) return false;

	if(r-l+1 >= PARAM_MAX_LENGTH) return false;

	strncpy(pvalue, l, r-l+1);
	pvalue[r-l+1] = '\0';

	return true;
}


//////////////////////////
// loadParameterFile
//////////////////////////
// Description:
//   loads a parameter file and creates an array of parameters declared therein
//
// Arguments:
//   filename   string containing the path to the parameter file
//   params     will contain the array of parameters (memory will be allocated)
//
// Returns:
//   number of parameters defined in the parameter file (= length of parameter array)
//
//////////////////////////

int loadParameterFile(const char * filename, parameter * & params)
{
	int numparam = 0;
	int i = 0;

	if(parallel.grid_rank()[0] == 0) // read file
	{
		FILE * paramfile;
		char line[PARAM_MAX_LINESIZE];
		char pname[PARAM_MAX_LENGTH];
		char pvalue[PARAM_MAX_LENGTH];

		paramfile = fopen(filename, "r");

		if(paramfile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Unable to open parameter file " << filename << "." << endl;
			parallel.abortForce();
		}

		while(!feof(paramfile) && !ferror(paramfile))
		{
			fgets(line, PARAM_MAX_LINESIZE, paramfile);
			if(readline(line, pname, pvalue) == true)
			{
				numparam++;
			}
		}

		if(numparam == 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! No valid data found in file " << filename << "." << endl;
			fclose(paramfile);
			parallel.abortForce();
		}

		params = (parameter *) malloc(sizeof(parameter) * numparam);

		if(params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			fclose(paramfile);
			parallel.abortForce();
		}

		rewind(paramfile);

		while(!feof(paramfile) && !ferror(paramfile) && i<numparam)
		{
			fgets(line, PARAM_MAX_LINESIZE, paramfile);
			if(readline(line, params[i].name, params[i].value) == true)
			{
				params[i].used = false;
				i++;
			}
		}

		fclose(paramfile);

		if(i<numparam)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! File may have changed or file pointer corrupted." << endl;
			free(params);
			parallel.abortForce();
		}

		parallel.broadcast_dim0<int>(numparam, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numparam, 0);

		params = (parameter *) malloc(sizeof(parameter) * numparam);

		if(params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			parallel.abortForce();
		}
	}
	parallel.broadcast_dim0<parameter>(params, numparam, 0);

	return numparam;
}


//////////////////////////
// saveParameterFile
//////////////////////////
// Description:
//   saves a parameter file
//
// Arguments:
//   filename   string containing the path to the parameter file
//   params     array of parameters
//   numparam   length of parameter array
//   used_only  if 'true', only the used parameters will be written (default)
//
// Returns:
//
//////////////////////////

void saveParameterFile(const char * filename, parameter * params, const int numparam, bool used_only = true)
{
	if(parallel.isRoot())
	{
		FILE * paramfile;

		paramfile = fopen(filename, "w");

		if(paramfile == NULL)
		{
			cout << " error in saveParameterFile! Unable to open file " << filename << "." << endl;
			parallel.abortForce();
		}
		else
		{
			for(int i=0; i<numparam; i++)
			{
				if(!used_only || params[i].used)
				{
					fprintf(paramfile, "%s = %s\n", params[i].name, params[i].value);
				}
			}

			fclose(paramfile);
		}
	}
}


//////////////////////////
// parseParameter (int)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, int & pvalue)
{
	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			if(sscanf(params[i].value, "%d", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}

	return false;
}


//////////////////////////
// parseParameter (long)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, long & pvalue)
{
	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			if(sscanf(params[i].value, "%ld", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}

	return false;
}


//////////////////////////
// parseParameter (double)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as double
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to double which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, double & pvalue)
{
	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			if(sscanf(params[i].value, "%lf", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}

	return false;
}


//////////////////////////
// parseParameter (char *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and retrieves its value as string
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     character string which will contain a copy of the parameter value (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, char * pvalue)
{
	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			strcpy(pvalue, params[i].value);
			params[i].used = true;
			return true;
		}
	}

	return false;
}


//////////////////////////
// parseParameter (double *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated double values
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of double values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, double * pvalue, int & nmax)
{
	char * start;
	char * comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;

	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if(nmax > 1)
			{
				while((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma-start);
					item[comma-start] = '\0';
					if(sscanf(item, " %lf ", pvalue+n) != 1)
					{
						nmax = n;
						return false;
					}
					start = comma+1;
					if(++n > nmax-2) break;
				}
			}
			if(sscanf(start, " %lf ", pvalue+n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}

	nmax = 0;
	return false;
}


//////////////////////////
// parseParameter (int *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated integer values
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of integer values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, int * pvalue, int & nmax)
{
	char * start;
	char * comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;

	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if(nmax > 1)
			{
				while((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma-start);
					item[comma-start] = '\0';
					if(sscanf(item, " %d ", pvalue+n) != 1)
					{
						nmax = n;
						return false;
					}
					start = comma+1;
					if(++n > nmax-2) break;
				}
			}
			if(sscanf(start, " %d ", pvalue+n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}

	nmax = 0;
	return false;
}


//////////////////////////
// parseParameter (char **)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated strings
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of character strings which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, char ** pvalue, int & nmax)
{
	char * start;
	char * comma;
	char * l;
	char * r;
	int n = 0;

	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if(nmax > 1)
			{
				while((comma = strchr(start, ',')) != NULL)
				{
					l = start;
					while(*l == ' ' || *l == '\t') l++;
					r = comma-1;
					while((*r == ' ' || *r == '\t') && r > start) r--;

					if(r < l)
					{
						nmax = n;
						return false;
					}

					strncpy(pvalue[n], l, r-l+1);
					pvalue[n][r-l+1] = '\0';

					start = comma+1;
					if(++n > nmax-2) break;
				}
			}
			l = start;
			while(*l == ' ' || *l == '\t') l++;
			r = l;
			while(*r != ' ' && *r != '\t' && *r != '\0') r++;
			r--;

			if(r < l)
			{
				nmax = n;
				return false;
			}

			strncpy(pvalue[n], l, r-l+1);
			pvalue[n][r-l+1] = '\0';

			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}

	nmax = 0;
	return false;
}


//////////////////////////
// parseFieldSpecifiers
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated field specifiers
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     integer which will contain the binary-encoded list of parsed specifiers (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseFieldSpecifiers(parameter * & params, const int numparam, const char * pname, int & pvalue)
{
	char * start;
	char * comma;
	int pos;
	char item[PARAM_MAX_LENGTH];

	for(int i=0; i<numparam; i++)
	{
		if(strcmp(params[i].name, pname) == 0)
		{
			pvalue = 0;
			start = params[i].value;
			while((comma = strchr(start, ',')) != NULL)
			{
				strncpy(item, start, comma-start);
				for(pos = comma-start; pos > 0; pos--)
				{
					if(item[pos-1] != ' ' && item[pos-1] != '\t') break;
				}
				item[pos] = '\0';
				if(strcmp(item, "Phi") == 0 || strcmp(item, "phi") == 0) pvalue |= MASK_PHI;
				else if( (strcmp(item, "xi") == 0 || strcmp(item, "Xi") == 0) ) pvalue |= MASK_XI;
				else if( (strcmp(item, "zeta") == 0 || strcmp(item, "Zeta") == 0) ) pvalue |= MASK_ZETA;
				else if( (strcmp(item, "deltaR") == 0 || strcmp(item, "DeltaR") == 0 || strcmp(item, "delta_R") == 0 || strcmp(item, "Delta_R") == 0) ) pvalue |= MASK_DELTAR;
				else if( (strcmp(item, "deltaT") == 0 || strcmp(item, "DeltaT") == 0 || strcmp(item, "delta_T") == 0 || strcmp(item, "Delta_T") == 0) ) pvalue |= MASK_DELTAT;
				else if(strcmp(item, "Chi") == 0 || strcmp(item, "chi") == 0) pvalue |= MASK_CHI;
				else if(strcmp(item, "Pot") == 0 || strcmp(item, "pot") == 0 || strcmp(item, "Psi_N") == 0 || strcmp(item, "psi_N") == 0 || strcmp(item, "PsiN") == 0 || strcmp(item, "psiN") == 0) pvalue |= MASK_POT;
				else if(strcmp(item, "B") == 0 || strcmp(item, "Bi") == 0) pvalue |= MASK_B;
				else if(strcmp(item, "P") == 0 || strcmp(item, "p") == 0 || strcmp(item, "v") == 0) pvalue |= MASK_P;
				else if(strcmp(item, "T00") == 0 || strcmp(item, "rho") == 0) pvalue |= MASK_T00;
				else if(strcmp(item, "Tij") == 0) pvalue |= MASK_TIJ;
				else if(strcmp(item, "rho_N") == 0 || strcmp(item, "rhoN") == 0) pvalue |= MASK_RBARE;
				else if(strcmp(item, "hij") == 0 || strcmp(item, "GW") == 0) pvalue |= MASK_HIJ;
				else if(strcmp(item, "Gadget") == 0 || strcmp(item, "Gadget2") == 0 || strcmp(item, "gadget") == 0 || strcmp(item, "gadget2") == 0) pvalue |= MASK_GADGET;
				else if(strcmp(item, "Particles") == 0 || strcmp(item, "particles") == 0 || strcmp(item, "pcls") == 0 || strcmp(item, "part") == 0) pvalue |= MASK_PCLS;
				else if(strcmp(item, "cross") == 0 || strcmp(item, "X-spectra") == 0 || strcmp(item, "x-spectra") == 0) pvalue |= MASK_XSPEC;
				else if(strcmp(item, "delta") == 0 || strcmp(item, "Ds") == 0 || strcmp(item, "D_s") == 0) pvalue |= MASK_DELTA;
				else if(strcmp(item, "delta_N") == 0 || strcmp(item, "deltaN") == 0) pvalue |= MASK_DBARE;
				else if(strcmp(item, "laplace_xi") == 0 || strcmp(item, "laplace xi") == 0) pvalue |= MASK_LAPLACE_XI;
				else if(strcmp(item, "phi_effective") == 0 || strcmp(item, "Phi_effective") == 0) pvalue |= MASK_PHI_EFFECTIVE;
				else if(strcmp(item, "phi_ddot") == 0 || strcmp(item, "Phi_ddot") == 0) pvalue |= MASK_PHI_DDOT;

				start = comma+1;
				while(*start == ' ' || *start == '\t') start++;
			}

			if(strcmp(start, "Phi") == 0 || strcmp(start, "phi") == 0) pvalue |= MASK_PHI;
			else if( (strcmp(start, "xi") == 0 || strcmp(start, "xi") == 0) ) pvalue |= MASK_XI;
			else if( (strcmp(start, "zeta") == 0 || strcmp(start, "Zeta") == 0) ) pvalue |= MASK_ZETA;
			else if( (strcmp(start, "deltaR") == 0 || strcmp(start, "DeltaR") == 0 || strcmp(start, "delta_R") == 0 || strcmp(start, "Delta_R") == 0) ) pvalue |= MASK_DELTAR;
			else if( (strcmp(start, "deltaT") == 0 || strcmp(start, "DeltaT") == 0 || strcmp(start, "delta_T") == 0 || strcmp(start, "Delta_T") == 0) ) pvalue |= MASK_DELTAT;
			else if(strcmp(start, "Chi") == 0 || strcmp(start, "chi") == 0) pvalue |= MASK_CHI;
			else if(strcmp(start, "Pot") == 0 || strcmp(start, "pot") == 0 || strcmp(start, "Psi_N") == 0 || strcmp(start, "psi_N") == 0 || strcmp(start, "PsiN") == 0 || strcmp(start, "psiN") == 0) pvalue |= MASK_POT;
			else if(strcmp(start, "B") == 0 || strcmp(start, "Bi") == 0) pvalue |= MASK_B;
			else if(strcmp(start, "P") == 0 || strcmp(start, "p") == 0 || strcmp(start, "v") == 0) pvalue |= MASK_P;
			else if(strcmp(start, "T00") == 0 || strcmp(start, "rho") == 0) pvalue |= MASK_T00;
			else if(strcmp(start, "Tij") == 0) pvalue |= MASK_TIJ;
			else if(strcmp(start, "rho_N") == 0 || strcmp(start, "rhoN") == 0) pvalue |= MASK_RBARE;
			else if(strcmp(start, "hij") == 0 || strcmp(start, "GW") == 0) pvalue |= MASK_HIJ;
			else if(strcmp(start, "Gadget") == 0 || strcmp(start, "Gadget2") == 0 || strcmp(start, "gadget") == 0 || strcmp(start, "gadget2") == 0) pvalue |= MASK_GADGET;
			else if(strcmp(start, "Particles") == 0 || strcmp(start, "particles") == 0 || strcmp(start, "pcls") == 0 || strcmp(start, "part") == 0) pvalue |= MASK_PCLS;
			else if(strcmp(start, "cross") == 0 || strcmp(start, "X-spectra") == 0 || strcmp(start, "x-spectra") == 0) pvalue |= MASK_XSPEC;
			else if(strcmp(start, "delta") == 0 || strcmp(start, "Ds") == 0 || strcmp(start, "D_s") == 0) pvalue |= MASK_DELTA;
			else if(strcmp(start, "delta_N") == 0 || strcmp(start, "deltaN") == 0) pvalue |= MASK_DBARE;
			else if(strcmp(start, "laplace_xi") == 0 || strcmp(start, "laplace xi") == 0) pvalue |= MASK_LAPLACE_XI;
			else if(strcmp(start, "phi_effective") == 0 || strcmp(start, "Phi_effective") == 0) pvalue |= MASK_PHI_EFFECTIVE;
			else if(strcmp(start, "phi_ddot") == 0 || strcmp(start, "Phi_ddot") == 0) pvalue |= MASK_PHI_DDOT;

			params[i].used = true;
			return true;
		}
	}

	return false;
}


//////////////////////////
// parseMetadata
//////////////////////////
// Description:
//   parses all metadata from the parameter array
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   sim        reference to metadata stucture (holds simulation parameters)
//   cosmo      reference to cosmology structure (holds cosmological parameters)
//   ic         reference to icsettings structure (holds settings for IC generation)
//
// Returns:
//   number of parameters parsed
//
//////////////////////////

int parseMetadata(parameter * & params, const int numparam, metadata & sim, cosmology & cosmo, icsettings & ic)
{
	char par_string[PARAM_MAX_LENGTH];
	char * pptr[MAX_PCL_SPECIES + METRICFILE_LENGTH];
	int usedparams = 0;
	int i;

	// parse settings for IC generator

	ic.pkfile[0] = '\0';
	ic.tkfile[0] = '\0';
	for(i=0; i<METRICFILE_LENGTH; i++)
	{
		ic.metricfile[i][0] = '\0';
	}
	ic.seed = 0;
	ic.flags = 0;
	ic.z_ic = -2.;
	ic.z_relax = -2.;
	ic.Cf = -1.0;
	ic.A_s = P_SPECTRAL_AMP;
	ic.n_s = P_SPECTRAL_INDEX;
	ic.k_pivot = P_PIVOT_SCALE;
	ic.restart_cycle = -1;
	ic.restart_tau = 0.;
	ic.restart_dtau = 0.;
	ic.restart_version = -1.;

	ic.restart_dtau_old = -1.;
	ic.restart_dtau_old_2 = -1.;
	ic.restart_dtau_osci = -1.;
	ic.restart_dtau_bg = -1.;
	ic.restart_a = -1.;
	ic.restart_Hubble = -1.;
	ic.restart_Rbar = -1.;
	ic.restart_dot_Rbar = -1.;

	// parse metadata
	sim.numpts = 0;
	sim.downgrade_factor = 1;
	for(i=0; i<MAX_PCL_SPECIES; i++)
	{
		sim.numpcl[i] = 0;
	}
	sim.vector_flag = VECTOR_PARABOLIC;
	sim.relativistic_flag = 0;
	sim.modified_gravity_flag = 0;
	sim.out_pk = 0;
	sim.out_snapshot = 0;
	sim.out_check = 0;
	sim.num_pk = MAX_OUTPUTS;
	sim.numbins = 0;
	sim.num_snapshot = MAX_OUTPUTS;
	sim.num_restart = MAX_OUTPUTS;
	for(i=0; i<MAX_PCL_SPECIES; i++) sim.tracer_factor[i] = 1;
	sim.Cf = 1.;
	sim.steplimit = 1.;
	sim.boxsize = -1.;
	sim.wallclocklimit = -1.;
	sim.z_in = 0.;
	sim.fR_model = 0;
	for(i=0; i<MAX_FR_PARAMS; i++) sim.fR_params[i] = 0;
	sim.num_fR_params = MAX_FR_PARAMS;

	// parse cosmological parameters
	if(!parseParameter(params, numparam, "h", cosmo.h))
	{
		cosmo.h = P_HUBBLE;
	}

	cosmo.num_ncdm = MAX_PCL_SPECIES-2;
	if(!parseParameter(params, numparam, "m_ncdm", cosmo.m_ncdm, cosmo.num_ncdm))
	{
		for(i=0; i<MAX_PCL_SPECIES-2; i++)
		{
			cosmo.m_ncdm[i] = 0.;
		}
		cosmo.num_ncdm = 0;
	}

	if(parseParameter(params, numparam, "N_ncdm", i))
	{
		if(i<0 || !isfinite(i))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of ncdm species not set properly!" << endl;
			parallel.abortForce();
		}
		if(i>cosmo.num_ncdm)
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": N_ncdm = " << i << " is larger than the number of mass parameters specified (" << cosmo.num_ncdm << ")!" << endl;
			parallel.abortForce();
		}
		cosmo.num_ncdm = i;
	}
	else if(cosmo.num_ncdm > 0)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": N_ncdm not specified, inferring from number of mass parameters in m_ncdm (" << cosmo.num_ncdm << ")!" << endl;
	}

	for(i=0; i<MAX_PCL_SPECIES-2; i++)
	{
		cosmo.T_ncdm[i] = P_T_NCDM;
		cosmo.deg_ncdm[i] = 1.0;
	}

	parseParameter(params, numparam, "T_ncdm", cosmo.T_ncdm, i);
	i = MAX_PCL_SPECIES-2;
	parseParameter(params, numparam, "deg_ncdm", cosmo.deg_ncdm, i);

	for(i=0; i<cosmo.num_ncdm; i++)
	{
		cosmo.Omega_ncdm[i] = cosmo.m_ncdm[i] * cosmo.deg_ncdm[i] / P_NCDM_MASS_OMEGA / cosmo.h / cosmo.h;
	}

	if(parseParameter(params, numparam, "T_cmb", cosmo.Omega_g))
	{
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g / cosmo.h;
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g * C_PLANCK_LAW; // Planck's law
	}
	else if(parseParameter(params, numparam, "omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g /= cosmo.h * cosmo.h;
	}
	else if(!parseParameter(params, numparam, "Omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g = 0.;
	}

	if(parseParameter(params, numparam, "N_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	else if(parseParameter(params, numparam, "N_eff", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	else if(parseParameter(params, numparam, "omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur /= cosmo.h * cosmo.h;
	}
	else if(!parseParameter(params, numparam, "Omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur = P_N_UR * (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}

	cosmo.Omega_rad = cosmo.Omega_g + cosmo.Omega_ur;
	if(parseParameter(params, numparam, "omega_b", cosmo.Omega_b))
	{
		cosmo.Omega_b /= cosmo.h * cosmo.h;
	}
	else if(!parseParameter(params, numparam, "Omega_b", cosmo.Omega_b))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": Omega_b not found in settings file, setting to default (0)." << endl;
		cosmo.Omega_b = 0.;
	}

	if(parseParameter(params, numparam, "omega_cdm", cosmo.Omega_cdm))
	{
		cosmo.Omega_cdm /= cosmo.h * cosmo.h;
	}
	else if(!parseParameter(params, numparam, "Omega_cdm", cosmo.Omega_cdm))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": Omega_cdm not found in settings file, setting to default (1)." << endl;
		cosmo.Omega_cdm = 1.;
	}

	cosmo.Omega_m = cosmo.Omega_cdm + cosmo.Omega_b;

	for(i=0; i<cosmo.num_ncdm; i++)
	{
		cosmo.Omega_m += cosmo.Omega_ncdm[i];
	}

	if(cosmo.Omega_m <= 0. || cosmo.Omega_m > 1.)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": total matter density out of range!" << endl;
		COUT << " cosmo.Omega_m = " << cosmo.Omega_m << endl;
		parallel.abortForce();
	}
	else if(cosmo.Omega_rad < 0. || cosmo.Omega_rad > 1. - cosmo.Omega_m)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": total radiation energy density out of range!" << endl;
		COUT << " cosmo.Omega_rad = " << cosmo.Omega_rad << endl;
		COUT << " 1 - cosmo.Omega_m = " << 1 - cosmo.Omega_m << endl;
		parallel.abortForce();
	}

	cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad;

	if(!parseParameter(params, numparam, "initial redshift", sim.z_in))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": initial redshift not specified!" << endl;
		parallel.abortForce();
	}

	parseParameter(params, numparam, "boxsize", sim.boxsize);
	if(sim.boxsize <= 0. || !isfinite(sim.boxsize))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": simulation box size not set properly!" << endl;
		parallel.abortForce();
	}

	// Parse Gravity Theory
	if(parseParameter(params, numparam, "gravity theory", par_string))
	{
		if(par_string[0] == 'N' || par_string[0] == 'n')
		{
			COUT << " gravity theory set to: " << COLORTEXT_CYAN << "Newtonian" << COLORTEXT_RESET << endl;
			sim.relativistic_flag = 0;
			if(ic.pkfile[0] == '\0' && ic.tkfile[0] != '\0'
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
				)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gauge transformation to N-body gauge can only be performed for the positions; the transformation for" << endl;
				COUT << "              the velocities requires time derivatives of transfer functions. Call CLASS directly to avoid this issue." << endl;
			}
		}
		else if(par_string[0] == 'G' || par_string[0] == 'g')
		{
			COUT << " Gravity theory set to: " << COLORTEXT_CYAN << "General Relativity" << COLORTEXT_RESET << endl;
			sim.relativistic_flag = 1;
		}
		else if(par_string[0] == 'F' || par_string[0] == 'f')
		{
			sim.modified_gravity_flag = MODIFIED_GRAVITY_FLAG_FR;

			if(!parseParameter(params, numparam, "lcdm background", sim.lcdm_background))
			{
				sim.lcdm_background = 0;
			}

			//TODO: read f(R) params and model
			if(!parseParameter(params, numparam, "f(R) parameters", sim.fR_params, sim.num_fR_params))
			{
				COUT << " /!\\ No f(R) parameters specified. Closing..." << endl;
				parallel.abortForce();
			}

			if(parseParameter(params, numparam, "f(R) model", par_string))
			{
				if(par_string[0] == 'R' || par_string[0] == 'r')
				{
					if(sim.fR_params[0] <= 0.)
					{
						COUT << " The coefficient a of f(R) = a * R^n is <= 0. This leads to a tachyonic instability for the scalaron. Closing...\n";
						parallel.abortForce();
					}
					else if(sim.fR_params[1] < 0.)
					{
						COUT << " The coefficient n of f(R) = a * R^n is <= 0. Closing...\n";
						parallel.abortForce();
					}
					else if(sim.fR_params[1] == 0)
					{
						COUT << " The exponent n of f(R) = a*R^n is set to 0, which corresponds to a pure Lambda term. Closing...\n";
						parallel.abortForce();
					}
					else if(sim.fR_params[1] == 1. || sim.fR_params[1] == 1)
					{
						COUT << " The exponent n of f(R) = a*R^n is set to 1. This is just a redefinition of Newton's constant. Closing...\n";
						parallel.abortForce();
					}
					else if(sim.fR_params[1] == 2. || sim.fR_params[1] == 2)
					{
						sim.fR_model = FR_MODEL_R2;
					}
					else
					{
						sim.fR_model = FR_MODEL_RN;
					}
				}
				else if(par_string[0] == 'H' || par_string[0] == 'h')
				{
					if(sim.num_fR_params < 3)
					{
						COUT << " Not enough parameters for Hu-Sawicki model! Closing...\n";
						parallel.abortForce();
					}
					else if(sim.fR_params[0] <= 0)
					{
						COUT << " The parameter m2 of Hu-Sawicki model should be strictly positive. Closing...\n";
						parallel.abortForce();
					}
					else if(sim.fR_params[2] < 1.)
					{
						COUT << " The parameter n of Hu-Sawicki model should be >= 1. Closing...\n";
						parallel.abortForce();
					}

					sim.fR_model = FR_MODEL_HU_SAWICKI;

					if(sim.lcdm_background)
					{
						cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad;
						COUT << " Hu-Sawicki Model selected, but LCDM background, so Omega_Lambda = 1 - Omega_m - Omega_r" << endl;
					}
					else
					{
						cosmo.Omega_Lambda = 0.;
						COUT << " Hu-Sawicki model, Omega_Lambda automatically set to zero.\n";
					}
				}
				else if(par_string[0] == 'D' || par_string[0] == 'd')
				{
					// if(sim.num_fR_params < 2)
					// {
					// 	COUT << " Not enough parameters for R^(1+delta) model! Closing...\n";
					// 	parallel.abortForce();
					// }
					// else if(sim.fR_params[0] <= 0)
					// {
					// 	COUT << " The parameter a of R * (R/a)^delta model must be positive. Closing...\n";
					// 	parallel.abortForce();
					// }
					// else if(sim.fR_params[1] <= 0 || sim.fR_params[1] > 1.)
					// {
					// 	COUT << " The parameter delta of R * (R/a)^delta model should be 0 < a <= 1. The case a = 0 is simply GR. Closing...\n";
					// 	parallel.abortForce();
					// }

					sim.fR_model = FR_MODEL_DELTA;
				}
				else
				{
					COUT << " /!\\ f(R) model not recognised. Closing..." << endl;
					parallel.abortForce();
				}

				// Added parser for Omega_Lambda in f(R) gravity -- Lambda will typically be zero, but can be nonzero
				if(sim.fR_model != FR_MODEL_HU_SAWICKI && !sim.lcdm_background)
				{
					COUT << " Using f(R) with" << COLORTEXT_YELLOW << " EXPLICIT " << COLORTEXT_RESET << "Lambda term, Omega_Lambda = ";

					if(parseParameter(params, numparam, "Omega_Lambda", cosmo.Omega_Lambda))
					{
						if(cosmo.Omega_Lambda > 0 && cosmo.Omega_Lambda <= 1 - cosmo.Omega_m - cosmo.Omega_rad)
						{
							COUT << cosmo.Omega_Lambda << endl; // TODO: might not be needed
						}
						else if(cosmo.Omega_Lambda > 1. - cosmo.Omega_m - cosmo.Omega_rad)
						{
							cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad;
							COUT << "1. - Omega_m - Omega_rad." << endl;
						}
						else
						{
							cosmo.Omega_Lambda = 0.;
							COUT << " 0." << endl;
						}
					}
					else
					{
						cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad;
						COUT << "1. - Omega_m - Omega_rad. (default)" << endl;
					}
				}
			}
			else
			{
				COUT << " /!\\ No f(R) model specified. Closing..." << endl;
				parallel.abortForce();
			}

			if(!parseParameter(params, numparam, "f(R) epsilon background", sim.fR_epsilon_bg))
			{
				sim.fR_epsilon_bg = FR_EPSILON_BACKGROUND_DEFAULT;
			}
			else if(sim.fR_epsilon_bg <= 0)
			{
				sim.fR_epsilon_bg = FR_EPSILON_BACKGROUND_DEFAULT;
				COUT << " /!\\ Wrong f(R) epsilon background specified. Using default: " << sim.fR_epsilon_bg << endl;
			}
		}
		else
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gravity theory unknown, using default (General Relativity)" << endl;
			sim.relativistic_flag = 1;
		}
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": gravity theory not selected, using default (General Relativity)" << endl;
		sim.relativistic_flag = 1;
	}

	COUT
	<< " cosmological parameters:" << endl
	<< "               Omega_m0 = " << cosmo.Omega_m << endl
	<< "             Omega_rad0 = " << cosmo.Omega_rad << endl
	<< "           Omega_Lambda = " << cosmo.Omega_Lambda << endl
	<< "                      h = " << cosmo.h << endl;

	parseParameter(params, numparam, "time step limit", sim.steplimit);

	if(!parseParameter(params, numparam, "output path", sim.output_path))
	{
		sim.output_path[0] = '\0';
	}
	if(!parseParameter(params, numparam, "hibernation path", sim.restart_path))
	{
		strcpy(sim.restart_path, sim.output_path);
	}
	if(!parseParameter(params, numparam, "hibernation file base", sim.basename_restart))
	{
		strcpy(sim.basename_restart, "restart");
	}
	if(!parseParameter(params, numparam, "BACKGROUND_NUMPTS", sim.BACKGROUND_NUMPTS))
	{
		sim.BACKGROUND_NUMPTS = -1;
	}

	// FIRST: BACKGROUND ONLY MODE -- Most parameters unused
	if(!parseParameter(params, numparam, "background only", sim.background_only))
	{
		sim.background_only = 0;
	}

	if(sim.background_only)
	{
		if(!parseParameter(params, numparam, "background final redshift", sim.z_fin))
		{
			sim.z_fin = 0.;
		}

		COUT
		<< " Background only mode. Initial redshift = " << sim.z_in << endl
		<< "                         Final redshift = " << sim.z_fin << endl;

		COUT << " ===================== UNUSED PARAMETERS =====================" << endl;
	}

	// FULL EVOLUTION
	if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
	{
		if(!parseParameter(params, numparam, "Newtonian f(R)", sim.relativistic_flag))
		{
			COUT << " Gravity theory set to: " << COLORTEXT_CYAN << "relativistic f(R)" << COLORTEXT_RESET << endl;
			sim.relativistic_flag = 1;
		}
		else if(sim.relativistic_flag)
		{
			sim.relativistic_flag = 0;
			COUT << " Gravity theory set to: " << COLORTEXT_CYAN << "Newtonian f(R)" << COLORTEXT_RESET << endl;
		}
		else
		{
			sim.relativistic_flag = 1;
			COUT << " Gravity theory set to: " << COLORTEXT_CYAN << "relativistic f(R)" << COLORTEXT_RESET << endl;
		}

		if(sim.relativistic_flag)
		{
			if(!parseParameter(params, numparam, "xi_Hubble", sim.xi_Hubble))
			{
				sim.xi_Hubble = 0;
			}
			else if(sim.xi_Hubble)
			{
				COUT << " xi_Hubble mode active." << endl;
			}
		}

		if(sim.lcdm_background)
		{
			COUT << " Background is LCDM.\n";
		}
		else // Only switch between LCDM and f(R) background if sim.lcdm_background is not flagged
		{
			if(!parseParameter(params, numparam, "switch to f(R) redshift", sim.z_switch_fR_background))
			{
				sim.z_switch_fR_background = -100.;
			}

			if(sim.z_switch_fR_background > 0. && sim.z_switch_fR_background < sim.z_in)
			{
				COUT << " Full f(R) background expansion will start after z = " << sim.z_switch_fR_background << endl;
				sim.lcdm_background = 1;
			}
		}
	}

	if(!parseParameter(params, numparam, "generic file base", sim.basename_generic))
	{
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			if(sim.relativistic_flag)
			{
				strcpy(sim.basename_generic, "fR");
			}
			else
			{
				strcpy(sim.basename_generic, "fR_Newtonian");
			}
		}
		else if(sim.relativistic_flag)
		{
			strcpy(sim.basename_generic, "lcdm");
		}
		else
		{
			strcpy(sim.basename_generic, "Newton");
		}
	}

	if(parseParameter(params, numparam, "IC generator", par_string))
	{
		if(par_string[0] == 'B' || par_string[0] == 'b')
		{
			ic.generator = ICGEN_BASIC;
		}
		else if(par_string[0] == 'R' || par_string[0] == 'r')
		{
			ic.generator = ICGEN_READ_FROM_DISK;
		}
#ifdef ICGEN_PREVOLUTION
		else if(par_string[0] == 'P' || par_string[0] == 'p')
		{
			ic.generator = ICGEN_PREVOLUTION;
		}
#endif
#ifdef ICGEN_SONG
		else if(par_string[0] == 'S' || par_string[0] == 's')
		{
			ic.generator = ICGEN_SONG;
		}
#endif
#ifdef ICGEN_FALCONIC
		else if(par_string[0] == 'F' || par_string[0] == 'f')
		{
			ic.generator = ICGEN_FALCONIC;
		}
#endif
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": IC generator not recognized!" << endl;
			parallel.abortForce();
		}
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": IC generator not specified, selecting default (basic)" << endl;
		ic.generator = ICGEN_BASIC;
	}

	parseParameter(params, numparam, "seed", ic.seed);

	for(i=0; i<MAX_PCL_SPECIES; i++)
	{
		pptr[i] = ic.pclfile[i];
	}

	if(ic.generator == ICGEN_READ_FROM_DISK)
	{
		if(!parseParameter(params, numparam, "particle file", pptr, i))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no particle file specified!" << endl;
			parallel.abortForce();
		}
	}
	else if(!parseParameter(params, numparam, "template file", pptr, i))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no template file specified!" << endl;
		parallel.abortForce();
	}

	for(; i<MAX_PCL_SPECIES; i++)
	{
		strcpy(ic.pclfile[i], ic.pclfile[i-1]);
	}

#ifdef ICGEN_FALCONIC
		if((!parseParameter(params, numparam, "mPk file", ic.pkfile) && !parseParameter(params, numparam, "Tk file", ic.tkfile) && ic.generator != ICGEN_READ_FROM_DISK && ic.generator != ICGEN_FALCONIC)
#else
		if((!parseParameter(params, numparam, "mPk file", ic.pkfile) && !parseParameter(params, numparam, "Tk file", ic.tkfile) && ic.generator != ICGEN_READ_FROM_DISK)
#endif
#ifdef ICGEN_PREVOLUTION
	  || ic.generator == ICGEN_PREVOLUTION)
#else
		)
#endif
	{
#ifdef HAVE_CLASS
			COUT << " initial transfer functions will be computed by calling CLASS" << endl;
#else
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no power spectrum file nor transfer function file specified!" << endl;
			parallel.abortForce();
#endif
	}

	if(parseParameter(params, numparam, "correct displacement", par_string))
	{
		if(par_string[0] == 'Y' || par_string[0] == 'y')
		{
			ic.flags |= ICFLAG_CORRECT_DISPLACEMENT;
		}
		else if(par_string[0] != 'N' && par_string[0] != 'n')
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": setting chosen for deconvolve displacement option not recognized, using default (no)" << endl;
		}
	}

	if(parseParameter(params, numparam, "k-domain", par_string))
	{
		if(par_string[0] == 'S' || par_string[0] == 's')
		{
			ic.flags |= ICFLAG_KSPHERE;
		}
		else if(par_string[0] != 'C' && par_string[0] != 'c')
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": setting chosen for k-domain option not recognized, using default (cube)" << endl;
		}
	}

	for(i=0; i<MAX_PCL_SPECIES; i++)
	{
		ic.numtile[i] = 0;
	}

	// TODO Commented this part out to make tiling factor automatic -- see later
	// if(!parseParameter(params, numparam, "tiling factor", ic.numtile, i) && ic.generator != ICGEN_READ_FROM_DISK)
	// {
	// 	COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling factor not specified, using default value for all species (1)" << endl;
	// 	ic.numtile[0] = 1;
	// 	i = 1;
	// }
	//
	// for(; i<MAX_PCL_SPECIES; i++)
	// {
	// 	ic.numtile[i] = ic.numtile[i-1];
	// }
	//
	// for(i=0; i<MAX_PCL_SPECIES; i++)
	// {
	// 	if(ic.numtile[i] <= 0 && ic.generator != ICGEN_READ_FROM_DISK)
	// 	{
	// 		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tiling number for particle template not set properly; using default value (1)" << endl;
	// 		ic.numtile[i] = 1;
	// 	}
	// }
	// END of commented part

	if(ic.pkfile[0] != '\0')
	{
		sim.baryon_flag = 0;
		if(parseParameter(params, numparam, "baryon treatment", par_string))
		{
			if(par_string[0] != 'i' && par_string[0] != 'I')
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using mPk file forces baryon treatment = ignore" << endl;
			}
		}
	}
	else if(parseParameter(params, numparam, "baryon treatment", par_string))
	{
		if(par_string[0] == 'i' || par_string[0] == 'I')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "ignore" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 0;
		}
		else if(par_string[0] == 's' || par_string[0] == 'S')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "sample" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 1;
		}
		else if(par_string[0] == 'b' || par_string[0] == 'B')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "blend" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 2;
		}
		else if(par_string[0] == 'h' || par_string[0] == 'H')
		{
			COUT << " baryon treatment set to: " << COLORTEXT_CYAN << "hybrid" << COLORTEXT_RESET << endl;
			sim.baryon_flag = 3;
		}
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": baryon treatment not supported!" << endl;
			parallel.abortForce();
		}
	}
	else if(ic.generator == ICGEN_READ_FROM_DISK)
	{
		sim.baryon_flag = 0;
	}
	else
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": baryon treatment not specified, using default (blend)" << endl;
		sim.baryon_flag = 2;
	}

	if(parseParameter(params, numparam, "radiation treatment", par_string))
	{
		if(par_string[0] == 'i' || par_string[0] == 'I')
		{
			sim.radiation_flag = 0;
			COUT << " radiation treatment set to: " << COLORTEXT_CYAN << "ignore" << COLORTEXT_RESET << endl;
		}
#ifdef HAVE_CLASS
		else if(par_string[0] == 'c' || par_string[0] == 'C')
		{
			sim.radiation_flag = 1;
			COUT << " radiation treatment set to: " << COLORTEXT_CYAN << "CLASS" << COLORTEXT_RESET << endl;
			if((ic.pkfile[0] != '\0' || ic.tkfile[0] != '\0')
#ifdef ICGEN_PREVOLUTION
				&& ic.generator != ICGEN_PREVOLUTION
#endif
			)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using radiation treatment = CLASS and providing initial power spectra / transfer functions independently" << endl;
				COUT << "              is dangerous! In order to ensure consistency, it is recommended to call CLASS directly." << endl;
			}
		}
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": radiation treatment not supported!" << endl;
			parallel.abortForce();
		}
#else
		else
		{
			sim.radiation_flag = 0;
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": CLASS is not available, setting radiation treatment = ignore" << endl;
#endif
		}
	}
	else
	{
		sim.radiation_flag = 0;
	}

	parseParameter(params, numparam, "relaxation redshift", ic.z_relax);

	if(ic.generator == ICGEN_READ_FROM_DISK)
	{
		parseParameter(params, numparam, "restart redshift", ic.z_ic);
		if(!parseParameter(params, numparam, "cycle", ic.restart_cycle))
		{
			if(sim.radiation_flag)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": using radiation treatment = CLASS and IC generator = read from disk does not guarantee" << endl;
				COUT << "              that the realization of the radiation perturbations is consistent!" << endl;
			}
		}
		parseParameter(params, numparam, "tau", ic.restart_tau);
		parseParameter(params, numparam, "dtau", ic.restart_dtau);
		for(i=0; i<METRICFILE_LENGTH; i++)
		{
			pptr[i] = ic.metricfile[i];
		}

		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			parseParameter(params, numparam, "dtau_old", ic.restart_dtau_old);
			parseParameter(params, numparam, "dtau_old_2", ic.restart_dtau_old_2);
			parseParameter(params, numparam, "dtau_osci", ic.restart_dtau_osci);
			parseParameter(params, numparam, "dtau_bg", ic.restart_dtau_bg);
			parseParameter(params, numparam, "scale_factor", ic.restart_a);
			parseParameter(params, numparam, "Hubble", ic.restart_Hubble);
			parseParameter(params, numparam, "Rbar", ic.restart_Rbar);
			parseParameter(params, numparam, "dot_Rbar", ic.restart_dot_Rbar);

			// TODO: Maybe check that all of these are correcty parsed?
		}

		parseParameter(params, numparam, "metric file", pptr, i);
		if(parseParameter(params, numparam, "gevolution version", ic.restart_version))
		{
			if(ic.restart_version - FREVOLUTION_VERSION > 0.0001)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": version number of settings file (" << ic.restart_version << ") is higher than version of executable (" << FREVOLUTION_VERSION << ")!" << endl;
			}
		}
	}
#ifdef ICGEN_PREVOLUTION
		else if(ic.generator == ICGEN_PREVOLUTION)
		{
			if(!parseParameter(params, numparam, "prevolution redshift", ic.z_ic))
			{
				COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no starting redshift specified for IC generator = prevolution" << endl;
				parallel.abortForce();
			}

			parseParameter(params, numparam, "prevolution Courant factor", ic.Cf);
		}
#endif

	if(!parseParameter(params, numparam, "A_s", ic.A_s) && (
#ifdef ICGEN_FALCONIC
		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
		ic.generator == ICGEN_PREVOLUTION ||
#endif
		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": power spectrum normalization not specified, using default value (2.215e-9)" << endl;
	}

	if(!parseParameter(params, numparam, "n_s", ic.n_s) && (
#ifdef ICGEN_FALCONIC
		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
		ic.generator == ICGEN_PREVOLUTION ||
#endif
		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": scalar spectral index not specified, using default value (0.9619)" << endl;
	}

	if(!parseParameter(params, numparam, "k_pivot", ic.k_pivot) && (
#ifdef ICGEN_FALCONIC
		ic.generator == ICGEN_FALCONIC ||
#endif
#ifdef ICGEN_PREVOLUTION
		ic.generator == ICGEN_PREVOLUTION ||
#endif
		sim.radiation_flag > 0 || (ic.pkfile[0] == '\0' && ic.generator != ICGEN_READ_FROM_DISK)))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": pivot scale not specified, using default value (0.05 / Mpc)" << endl;
	}

	if(!parseParameter(params, numparam, "CYCLE_INFO_INTERVAL", sim.CYCLE_INFO_INTERVAL))
	{
		sim.CYCLE_INFO_INTERVAL = 10;
	}
	else if(sim.CYCLE_INFO_INTERVAL <= 0)
	{
		sim.CYCLE_INFO_INTERVAL = 10;
		COUT << " /!\\ Wrong CYCLE_INFO_INTERVAL specified. Using default: " << sim.CYCLE_INFO_INTERVAL << endl;
	}

	if(parseParameter(params, numparam, "vector method", par_string))
	{
		if(par_string[0] == 'p' || par_string[0] == 'P')
		{
			COUT << " vector method set to: " << COLORTEXT_CYAN << "parabolic" << COLORTEXT_RESET << endl;
			sim.vector_flag = VECTOR_PARABOLIC;
		}
		else if(par_string[0] == 'e' || par_string[0] == 'E')
		{
			COUT << " vector method set to: " << COLORTEXT_CYAN << "elliptic" << COLORTEXT_RESET << endl;
			sim.vector_flag = VECTOR_ELLIPTIC;
		}
		else
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": vector method not supported!" << endl;
			parallel.abortForce();
		}
	}

	parseParameter(params, numparam, "Ngrid", sim.numpts);
	if(sim.numpts < 2 || !isfinite(sim.numpts))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of grid points not set properly!" << endl;
		parallel.abortForce();
	}

	// TODO Check if this is correct
	if(ic.generator != ICGEN_READ_FROM_DISK)
	{
		COUT << " Setting tiling factor to Ngrid/4" << endl;
		for(i=0; i<MAX_PCL_SPECIES; i++)
		{
			ic.numtile[i] = sim.numpts / 4;
		}
		COUT << "sim.numpts = " << sim.numpts << endl;
		COUT << "ic.numtile[0] = " << ic.numtile[0] << endl;
	}

	if(parseParameter(params, numparam, "downgrade factor", sim.downgrade_factor))
	{
		if(sim.downgrade_factor < 1 || sim.downgrade_factor >= sim.numpts || !isfinite(sim.downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor makes no sense!" << endl;
			parallel.abortForce();
		}
		if(sim.downgrade_factor > 1 && (sim.numpts % parallel.grid_size()[0] || sim.numpts % parallel.grid_size()[1] || (sim.numpts / parallel.grid_size()[0]) % sim.downgrade_factor || (sim.numpts / parallel.grid_size()[1]) % sim.downgrade_factor))
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": downgrade factor does not appear to be compatible with process layout at given Ngrid!" << endl;
			parallel.abortForce();
		}
	}

	if(!parseParameter(params, numparam, "move limit", sim.movelimit))
	{
		sim.movelimit = (double) sim.numpts;
	}

	if(ic.z_relax < -1.)
	{
		ic.z_relax = sim.z_in;
	}
#ifdef ICGEN_PREVOLUTION
		else if(ic.generator == ICGEN_PREVOLUTION && ic.z_relax < sim.z_in)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": relaxation redshift cannot be below initial redshift for IC generator = prevolution; reset to initial redshift!" << endl;
			ic.z_relax = sim.z_in;
		}
#endif

	if(ic.z_ic < sim.z_in && ic.generator != ICGEN_READ_FROM_DISK)
	{
		ic.z_ic = sim.z_in;
	}

	parseParameter(params, numparam, "snapshot redshifts", sim.z_snapshot, sim.num_snapshot);
	if(sim.num_snapshot > 0)
	{
		qsort((void *) sim.z_snapshot, (size_t) sim.num_snapshot, sizeof(double), sort_descending);
	}

	parseParameter(params, numparam, "Pk redshifts", sim.z_pk, sim.num_pk);
	if(sim.num_pk > 0)
	{
		qsort((void *) sim.z_pk, (size_t) sim.num_pk, sizeof(double), sort_descending);
	}

	parseParameter(params, numparam, "hibernation redshifts", sim.z_restart, sim.num_restart);
	if(sim.num_restart > 0)
	{
		qsort((void *) sim.z_restart, (size_t) sim.num_restart, sizeof(double), sort_descending);
	}

	parseParameter(params, numparam, "hibernation wallclock limit", sim.wallclocklimit);
	i = MAX_PCL_SPECIES;

	parseParameter(params, numparam, "tracer factor", sim.tracer_factor, i);
	for(; i>0; i--)
	{
		if(sim.tracer_factor[i-1] < 1)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": tracer factor not set properly; using default value (1)" << endl;
			sim.tracer_factor[i-1] = 1;
		}
	}

	if(!parseParameter(params, numparam, "Pk bins", sim.numbins))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of Pk bins not set properly; using default value (64)" << endl;
		sim.numbins = 64;
	}

	parseParameter(params, numparam, "Courant factor", sim.Cf);
	if(ic.Cf < 0.)
	{
		ic.Cf = sim.Cf;
	}

	if(!parseParameter(params, numparam, "multigrid n-grids", sim.multigrid_n_grids))
	{
		sim.multigrid_n_grids = 1; //  TODO: Must default to 1 otherwise it gives malloc() errors
	}
	else if(sim.multigrid_n_grids < 1)
	{
		sim.multigrid_n_grids = 1; //  TODO: Must default to 1 otherwise it gives malloc() errors
		COUT << " /!\\ Wrong multigrid_n_grids specified. Using default: " << sim.multigrid_n_grids << endl;
	}

	if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
	{
		if(!parseParameter(params, numparam, "quasi-static", sim.quasi_static))
		{
			sim.quasi_static = 0;
		}
		else if(sim.quasi_static)
		{
			COUT << " Quasi-static mode active.\n";
		}

		if(!parseParameter(params, numparam, "GR curvature", sim.gr_curvature))
		{
			sim.gr_curvature = 0;
		}
		else if(sim.gr_curvature)
		{
			COUT << " Curvature is as in GR: deltaR = 8piG_deltaT.\n";
		}

		if(!parseParameter(params, numparam, "f(R) epsilon fields", sim.fR_epsilon_fields))
		{
			sim.fR_epsilon_fields = FR_EPSILON_FIELDS_DEFAULT;
			COUT << " /!\\ f(R) epsilon fields not specified. Using default: " << sim.fR_epsilon_fields << endl;
		}
		else if(sim.fR_epsilon_fields <= 0.)
		{
			sim.fR_epsilon_fields = FR_EPSILON_FIELDS_DEFAULT;
			COUT << " /!\\ Wrong f(R) epsilon fields specified. Using default: " << sim.fR_epsilon_fields << endl;
		}

		if(!parseParameter(params, numparam, "f(R) target precision", sim.fR_target_precision))
		{
			sim.fR_target_precision = FR_TARGET_PRECISION_DEFAULT;
			COUT << " /!\\ f(R) target precision not specified. Using default: " << sim.fR_target_precision << endl;
		}

		if(sim.fR_target_precision <= 0.)
		{
			sim.fR_target_precision = FR_TARGET_PRECISION_DEFAULT;
			COUT << " /!\\ Wrong f(R) target precision specified. Using default: " << sim.fR_target_precision << endl;
		}

		// TODO: Find a way to remove fR_count_max
		if(!parseParameter(params, numparam, "f(R) count max", sim.fR_count_max))
		{
			sim.fR_count_max = FR_COUNT_MAX_DEFAULT;
			COUT << "/!\\ f(R) count max not specified. Using default: " << sim.fR_count_max << endl;
		}
		else if(sim.fR_count_max < 1)
		{
			sim.fR_count_max = FR_COUNT_MAX_DEFAULT;
			COUT << "/!\\ Wrong f(R) count max specified. Using default: " << sim.fR_count_max << endl;
		}

		if(sim.fR_model != FR_MODEL_R2)
		{
			if(parseParameter(params, numparam, "relaxation method", sim.relaxation_method))
			{
				sim.multigrid_or_not = 1;
				if(sim.relaxation_method == METHOD_FMG)
				{
					COUT << " Relaxation: multigrid FMG mode" << endl;
				}
				else if(sim.relaxation_method == METHOD_RELAX)
				{
					sim.multigrid_or_not = 0;
					COUT << " Relaxation: single layer (no multigrid)" << endl;
				}
				else if(sim.relaxation_method == METHOD_MULTIGRID)
				{
					COUT << " Relaxation: multigrid" << endl;
					if(!parseParameter(params, numparam, "multigrid shape", par_string))
					{
						sim.multigrid_shape = MULTIGRID_SHAPE_V;
						COUT << " /!\\ Multigrid shape not specified: Using default: V" << endl;
					}
					else if(par_string[0] == 'V' || par_string[0] == 'v')
					{
						sim.multigrid_shape = MULTIGRID_SHAPE_V;
						COUT << " Multigrid shape = V" << endl;
					}
					else if(par_string[0] == 'W' || par_string[0] == 'w')
					{
						sim.multigrid_shape = MULTIGRID_SHAPE_W;
						COUT << " Multigrid shape = W" << endl;
					}
					else
					{
						sim.multigrid_shape = MULTIGRID_SHAPE_V;
						COUT << " /!\\ Multigrid shape not recognised. Using default: V" << endl;
					}
				}
				else
				{
					sim.relaxation_method = METHOD_RELAX;
					COUT << " /!\\ Wrong relaxation method specified. Using default: single layer (no multigrid)" << endl;
				}
			}
			else
			{
				sim.relaxation_method = METHOD_RELAX;
				COUT << " /!\\ Relaxation method not specified. Using default: single layer (no multigrid)" << endl;
			}

			if(!parseParameter(params, numparam, "truncate relaxation", sim.truncate_relaxation))
			{
				COUT << " Relaxation solver never truncated, even if convergence slows down." << endl;
				sim.truncate_relaxation = 0;
			}
			else if(sim.truncate_relaxation <= 0)
			{
				COUT << " Relaxation solver never truncated, even if convergence slows down." << endl;
			}
			else
			{
				COUT << " Relaxation truncated after " << sim.truncate_relaxation << " steps with (almost) no improvement." << endl;

				if(!parseParameter(params, numparam, "truncation threshold", sim.relaxation_truncation_threshold))
				{
					COUT << " Setting truncation threshold at default = 0.01." << endl;
					sim.relaxation_truncation_threshold = 0.01;
				}
				else
				{
					COUT << " Truncation threshold set at: " << sim.relaxation_truncation_threshold << endl;
				}
			}

			if(!parseParameter(params, numparam, "red black", sim.red_black))
			{
				sim.red_black = 0;
				COUT << " Using simple sequential sweep of sites." << endl;
			}
			else if(sim.red_black)
			{
				COUT << " Using red-black sweep scheme." << endl;
			}
			else
			{
				sim.red_black = 0; // TODO: This is most probably redundant. Just making sure that it's really zero :)
				COUT << " Using simple sequential sweep of sites." << endl;
			}

			if(!parseParameter(params, numparam, "relaxation error", sim.relaxation_error))
			{
				sim.relaxation_error = RELAXATION_ERROR_DEFAULT;
				COUT << "/!\\ Relaxation error not specified. Using default: " << sim.relaxation_error << endl;
			}
			else if(sim.relaxation_error <= 0.)
			{
				sim.relaxation_error = RELAXATION_ERROR_DEFAULT;
				COUT << "/!\\ Wrong relaxation error specified. Using default: " << sim.relaxation_error << endl;
			}

			if(!parseParameter(params, numparam, "overrelaxation factor", sim.overrelaxation_factor))
			{
				sim.overrelaxation_factor = 1.;
				COUT << " No overrelaxation fator set. Using default value = " << sim.overrelaxation_factor << endl;
			}
			else if(sim.overrelaxation_factor >= 2. || sim.overrelaxation_factor <= 0.)
			{
				sim.overrelaxation_factor = 1.;
				COUT << " Wrong value for overrelaxation fator. Using default value = " << sim.overrelaxation_factor << endl;
			}

			if(sim.multigrid_or_not)
			{
				if(!parseParameter(params, numparam, "pre-smoothing", sim.pre_smoothing))
				{
					sim.pre_smoothing = PRE_SMOOTHING_DEFAULT;
					COUT << "/!\\ Pre-smoothing not specified. Using default: " << sim.pre_smoothing << endl;
				}
				else if(sim.pre_smoothing < 1)
				{
					sim.pre_smoothing = PRE_SMOOTHING_DEFAULT;
					COUT << "/!\\ Wrong pre-smoothing specified. Using default: " << sim.pre_smoothing << endl;
				}

				if(!parseParameter(params, numparam, "post-smoothing", sim.post_smoothing))
				{
					sim.post_smoothing = POST_SMOOTHING_DEFAULT;
					COUT << "/!\\ Post-smoothing not specified. Using default: " << sim.post_smoothing << endl;
				}
				else if(sim.post_smoothing < 1)
				{
					sim.post_smoothing = POST_SMOOTHING_DEFAULT;
					COUT << "/!\\ Wrong post-smoothing specified. Using default: " << sim.post_smoothing << endl;
				}

				if(!parseParameter(params, numparam, "multigrid n-cycles", sim.multigrid_n_cycles))
				{
					sim.multigrid_n_cycles = 1;
					COUT << " /!\\ Multigrid n_cycles not specified. Using default: " << sim.multigrid_n_cycles << endl;
				}
				else if(sim.multigrid_n_cycles < 1)
				{
					sim.multigrid_n_cycles = 1;
					COUT << " /!\\ Wrong multigrid n_cycles specifiied. Using default: " << sim.multigrid_n_cycles << endl;
				}

				if(!parseParameter(params, numparam, "restrict mode", par_string))
				{
					sim.multigrid_restrict_mode = RESTRICT_DELTAR;
					COUT << " /!\\ Restrict mode not specified. Using default: deltaR." << endl;
				}
				else if(par_string[0] == 's' || par_string[0] == 'S')
				{
					sim.multigrid_restrict_mode = RESTRICT_XI;
					COUT << " Restriction variable: xi." << endl;
				}
				else if(par_string[0] == 'r' || par_string[0] == 'R' || par_string[0] == 'd' || par_string[0] == 'D')
				{
					sim.multigrid_restrict_mode = RESTRICT_DELTAR;
					COUT << " Restriction variable: deltaR." << endl;
				}
				else
				{
					sim.multigrid_restrict_mode = RESTRICT_DELTAR;
					COUT << " /!\\ Wrong restrict mode specified. Using default: deltaR." << endl;
				}

				if(!parseParameter(params, numparam, "check shape", sim.multigrid_check_shape))
				{
					sim.multigrid_check_shape = 0;
				}
			}
		}
	}

	if(!parseParameter(params, numparam, "check fields", sim.check_fields))
	{
		sim.check_fields = 0;
	}
	else if(sim.check_fields)
	{
		if(!parseParameter(params, numparam, "check pause", sim.check_pause))
		{
			sim.check_pause = 0;
		}

		if(!parseParameter(params, numparam, "check precision", sim.check_fields_precision))
		{
			sim.check_fields_precision = 6;
		}
	}

	parseParameter(params, numparam, "check redshift", sim.z_check);
	parseFieldSpecifiers(params, numparam, "snapshot outputs", sim.out_snapshot);
	parseFieldSpecifiers(params, numparam, "Pk outputs", sim.out_pk);
	parseFieldSpecifiers(params, numparam, "check outputs", sim.out_check);

	if((sim.num_snapshot <= 0 || sim.out_snapshot == 0) && (sim.num_pk <= 0 || sim.out_pk == 0))
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": no output specified!" << endl;
	}

	if(!parseParameter(params, numparam, "snapshot file base", sim.basename_snapshot))
	{
		strcpy(sim.basename_snapshot, "snap");
	}

	if(!parseParameter(params, numparam, "Pk file base", sim.basename_pk))
	{
		strcpy(sim.basename_pk, "pk");
	}

	if(parseParameter(params, numparam, "particle save mode", par_string))
	{
		if(par_string[0] == 'h' || par_string[0] == 'H')
		{
			sim.hibernation_save_mode = HIB_SAVE_HDF5;
		}
		else if(par_string[0] == 'g' || par_string[0] == 'G')
		{
			sim.hibernation_save_mode = HIB_SAVE_GADGET2;
		}
		else
		{
			sim.hibernation_save_mode = HIB_SAVE_GADGET2;
			COUT << " Wrong particle save mode specified for hibernation. Defaults to HDF5." << endl;
		}
	}

	if(!parseParameter(params, numparam, "switch delta_rad", sim.z_switch_deltarad))
	{
		sim.z_switch_deltarad = 0.;
	}

	i = MAX_PCL_SPECIES-2;
	if(!parseParameter(params, numparam, "switch delta_ncdm", sim.z_switch_deltancdm, i))
	{
		sim.z_switch_deltancdm[0] = sim.z_in;
		i = 1;
	}
	for(; i<MAX_PCL_SPECIES-2; i++)
	{
		sim.z_switch_deltancdm[i] = sim.z_switch_deltancdm[i-1];
	}

	if(!parseParameter(params, numparam, "switch linear chi", sim.z_switch_linearchi))
	{
		if(sim.relativistic_flag > 0)
		{
			sim.z_switch_linearchi = 0.;
		}
		else
		{
			sim.z_switch_linearchi = 0.011;
		}

		for(i=0; i<cosmo.num_ncdm; i++)
		{
			if(sim.z_switch_linearchi < sim.z_switch_deltancdm[i])
			{
				sim.z_switch_linearchi = sim.z_switch_deltancdm[i];
			}
		}
	}
	else if(sim.relativistic_flag == 0 && sim.z_switch_linearchi <= 0.01)
	{
		COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": with gravity theory = Newton the switch linear chi redshift must be larger than 0.01." << endl;
		COUT << "              setting switch linear chi = 0.011" << endl;
		sim.z_switch_linearchi = 0.011;
	}

	i = MAX_PCL_SPECIES-2;
	if(!parseParameter(params, numparam, "switch B ncdm", sim.z_switch_Bncdm, i))
	{
		for(i=0; i<MAX_PCL_SPECIES-2; i++)
		{
			sim.z_switch_Bncdm[i] = sim.z_switch_deltancdm[i];
		}
	}
	for(; i<MAX_PCL_SPECIES-2; i++)
	{
		sim.z_switch_Bncdm[i] = sim.z_switch_Bncdm[i-1];
	}

	if(sim.background_only)
	{
		COUT << " =================== END UNUSED PARAMETERS ===================" << endl;
	}

	for(i=0; i<numparam; i++)
	{
		if(params[i].used)
		{
			usedparams++;
		}
	}

	return usedparams;
}

#endif
