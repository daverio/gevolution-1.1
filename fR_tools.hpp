//////////////////////////
// FR_tools.hpp
//////////////////////////
//
// code components related to f(R)
//
// Authors: David Daverio (Cambridge University)
//          Lorenzo Reverberi (Cape Town University)
//
// Last modified: February 2016
//
//////////////////////////

#ifndef FR_TOOLS_HEADER
#define FR_TOOLS_HEADER

void hold()
{
	COUT << "Press ENTER to continue..." << std::flush;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}


///////////////////// f(R) FUNCTIONS /////////////////////

Real f(const double R, const metadata & sim, int code)
{
	double output,
				 Rpow;

	if(sim.fR_type == FR_TYPE_RN)
	{
		Rpow = pow(R, sim.fR_params[1]);
		output = sim.fR_params[0] * Rpow;
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		output = sim.fR_params[0] * R * R;
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = sim.fR_params[0];
		c2 = sim.fR_params[1];
		n  = sim.fR_params[2];
		c1 = sim.fR_params[3];
		Rpow = pow(R/m2, n);
		output = - m2 * c1 * Rpow / (1. + c2 * Rpow);
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];
		output = R * (pow(R/a, delta) - 1.);
	}
	else
	{
		cout << " Something went wrong when computing f(R) (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(!std::isnan(output))
	{
		return output;
	}
	else
	{
		cout << " f(R) Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", f(R) = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}

// fR(R) -- First derivative with respect to R
Real fR(const double R, const metadata & sim, int code)
{
	double output,
				 Rpow;

	if(sim.fR_type == FR_TYPE_RN)
	{
		Rpow = pow(R,sim.fR_params[1] - 1.);
		output = sim.fR_params[0] * sim.fR_params[1] * Rpow;
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		output = 2. * sim.fR_params[0] * R;
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = sim.fR_params[0];
		c2 = sim.fR_params[1];
		n  = sim.fR_params[2];
		c1 = sim.fR_params[3];
		if(n == 1)
		{
			output = -c1 / (1. + c2 * R / m2) / (1. + c2 * R / m2);
		}
		else
		{
			Rpow = pow(R/m2, n);
			output = - c1 * n * pow(R/m2, n-1) / (1. + c2*Rpow) / (1. + c2*Rpow);
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];
		output = (1. + delta) * pow(R/a, delta) - 1.;
	}
	else
	{
		cout << " Something went wrong when computing FR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(!std::isnan(output))
	{
		return output;
	}
	else
	{
		cout << " fR Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", fR = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}

// F_{RR}(R) -- Second derivative with respect to R
Real fRR(const double R, const metadata & sim, int code)
{
	double output,
				 Roverm2,
				 Rpow;

	if(sim.fR_type == FR_TYPE_RN)
	{
		Rpow = pow(R, sim.fR_params[1] - 2.);
		output = sim.fR_params[0] * sim.fR_params[1] * (sim.fR_params[1] - 1.) * Rpow;
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		output = 2. * sim.fR_params[0];
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = sim.fR_params[0];
		c2 = sim.fR_params[1];
		n  = sim.fR_params[2];
		c1 = sim.fR_params[3];
		Roverm2 = R/m2;
		Rpow = pow(Roverm2, n);
		if(n == 1.)
		{
			output = 2. * c1 * c2 / ( m2 * (1. + c2*Roverm2) * (1. + c2*Roverm2) * (1. + c2*Roverm2) );
		}
		else
		{
			output = c1 * n * Rpow * (1. - n + c2*(1. + n)*Rpow ) / ( Roverm2 * Roverm2 * m2 * (1. + c2*Rpow) * (1. + c2*Rpow) * (1. + c2*Rpow) );
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];
		output = delta * (delta + 1.) * pow(R/a, delta) / R;
	}
	else
	{
		cout << " Something went wrong when computing fRR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(output && !std::isnan(output))
	{
		return output;
	}
	else if(!output)
	{
		cout << " fRR Evaluated to 0 (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", fRR = " << output << endl;
		return output;
	}
	else if(std::isnan(output))
	{
		cout << " fRR Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", fRR = " << output << endl;
		return output;
	}
}

// F_{RRR}(R) -- Third derivative with respect to R
Real fRRR(const double R, const metadata & sim, int code)
{
	double output,
				 Roverm2,
				 Rpow;

	if(sim.fR_type == FR_TYPE_RN)
	{
		Rpow = pow(R, sim.fR_params[1] - 3.);
		output = sim.fR_params[0] * sim.fR_params[1] * (sim.fR_params[1] - 1.) * (sim.fR_params[1] - 2.) * Rpow;
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		output = 0.;
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = sim.fR_params[0];
		c2 = sim.fR_params[1];
		n  = sim.fR_params[2];
		c1 = sim.fR_params[3];
		Roverm2 = R/m2;
		Rpow = pow(Roverm2, n);
		if(n == 1)
		{
			output = -6. * c1 * c2 * c2 / ( m2 * m2 * (1. + c2 * Roverm2) * (1. + c2 * Roverm2) * (1. + c2 * Roverm2) * (1. + c2 * Roverm2) );
		}
		else if(n == 2)
		{
			output = -24. * c1 * c2 * Roverm2 * (-1. + c2 * Rpow) / ( m2 * m2 * (1. + c2 * Rpow) * (1. + c2 * Rpow) * (1. + c2 * Rpow) * (1. + c2 * Rpow) );
		}
		else
		{
			output = c1 * n * Rpow * ( -2. + 3. * n - n * n + 4. * c2 * (n * n - 1.) * Rpow - c2 * c2 * (1. + n) * (2. + n) * Rpow * Rpow ) / (Roverm2 * Roverm2 * Roverm2 * m2 * m2 * (1. + c2 * Rpow) * (1. + c2 * Rpow) * (1. + c2 * Rpow) * (1. + c2 * Rpow));
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];
		output = delta * (delta * delta - 1.) * pow(R/a, delta) / R / R;
	}
	else
	{
		cout << " Something went wrong when computing fRRR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(!std::isnan(output))
	{
		return output;
	}
	else if(std::isnan(output))
	{
		cout << " fRRR Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", fRRR = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}

// Print out f(R) model details
void fR_details(const cosmology cosmo, metadata * sim, const double fourpiG)
{
	if(sim->fR_type == FR_TYPE_HU_SAWICKI)
	{
		//TODO: Check whether we need to rescale something by some power of the boxsize
		// Hu-Sawicki Model: arXiv:0705.1158 -- See paper for details
		// m2 := m^2 (comparison between these values and the nomenclature of the original paper)
		// f(R) = -m2 * c1 * (R/m2)^n / (1 + c2*(R/m2)^n)
		// In settings.ini:
		// params[0] = m^2 is given in units of (matter energy density today) * 8piG / 3.
		// params[1] = |fR0| (<< 1)
		// params[2] = n
		// Constructed parameters:
		// params[1] = c2 = 16piG (1 - Omega_m - Omega_rad) / |fR0| / m2 / (12/Omega_m - 9)^(n+1)
		// params[3] = c1, calibrated to give c1*m2/c2 = 2 * 8piG * rho_Lambda [see 0705.1158 equations (25,26)]
		COUT << " Model: Hu-Sawicki\n f(R) = - m2 * c1 * pow(R/m2, n) / (1. + c2 * pow(R/m2, n))" << endl
				 << " with m2 = " << sim->fR_params[0] << " * 8piG * rho_{m0} / 3. = ";
		sim->fR_params[0] *= fourpiG * cosmo.Omega_m / 1.5;
		sim->fR_params[1] = 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[1] / pow(12./cosmo.Omega_m - 9., sim->fR_params[2] + 1.) / sim->fR_params[0];
		sim->fR_params[3] = sim->fR_params[1] * 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[0];
		COUT << sim->fR_params[0] << endl;
		COUT << "      c1 = " << sim->fR_params[3] << endl
				 << "      c2 = " << sim->fR_params[1] << endl
				 << "      n  = " << sim->fR_params[2] << endl
				 << "so |fR0| ~ -n*c1/c2^2/(12/Omega_m - 9)^(n+1) = " << sim->fR_params[2] * sim->fR_params[3] / sim->fR_params[1]/sim->fR_params[1] / pow( 12. / cosmo.Omega_m - 9., sim->fR_params[2] + 1.) << " (should be << 1)\n";
	}
	else if(sim->fR_type == FR_TYPE_RN)
	{
		COUT << " f(R) model: a*R^n, with a = " << sim->fR_params[0] << " and n = " << sim->fR_params[1] << endl;
	}
	else if(sim->fR_type == FR_TYPE_R2)
	{
		// In R + alpha * R^2, alpha is expressed in units of (inverse) 8piG * density today (= 1.)
		sim->fR_params[0] /= 2. * fourpiG;
		COUT << " f(R) model: a*R^2, with a = " << sim->fR_params[0] << endl;
	}
	else if(sim->fR_type == FR_TYPE_DELTA) //  TODO: Fix for delta < 0: fRR < 0, so everything blows up!
	{
		double temp = sim->fR_params[0];
		sim->fR_params[0] *= 2 * fourpiG * cosmo.Omega_m;
		COUT << " f(R) model: a * (R/a)^(1+delta) - R" << endl
				 << " with    a = " << temp << " * 8piG * rho_{m0} = " << sim->fR_params[0] << endl
				 << "         delta = " << sim->fR_params[1] << endl;
	}
}


/////////////////////////////////////////////////
// Check (scalar) field:
// Max = max(fabs(field))
// |hom| = sum(field)/number_of_points
// |avg| average of fabs(field), not fabs(average)!
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
double check_field(Field<FieldType> & field, string field_name, long n3, string message = "") // TODO: correct lattice size
{
  Site x(field.lattice());
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);

	double max = 0.,
				 hom = 0.,
				 sum = 0.,
				 temp;

  for(x.first(); x.test(); x.next())
  {
    temp = field(x);
    hom += temp;
    sum += fabs(temp);
    if(fabs(temp) >= max)
    {
      max = fabs(temp);
    }
  }
  parallel.max(max);
  parallel.sum(sum);
  parallel.sum(hom);
  sum /= n3;
  hom /= n3;

	int prec = 16;
	COUT << scientific << setprecision(prec);
  COUT << message
       << setw(17)   << field_name
       << "  Max = " << setw(prec + 3) << max
       << "  hom = " << setw(prec + 3) << hom
       << endl;

	std::cout.copyfmt(oldState);

	return max;
}


/////////////////////////////////////////////////
// Checks vector field, each component separately
// Quantities are the same as check_field
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void check_vector_field(Field<FieldType> & field, string field_name, long n3, string message = "") // TODO: correct lattice size
{
  Site x(field.lattice());
	std::ios oldState(nullptr);
	oldState.copyfmt(std::cout);
  double max[3] = {0.,0.,0.},
				 hom[3] = {0.,0.,0.},
				 temp;
  int i;

  COUT << message;
	COUT << scientific << setprecision(6);

  for(i=0; i<3; i++)
  {
    for(x.first(); x.test(); x.next())
    {
      temp = field(x,i);
      hom[i] += temp;
      if(fabs(temp) > max[i])
      {
        max[i] = fabs(temp);
      }
    }
  parallel.max(max[i]);
	hom[i] /= n3;

  COUT << setw(14)    << field_name << "[" << i << "]"
       << "  Max = " << setw(9) << max[i]
       << "  hom = " << setw(9) << hom[i]
       << endl;
  }
	std::cout.copyfmt(oldState);

	return;
}


// Compute max of fRR
// TODO Add description
template <class FieldType>
double compute_max_fRR(Field<FieldType> &deltaR, double Rbar, const metadata & sim)
{
	if(sim.fR_type == FR_TYPE_R2)// In R + R^2 gravity, fRR is a constant
	{
		return 2. * sim.fR_params[0];
	}

	Site x(deltaR.lattice());
	double temp,
				 max = 0.;

	for(x.first(); x.test(); x.next())
	{
		temp = fabs(fRR(deltaR(x) + Rbar, sim, 555));
		if(temp > max) max = temp;
	}
	if(!max || std::isnan(max))
	{
		cout << "Something went wrong in computing max_fRR. Closing...\n";
		parallel.abortForce();
	}
	else
	{
		return max;
	}
}

/////////////////////////////////////////////////
// Flips (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void flip_fields(Field<FieldType> & field1, Field<FieldType> & field2)
{
	Site x(field1.lattice());
	FieldType temp;

	for(x.first(); x.test(); x.next())
	{
		temp = field1(x);
		field1(x) = field2(x);
		field2(x) = temp;
	}
	return;
}

// Writes field1 --> field2, field2 --> field3, field3 --> field1
template <class FieldType>
void flip_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & field3)
{
	Site x(field1.lattice());
	FieldType temp;

	for(x.first(); x.test(); x.next())
	{
		temp = field3(x);
		field3(x) = field2(x);
		field2(x) = field1(x);
		field1(x) = temp;
	}
	return;
}


/////////////////////////////////////////////////
// Sets whole field to zero
/////////////////////////////////////////////////
template <class FieldType>
void zero_field(Field<FieldType> & field)
{
	Site x(field.lattice());
	for(x.first(); x.test(); x.next())
	{
		field(x) = 0.;
	}
}

/////////////////////////////////////////////////
// Copies (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void copy_field(Field<FieldType> & source, Field<FieldType> & destination, double coeff = 1.)
{
	Site x(source.lattice());

	if(coeff == 1.)
	{
		for(x.first(); x.test(); x.next())
		{
			destination(x) = source(x);
		}
	}
	else
	{
		for(x.first(); x.test(); x.next())
		{
			destination(x) = coeff * source(x);
		}
	}
	return;
}

/////////////////////////////////////////////////
// Sums (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void add_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) + field2(x);
	}
	return;
}

template <class FieldType>
void subtract_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) - field2(x);
	}
	return;
}

template <class FieldType>
void add_fields(Field<FieldType> & field1, double coeff1, Field<FieldType> & field2, double coeff2, Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = coeff1 * field1(x) + coeff2 * field2(x);
	}
	return;
}

template <class FieldType>
void scatter_field(Field<FieldType> & field1, double c1, Field<FieldType> & field2, double c2, Field<FieldType> & result, double k)
{
	Site x(field1.lattice());
	double r;

	for(x.first(); x.test(); x.next())
	{
		r = (double) rand()/RAND_MAX;
		r *= k;
		result(x) = c1 * (1. - r) * field1(x) + c2 * r * field2(x);
	}
	return;
}

/////////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void leapfrog_dotR(Field<FieldType> & deltaR, Field<FieldType> & deltaT, Field<FieldType> & dot_deltaR_old, Field<FieldType> & dot_deltaR_new, double Hubble, double coeff, double dtau_old, double dtau, double dx2)
{
	Site x(deltaR.lattice());
	double temp,
				 coeff1 = (dtau_old + dtau) / 2.;
	coeff1 *= 1. - coeff1 * Hubble;

	for(x.first(); x.test(); x.next())
	{
		temp = deltaR(x+0) + deltaR(x-0) + deltaR(x+1) + deltaR(x-1) + deltaR(x+2) + deltaR(x-2) - 6. * deltaR(x);
		temp /= dx2;
		temp -= coeff * (deltaR(x) + deltaT(x)) + 2. * Hubble * dot_deltaR_old(x);
		dot_deltaR_new(x) = dot_deltaR_old(x) + temp * coeff1;
	}
	return;
}

/////////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void leapfrog_dotxi(Field<FieldType> & xi, Field<FieldType> & zeta, Field<FieldType> & dot_xi_old, Field<FieldType> & dot_xi_new, double Hubble, double a2, double dtau_old, double dtau, double dx2)
{
	Site x(xi.lattice());
	double temp,
				 coeff1 = (dtau_old + dtau) / 2.;
	coeff1 *= 1. - coeff1 * Hubble;

	for(x.first(); x.test(); x.next())
	{
		temp = xi(x+0) + xi(x-0) + xi(x+1) + xi(x-1) + xi(x+2) + xi(x-2) - 6. * xi(x);
		temp /= dx2;
		temp -= a2 * zeta(x) / 3. + 2. * Hubble + dot_xi_old(x);
		dot_xi_new(x) = dot_xi_old(x) + temp * coeff1;
	}
	return;
}

/////////////////////////////////////////////////
// Converts deltaR to xi
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void convert_deltaR_to_xi(Field<FieldType> & xi, Field<FieldType> & deltaR, double const Rbar, double const fRbar, const metadata & sim)
{
	Site x(xi.lattice());

	if(sim.fR_type == FR_TYPE_R2)
	{
		for(x.first(); x.test(); x.next())
		{
			xi(x) = 2. * sim.fR_params[0] * deltaR(x);
		}
	}
	else
	{
		for(x.first(); x.test(); x.next())
		{
			xi(x) = fR(Rbar + deltaR(x), sim, 832) - fRbar;
		}
	}
	return;
}


/////////////////////////////////////////////////
// Converts deltaR to u
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void convert_deltaR_to_u(Field<FieldType> & u, Field<FieldType> & deltaR, double const Rbar, double const fRbar, const metadata & sim)
{
	Site x(u.lattice());

	for(x.first(); x.test(); x.next())
	{
		u(x) = log(fR(Rbar + deltaR(x), sim, 833)/fRbar);
	}
	return;
}


/////////////////////////////////////////////////
// Converts u to xi
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void convert_u_to_xi(Field<FieldType> & u, Field<FieldType> & xi, double const fRbar)
{
	Site x(u.lattice());

	for(x.first(); x.test(); x.next())
	{
		xi(x) = fRbar * (exp(u(x)) - 1.);
	}
	return;
}


/////////////////////////////////////////////////
// Converts xi to u
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void convert_xi_to_u(Field<FieldType> & xi, Field<FieldType> & u, double const fRbar)
{
	Site x(u.lattice());

	for(x.first(); x.test(); x.next())
	{
		u(x) = log(xi(x)/fRbar + 1.);
	}
	return;
}


/////////////////////////////////////////////////
// Initial conditions for xi_prev
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void xi_prev_initial_conditions(Field<FieldType> & xi_prev, Field<FieldType> & xi)
{
  Site x(xi_prev.lattice());

	for(x.first(); x.test(); x.next())
	{
		xi_prev(x) = xi(x); // TODO: first approximation, try different initial conditions
	}
}

/////////////////////////////////////////////////
// Computes zeta and deltaR
//
// Returns maximum value of fRR over the entire lattice.
// the typical oscillation frenquency will be of order a/(sqrt(3*fRR_max))
/////////////////////////////////////////////////
template <class FieldType>
double convert_xi_to_deltaR(Field<FieldType> & eightpiG_deltaT,
                      	 		Field<FieldType> & deltaR,
									    	 		Field<FieldType> & xi,
									    	 		const double Rbar,
								 	    	 		const double fRbar,
								 	    	 		const metadata & sim)
{
  // Using Newton-Raphson Method
  Site x(xi.lattice());
  double R_temp,
				 temp,
				 max_fRR = 0.,
				 fRR_temp;
  int count,
			count_n = 1,
			my_check = 0;

	if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double correction = 0.99;
		double m2 = sim.fR_params[0],
					 c2 = sim.fR_params[1],
					  n = sim.fR_params[2],
					 c1 = sim.fR_params[3];

		if(n == 1.) // xi <-> R relation is invertible algebraically
		{
			for(x.first(); x.test(); x.next())
			{
				R_temp = -c1 / (xi(x) + fRbar);
				if(R_temp <= 1.)
				{
					xi(x) = - correction * fRbar;
					R_temp = -c1 / (1. - correction) / fRbar;
				}
				R_temp = m2 / c2 * ( sqrt(R_temp) - 1.);
				fRR_temp = fabs(fRR(R_temp, sim, 5911));
				if(max_fRR < fRR_temp)
				{
					max_fRR = fRR_temp;
				}
				deltaR(x) = R_temp - Rbar;
			}
		}
		else
		{
			for(x.first(); x.test(); x.next())
			{
				// The initial guess for R is:
				// 1) found from a large-R expansion of FR for R >> m2, or
				// 2) equal to the new GR solution (Rbar - eightpiG_deltaT)
				R_temp = Rbar - eightpiG_deltaT(x);
				if(R_temp > 10. * m2)
				{
					R_temp = -c1 * n / c2 / c2 / (xi(x) + fRbar);
					R_temp = m2 * pow( R_temp, 1./(n+1.));
				}

				count = 0;
				while(true)
				{
					temp = fabs(xi(x)/(fR(R_temp, sim, 56) - fRbar) - 1.);
					if(temp < sim.fR_target_precision) break;
					fRR_temp = fRR(R_temp, sim, 57);

					if(fRR_temp > 0 && fabs(R_temp) <= 1.E+20)
					{
						R_temp += (xi(x) - fR(R_temp, sim, 58) + fRbar) / fRR_temp;
					}
					else
					{
						R_temp = Rbar - eightpiG_deltaT(x);
						R_temp += ((double) rand() / RAND_MAX) * 0.1 * R_temp;
					}

					count ++;
					if(count > count_n * sim.fR_count_max)
					{
						cout << "Could only reach a precision of " << temp << " in " << count_n * sim.fR_count_max << " steps.\n";
						count_n++;
					}
				}

				fRR_temp = fabs(fRR(R_temp, sim, 5912));
				if(max_fRR < fRR_temp)
				{
					max_fRR = fRR_temp;
				}
				deltaR(x) = R_temp - Rbar;
			}
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		for(x.first(); x.test(); x.next())
		{
			R_temp = pow(Rbar, sim.fR_params[1] - 1.) + xi(x)/sim.fR_params[0]/sim.fR_params[1];
			R_temp = pow(R_temp, 1./(sim.fR_params[1] - 1.));
			fRR_temp = fabs(fRR(R_temp, sim, 60));
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;
			deltaR(x) = R_temp - Rbar;
		}
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		cout << "Something is wrong. R + R^2 shouldn not call this function. Closing...\n";
		parallel.abortForce();
	}
	else // For other f(R) models TODO check!!
	{
		for(x.first(); x.test(); x.next())
		{
			R_temp = Rbar - eightpiG_deltaT(x);
			count = 0;
			while(true)
			{
				temp = fabs(xi(x)/(fR(R_temp, sim, 61) - fRbar) - 1.);
				if(temp < sim.fR_target_precision) break;
				fRR_temp = fRR(R_temp, sim, 62);
				R_temp += fRR_temp > 0 ? (xi(x) - fR(R_temp, sim, 63) + fRbar) / fRR_temp : .01 * R_temp; // Displace R_temp slightly, if fRR(R_temp) is "bad"
				count ++;
				if(count > sim.fR_count_max)
				{
					cout << "Could only reach a precision of " << temp << " in " << sim.fR_count_max << " steps.\n";
					break;
				}
			}
			fRR_temp = fabs(fRR(R_temp, sim, 64));
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;
			deltaR(x) = R_temp - Rbar;
		}
	}

	parallel.max(max_fRR);
  if(max_fRR and !std::isnan(max_fRR))
	{
		return max_fRR;
	}
  else
  {
    COUT << " Warning: returning fRRbar\n";
    return fRR(Rbar, sim, 65);
  }
}



// Convert u = log(xi/fRbar + 1.)
template <class FieldType>
double convert_u_to_deltaR(Field<FieldType> & eightpiG_deltaT,
                      	 	 Field<FieldType> & deltaR,
									    	 	 Field<FieldType> & u,
									    	 	 const double Rbar,
								 	    	 	 const double fRbar,
								 	    	 	 const metadata & sim)
{
  // Using Newton-Raphson Method
  Site x(u.lattice());
  double R_temp,
			   temp,
				 max_fRR = 0.,
				 fRR_temp;
  int count,
			count_n = 1,
			my_check = 0;

	if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2 = sim.fR_params[0],
					 c2 = sim.fR_params[1],
					  n = sim.fR_params[2],
					 c1 = sim.fR_params[3];

		if(n == 1.) // u <-> R relation is invertible algebraically
		{
			for(x.first(); x.test(); x.next())
			{
				R_temp = m2 / c2 * ( sqrt(-c1 * exp(-u(x)) / fRbar ) - 1.);
				fRR_temp = fRR(R_temp, sim, 5913);
				if(std::isnan(fRR_temp) || !fRR_temp)
				{
					cout << " u(x) = " << u(x) << "\n R_temp = " << R_temp << "\n fRR_temp = " << fRR_temp << endl;
					hold();
				}

				fRR_temp = fabs(fRR_temp);
				if(max_fRR < fRR_temp)
				{
					max_fRR = fRR_temp;
				}
				deltaR(x) = R_temp - Rbar;
			}
		}
		else
		{
			for(x.first(); x.test(); x.next())
			{
				u(x) = fRbar * (exp(u(x)) - 1.); // Convert xi -> u
				// The initial guess for R is:
				// 1) found from a large-R expansion of FR for R >> m2, or
				// 2) equal to the new GR solution (Rbar - eightpiG_deltaT)
				R_temp = Rbar - eightpiG_deltaT(x);
				if(R_temp > 10. * m2)
				{
					R_temp = -c1 * n / c2 / c2 / (u(x) + fRbar);
					R_temp = m2 * pow( R_temp, 1./(n+1.));
				}

				count = 0;
				while(true)
				{
					temp = fabs(u(x)/(fR(R_temp, sim, 56) - fRbar) - 1.);
					if(temp < sim.fR_target_precision) break;
					fRR_temp = fRR(R_temp, sim, 57);

					if(fRR_temp > 0 && fabs(R_temp) <= 1.E+20)
					{
						R_temp += (u(x) - fR(R_temp, sim, 58) + fRbar) / fRR_temp;
					}
					else
					{
						R_temp = Rbar - eightpiG_deltaT(x);
						R_temp += ((double) rand() / RAND_MAX) * 0.1 * R_temp;
					}

					count ++;
					if(count > count_n * sim.fR_count_max)
					{
						cout << "Could only reach a precision of " << temp << " in " << count_n * sim.fR_count_max << " steps.\n";
						count_n++;
					}
				}

				fRR_temp = fabs(fRR(R_temp, sim, 5914));
				if(max_fRR < fRR_temp)
				{
					max_fRR = fRR_temp;
				}
				deltaR(x) = R_temp - Rbar;
				u(x) = log(u(x)/fRbar + 1.); // Convert back u -> xi
			}
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		for(x.first(); x.test(); x.next())
		{
			u(x) = fRbar * (exp(u(x)) - 1.); // Convert xi -> u
			R_temp = pow(Rbar, sim.fR_params[1] - 1.) + u(x)/sim.fR_params[0]/sim.fR_params[1];
			R_temp = pow(R_temp, 1./(sim.fR_params[1] - 1.));
			fRR_temp = fabs(fRR(R_temp, sim, 60));
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;
			deltaR(x) = R_temp - Rbar;
			u(x) = log(u(x)/fRbar + 1.); // Convert back u -> xi
		}
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		cout << "Something is wrong. R + R^2 shouldn not call this function. Closing...\n";
		parallel.abortForce();
	}
	else // For other f(R) models TODO check!!
	{
		for(x.first(); x.test(); x.next())
		{
			u(x) = fRbar * (exp(u(x)) - 1.); // Convert xi -> u
			R_temp = Rbar - eightpiG_deltaT(x);
			count = 0;
			while(true)
			{
				temp = fabs(u(x)/(fR(R_temp, sim, 61) - fRbar) - 1.);
				if(temp < sim.fR_target_precision) break;
				fRR_temp = fRR(R_temp, sim, 62);
				R_temp += fRR_temp > 0 ? (u(x) - fR(R_temp, sim, 63) + fRbar) / fRR_temp : .01 * R_temp; // Displace R_temp slightly, if fRR(R_temp) is "bad"
				count ++;
				if(count > sim.fR_count_max)
				{
					cout << "Could only reach a precision of " << temp << " in " << sim.fR_count_max << " steps.\n";
					break;
				}
			}
			fRR_temp = fabs(fRR(R_temp, sim, 64));
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;
			deltaR(x) = R_temp - Rbar;
			u(x) = log(u(x)/fRbar + 1.); // Convert back u -> xi
		}
	}

	parallel.max(max_fRR);
  if(max_fRR and !std::isnan(max_fRR))
	{
		return max_fRR;
	}
  else
  {
    COUT << " Warning: returning fRRbar\n";
    return fRR(Rbar, sim, 65);
  }
}


// Builds laplacian from field
template <class FieldType>
void build_laplacian(Field<FieldType> & field, Field<FieldType> & laplace, double dx)
{
	double dx2 = dx*dx;
	Site x(field.lattice());

	for(x.first(); x.test(); x.next())
	{
		laplace(x) = field(x+0) + field(x-0) + field(x+1) + field(x-1) + field(x+2) + field(x-2) - 6. * field(x);
		laplace(x) /= dx2;
	}
	return;
}


// TODO: Maybe put this in gevolution.hpp?
template <class FieldType>
void prepareFTsource_leapfrog_R2(Field<FieldType> & eightpiG_deltaT, Field<FieldType> & source, const double dx)
{
	double dx2 = dx*dx;
	Site x(source.lattice());

	for(x.first(); x.test(); x.next())
	{
		source(x) = eightpiG_deltaT(x+0) + eightpiG_deltaT(x-0) + eightpiG_deltaT(x+1) + eightpiG_deltaT(x-1) + eightpiG_deltaT(x+2) + eightpiG_deltaT(x-2) - 6. * eightpiG_deltaT(x);
		source(x) /= dx2;
	}
	return;
}



#endif
