//////////////////////////
// fofR_tools.hpp
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
	sim->fofR_params[3] = sim->fofR_params[1] * 6. * (1 - cosmo.Omega_m - cosmo.Omega_rad) / cosmo.Omega_m / sim->fofR_params[0]; //Omega_m > 0, checked in parseMetadata()
	sim->fofR_params[0] *= fourpiG * cosmo.Omega_m / 1.5;
}


///////////////////// F(R) FUNCTIONS /////////////////////

// F(R) = f(R) - R
Real F(double R, double * params, const int fofR_type, int code = -1)
{
	double output, Rpow;

	if(fofR_type == FOFR_TYPE_RN)
	{
		Rpow = pow(R,params[1]);
		output = params[0] * Rpow;
	}
	else if(fofR_type == FOFR_TYPE_R2)
	{
		output = params[0] * R * R;
	}
	else if(fofR_type == FOFR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = params[0];
		c2 = params[1];
		n  = params[2];
		c1 = params[3];
		Rpow = pow(R/m2, n);
		output = - m2 * c1 * Rpow / (1. + c2 * Rpow);
	}
	else
	{
		COUT << " Something went wrong when computing F(R) (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(!std::isnan(output)) return output;
	else
	{
		COUT << " F(R) Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", F(R) = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}



// F_R(R) -- First derivative with respect to R
Real FR(double R, double * params, const int fofR_type, int code = -1)
{
	double output, Rpow;

	if(fofR_type == FOFR_TYPE_RN)
	{
		Rpow = pow(R,params[1] - 1.);
		output = params[0] * params[1] * Rpow;
	}
	else if(fofR_type == FOFR_TYPE_R2)
	{
		output = 2. * params[0] * R;
	}
	else if(fofR_type == FOFR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = params[0];
		c2 = params[1];
		n  = params[2];
		c1 = params[3];
		Rpow = pow(R/m2, n);
		output = - c1 * n * pow(R/m2, n-1) / ( (1. + c2*Rpow) * (1. + c2*Rpow) );
	}
	else
	{
		COUT << " Something went wrong when computing FR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(!std::isnan(output))
	{
		return output;
	}
	else
	{
		COUT << " FR Evaluated to NaN (code " << code << ").\n"
				 << " params = " << params[0] << " " << params[1] << " " << params[2] << " " << params[3] << endl
		     << " R = " << R << ", Rpow = " << Rpow << ", FR = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}


// F_{RR}(R) -- Second derivative with respect to R
Real FRR(double R, double * params, const int fofR_type, int code = -1)
{
	double output, Roverm2, Rpow;

	if(fofR_type == FOFR_TYPE_RN)
	{
		Rpow = pow(R, params[1] - 2.);
		output = params[0] * params[1] * (params[1] - 1.) * Rpow;
	}
	else if(fofR_type == FOFR_TYPE_R2)
	{
		output = 2. * params[0];
	}
	else if(fofR_type == FOFR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = params[0];
		c2 = params[1];
		n  = params[2];
		c1 = params[3];
		Roverm2 = R/m2;
		Rpow = pow(Roverm2, n);
		if(n == 1) output = 2. * c1 * c2 / ( m2 * (1. + c2*Roverm2) * (1. + c2*Roverm2) * (1. + c2*Roverm2) );
		else
		{
			output = c1 * n * Rpow * (1. - n + c2*(1. + n)*Rpow ) / ( Roverm2 * Roverm2 * m2 * (1. + c2*Rpow) * (1. + c2*Rpow) * (1. + c2*Rpow) );
		}
	}
	else
	{
		COUT << " Something went wrong when computing FRR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(output && !std::isnan(output))
	{
		return output;
	}
	else if(!output)
	{
		COUT << " FRR Evaluated to 0 (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", FRR = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
	else if(std::isnan(output))
	{
		COUT << " FRR Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", FRR = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}



// F_{RRR}(R) -- Third derivative with respect to R
Real FRRR(double R, double * params, const int fofR_type, int code = -1)
{
	double output, Roverm2, Rpow;

	if(fofR_type == FOFR_TYPE_RN)
	{
		Rpow = pow(R, params[1] - 3.);
		output = params[0] * params[1] * (params[1] - 1.) * (params[1] - 2.) * Rpow;
	}
	else if(fofR_type == FOFR_TYPE_R2)
	{
		output = 0.;
	}
	else if(fofR_type == FOFR_TYPE_HU_SAWICKI)
	{
		double m2, c1, c2, n;
		m2 = params[0];
		c2 = params[1];
		n  = params[2];
		c1 = params[3];
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
	else
	{
		COUT << " Something went wrong when computing FRRR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

	if(!std::isnan(output))
	{
		return output;
	}
	else if(std::isnan(output))
	{
		COUT << " FRRR Evaluated to NaN (code " << code << ").\n"
		     << " R = " << R << ", Rpow = " << Rpow << ", FRRR = " << output << "\n Closing...\n";
		parallel.abortForce();
	}
}



// Print out F(R) model details
void fR_details(const cosmology cosmo, metadata * sim, const double fourpiG)
{
	if(sim->fofR_type == FOFR_TYPE_HU_SAWICKI)
	{
		double temp = sim->fofR_params[0];
		rescale_params_Hu_Sawicki(cosmo, fourpiG, sim);
		COUT << " f(R) model: Hu-Sawicki    F(R) = - m2 * c1 * pow(R/m2, n) / (1. + c2 * pow(R/m2, n))" << endl
				 << "  with m2 = " << temp << " * 8piG * rho_{m0} / 3. = " << sim->fofR_params[0] << endl;
		COUT << "       c1 = " << sim->fofR_params[3] << endl
				 << "       c2 = " << sim->fofR_params[1] << endl
				 << "       n  = " << sim->fofR_params[2] << endl
				 << " so |fR0| ~ -n*c1/c2^2/(12/Omega_m - 9)^(n+1) = " << sim->fofR_params[2] * sim->fofR_params[3] / sim->fofR_params[1]/sim->fofR_params[1] / pow( 12. / cosmo.Omega_m - 9., sim->fofR_params[2] + 1.) << "   (should be << 1)\n";
	}
	else if(sim->fofR_type == FOFR_TYPE_RN || sim->fofR_type == FOFR_TYPE_R2)
	{
		COUT << " f(R) model: a*R^n, with a = " << sim->fofR_params[0] << " and n = " << sim->fofR_params[1] << endl;
	}
}

void compute_dtrace_eq(double a, double H, double R, double T, double fourpiG, cosmology * cosmo, metadata * sim)
{
	COUT << " T^mu_mu (hom.) = " << T << "  model = " << - cosmo->Omega_m/a/a/a << endl;
}


/////////////////////////////////////////////////
// Check (scalar) field:
// Max = max(fabs(field))
// |hom| = sum(field)/number_of_points
// |avg| average of fabs(field), not fabs(average)!
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
Site check_field(Field<FieldType> & field, string field_name, string message = "", int n = 64) // TODO: correct lattice size
{
  Site x(field.lattice());
  Site y, z[4];
  double max = 0., hom = 0., sum = 0., temp;
  for(x.first(); x.test(); x.next())
  {
    temp = field(x);
    hom += temp;
    sum += fabs(temp);
    if(fabs(temp) >= max)
    {
      max = fabs(temp);
      z[parallel.rank()] = x;
    }
  }
  parallel.max(max);
  parallel.sum(sum);
  parallel.sum(hom);
  sum /= n * n * n;
  hom /= n * n * n;

  y = z[parallel.rank()];

  COUT << message
  // MPI_Barrier(MPI_COMM_WORLD);
       << " " << field_name
       << ",  Max = " << max
       << ",  hom = " << hom
      //  << ",  |avg| = " << sum
      //  << ",  val(x_max) = " << field(y)
      //  << ", rank = " << parallel.rank()
       << endl;
  return y;
}


/////////////////////////////////////////////////
// Checks vector field, each component separately
// Quantities are the same as check_field
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
Site check_vector_field(Field<FieldType> & field, string field_name, string message = "", int n = 64) // TODO: correct lattice size
{
  Site x(field.lattice()), y;
  double max[3] = {0.,0.,0.}, hom[3] = {0.,0.,0.}, sum[3] = {0.,0.,0.}, temp;
  int i;

  COUT << message;

  for(i=0; i<3; i++)
  {
    for(x.first(); x.test(); x.next())
    {
      temp = field(x,i);
      hom[i] += temp;
      sum[i] += fabs(temp);
      if(fabs(temp) > max[i])
      {
        max[i] = fabs(temp);
        y = x;
      }
    }
  parallel.max(max[i]);
  parallel.sum(sum[i]);
  parallel.sum(hom[i]);
  sum[i] = sum[i] / n / n / n;
  // MPI_Barrier(MPI_COMM_WORLD);
  COUT << "  Field: " << field_name << "[" << i << "]"
       << ",  Max = " << max[i]
       << ",  hom = " << hom[i]
       << ",  |avg| = " << sum[i]
       << ",  val(x_max) = " << field(y)
      //  << ", rank = " << parallel.rank()
       << endl;
  }
  return y;
}


// Compute max of FRR
// TODO Add description

template <class FieldType>
double compute_max_FRR(Field<FieldType> &deltaR, double Rbar, double * params, int type)
{
	Site x(deltaR.lattice());
	double temp, max = 0.;
	for(x.first(); x.test(); x.next())
	{
		temp = fabs(FRR(deltaR(x) + Rbar, params, type, 555));
		if(temp > max) max = temp;
	}
	if(!max || std::isnan(max))
	{
		COUT << "Something went wrong in computing max_FRR. Closing...\n";
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
// Copies (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void copy_field(Field<FieldType> & source, Field<FieldType> & destination)
{
	Site x(source.lattice());
	for(x.first(); x.test(); x.next())
	{
		destination(x) = source(x);
	}
	return;
}

template <class FieldType>
void copy_field(Field<FieldType> & source, Field<FieldType> & destination, double coeff)
{
	Site x(source.lattice());
	for(x.first(); x.test(); x.next())
	{
		destination(x) = coeff * source(x);
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
void add_fields(Field<FieldType> & field1, Field<FieldType> & field2, Field<FieldType> & result, double coeff1, double coeff2)
{
	Site x(field1.lattice());
	for(x.first(); x.test(); x.next())
	{
		result(x) = coeff1 * field1(x) + coeff2 * field2(x);
	}
	return;
}


/////////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void deltaR_alternative(Field<FieldType> & R_old, Field<FieldType> & R, Field<FieldType> & R_new, double Hubble, double m2, double dtau_old, double dtau, double dx2)
{
	Site x(R.lattice());
	double temp;
	for(x.first(); x.test(); x.next())
	{
		R_new(x) = R(x+0) + R(x-0) + R(x+1) + R(x-1) + R(x+2) + R(x-2) - 6. * R(x);
		R_new(x) /= dx2;
		R_new(x) += ( 2./dtau_old/dtau - m2 + 2. * Hubble * (1./dtau - 1./dtau_old) ) * R(x);
		R_new(x) *= (dtau_old + dtau) * dtau / 2.;
		R_new(x) -= dtau * (1. - Hubble * dtau) * R_old(x) / dtau_old;
		R_new(x) /= 1. + Hubble * dtau_old;
	}
	return;
}


/////////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void lp_update_dotR(Field<FieldType> & deltaR, Field<FieldType> & deltaT, Field<FieldType> & dot_deltaR_old, Field<FieldType> & dot_deltaR_new, double Hubble, double coeff, double dtau_old, double dtau, double dx2)
{
	Site x(deltaR.lattice());
	double temp, coeff1, coeff2, coeff3;
	coeff1 = 1. + Hubble * dtau_old;
	coeff2 = 1. - Hubble * dtau;
	coeff3 = (dtau + dtau_old) / 2.;
	for(x.first(); x.test(); x.next())
	{
		temp = deltaR(x+0) + deltaR(x-0) + deltaR(x+1) + deltaR(x-1) + deltaR(x+2) + deltaR(x-2) - 6. * deltaR(x);
		temp /= dx2;
		temp -= coeff * (deltaR(x) + deltaT(x));
		temp *= coeff3;
		dot_deltaR_new(x) = (coeff2 * dot_deltaR_old(x) + temp) / coeff1;
	}
	return;
}


/////////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void lp_update_dotxi(Field<FieldType> & xi, Field<FieldType> & zeta, Field<FieldType> & dot_xi_old, Field<FieldType> & dot_xi_new, double Hubble, double a2, double dtau_old, double dtau, double dx2)
{
	Site x(xi.lattice());
	double temp, coeff1, coeff2, coeff3;
	coeff1 = 1. + Hubble * dtau_old;
	coeff2 = 1. - Hubble * dtau;
	coeff3 = (dtau + dtau_old)/6.;
	for(x.first(); x.test(); x.next())
	{
		temp = xi(x+0) + xi(x-0) + xi(x+1) + xi(x-1) + xi(x+2) + xi(x-2) - 6. * xi(x);
		temp *= 3. / dx2;
		temp -= a2 * zeta(x);
		temp *= coeff3;
		temp += coeff2 * dot_xi_old(x);
		dot_xi_new(x) = temp / coeff1;
	}
	return;
}


/////////////////////////////////////////////////
// Computes ddot_deltaR
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void laplace_ddot_deltaR_set(Field<FieldType> & R_old, Field<FieldType> & R, Field<FieldType> & R_new, Field<FieldType> & laplace_deltaR, Field<FieldType> & ddot_R, double dt_old, double dt_new, double dx2)
{
	Site x(R.lattice());
	for(x.first(); x.test(); x.next())
	{
		laplace_deltaR(x) = (R_new(x+0) + R_new(x-0) + R_new(x+1) + R_new(x-1) + R_new(x+2) + R_new(x-2) - 6. * R_new(x)) / dx2;
	}
	if(dt_old * dt_new)
		for(x.first(); x.test(); x.next())
		{
			ddot_R(x) = 2./(dt_old + dt_new) * (R_new(x)/dt_new + R_old(x)/dt_old) - 2. * R(x)/(dt_old*dt_new);
		}
	return;
}


template <class FieldType>
void laplace_ddot_deltaR_set(Field<FieldType> & R, Field<FieldType> & laplace_deltaR, double dx2)
{
	Site x(R.lattice());
	for(x.first(); x.test(); x.next())
	{
		laplace_deltaR(x) = (R(x+0) + R(x-0) + R(x+1) + R(x-1) + R(x+2) + R(x-2) - 6. * R(x)) / dx2;
	}
	return;
}


/////////////////////////////////////////////////
// Initial conditions for xi
// TODO: Add comments here
/////////////////////////////////////////////////
template <class FieldType>
void xi_set(Field<FieldType> & xi, Field<FieldType> & deltaR, double Rbar, double FRbar, double * fofR_params, int fofR_type)
{
	Site x(xi.lattice());
	for(x.first(); x.test(); x.next())
	{
		xi(x) = FR(Rbar + deltaR(x), fofR_params, fofR_type, 832) - FRbar;
	}
	xi.updateHalo();
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
// Returns maximum value of FRR over the entire lattice.
// the typical oscillation frenquency will be of order a/(sqrt(3*FRR_max))
/////////////////////////////////////////////////
template <class FieldType>
double computeDRzeta(Field<FieldType> & eightpiG_deltaT,
                     Field<FieldType> & deltaR,
									   Field<FieldType> & zeta,
									   Field<FieldType> & xi,
									   double Rbar,
								 	   double FRbar,
								 	   double * params,
									   const int fofR_type,
                     double prec)
{
  // Using Newton-Raphson Method
  Site x(xi.lattice());
  double R_temp, temp, max_FRR = 0., FRR_temp;
  int count, count_max = 10000;
  for(x.first(); x.test(); x.next())
  {
    R_temp = Rbar - eightpiG_deltaT(x);
    count = 0;
    while(true)
    {
      temp = fabs(xi(x)/(FR(R_temp, params, fofR_type, 56) - FRbar) - 1.);
      if(temp < prec) break;
      FRR_temp = FRR(R_temp, params, fofR_type, 159);
      R_temp = FRR_temp > 0 ? R_temp + (xi(x) - FR(R_temp, params, fofR_type, 373) + FRbar) / FRR_temp : 1.1 * R_temp; // Displace R_temp slightly, if FRR(R_temp) is "bad"
      count ++;
      if(count > count_max)
      {
        // cout << "Could only reach a precision of " << temp << " in " << count_max << " steps.\n";
        break;
      }
    }
    FRR_temp = fabs(FRR(R_temp, params, fofR_type, 157));
    if(max_FRR < FRR_temp) max_FRR = FRR_temp;
    deltaR(x) = R_temp - Rbar;
    zeta(x) = deltaR(x) + eightpiG_deltaT(x);
  }

  parallel.max(max_FRR);
  if(max_FRR and !std::isnan(max_FRR)) return max_FRR;
  else
  {
    COUT << " Warning: returning FRRbar\n";
    return FRR(Rbar, params, fofR_type, 45);
  }
}



#endif
