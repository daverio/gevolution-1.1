//////////////////////////
// FR_tools.hpp
//////////////////////////
//
// code components related to f(R)
//
// Authors: David Daverio (Cambridge University)
//          Lorenzo Reverberi (CEICO - Czech Academy of Sciences, Prague)
//
// Last modified: February 2016
//
//////////////////////////

using namespace std;

#ifndef FR_TOOLS_HEADER
#define FR_TOOLS_HEADER

void hold(int n = 1)
{
	COUT << n << " - Press ENTER to continue..." << flush;
	MPI_Barrier(MPI_COMM_WORLD);
	cin.ignore(numeric_limits<streamsize>::max(), '\n');
}


////////////////////////// f(R) FUNCTIONS //////////////////////////
//
// Full gravitational Lagrangian is
// L_grav = R + f(R)
//
////////////////////////////// MODELS //////////////////////////////
/////////////// FR_MODEL_RN
// f(R) = a * R^n = R * (R/(b*Ra))^(n-1)
// For convenience pass as arguments: b, n
// where Ra := 8piG * rho_0 (= 8piG in code units)
// then compute a = (b*Ra)^(1-n)
//
/////////////// FR_MODEL_R2
// Sub-class of FR_MODEL_RN with n=2 -- "special" because F'' = const.
// non-linearities disappear from the evolution equation for deltaR
//
/////////////// FR_MODEL_DELTA
// f(R) = R * (R/(a * Ra))^delta - R
// where Ra := 8piG * rho_0 (= 8piG in code units), a = dimensionless constant
// Pass a, delta

/////////////// FR_MODEL_HU_SAWICKI
// Hu-Sawicki Model: arXiv:0705.1158 -- See paper for 2
// m2 := m^2 (comparison between these values and the nomenclature of the original paper)
// f(R) = -m2 * c1 * (R/m2)^n / (1 + c2*(R/m2)^n)
// In settings.ini:
// params[0] = m^2 is given in units of (matter energy density today) * 8piG / 3.
// params[1] = |fR0| (<< 1)
// params[2] = n
// Constructed parameters:
// params[1] = c2 = 16piG (1 - Omega_m - Omega_rad) / |fR0| / m2 / (12/Omega_m - 9)^(n+1)
// params[3] = c1, calibrated to give c1*m2/c2 = 2 * 8piG * rho_Lambda [see 0705.1158 equations (25,26)]
//



/////////////////////////////////////////////
// f(R)
/////////////////////////////////////////////
Real f(
  double R,
  const metadata & sim,
  int code)
{
	double output = 0.,
				 Rpow;


  if(R <= 0.)
  {
    cout << " R = " << R << " in f (code " << code << ")." << endl;
    return FR_WRONG_RETURN;
  }

  if(sim.fR_model == FR_MODEL_R2)
	{
		output = sim.fR_params[0] * R * R;
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];
		output = a * pow(R, n);
	}
	else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
	{
		double
		m2 = sim.fR_params[0],
		c2 = sim.fR_params[1],
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3];

    if(n == 1. || n == 1)
    {
      output = - m2 * c1 / (c2 + m2 / R);
    }
    else
    {
      output = - m2 * c1 / (c2 + 1. / pow(R/m2, n));
    }

    if(isnan(output))
    {
      cout << " f(R) Evaluated to NaN (code " << code << ")" << endl
      << " R = " << R << endl;
      return FR_WRONG_RETURN;
    }
	}
	else if(sim.fR_model == FR_MODEL_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];

		Rpow = pow(R, delta);
		output = R * (a * Rpow - 1.);

		if(isnan(output))
		{
			cout
			<< " f(R) Evaluated to NaN (code " << code << ")" << endl
			<< " R = " << R << ", R^delta = " << Rpow << ", f(R) = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else
	{
		cout << " Something went wrong when computing f(R) (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

  return output;
}

/////////////////////////////////////////////
// fR(R) -- First derivative with respect to R
/////////////////////////////////////////////
Real fR(
  double R,
  const metadata & sim,
  int code)
{
	double output = 0.,
				 Rpow;

  if(R <= 0.)
  {
    cout << " R = " << R << " in fR (code " << code << ")." << endl;
    return FR_WRONG_RETURN;
  }

  if(sim.fR_model == FR_MODEL_R2)
	{
    output = 2. * sim.fR_params[0] * R;

		if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ")" << endl
			<< " R = " << R << ", fR = " << output << endl;
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];

    // Rpow = n > 1. ? pow(R, n - 1.) : pow(R, 1. - n);
		Rpow = pow(R, n - 1.);

		output = a * n * Rpow;

		if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ")" << endl
			     << " R = " << R << ", R^(n-1) = " << Rpow << ", fR = " << output << endl;
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
	{
		double
		m2 = sim.fR_params[0],
		c2 = sim.fR_params[1],
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3];

		if(n == 1 || n == 1.)
		{
			output = (1. + c2 * R / m2);
			output = - c1 / output / output;
		}
		else
		{
      Rpow = pow(R/m2, n);

			output = (1. + c2 * Rpow);
			output = - c1 * n * pow(R/m2, n-1.) / output / output;
		}

		if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ")" << endl;
			cout << " R = " << R << ", (R/m2)^n = " << Rpow << ", fR = " << output << endl;
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];

		Rpow = pow(R, delta);
		output = a * (1. + delta) * Rpow - 1.;

    if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ")" << endl
			     << " R = " << R << ", R^delta = " << Rpow << ", fR = " << output << endl;
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else
	{
		cout << " Something went wrong when computing FR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

  return output;
}

/////////////////////////////////////////////
// F_{RR}(R) -- Second derivative with respect to R
/////////////////////////////////////////////
Real fRR(
  double R,
  const metadata & sim,
  int code)
{
	double output = 0.,
	       R_over_m2,
         Rpow;

  if(R <= 0.)
  {
    cout << " R = " << R << " in fRR (code " << code << ")." << endl;
    return FR_WRONG_RETURN;
  }

  if(sim.fR_model == FR_MODEL_R2)
  {
		output = 2. * sim.fR_params[0];
		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ")" << endl
			<< " R = " << R << ", fRR = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];

    Rpow = pow(R, n) / R / R;

		output = a * n * (n - 1.) * Rpow;

		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ")" << endl
					 << " R = " << R << ", R^(n-2) = " << Rpow << ", fRR = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
	{
		double m2 = sim.fR_params[0],
		       c2 = sim.fR_params[1],
           n  = sim.fR_params[2],
           c1 = sim.fR_params[3];

    R_over_m2 = R/m2;
    Rpow = pow(R_over_m2, n);
		output = (1. + c2 * Rpow);

		if(n == 1. || n == 1)
		{
			output = 2. * c1 * c2 / (m2 * output * output * output);
		}
		else
		{
			if(R_over_m2 * output) output = c1 * n * Rpow * (1. - n + (1. + n) * c2 * Rpow ) / (R_over_m2 * R_over_m2 * m2 * output * output * output);
      else
      {
        cout << " R_over_m2 * output = 0 in fRR [FR_MODEL_HU_SAWICKI] (code " << code << ")" << endl
        << " R = " << R << ", R/m2 = " << R_over_m2 << ", output = " << output << ", n = " << n << endl;
        parallel.abortForce();
      }
		}

		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ")" << endl
			     << " R = " << R << ", (R/m2)^n = " << Rpow << ", fRR = " << output << ", n = " << n << endl;
           parallel.abortForce();
		}
	}
	else if(sim.fR_model == FR_MODEL_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];

		Rpow = pow(R, delta - 1.);
		output = a * delta * (delta + 1.) * Rpow;

		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ")" << endl
			     << " R = " << R << ", R^(delta-1) = " << Rpow << ", fRR = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else
	{
		cout << " Something went wrong when computing fRR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

  return output;
}

/////////////////////////////////////////////
// F_{RRR}(R) -- Third derivative with respect to R
/////////////////////////////////////////////
Real fRRR(
  double R,
  const metadata & sim,
  int code
)
{
	double output = 0.,
				 R_over_m2,
				 Rpow;

  if(R <= 0.)
  {
    cout << " R = " << R << " in fRRR (code " << code << ")." << endl;
    parallel.abortForce();
    return FR_WRONG_RETURN;
  }

  if(sim.fR_model == FR_MODEL_R2)
  {
		return 0.;
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];

    Rpow = pow(R, n) / R / R / R;
		output = a * n * (n - 1.) * (n - 2.) * Rpow;

		if(isnan(output))
		{
			cout << " fRRR Evaluated to NaN (code " << code << ")" << endl
			<< " R = " << R << ", R^(n-3) = " << Rpow << ", fRRR = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
	{
		double
		m2 = sim.fR_params[0],
		c2 = sim.fR_params[1],
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3];

		R_over_m2 = R/m2;
    Rpow = pow(R_over_m2, n);
    output = 1. + c2 * Rpow;

		if(n == 1 || n == 1.)
		{
      output = c2 / m2 / output / output;
			output = -6. * c1 * output * output;
		}
		else if(n == 2 || n == 2.)
		{
      output = 1. / m2 / output / output;
			output = -24. * c1 * c2 * R_over_m2 * (-1. + c2 * Rpow) * output * output;
		}
		else
		{
			output = 1. / m2 / output / output / R_over_m2;
      output = - c1 * n * Rpow * (2. - 3. * n + n * n + (4. - 4. * n * n) * Rpow + (2. + 3. * n + n * n) * Rpow * Rpow) * output * output / R_over_m2;
		}

		if(isnan(output))
		{
			cout << " fRRR Evaluated to NaN (code " << code << ")" << endl
			     << " R = " << R << ", (R/m2)^n = " << Rpow << ", fRRR = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_model == FR_MODEL_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];

		Rpow = pow(R, delta - 2.);
		output = a * delta * (delta * delta - 1.) * Rpow;

		if(isnan(output))
		{
			cout << " fRRR Evaluated to NaN (code " << code << ")" << endl
			<< " R = " << R << ", R(delta-2) = " << Rpow << ", fRRR = " << output << endl;
			parallel.abortForce();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else
	{
		cout << " Something went wrong when computing fRRR (code " << code << "). Closing...\n";
		parallel.abortForce();
	}

  return output;
}


/////////////////////////////////////////////
// allowed values for xi
/////////////////////////////////////////////
bool xi_allowed(
	double xi,
	double Rbar,
	double fRbar,
	const metadata & sim
)
{
	if(isnan(xi) || fabs(xi) >= 1.E20)
	{
		return false;
	}

	if(sim.fR_model == FR_MODEL_RN)
	{
		double xi_min = - fRbar;
		if(xi <= xi_min)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	else if(sim.fR_model == FR_MODEL_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];
		double xi_min = - a * (1. + delta) * pow(Rbar, delta);
		if(xi <= xi_min)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
	{
		double
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3],
		c2 = sim.fR_params[1],
		xi_min;

		if(n == 1 || n == 1.)
		{
			xi_min = - c1 - fRbar;
		}
		else
		{
			xi_min = - (1.+n)*(1.+n) * pow((n-1.) / c2 / (1.+n), 1.-1./n) / n / 4.- fRbar;
		}

		if(xi <= xi_min || xi + fRbar >= 0.)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
}



/////////////////////////////////////////////
// Print out f(R) model details
/////////////////////////////////////////////
void fR_details(
  const cosmology cosmo,
  metadata * sim,
  double fourpiG
)
{
	if(sim->fR_model == FR_MODEL_HU_SAWICKI)
	{
		//TODO: Check whether we need to rescale something by some power of the boxsize
		COUT << " Model: Hu-Sawicki\n f(R) = - m2 * c1 * pow(R/m2, n) / (1. + c2 * pow(R/m2, n))" << endl;
		COUT << " with m2 = " << sim->fR_params[0] << " * 8piG * rho_{m0} / 3. = ";

    sim->fR_params[0] *= fourpiG * cosmo.Omega_m / 1.5;
    sim->fR_params[1] = 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[1] / pow(12./cosmo.Omega_m - 9., sim->fR_params[2] + 1.) / sim->fR_params[0];
    sim->fR_params[3] = sim->fR_params[1] * 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[0];

    COUT << sim->fR_params[0] << endl;
		COUT << "      c1 = " << sim->fR_params[3] << endl;
		COUT << "      c2 = " << sim->fR_params[1] << endl;
		COUT << "       n = " << sim->fR_params[2] << endl;
		COUT << "so |fR0| ~ -n*c1/c2^2/(12/Omega_m - 9)^(n+1) = " << sim->fR_params[2] * sim->fR_params[3] / sim->fR_params[1]/sim->fR_params[1] / pow(12. / cosmo.Omega_m - 9.,
			sim->fR_params[2] + 1.) << " (should be << 1)\n";
	}
	else if(sim->fR_model == FR_MODEL_RN)
	{
		COUT
		<< " f(R) model:" << endl;

    COUT
    << " f(R) = R * (R / (b*8piG*rho0))^(n-1), with b = " << sim->fR_params[0] << endl
    << "                                            n = " << sim->fR_params[1] << endl;

		sim->fR_params[0] = 1. / pow(2. * fourpiG * sim->fR_params[0], sim->fR_params[1] - 1.);

		// if(sim->fR_params[1] > 1.) sim->fR_params[0] = 1. / pow(2. * fourpiG * sim->fR_params[0], sim->fR_params[1] - 1.);
    // else sim->fR_params[0] = pow(2. * fourpiG * sim->fR_params[0], 1. - sim->fR_params[1]);

    COUT << " rescaled as: f(R) = a * R^n, with a = " << sim->fR_params[0] << endl;
	}
	else if(sim->fR_model == FR_MODEL_R2)
	{
		COUT
		<< " f(R) model:" << endl
		<< " f(R) = R^2 / (b*8piG*rho0), with b = " << sim->fR_params[0] << endl;

		sim->fR_params[0] = 0.5 / fourpiG / sim->fR_params[0];

		COUT << " rescaled as: f(R) = a * R^2, with a = " << sim->fR_params[0] << endl;
	}
	else if(sim->fR_model == FR_MODEL_DELTA)
	{
		// Specify fR0 and delta
		double delta = sim->fR_params[1];

		sim->fR_params[0] = (1. + sim->fR_params[0]) / (1. + delta);
		sim->fR_params[0] /= pow(2.*fourpiG*(4.*cosmo.Omega_Lambda + cosmo.Omega_m), delta); // if you initially pass fR0

		double a = sim->fR_params[0];

		COUT << " f(R) model: f(R) = a * (R / (8piG * rho0))^(1+delta) - R" << endl;
		COUT << "              a-1 = " << a-1. << endl;
		COUT << "            delta = " << delta << endl;
	}
}


// TODO CHECK AND REWRITE
void fR_rescale_before_hibernate(
  const cosmology cosmo,
  metadata * sim,
  double fourpiG
)
{
	if(sim->fR_model == FR_MODEL_HU_SAWICKI)
	{
    sim->fR_params[1] = 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[1] / pow(12./cosmo.Omega_m - 9., sim->fR_params[2] + 1.) / sim->fR_params[0];
    sim->fR_params[0] /= fourpiG * cosmo.Omega_m / 1.5;
	}
	else if(sim->fR_model == FR_MODEL_RN)
	{
    if(sim->fR_params[1] > 1.) sim->fR_params[0] = 1. / pow(sim->fR_params[0], 1. / (sim->fR_params[1] - 1.)) / 2. / fourpiG;
    else sim->fR_params[0] = pow(sim->fR_params[0], 1. / (1. - sim->fR_params[1])) / 2. / fourpiG;
	}
	else if(sim->fR_model == FR_MODEL_R2)
	{
    sim->fR_params[0] = 0.5 / fourpiG / sim->fR_params[0];
	}
	else if(sim->fR_model == FR_MODEL_DELTA)
	{
	}
}




/////////////////////////////////////////////
// Output for check_field:
// TODO: Add comments here
/////////////////////////////////////////////
void output_check_field(
  double max,
  double absmax,
  double min,
  double absmin,
  double avg,
  double absavg,
  int prec,
  string field_name,
  string message = "")
{
  int spaces = 20;

  COUT << scientific << setprecision(prec) << setw(spaces) << field_name;
  COUT << "  Max|.| =  " << absmax;
  COUT << "  Min|.| =  " << absmin;
  COUT << "  Avg|.| =  " << absavg;

  COUT << endl << setw(spaces) << " ";

  if(max < 0)
  {
    COUT << "  Max    = " << max;
  }
  else
  {
    COUT << "  Max    =  " << max;
  }

  if(min < 0)
  {
    COUT << "  Min    = " << min;
  }
  else
  {
    COUT << "  Min    =  " << min;
  }

  if(avg < 0)
  {
    COUT << "  Avg    = " << avg;
  }
  else
  {
    COUT << "  Avg    =  " << avg;
  }

  COUT << " " << message << endl;

  return;
}



/////////////////////////////////////////////
// Check (scalar) field:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double check_field(
  Field<FieldType> & field,
  string field_name,
  long n3,
  const metadata & sim,
  string message = "")
{
  Site x(field.lattice());
	ios oldState(nullptr);
	oldState.copyfmt(cout);

	double max = -1.E+30, absmax = 0., min = 1.E+30, absmin = 1.E+30, avg = 0., absavg = 0., temp;

  for(x.first(); x.test(); x.next())
  {
    temp = field(x);
    avg += temp;

    if(temp >= max)
    {
      max = temp;
    }

		if(temp <= min)
		{
			min = temp;
		}

		temp = fabs(temp);
    absavg += temp;

		if(temp <= absmin)
		{
			absmin = temp;
		}

    if(temp >= absmax)
    {
      absmax = temp;
    }
  }

  parallel.max(max);
  parallel.max(absmax);
  parallel.min(min);
  parallel.min(absmin);
  parallel.sum(avg);
  parallel.sum(absavg);
  avg /= n3;
  absavg /= n3;

  output_check_field(max, absmax, min, absmin, avg, absavg, sim.check_fields_precision, field_name, message);
	cout.copyfmt(oldState);
	return absmax;
}


/////////////////////////////////////////////
template <class FieldType>
double check_field(
  Field<FieldType> & field,
  string field_name,
  long n3,
	string message = "",
  int prec = 6
)
{
  Site x(field.lattice());
	ios oldState(nullptr);
	oldState.copyfmt(cout);

	double max = -1.E+30, absmax = 0., min = 1.E+30, absmin = 1.E+30, avg = 0., absavg = 0., temp;

  for(x.first(); x.test(); x.next())
  {
    temp = field(x);
    avg += temp;

    if(temp >= max)
    {
      max = temp;
    }

		if(temp <= min)
		{
			min = temp;
		}

		temp = fabs(temp);
    absavg += temp;

		if(temp <= absmin)
		{
			absmin = temp;
		}

    if(temp >= absmax)
    {
      absmax = temp;
    }
  }

  parallel.max(max);
  parallel.max(absmax);
  parallel.min(min);
  parallel.min(absmin);
  parallel.sum(avg);
  parallel.sum(absavg);
  avg /= n3;
  absavg /= n3;

  output_check_field(max, absmax, min, absmin, avg, absavg, prec, field_name, message);
	cout.copyfmt(oldState);
	return absmax;
}

/////////////////////////////////////////////
// Check correlation of two fields
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double check_correl(
  Field<FieldType> & field1,
  string field_name1,
  Field<FieldType> & field2,
  string field_name2,
  long n3,
  const metadata & sim)
{
	Site x(field1.lattice());
	ios oldState(nullptr);
	oldState.copyfmt(cout);

  double sum1 = 0., sum2 = 0., sq1 = 0., sq2 = 0., sum12 = 0., temp1, temp2, correl_num, correl_den;

  for(x.first(); x.test(); x.next())
  {
    temp1 = field1(x);
    temp2 = field2(x);

    sum1 += temp1;
    sq1 += temp1 * temp1;
    sum2 += temp2;
    sq2 += temp2 * temp2;
    sum12 += temp1 * temp2;
  }

  parallel.sum(sum1);
  parallel.sum(sum2);
  parallel.sum(sum12);
  parallel.sum(sq1);
  parallel.sum(sq2);

  correl_num = sum12 - (2. - 1./n3) * sum1 * sum2 / n3;
  correl_den = sqrt( sq1 - (2. - 1./n3) * sum1 * sum1 / n3 ) * sqrt( sq2 - (2. - 1./n3) * sum2 * sum2 / n3 );

  COUT << " Correlator for " + field_name1 + " and " + field_name2 + " = " << correl_num / correl_den << endl;

	cout.copyfmt(oldState);

	return sum12;
}


/////////////////////////////////////////////
// Check linear combination of two scalar fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double check_field(
  Field<FieldType> & field1,
  double c1,
  Field<FieldType> & field2,
  double c2,
  string field_name,
  long n3,
  const metadata & sim,
  string message = "") // TODO: correct lattice size
{
  Site x(field1.lattice());
	ios oldState(nullptr);
	oldState.copyfmt(cout);

  double max = -1.E+30, absmax = 0., min = 1.E+30, absmin = 1.E+30, avg = 0., absavg = 0., temp;

  for(x.first(); x.test(); x.next())
  {
    temp = c1 * field1(x) + c2 * field2(x);
    avg += temp;

    if(temp >= max)
    {
      max = temp;
    }

		if(temp <= min)
		{
			min = temp;
		}

		temp = fabs(temp);
    absavg += temp;

		if(temp <= absmin)
		{
			absmin = temp;
		}

    if(temp >= absmax)
    {
      absmax = temp;
    }
  }

  parallel.max(max);
  parallel.max(absmax);
  parallel.min(min);
  parallel.min(absmin);
  parallel.sum(avg);
  parallel.sum(absavg);
  avg /= n3;
  absavg /= n3;

  output_check_field(max, absmax, min, absmin, avg, absavg, sim.check_fields_precision, field_name, message);
	cout.copyfmt(oldState);

	return absmax;
}


/////////////////////////////////////////////
// Check linear combination of three scalar fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double check_field(
  Field<FieldType> & field1,
  Field<FieldType> & field2,
  Field<FieldType> & field3,
  string field_name,
  double c1,
  double c2,
  long n3,
  const metadata & sim,
  string message = "") // TODO: correct lattice size
{
  Site x(field1.lattice());
	ios oldState(nullptr);
	oldState.copyfmt(cout);

  double max = -1.E+30, absmax = 0., min = 1.E+30, absmin = 1.E+30, avg = 0., absavg = 0., temp;

  for(x.first(); x.test(); x.next())
  {
    temp = field1(x) + c1 * field2(x) + c2 * field3(x);
    avg += temp;

    if(temp >= max)
    {
      max = temp;
    }

		if(temp <= min)
		{
			min = temp;
		}

		temp = fabs(temp);
    absavg += temp;

		if(temp <= absmin)
		{
			absmin = temp;
		}

    if(temp >= absmax)
    {
      absmax = temp;
    }
  }

  parallel.max(max);
  parallel.max(absmax);
  parallel.min(min);
  parallel.min(absmin);
  parallel.sum(avg);
  parallel.sum(absavg);
  avg /= n3;
  absavg /= n3;

  output_check_field(max, absmax, min, absmin, avg, absavg, sim.check_fields_precision, field_name, message);
	cout.copyfmt(oldState);

	return absmax;
}

/////////////////////////////////////////////
// Checks vector field, each component separately
// Quantities are the same as check_field
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void check_vector_field(
  Field<FieldType> & field,
  string field_name,
  long n3,
  const metadata & sim,
  string message = "") // TODO: correct lattice size
{
  Site x(field.lattice());
	ios oldState(nullptr);
	oldState.copyfmt(cout);

  double max = -1.E+30, absmax = 0., min = 1.E+30, absmin = 1.E+30, avg = 0., absavg = 0., temp;

	for(int i=0; i<3; i++)
  {
    max = -1.E+30;
    absmax = 0.;
    min = 1.E+30;
    absmin = 1.E+30;
    avg = 0.;
    absavg = 0.;

    for(x.first(); x.test(); x.next())
    {
      temp = field(x,i);
      avg += temp;

      if(temp >= max)
      {
        max = temp;
      }

  		if(temp <= min)
  		{
  			min = temp;
  		}

  		temp = fabs(temp);
      absavg += temp;

  		if(temp <= absmin)
  		{
  			absmin = temp;
  		}

      if(temp >= absmax)
      {
        absmax = temp;
      }
    }

    parallel.max(max);
    parallel.max(absmax);
    parallel.min(min);
    parallel.min(absmin);
    parallel.sum(avg);
    parallel.sum(absavg);
    avg /= n3;
    absavg /= n3;

    output_check_field(max, absmax, min, absmin, avg, absavg, sim.check_fields_precision, field_name + "[" + std::to_string(i) + "]", message);
  }

	cout.copyfmt(oldState);

	return;
}

/////////////////////////////////////////////
// Check all scalar fields
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void check_all_fields(
	Field<FieldType> &phi,
	Field<FieldType> &xi,
	Field<FieldType> &laplace_xi,
	Field<FieldType> &chi,
	Field<FieldType> &rho,
	Field<FieldType> &deltaR,
	Field<FieldType> &eightpiG_deltaT,
	Field<FieldType> &zeta,
  Field<FieldType> &phi_effective,
  Field<FieldType> &phi_ddot,
	Field<FieldType> &Bi,
	long n3,
  const metadata & sim)
{
  int oc = sim.out_check;

	if(oc & MASK_PHI) check_field(phi, "phi", n3, sim);
	if(oc & MASK_CHI) check_field(chi, "chi", n3, sim);
	if(oc & MASK_CHI) check_field(rho, "rho", n3, sim);
	if(oc & MASK_DELTAR) check_field(deltaR, "deltaR", n3, sim);
	if(oc & MASK_DELTAT) check_field(eightpiG_deltaT, "eightpiG_deltaT", n3, sim);
	if(oc & MASK_B) check_vector_field(Bi, "Bi", n3, sim);
  if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
  {
    if(oc & MASK_XI) check_field(xi, "xi", n3, sim);
    if(oc & MASK_LAPLACE_XI) check_field(laplace_xi, "laplace_xi", n3, sim);
    if(oc & MASK_PHI_EFFECTIVE) check_field(phi_effective, "phi_effective", n3, sim);
    if(oc & MASK_ZETA) check_field(zeta, "zeta", n3, sim);
  }
	if(oc & MASK_PHI_DDOT) check_field(phi_ddot, "phi_ddot", n3, sim);

	return;
}



/////////////////////////////////////////////
// Check particles (with some ID#, TODO: make this more universal, not with random ID
/////////////////////////////////////////////
void check_particles(
  const metadata & sim,
  Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm)
{
	int i, num[8] = {0, 262143, 1, 26, 30000, 180024, 220023, 262142};

	for(i=0; i<2; i++)
	{
		pcls_cdm->coutPart(num[i]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

/////////////////////////////////////////////
// Compute max of fRR
// TODO Add description
/////////////////////////////////////////////
template <class FieldType>
double compute_max_fRR(
  Field<FieldType> &deltaR,
  double Rbar,
  const metadata & sim)
{
	if(sim.fR_model == FR_MODEL_R2)// In R + R^2 gravity, fRR is a constant
	{
		return 2. * sim.fR_params[0];
	}

	Site x(deltaR.lattice());
	double temp,
				 max = -1.E+30;

	for(x.first(); x.test(); x.next())
	{
		temp = fabs(fRR(deltaR(x) + Rbar, sim, 555));
		if(temp > max) max = temp;
	}

	parallel.max(max);

	if(!max || isnan(max))
	{
		cout << "Something went wrong in computing max_fRR. Closing...\n";
		parallel.abortForce();
	}
  return max;
}
/////////////////////////////////////////////
// Flips (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void flip_fields(
  Field<FieldType> & field1,
  Field<FieldType> & field2)
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

/////////////////////////////////////////////
// Writes field1 --> field2, field2 --> field3, field3 --> field1
/////////////////////////////////////////////
template <class FieldType>
void flip_fields(
  Field<FieldType> & field1,
  Field<FieldType> & field2,
  Field<FieldType> & field3)
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


/////////////////////////////////////////////
// Flips vector fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void flip_vector_fields(
  Field<FieldType> & field1,
  Field<FieldType> & field2
)
{
	Site x(field1.lattice());
	int i;
	double t;

	for(x.first(); x.test(); x.next())
	{
		for(i=0; i<3; ++i)
		{
			t = field1(x,i);
			field1(x,i) = field2(x,i);
			field2(x,i) = t;
		}
	}

	field1.updateHalo();
	field2.updateHalo();

	return;
}

/////////////////////////////////////////////
// Sets whole field to zero
/////////////////////////////////////////////
template <class FieldType>
void erase_field(Field<FieldType> & field)
{
	Site x(field.lattice());
	for(x.first(); x.test(); x.next())
	{
		field(x) = 0.;
	}
}

/////////////////////////////////////////////
// Copies (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void copy_field(
  Field<FieldType> & source,
  Field<FieldType> & destination,
  double coeff = 1.)
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


/////////////////////////////////////////////
// Copies vector fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void copy_vector_field(
  Field<FieldType> & source,
  Field<FieldType> & destination,
  double coeff = 1.)
{
	Site x(source.lattice());
	int i;

	if(coeff == 1.)
	{
		for(x.first(); x.test(); x.next())
		{
			for(i=0; i<3; ++i)
			{
				destination(x,i) = source(x,i);
			}
		}
	}
	else
	{
		for(i=0; i<3; ++i)
		{
			destination(x,i) = source(x,i);
		}
	}
	return;
}


/////////////////////////////////////////////
// Sums (scalar) fields:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void add_fields(
  Field<FieldType> & field1,
  Field<FieldType> & field2,
  Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) + field2(x);
	}
	return;
}

/////////////////////////////////////////////
// Adds a constant to scalar field:
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void add_fields(
  Field<FieldType> & field1,
  double c,
  Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) + c;
	}
	return;
}

template <class FieldType>
void subtract_fields(
  Field<FieldType> & field1,
  Field<FieldType> & field2,
  Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = field1(x) - field2(x);
	}
	return;
}

template <class FieldType>
void add_fields(
  Field<FieldType> & field1,
  double coeff1,
  Field<FieldType> & field2,
  double coeff2,
  Field<FieldType> & result)
{
	Site x(field1.lattice());

	for(x.first(); x.test(); x.next())
	{
		result(x) = coeff1 * field1(x) + coeff2 * field2(x);
	}
	return;
}

template <class FieldType>
void scatter_field(
  Field<FieldType> & field1,
  double c1,
  Field<FieldType> & field2,
  double c2,
  Field<FieldType> & result,
  double k)
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

/////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void leapfrog_dotR(
  Field<FieldType> & deltaR,
  Field<FieldType> & deltaT,
  Field<FieldType> & dot_deltaR_old,
  Field<FieldType> & dot_deltaR_new,
  double Hubble,
  double coeff,
  double dtau_old,
  double dtau,
  double dx2)
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

/////////////////////////////////////////////
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
void leapfrog_dot_xi(
  Field<FieldType> & xi,
  Field<FieldType> & zeta,
  Field<FieldType> & dot_xi_old,
  Field<FieldType> & dot_xi_new,
  double Hubble,
  double a2,
  double dtau_old,
  double dtau,
  double dx2)
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

///// TODO: FOR EVERY INSTANCE OF convert_*_to_*, CHECK WITH
// if(convert_*_to_*() > FR_WRONG) return FR_WRONG_RETURN;
// (this could be important)

/////////////////////////////////////////////
// Converts deltaR to xi
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double convert_deltaR_to_xi(
  Field<FieldType> & deltaR,
  Field<FieldType> & xi,
  double Rbar,
  double fRbar,
  const metadata & sim)
{
	Site x(xi.lattice());
	double temp;

	if(sim.fR_model == FR_MODEL_R2)
	{
		for(x.first(); x.test(); x.next())
		{
			xi(x) = 2. * sim.fR_params[0] * deltaR(x);
		}
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1], n_minus_1 = n - 1.;
		for(x.first(); x.test(); x.next())
		{
			xi(x) = a * n * pow(Rbar + deltaR(x), n_minus_1) - fRbar;
		}
	}
	else
	{
		for(x.first(); x.test(); x.next())
		{
			temp = fR(Rbar + deltaR(x), sim, 832);
			//TODO REMOVE after debugging
			if(temp > FR_WRONG)
			{
	      cout << parallel.rank() << " 3 temp = FR_WRONG in convert_deltaR_to_xi -- xi, Rbar, temp = " << xi(x) << ", " << Rbar << ", " << temp << endl;
				parallel.abortForce();
				return FR_WRONG_RETURN; // Returns a huge number to throw some exception
			}
			//END REMOVE
			xi(x) =  temp - fRbar;
		}
	}
return 1.;
}


/////////////////////////////////////////////
// For single point
/////////////////////////////////////////////
double convert_deltaR_to_xi_single(
  double deltaR,
  double Rbar,
  double fRbar,
  const metadata & sim)
{
	if(sim.fR_model == FR_MODEL_R2)
	{
    return 2. * sim.fR_params[0] * deltaR;
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
    return sim.fR_params[0] * sim.fR_params[1] * pow(Rbar + deltaR, sim.fR_params[1] - 1.) - fRbar;
	}
	else
	{
    double fR_temp;
    fR_temp = fR(Rbar + deltaR, sim, 832);
    //TODO REMOVE after debugging
    if(fR_temp > FR_WRONG)
    {
      cout << parallel.rank() << " 3 fR_temp = FR_WRONG in convert_deltaR_to_xi -- Rbar, fR_temp = " << Rbar << ", " << fR_temp << endl;
      parallel.abortForce();
      return FR_WRONG_RETURN; // Returns a huge number to throw some exception
    }
    //END REMOVE
    return fR_temp - fRbar;
	}
}


/////////////////////////////////////////////
// For single point
// TODO This can probably be greatly optimised.
//      All the model-dependent tests are done for each point instead of just once.
/////////////////////////////////////////////
double convert_xi_to_deltaR_single(
	double & deltaR,
	double xi,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  double R_temp,
				 temp,
         fr,
				 max_fRR = 0.,
				 fRR_temp;
  int count,
			count_n = 1;

	if(fabs(xi) >= 1.E20 || isnan(xi))
	{
		cout << endl;
		cout << " convert_xi_to_deltaR_single before deltaR = " << deltaR << endl;
		cout << " convert_xi_to_deltaR_single before Rbar = " << Rbar << endl;
		cout << " convert_xi_to_deltaR_single before deltaR + Rbar = " << deltaR + Rbar << endl;
		cout << " convert_xi_to_deltaR_single before xi = " << xi << endl;
		cout << endl;
		return FR_WRONG_RETURN;
	}

	if(sim.fR_model == FR_MODEL_HU_SAWICKI)
	{
		double m2 = sim.fR_params[0],
					 c2 = sim.fR_params[1],
					 n  = sim.fR_params[2],
					 c1 = sim.fR_params[3];

		if(n == 1. || n == 1) // xi <-> R relation is invertible algebraically
		{
      fr = xi + fRbar;
			R_temp = - c1 / fr;
			R_temp = m2 * (sqrt(R_temp) - 1.) / c2;

      if(isnan(R_temp) || R_temp < 0.)
			{
				cout << endl;
				cout << " convert_xi_to_deltaR_single HS = " << endl;
				cout << "     fr = " << fr << endl;
				cout << "     xi = " << xi << endl;
				cout << "  fRbar = " << fRbar << endl;
				cout << " R_temp = " << R_temp << endl;
				cout << endl;
				R_temp = 0.;
			}
		}
		else
		{
      fr = xi + fRbar;
      // The initial guess for R is:
      // 1) found from a large-R expansion of FR or
      if(fr > - c1 * n / pow(c2, 1. - 1./n) / 4.)
      {
        R_temp = - c1 * n / c2 / c2 / fr;
        R_temp = m2 * pow(fabs(R_temp), 1./(n + 1.));
      }
      else
      {
        R_temp = m2 * pow(- fr / c1 / n, 1./(n - 1.));
      }

      count = 0;
      while(true)
      {
        temp = fR(R_temp, sim, 56);
        fRR_temp = fRR(R_temp, sim, 57);

        if(fabs(temp / fr - 1.) < sim.fR_target_precision) break;

        R_temp += (fr - temp) / fRR_temp;

        count ++;
        if(count > count_n * sim.fR_count_max)
        {
          cout << "Could only reach a precision of " << temp << " in " << count_n * sim.fR_count_max << " steps.\n";
          count_n++;
        }
      }
		}
	}
	else if(sim.fR_model == FR_MODEL_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1], n_minus_1 = n - 1.;

    R_temp = 1. + xi / (a * n * pow(Rbar, n_minus_1));
		if(R_temp > 0.)
		{
			R_temp = Rbar * pow(R_temp, 1./n_minus_1);
		}
	}
	else if(sim.fR_model == FR_MODEL_DELTA)
	{
		double a = sim.fR_params[0], delta = sim.fR_params[1];

		// METHOD 1
		R_temp = xi / a / (1. + delta) / pow(Rbar, delta) + 1.;
		R_temp = Rbar * pow(R_temp, 1./delta);

		// METHOD 2
		// R_temp = pow( (1. + xi + fRbar) / (a * (1. + delta)) , 1./delta);
	}

	deltaR = R_temp - Rbar;

	if(isnan(deltaR) || R_temp <= 0.)
	{
		cout << " convert_xi_to_deltaR_single END" << endl;
		cout << " deltaR = " << deltaR << endl;
		cout << " R_temp = " << R_temp << endl;
		cout << "   Rbar = " << Rbar << endl;
		return FR_WRONG_RETURN;
	}

  return deltaR;
}

/////////////////////////////////////////////
// TODO Details
// Returns maximum value of fRR over the entire lattice.
// the typical oscillation frenquency will be of order a/(sqrt(3*fRR_max))
/////////////////////////////////////////////
template <class FieldType>
double convert_xi_to_deltaR(
	Field<FieldType> & deltaR,
	Field<FieldType> & xi,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  // Using Newton-Raphson Method
  Site x(xi.lattice());

  double fRR_temp,
				 max_fRR = 0.;

  for(x.first(); x.test(); x.next())
  {
    convert_xi_to_deltaR_single(deltaR(x), xi(x), Rbar, fRbar, sim);
    fRR_temp = fabs(fRR(deltaR(x) + Rbar, sim, 656));
    if(fRR_temp > max_fRR) max_fRR = fRR_temp;
  }

  parallel.max(max_fRR);

	if(max_fRR and !isnan(max_fRR))
	{
		return max_fRR;
	}
  else
  {
    COUT << " Warning 19: returning fRRbar\n";
    return fRR(Rbar, sim, 65);
  }
}

/////////////////////////////////////////////////
// Builds laplacian from field
/////////////////////////////////////////////////
template <class FieldType>
void build_laplacian(
  Field<FieldType> & field,
  Field<FieldType> & laplacian,
  double dx)
{
	double dx2 = dx*dx;
	Site x(field.lattice());

  field.updateHalo(); // TODO Makes sure that halo is updated before computing derivatives

	for(x.first(); x.test(); x.next())
	{
		laplacian(x) = field(x+0) + field(x-0) + field(x+1) + field(x-1) + field(x+2) + field(x-2) - 6. * field(x);
		laplacian(x) /= dx2;
	}

	return;
}

/////////////////////////////////////////////////
// TODO: Maybe put this in gevolution.hpp?
/////////////////////////////////////////////////
template <class FieldType>
void prepareFTsource_leapfrog_R2(
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & source,
  double dx)
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


/////////////////////////////////////////////////
/////////////////////////////////////////////////
template <class FieldType>
double remove_homogeneous_part(
	Field<FieldType> & original_field,
	Field<FieldType> & pert_field,
	long n3
)
{
	Site x(original_field.lattice());
	double hom = 0;

	for(x.first(); x.test(); x.next())
	{
		hom += original_field(x);
	}

	parallel.sum(hom);
	hom /= n3;

	for(x.first(); x.test(); x.next())
	{
		pert_field(x) = original_field(x) - hom;
	}

	return hom;
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
template <class FieldType>
double build_phi_ddot(
	Field<FieldType> & phi,
	Field<FieldType> & chi,
	Field<FieldType> & phi_dot,
	Field<FieldType> & chi_dot,
	Field<FieldType> & phi_ddot,
	Field<FieldType> & deltaR,
	double dx,
	double a,
	double Hubble,
	double Rbar
)
{
	Site x(phi.lattice());
	int i;
	double
		temp1,
		temp2,
	  a2 = a*a,
		dx2 = dx*dx,
		lap_phi,
		lap_chi,
		grad_squared_phi,
		grad_squared_chi,
		grad_squared_mixed;

	for(x.first(); x.test(); x.next())
	{
		lap_phi = phi(x+0) + phi(x-0) + phi(x+1) + phi(x-1) + phi(x+2) + phi(x-2) - 6.*phi(x);
		lap_chi = chi(x+0) + chi(x-0) + chi(x+1) + chi(x-1) + chi(x+2) + chi(x-2) - 6.*chi(x);
		lap_phi /= dx2;
		lap_chi /= dx2;

		grad_squared_phi = 0.;
		grad_squared_chi = 0.;
		grad_squared_mixed = 0.;
		for(i=0; i<3; ++i)
		{
			temp1 = phi(x+i) - phi(x-i);
			grad_squared_phi += temp1 * temp1;
			temp2 = chi(x+i) - chi(x-i);
			grad_squared_chi += temp2 * temp2;
			grad_squared_mixed += temp1 * temp2;
		}
		grad_squared_phi /= 4. * dx2;
		grad_squared_chi /= 4. * dx2;
		grad_squared_mixed /= 4. * dx2;

		temp1 = -a2 * deltaR(x);
		temp1 += 2. * (1. + 8.*phi(x) - 2.*chi(x)) * lap_phi;
		temp1 += 2. * lap_chi;
		temp1 += 2. * a2 * Rbar * (chi(x) - phi(x));
		temp1 += 10. * grad_squared_phi - 6. * grad_squared_mixed + 2. * grad_squared_chi;
		temp1 += 6. * Hubble * (chi_dot(x) - 4. * phi_dot(x));

		phi_ddot(x) = temp1 / 6.;
	}

	return 1.;
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
template <class FieldType>
void generate_density_field_BH(
	Field<FieldType> & source,
	long numpts3d,
	const cosmology cosmo,
	const metadata & sim
)
{
	Site x(source.lattice());
	int xPointMass, yPointMass, zPointMass;
	double rho_point_mass, rho_elsewhere;
	xPointMass = sim.numpts/2;
	yPointMass = sim.numpts/2;
	zPointMass = sim.numpts/2;
	rho_point_mass = (1. + 1.E-4 * (numpts3d - 1)) * cosmo.Omega_m;
	rho_elsewhere = (1. - 1.E-4) * cosmo.Omega_m;

	for(x.first(); x.test(); x.next())
	{
		if(x.coord(0) == xPointMass && x.coord(1) == yPointMass && x.coord(2) == zPointMass)
		{
			source(x) = rho_point_mass;
		}
		else
		{
			source(x) = rho_elsewhere;
		}
	}


	source.updateHalo();

	return;

}




#endif
