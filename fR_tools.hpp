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

template <class FieldType>
void check_xi_first(
  Field<FieldType> & xi,
  double fRbar);


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
/////////////// FR_TYPE_RN
// f(R) = a * R^n = R * (R/(b*Ra))^(n-1)
// For convenience pass as arguments: b, n
// where Ra := 8piG * rho_0 (= 8piG in code units)
// then compute a = (b*Ra)^(1-n)
//
/////////////// FR_TYPE_R2
// Sub-class of FR_TYPE_RN with n=2 -- "special" because F'' = const.
// non-linearities disappear from the evolution equation for deltaR
//
/////////////// FR_TYPE_DELTA
// f(R) = R * (R/(a * Ra))^delta - R
// where Ra := 8piG * rho_0 (= 8piG in code units), a dimensionless constant
// Pass a, delta

/////////////// FR_TYPE_HU_SAWICKI
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

	if(sim.fR_type == FR_TYPE_R2)
	{
		output = sim.fR_params[0] * R * R;

		if(isnan(output))
		{
      cout << " f(R) Evaluated to NaN (code " << code << ").\n"
      << " R = " << R << ", R^2 = " << R*R << ", f(R) = " << output << endl;
      cin.get();
      return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];
		if(R > 0.)
    {
      Rpow = pow(R, n);
    }
		else
		{
			cout << " /!\\ Warning 3! Rpow is negative in f [FR_TYPE_RN]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = a * Rpow;

    if(isnan(output))
		{
			cout << " f(R) Evaluated to NaN (code " << code << ").\n"
			<< " R = " << R << ", R^n = " << Rpow << ", f(R) = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double
		m2 = sim.fR_params[0],
		c2 = sim.fR_params[1],
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3];

    if(n == 1. || n == 1)
    {
      if(R > 0)
      {
        output = - m2 * c1 / (c2 + m2 / R);
      }
      else
      {
        cout << " /!\\ Warning 4! Rpow is negative in f [FR_TYPE_HU_SAWICKI]. R = " << R  << endl;
        return FR_WRONG_RETURN;
      }
    }
    else
    {
      if(R > 0.)
      {
        Rpow = pow(m2/R, n); // Notice the inverted power
      }
      else
      {
        cout << " /!\\ Warning 4 bis! Rpow is negative in f [FR_TYPE_HU_SAWICKI]. R = " << R  << endl;
        return FR_WRONG_RETURN;
      }

      output = - m2 * c1 / (c2 + Rpow);
    }
    if(isnan(output))
    {
      cout << " f(R) Evaluated to NaN (code " << code << ").\n"
      << " R = " << R << ", Rpow = " << Rpow << ", f(R) = " << output << endl;
      cin.get();
      return FR_WRONG_RETURN; // Returns a huge number to throw some exception
    }
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];

		if(R > 0.)
    {
      if(delta > 0)
      {
        Rpow = pow(R, delta);
      }
      else
      {
        Rpow = 1. / pow(R, -delta);
      }
    }
		else
		{
			cout << " /!\\ Warning 5! R is negative in f [FR_TYPE_DELTA]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = R * (b * Rpow - 1.);

		if(isnan(output))
		{
			cout
			<< " f(R) Evaluated to NaN (code " << code << ").\n"
			<< " R = " << R << ", R^delta = " << Rpow << ", f(R) = " << output << endl;
			cin.get();
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

	if(sim.fR_type == FR_TYPE_R2)
	{
    output = 2. * sim.fR_params[0] * R;
		if(R < 0.)
		{
			cout << " /!\\ Warning 6! R is negative in fR [FR_TYPE_R2]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ").\n"
			<< " R = " << R << ", fR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];

		if(R > 0.)
    {
      Rpow = n > 1. ? pow(R, n - 1.) : pow(R, 1. - n);
    }
		else
		{
			cout << " /!\\ Warning 7! R is negative in fR [FR_TYPE_RN]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = a * n * Rpow;

		if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ").\n"
			     << " R = " << R << ", R^(n-1) = " << Rpow << ", fR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
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
			if(R > 0.)
      {
        Rpow = pow(R/m2, n);
      }
			else
			{
				cout << " /!\\ Warning 8! R is negative in fR [FR_TYPE_HU_SAWICKI]. R = " << R  << endl;
				return FR_WRONG_RETURN;
			}

			output = (1. + c2 * Rpow);

      if(n > 1.)
      {
        output = - c1 * n * pow(R/m2, n-1.) / output / output;
      }
      else
      {
        output = - c1 * n / pow(R/m2, 1.-n) / output / output;
      }
		}

		if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ").\n"
			     << " R = " << R << ", (R/m2)^n = " << Rpow << ", fR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];

		if(R > 0.)
    {
      if(delta > 0)
      {
        Rpow = pow(R, delta);
      }
      else
      {
        Rpow = 1. / pow(R, -delta);
      }
    }
		else
		{
			cout << " /!\\ Warning 9! R is negative in fR [FR_TYPE_DELTA]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = b * (1. + delta) * Rpow - 1.;

    if(isnan(output))
		{
			cout << " fR Evaluated to NaN (code " << code << ").\n"
			     << " R = " << R << ", R^delta = " << Rpow << ", fR = " << output << endl;
			cin.get();
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

	if(sim.fR_type == FR_TYPE_R2)
	{
		output = 2. * sim.fR_params[0];
		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ").\n"
			<< " R = " << R << ", fRR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];

		if(R > 0)
    {
      Rpow = pow(R, n) / R / R;
    }
		else
		{
			cout << " /!\\ Warning 10! R is negative in fRR [FR_TYPE_RN]. R = " << R << endl;
			return FR_WRONG_RETURN;
		}

		output = a * n * (n - 1.) * Rpow;

		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ").\n"
					 << " R = " << R << ", R^(n-2) = " << Rpow << ", fRR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double
		m2 = sim.fR_params[0],
		c2 = sim.fR_params[1],
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3];

		R_over_m2 = R/m2;

		if(R > 0.)
    {
      Rpow = pow(R_over_m2, n);
    }
		else
		{
			cout << " /!\\ Warning 11! R is negative in fRR [FR_TYPE_HU_SAWICKI]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

    output = (1. + c2 * Rpow);

		if(n == 1. || n == 1)
		{
			output = 2. * c1 * c2 / ( m2 * output * output * output);
		}
		else
		{
			output = c1 * n * Rpow * (1. - n + (1. + n) * c2 * Rpow ) / (R_over_m2 * R_over_m2 * m2 * output * output * output);
		}

		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ").\n"
			     << " R = " << R << ", (R/m2)^n = " << Rpow << ", fRR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];

		if(R > 0.)
    {
      if(delta > 0.)
      {
        Rpow = pow(R, delta) / R;
      }
      else
      {
        Rpow = 1. / R / pow(R, -delta);
      }
    }
		else
		{
			cout << " /!\\ Warning 12! R is negative in fRR [FR_TYPE_DELTA]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = b * delta * (delta + 1.) * Rpow;

		if(isnan(output))
		{
			cout << " fRR Evaluated to NaN (code " << code << ").\n"
			     << " R = " << R << ", R^(delta-1) = " << Rpow << ", fRR = " << output << endl;
			cin.get();
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
  int code)
{
	double output = 0.,
				 R_over_m2,
				 Rpow;

	if(sim.fR_type == FR_TYPE_R2)
	{
		return 0.;
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1];

		if(R > 0.)
    {
      Rpow = pow(R, n) / R / R / R;
    }
		else
		{
			cout << " /!\\ Warning 13! R is negative in fRRR [FR_TYPE_RN]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = a * n * (n - 1.) * (n - 2.) * Rpow;

		if(isnan(output))
		{
			cout << " fRRR Evaluated to NaN (code " << code << ").\n"
			<< " R = " << R << ", R^(n-3) = " << Rpow << ", fRRR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double
		m2 = sim.fR_params[0],
		c2 = sim.fR_params[1],
		n  = sim.fR_params[2],
		c1 = sim.fR_params[3];

		R_over_m2 = R/m2;

		if(R > 0.)
    {
      Rpow = pow(R_over_m2, n);
    }
		else
		{
			cout << " /!\\ Warning 14! R is negative in fRRR [FR_TYPE_HU_SAWICKI]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

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
			cout << " fRRR Evaluated to NaN (code " << code << ").\n"
			     << " R = " << R << ", (R/m2)^n = " << Rpow << ", fRRR = " << output << endl;
			cin.get();
			return FR_WRONG_RETURN; // Returns a huge number to throw some exception
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];

		if(R > 0.)
    {
      if(delta > 0)
      {
        Rpow = pow(R, delta) / R / R;
      }
      else
      {
        Rpow = 1. / R / R / pow(R, -delta);
      }
    }
		else
		{
			cout << " /!\\ Warning 15! R is negative in fRR [FR_TYPE_DELTA]. R = " << R  << endl;
			return FR_WRONG_RETURN;
		}

		output = b * delta * (delta * delta - 1.) * Rpow;

		if(isnan(output))
		{
			cout << " fRRR Evaluated to NaN (code " << code << ").\n"
			<< " R = " << R << ", R(delta-2) = " << Rpow << ", fRRR = " << output << endl;
			cin.get();
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
// Print out f(R) model details
/////////////////////////////////////////////
void fR_details(
  const cosmology cosmo,
  metadata * sim,
  double fourpiG)
{
	if(sim->fR_type == FR_TYPE_HU_SAWICKI)
	{
		//TODO: Check whether we need to rescale something by some power of the boxsize
		COUT << " Model: Hu-Sawicki\n f(R) = - m2 * c1 * pow(R/m2, n) / (1. + c2 * pow(R/m2, n))" << endl
				 << " with m2 = " << sim->fR_params[0] << " * 8piG * rho_{m0} / 3. = ";
		sim->fR_params[0] *= fourpiG * cosmo.Omega_m / 1.5;
		sim->fR_params[1] = 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[1] / pow(12./cosmo.Omega_m - 9., sim->fR_params[2] + 1.) / sim->fR_params[0];
		sim->fR_params[3] = sim->fR_params[1] * 4. * fourpiG * (1. - cosmo.Omega_m - cosmo.Omega_rad) / sim->fR_params[0];
		COUT
		<< sim->fR_params[0] << endl
		<< "      c1 = " << sim->fR_params[3] << endl
		<< "      c2 = " << sim->fR_params[1] << endl
		<< "       n = " << sim->fR_params[2] << endl
		<< "so |fR0| ~ -n*c1/c2^2/(12/Omega_m - 9)^(n+1) = " << sim->fR_params[2] * sim->fR_params[3] / sim->fR_params[1]/sim->fR_params[1] / pow(12. / cosmo.Omega_m - 9., sim->fR_params[2] + 1.) << " (should be << 1)\n";
	}
	else if(sim->fR_type == FR_TYPE_RN)
	{
		COUT
		<< " f(R) model:" << endl
		<< " f(R) = R * (R / (b*8piG*rho0))^(n-1), with b = " << sim->fR_params[0] << endl
		<< "                                            n = " << sim->fR_params[1] << endl;
		if(sim->fR_params[1] > 1.) sim->fR_params[0] = 1. / pow(2. * fourpiG * sim->fR_params[0], sim->fR_params[1] - 1.);
		else sim->fR_params[0] = 0.5 / pow(2. * fourpiG * sim->fR_params[0], sim->fR_params[1]) / (fourpiG * sim->fR_params[0]);
		COUT << " rescaled as: f(R) = a * R^n, with a = " << sim->fR_params[0] << endl;
	}
	else if(sim->fR_type == FR_TYPE_R2)
	{
		COUT
		<< " f(R) model:" << endl
		<< " f(R) = R^2 / (b*8piG*rho0), with b = " << sim->fR_params[0] << endl;
		sim->fR_params[0] = 0.5 / fourpiG / sim->fR_params[0];
		COUT << " rescaled as: f(R) = a * R^2, with a = " << sim->fR_params[0] << endl;
	}
	else if(sim->fR_type == FR_TYPE_DELTA) //  TODO: Fix for delta < 0: fRR < 0, so everything blows up!
	{
		double temp = sim->fR_params[0];
		sim->fR_params[0] = 1. / pow(sim->fR_params[0] * 2. * fourpiG, sim->fR_params[1]);
		COUT
		<< " f(R) model: f(R) = (R / (a * 8piG * rho0))^(1+delta) - R" << endl
		<< "                a = " << temp << endl
		<< "            delta = " << sim->fR_params[1] << endl
		<< "         or: f(R) = b * R^(1+delta) - R" << endl
		<< "                b = " << sim->fR_params[0] << endl;
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
	return max;
}


/////////////////////////////////////////////
template <class FieldType>
double check_field(
  Field<FieldType> & field,
  string field_name,
  long n3,
  int prec = 6,
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

  output_check_field(max, absmax, min, absmin, avg, absavg, prec, field_name, message);
	cout.copyfmt(oldState);
	return max;
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

	return max;
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

	return max;
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
	Field<FieldType> &chi,
	Field<FieldType> &deltaR,
	Field<FieldType> &eightpiG_deltaT,
	Field<FieldType> &zeta,
	Field<FieldType> &Bi,
	long n3,
  const metadata & sim)
{
  int oc = sim.out_check;

	if(oc & MASK_PHI) check_field(phi, "phi", n3, sim);
	if(oc & MASK_CHI) check_field(chi, "chi", n3, sim);
	if(oc & MASK_DELTAR) check_field(deltaR, "deltaR", n3, sim);
	if(oc & MASK_DELTAT) check_field(eightpiG_deltaT, "eightpiG_deltaT", n3, sim);
	if(oc & MASK_B) check_vector_field(Bi, "Bi", n3, sim);
  if(sim.mg_flag == FLAG_FR)
  {
    if(oc & MASK_XI) check_field(xi, "xi", n3, sim);
    if(oc & MASK_ZETA) check_field(zeta, "zeta", n3, sim);
    if(oc & (MASK_XI * MASK_CHI)) check_correl(xi, "xi", chi, "chi", n3, sim);
  }

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
	if(sim.fR_type == FR_TYPE_R2)// In R + R^2 gravity, fRR is a constant
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
		if(result(x) <= 0.)
		{
			cout << " result(x) = " << field1(x) << " + " << c << " = " << result(x);
			cin.get();
		}
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
void leapfrog_dotxi(
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

	if(sim.fR_type == FR_TYPE_R2)
	{
		for(x.first(); x.test(); x.next())
		{
			xi(x) = 2. * sim.fR_params[0] * deltaR(x);
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
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
				cin.get();
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
	if(sim.fR_type == FR_TYPE_R2)
	{
    return 2. * sim.fR_params[0] * deltaR;
	}
	else if(sim.fR_type == FR_TYPE_RN)
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
      cin.get();
      return FR_WRONG_RETURN; // Returns a huge number to throw some exception
    }
    //END REMOVE
    return fR_temp - fRbar;
	}
}


/////////////////////////////////////////////
// Converts deltaR to u
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double convert_deltaR_to_u(
	Field<FieldType> & deltaR,
	Field<FieldType> & u,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
	Site x(u.lattice());
	double fR_temp;

	if(sim.fR_type == FR_TYPE_R2)
	{
		for(x.first(); x.test(); x.next())
		{
			u(x) = log(1. + deltaR(x)/Rbar);
		}
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double n_minus_1 = sim.fR_params[1] - 1.;
		for(x.first(); x.test(); x.next())
		{
			u(x) = n_minus_1 * log(1. + deltaR(x)/Rbar);
		}
	}
	else
	{
		for(x.first(); x.test(); x.next())
		{
			fR_temp = fR(Rbar + deltaR(x), sim, 833);
			if(fR_temp > FR_WRONG)
			{
				//TODO REMOVE after debugging
				cout << parallel.rank() << " 4 fR_temp = FR_WRONG in convert_deltaR_to_u -- u, Rbar, fR_temp = " << u(x) << ", " << Rbar << ", " << fR_temp << endl;
				//END REMOVE
				cin.get();
				return FR_WRONG_RETURN; // Returns a huge number to throw some exception
			}
			u(x) = log(fR_temp/fRbar);
		}
	}
	return 1.;
}

/////////////////////////////////////////////
// For single point
/////////////////////////////////////////////
double convert_deltaR_to_u_single(
	double deltaR,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
	if(sim.fR_type == FR_TYPE_R2)
	{
    return log(1. + deltaR/Rbar);
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
    return (sim.fR_params[1] - 1.) * log(1. + deltaR/Rbar);
	}
	else
	{
    double fR_temp;
    fR_temp = fR(Rbar + deltaR, sim, 833);
    if(fR_temp > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 4 fR_temp = FR_WRONG in convert_deltaR_to_u -- Rbar, fR_temp = " << Rbar << ", " << fR_temp << endl;
      //END REMOVE
      cin.get();
      return FR_WRONG_RETURN; // Returns a huge number to throw some exception
    }
    return log(fR_temp/fRbar);
	}
}

/////////////////////////////////////////////
// Converts deltaR to scalaron
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double convert_deltaR_to_scalaron(
	Field<FieldType> & deltaR,
	Field<FieldType> & scalaron,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  double t;

  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    t = convert_deltaR_to_u(deltaR, scalaron, Rbar, fRbar, sim);
  }
  else
  {
    t = convert_deltaR_to_xi(deltaR, scalaron, Rbar, fRbar, sim);
  }

  if(t > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << parallel.rank() << " 44 fR_temp = FR_WRONG in convert_deltaR_to_scalaron" << endl;
    //END REMOVE
    cin.get();
    return FR_WRONG_RETURN; // Returns a huge number to throw some exception
  }

  return t;
}

/////////////////////////////////////////////
// Converts u to xi
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double convert_u_to_xi(
  Field<FieldType> & u,
  Field<FieldType> & xi,
  double fRbar)
{
	Site x(u.lattice());

	for(x.first(); x.test(); x.next())
	{
		xi(x) = fRbar * (exp(u(x)) - 1.);
	}

	return 1.;
}

/////////////////////////////////////////////
// For single point
/////////////////////////////////////////////
inline double convert_u_to_xi_single(
  double u,
  double fRbar)
{
  return fRbar * (exp(u) - 1.);
}


/////////////////////////////////////////////
// Converts xi to u
// TODO: Add comments here
/////////////////////////////////////////////
template <class FieldType>
double convert_xi_to_u(
  Field<FieldType> & xi,
  Field<FieldType> & u,
  double fRbar)
{
	Site x(u.lattice());

	for(x.first(); x.test(); x.next())
	{
		u(x) = log(1. + xi(x)/fRbar);
	}
	return 1.;
}

/////////////////////////////////////////////
// For single point
/////////////////////////////////////////////
inline double convert_xi_to_u_single(
  double xi,
  double fRbar)
{
  return log(1. + xi/fRbar);
}

/////////////////////////////////////////////
// TODO Details
// Returns maximum value of fRR over the entire lattice.
// the typical oscillation frenquency will be of order a/(sqrt(3*fRR_max))
/////////////////////////////////////////////
template <class FieldType>
double convert_xi_to_deltaR(
	Field<FieldType> & eightpiG_deltaT,
	Field<FieldType> & deltaR,
	Field<FieldType> & xi,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  // Using Newton-Raphson Method
  Site x(xi.lattice());

  double R_temp,
				 temp,
				 max_fRR = 0.,
				 fRR_temp;
  int count,
			count_n = 1;

	if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
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
					COUT << " R has become negative at one point. Check what's going on..." << endl;
					cin.get();
					return FR_WRONG_RETURN;
				}

				R_temp = m2 / c2 * ( sqrt(R_temp) - 1.);
				fRR_temp = fabs(fRR(R_temp, sim, 5911));

        if(fRR_temp > FR_WRONG)
				{
					//TODO REMOVE after debugging
		      cout << parallel.rank() << " 5 fRR_temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, R_temp = " << xi(x) << ", " << Rbar << ", " << R_temp << endl;
		      //END REMOVE
					cin.get();
					return FR_WRONG_RETURN; // Returns a huge number to throw some exception
				}

				if(max_fRR < fRR_temp) max_fRR = fRR_temp;

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
					R_temp = - c1 * n / c2 / c2 / (xi(x) + fRbar);
					if(R_temp > 0.) R_temp = m2 * pow(R_temp, 1./(n+1.));
					else
					{
						cout << " /!\\ Warning 16! R_temp is negative in convert_xi_to_deltaR [FR_TYPE_HU_SAWICKI]. R_temp = " << R_temp  << ", x = " << x << endl;
						return FR_WRONG_RETURN;
					}
				}

				count = 0;
				while(true)
				{
					temp = fR(R_temp, sim, 56);
					if(temp > FR_WRONG)
					{
						//TODO REMOVE after debugging
						cout << parallel.rank() << " 6 temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, temp = " << xi(x) << ", " << Rbar << ", " << temp << endl;
			      //END REMOVE
						cin.get();
						return FR_WRONG_RETURN;
					}

					fRR_temp = fRR(R_temp, sim, 57);
					if(fRR_temp > FR_WRONG)
					{
						//TODO REMOVE after debugging
						cout << parallel.rank() << " 7 fRR_temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, temp = " << xi(x) << ", " << Rbar << ", " << temp << endl;
			      //END REMOVE
						cin.get();
						return FR_WRONG_RETURN;
					}

          if(fabs((temp - xi(x))/fRbar - 1.) < sim.fR_target_precision) break;

          R_temp += (xi(x) + fRbar - temp) / fRR_temp;

					count ++;
					if(count > count_n * sim.fR_count_max)
					{
						cout << "Could only reach a precision of " << temp << " in " << count_n * sim.fR_count_max << " steps.\n";
						count_n++;
					}
				}

				if(max_fRR < fRR_temp) max_fRR = fRR_temp;

				deltaR(x) = R_temp - Rbar;
			}
		}
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		cout << "Something is wrong. R + R^2 should not call this function. Closing...\n";
		parallel.abortForce();
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1], n_minus_1 = n - 1.;

		for(x.first(); x.test(); x.next())
		{
			R_temp = (xi(x) + fRbar) / a / n;
			if(R_temp > 0.)
			{
				if(n > 1.) R_temp = pow(R_temp, 1./n_minus_1);
				else R_temp = 1. / pow(R_temp, -1./n_minus_1);
			}
			else
			{
				cout << " /!\\ Warning 17! R is negative in convert_xi_to_deltaR. R_temp = " << R_temp << ", xi = " << xi(x) << ", fRbar = " << fRbar << ", xi + fRbar = " << xi(x) + fRbar << ", x =" << x << endl;
				return FR_WRONG_RETURN;
			}

			fRR_temp = fabs(fRR(R_temp, sim, 60));
			if(fRR_temp > FR_WRONG)
			{
				//TODO REMOVE after debugging
				cout << parallel.rank() << " 10 fRR_temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, R_temp = " << xi(x) << ", " << Rbar << ", " << R_temp << endl;
	      //END REMOVE
				cin.get();
				return FR_WRONG_RETURN;
			}
 			if(max_fRR < fRR_temp) max_fRR = fRR_temp;

			deltaR(x) = R_temp - Rbar;
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];
		for(x.first(); x.test(); x.next())
		{
			// TODO: Expect big error propagation here, maybe try to find a more reliable conversion?
			if(delta > 0.)
			{
				R_temp = xi(x) / b / (1. + delta) + pow(Rbar, delta);

				if(R_temp > 0.)
        {
          R_temp = pow(R_temp, 1./delta);
        }
				else
				{
					cout << " /!\\ Warning 18! R_temp is negative in convert_xi_to_deltaR. R_temp = " << R_temp << ", x = " << x << endl;
					return FR_WRONG_RETURN;
				}
			}
			else
			{
				R_temp = xi(x) / b / (1. + delta) + 1. / pow(Rbar, -delta);

        if(R_temp > 0.)
        {
          R_temp = 1. / pow(R_temp, -1./delta);
        }
				else
				{
					cout << " /!\\ Warning 18! R_temp is negative in convert_xi_to_deltaR. R_temp = " << R_temp << ", x = " << x << endl;
					return FR_WRONG_RETURN;
				}
			}

			temp = fRR(R_temp, sim, 61);
			if(temp > FR_WRONG)
			{
				//TODO REMOVE after debugging
				cout << parallel.rank() << " 11 temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, R_temp = " << xi(x) << ", " << Rbar << ", " << R_temp << endl;
	      //END REMOVE
				cin.get();
				return FR_WRONG_RETURN;
			}
			fRR_temp = fabs(temp);
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;

			deltaR(x) = R_temp - Rbar;
		}
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


/////////////////////////////////////////////
// For single point
/////////////////////////////////////////////
double convert_xi_to_deltaR_single(
	double eightpiG_deltaT,
	double xi,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  // Using Newton-Raphson Method
  double R_temp = 0.,
				 temp,
				 fRR_temp;
  int count,
			count_n = 1;

	if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2 = sim.fR_params[0],
					 c2 = sim.fR_params[1],
					  n = sim.fR_params[2],
					 c1 = sim.fR_params[3];

		if(n == 1. || n == 1) // xi <-> R relation is invertible algebraically
		{
      R_temp = -c1 / (xi + fRbar);

      if(R_temp <= 1.)
      {
        COUT << " R has become negative at one point. Check what's going on..." << endl;
        cin.get();
        return FR_WRONG_RETURN;
      }

      R_temp = m2 / c2 * (sqrt(R_temp) - 1.);
		}
		else
		{
      // The initial guess for R is:
      // 1) found from a large-R expansion of FR for R >> m2, TODO: Does this work??
      // 2) equal to the GR solution (Rbar - eightpiG_deltaT)
      R_temp = Rbar - eightpiG_deltaT;

      if(R_temp > 10. * m2)
      {
        R_temp = -c1 * n / c2 / c2 / (xi + fRbar);
        if(R_temp > 0.)
        {
          return m2 * pow(R_temp, 1./(n+1.));
        }
        else
        {
          cout << " /!\\ Warning 16! R_temp is negative in convert_xi_to_deltaR [FR_TYPE_HU_SAWICKI]. R_temp = " << R_temp << endl;
          return FR_WRONG_RETURN;
        }
      }

      count = 0;
      while(true)
      {
        temp = fR(R_temp, sim, 56);
        if(temp > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 6 temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, temp = " << xi << ", " << Rbar << ", " << temp << endl;
          //END REMOVE
          cin.get();
          return FR_WRONG_RETURN;
        }


        fRR_temp = fRR(R_temp, sim, 57);
        if(fRR_temp > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 6 temp = FRR_WRONG in convert_xi_to_deltaR -- xi, Rbar, temp = " << xi << ", " << Rbar << ", " << fRR_temp << endl;
          //END REMOVE
          cin.get();
          return FR_WRONG_RETURN;
        }

        if(fabs((temp - xi)/fRbar - 1.) < sim.fR_target_precision)
        {
          break;
        }

        R_temp += (xi + fRbar - temp) / fRR_temp;

        count ++;
        if(count > count_n * sim.fR_count_max)
        {
          cout << "Could only reach a precision of " << temp << " in " << count_n * sim.fR_count_max << " steps.\n";
        }
      }
		}
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		cout << "Something is wrong. R + R^2 should not call this function. Closing...\n";
		parallel.abortForce();
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double a = sim.fR_params[0], n = sim.fR_params[1], n_minus_1 = n - 1.;

    R_temp = (xi + fRbar) / a / n;

    if(R_temp > 0.)
    {
      if(n > 1.) R_temp = pow(R_temp, 1./n_minus_1);
      else R_temp = 1. / pow(R_temp, -1./n_minus_1);
    }
    else
    {
      cout << " /!\\ Warning 17! R is negative in convert_xi_to_deltaR. R_temp = " << R_temp << ", xi = " << xi << ", fRbar = " << fRbar << ", xi + fRbar = " << xi + fRbar << ", x =" << endl;
      return FR_WRONG_RETURN;
    }
	}
	else if(sim.fR_type == FR_TYPE_DELTA)
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];
    // TODO: Expect big error propagation here, maybe try to find a more reliable conversion?
    if(delta > 0.)
    {
      R_temp = xi / b / (1. + delta) + pow(Rbar, delta);

      if(R_temp > 0.)
      {
        R_temp = pow(R_temp, 1./delta);
      }
      else
      {
        cout << " /!\\ Warning 18! R_temp is negative in convert_xi_to_deltaR. R_temp = " << R_temp << endl;
        return FR_WRONG_RETURN;
      }
    }
    else
    {
      R_temp = xi / b / (1. + delta) + 1. / pow(Rbar, -delta);
      if(R_temp > 0.) R_temp = 1. / pow(R_temp, -1./delta);
      else
      {
        cout << " /!\\ Warning 18! R_temp is negative in convert_xi_to_deltaR. R_temp = " << R_temp << endl;
        return FR_WRONG_RETURN;
      }
    }

    temp = fRR(R_temp, sim, 61);
    if(temp > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 11 temp = FR_WRONG in convert_xi_to_deltaR -- xi, Rbar, R_temp = " << xi << ", " << Rbar << ", " << R_temp << endl;
      //END REMOVE
      cin.get();
      return FR_WRONG_RETURN;
    }
	}

  return R_temp - Rbar;
}

/////////////////////////////////////////////////
// Convert u = log(xi/fRbar + 1.)
/////////////////////////////////////////////////
template <class FieldType>
double convert_u_to_deltaR(
	Field<FieldType> & eightpiG_deltaT,
	Field<FieldType> & deltaR,
	Field<FieldType> & u,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  // Using Newton-Raphson Method
  Site x(u.lattice());
  double R_pow,
				 R_temp = 0.,
				 R1,
			   temp,
				 max_fRR = 0.,
				 fRR_temp;
  int count,
			count_n = 1;

	if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2 = sim.fR_params[0],
					 c2 = sim.fR_params[1],
					  n = sim.fR_params[2],
					 c1 = sim.fR_params[3];

		if(n == 1.) // u <-> R relation is invertible algebraically
		{
			c2 = c2 / m2; // they always appear in this combination for n = 1
			R_pow = fabs(1. + c2 * Rbar);
			for(x.first(); x.test(); x.next())
			{
				R_temp = (R_pow / exp(u(x)/2.) - 1.) / c2;
				R1 = - (R_pow / exp(u(x)/2.) + 1.) / c2;
				if(fabs(R_temp + eightpiG_deltaT(x) - Rbar) > fabs(R1 + eightpiG_deltaT(x) - Rbar)) R_temp = R1;
				fRR_temp = fRR(R_temp, sim, 5913);
				if(fRR_temp > FR_WRONG || isnan(fRR_temp) || !fRR_temp)
				{
					cout << x << endl;
					cout << parallel.rank() << " has a problem, u, R_temp, fRR_temp = ";
					cout << u(x) << ", " << R_temp << ", " << fRR_temp << endl;
					cin.get();
					return FR_WRONG_RETURN;
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
				//TODO CHECK THIS!
				// The initial guess for R is:
				// 1) found from a large-R expansion of FR for R >> m2, or
				// 2) equal to the new GR solution (Rbar - eightpiG_deltaT)
				R_temp = Rbar - eightpiG_deltaT(x);
				if(R_temp > 10. * m2)
				{
					R_temp = -c1 * n / c2 / c2 / fRbar / exp(u(x));
					R_temp = m2 * pow(R_temp, 1./(n+1.));
				}

				count = 0;
				while(true)
				{
					temp = fR(R_temp, sim, 56);
					if(temp > FR_WRONG)
					{
						//TODO REMOVE after debugging
						cout << parallel.rank() << " 12 temp = FR_WRONG in convert_u_to_deltaR -- x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
			      //END REMOVE
						cin.get();
						return FR_WRONG_RETURN;
					}

					temp = fabs(u(x)/(temp - fRbar) - 1.);
					if(temp < sim.fR_target_precision) break;

					fRR_temp = fRR(R_temp, sim, 57);
					if(fRR_temp > FR_WRONG)
					{
						//TODO REMOVE after debugging
						cout << parallel.rank() << " 13 fRR_temp = FR_WRONG in convert_u_to_deltaR -- x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
			      //END REMOVE
						cin.get();
						return FR_WRONG_RETURN;
					}

					if(fRR_temp > 0 && fabs(R_temp) <= 1.E+20)
					{
						temp = fR(R_temp, sim, 58);
						if(temp > FR_WRONG)
						{
							//TODO REMOVE after debugging
							cout << parallel.rank() << " 14 temp = FR_WRONG in convert_u_to_deltaR -- x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
				      //END REMOVE
							cin.get();
							return FR_WRONG_RETURN;
						}

						R_temp += (fRbar * u(x) - temp) / fRR_temp;
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

				temp = fRR(R_temp, sim, 5914);
				if(temp > FR_WRONG)
				{
					//TODO REMOVE after debugging
					cout << parallel.rank() << " 15 temp = FR_WRONG in convert_u_to_deltaR -- x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
		      //END REMOVE
					cin.get();
					return FR_WRONG_RETURN;
				}

				fRR_temp = fabs(temp);
				if(max_fRR < fRR_temp)
				{
					max_fRR = fRR_temp;
				}
				deltaR(x) = R_temp - Rbar;
			}
		}
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		cout << "Something is wrong. R + R^2 shouldn not call this function. Closing...\n";
		parallel.abortForce();
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double n = sim.fR_params[1];
		for(x.first(); x.test(); x.next())
		{
			R_temp = Rbar * exp(u(x)/(n-1.));

			if(R_temp > 0.) temp = fRR(R_temp, sim, 62);
			else
			{
				cout << " /!\\ Warning 1! R_temp is negative in convert_u_to_deltaR. x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
				return FR_WRONG_RETURN;
			}

			if(temp > FR_WRONG)
			{
				//TODO REMOVE after debugging
	      cout << parallel.rank() << " 1 temp = FR_WRONG in convert_u_to_deltaR -- x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
	      //END REMOVE
				cin.get();
				return FR_WRONG_RETURN;
			}
			fRR_temp = fabs(temp);
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;

			deltaR(x) = R_temp - Rbar;
		}
	}
	else if(sim.fR_type == FR_TYPE_DELTA) // For other f(R) models TODO check!!
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];
		for(x.first(); x.test(); x.next())
		{
			// TODO: Check if it's possible to find a more accurate way if this is not good enough
			if(delta > 0.)
			{
				R_temp = exp(u(x)) * (b * (1. + delta) * pow(Rbar, delta) - 1.) + 1.;
				R_temp /= b * (1. + delta);

				if(R_temp > 0.)
        {
          R_temp = pow(R_temp, 1./delta);
        }
				else
				{
					cout << " /!\\ Warning 2! R_temp is negative in convert_u_to_deltaR. R_temp = " << R_temp << ", x = " << x << endl;
					return FR_WRONG_RETURN;
				}
			}
			else
			{
				R_temp = exp(u(x)) * (b * (1. + delta) / pow(Rbar, -delta) - 1.) + 1.;
				R_temp /= b * (1. + delta);
				if(R_temp > 0.) R_temp = 1. / pow(R_temp, -1./delta);
				else
				{
					cout << " /!\\ Warning 2! R_temp is negative in convert_u_to_deltaR. R_temp = " << R_temp << ", x = " << x << endl;
					return FR_WRONG_RETURN;
				}
			}

			temp = fRR(R_temp, sim, 63);
			if(temp > FR_WRONG)
			{
				//TODO REMOVE after debugging
				cout << parallel.rank() << " 2 temp = FR_WRONG in convert_u_to_deltaR -- x, u, Rbar, R_temp = " << x << ", " << u(x) << ", " << Rbar << ", " << R_temp << endl;
	      //END REMOVE
				cin.get();
				return FR_WRONG_RETURN;
			}
			fRR_temp = fabs(temp);
			if(max_fRR < fRR_temp) max_fRR = fRR_temp;

			deltaR(x) = R_temp - Rbar;
		}
	}

	parallel.max(max_fRR);

  if(max_fRR and !isnan(max_fRR))
	{
		return max_fRR;
	}
  else
  {
    COUT << " Warning: returning fRRbar\n";
    return fRR(Rbar, sim, 65);
  }
}


/////////////////////////////////////////////////
// For single point
/////////////////////////////////////////////////
double convert_u_to_deltaR_single(
	double eightpiG_deltaT,
	double u,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  // Using Newton-Raphson Method
  double R_pow,
				 R_temp = 0.,
				 R1,
			   temp,
				 fRR_temp;
  int count,
			count_n = 1;

	if(sim.fR_type == FR_TYPE_HU_SAWICKI)
	{
		double m2 = sim.fR_params[0],
					 c2 = sim.fR_params[1],
					  n = sim.fR_params[2],
					 c1 = sim.fR_params[3];

		if(n == 1.) // u <-> R relation is invertible algebraically
		{
			c2 = c2 / m2; // they always appear in this combination for n = 1
			R_pow = fabs(1. + c2 * Rbar);
      R_temp = (R_pow / exp(u/2.) - 1.) / c2;
      R1 = - (R_pow / exp(u/2.) + 1.) / c2;

      if(fabs(R_temp + eightpiG_deltaT - Rbar) > fabs(R1 + eightpiG_deltaT - Rbar)) R_temp = R1;
		}
		else
		{
      //TODO CHECK THIS!
      // The initial guess for R is:
      // 1) found from a large-R expansion of FR for R >> m2, or
      // 2) equal to the new GR solution (Rbar - eightpiG_deltaT)
      R_temp = Rbar - eightpiG_deltaT;
      if(R_temp > 10. * m2)
      {
        R_temp = -c1 * n / c2 / c2 / fRbar / exp(u);
        R_temp = m2 * pow(R_temp, 1./(n+1.));
      }

      count = 0;
      while(true)
      {
        temp = fR(R_temp, sim, 56);
        if(temp > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 12 temp = FR_WRONG in convert_u_to_deltaR -- u, Rbar, R_temp = " << u << ", " << Rbar << ", " << R_temp << endl;
          //END REMOVE
          cin.get();
          return FR_WRONG_RETURN;
        }

        if(fabs(u/(temp - fRbar) - 1.) < sim.fR_target_precision) break;

        fRR_temp = fRR(R_temp, sim, 57);
        if(fRR_temp > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 13 fRR_temp = FR_WRONG in convert_u_to_deltaR -- u, Rbar, R_temp = " << u << ", " << Rbar << ", " << R_temp << endl;
          //END REMOVE
          cin.get();
          return FR_WRONG_RETURN;
        }

        R_temp += (fRbar * u - temp) / fRR_temp;

        count ++;
        if(count > count_n * sim.fR_count_max)
        {
          cout << "Could only reach a precision of " << temp << " in " << count_n * sim.fR_count_max << " steps.\n";
          count_n++;
        }
      }
		}
	}
	else if(sim.fR_type == FR_TYPE_R2)
	{
		cout << "Something is wrong. R + R^2 shouldn not call this function. Closing...\n";
		parallel.abortForce();
	}
	else if(sim.fR_type == FR_TYPE_RN)
	{
		double n = sim.fR_params[1];

    R_temp = Rbar * exp(u/(n-1.));

    if(R_temp < 0.)
    {
      cout << " /!\\ Warning 1! R_temp is negative in convert_u_to_deltaR. u, Rbar, R_temp = " << u << ", " << Rbar << ", " << R_temp << endl;
      return FR_WRONG_RETURN;
    }
	}
	else if(sim.fR_type == FR_TYPE_DELTA) // For other f(R) models TODO check!!
	{
		double b = sim.fR_params[0], delta = sim.fR_params[1];

    // TODO: Check if it's possible to find a more accurate way if this is not good enough
    if(delta > 0.)
    {
      R_temp = exp(u) * (b * (1. + delta) * pow(Rbar, delta) - 1.) + 1.;
      R_temp /= b * (1. + delta);

      if(R_temp > 0.)
      {
        R_temp = pow(R_temp, 1./delta);
      }
      else
      {
        cout << " /!\\ Warning 2! R_temp is negative in convert_u_to_deltaR. R_temp = " << R_temp << endl;
        return FR_WRONG_RETURN;
      }
    }
    else
    {
      R_temp = exp(u) * (b * (1. + delta) / pow(Rbar, -delta) - 1.) + 1.;
      R_temp /= b * (1. + delta);
      if(R_temp > 0.)
      {
        R_temp = 1. / pow(R_temp, -1./delta);
      }
      else
      {
        cout << " /!\\ Warning 2! R_temp is negative in convert_u_to_deltaR. R_temp = " << R_temp << endl;
        return FR_WRONG_RETURN;
      }
    }
	}

  return R_temp - Rbar;
}


/////////////////////////////////////////////////
// Convert scalaron to deltaR
/////////////////////////////////////////////////
template <class FieldType>
double convert_scalaron_to_deltaR(
	Field<FieldType> & eightpiG_deltaT,
	Field<FieldType> & deltaR,
	Field<FieldType> & scalaron,
	double Rbar,
	double fRbar,
	const metadata & sim)
{
  double t;

  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    t = convert_u_to_deltaR(eightpiG_deltaT, deltaR, scalaron, Rbar, fRbar, sim);
  }
  else
  {
    t = convert_xi_to_deltaR(eightpiG_deltaT, deltaR, scalaron, Rbar, fRbar, sim);
  }

  if(t > FR_WRONG)
  {
    // TODO: REMOVE
    cout << "Something went wrong in convert_scalaron_to_deltaR.";
    cin.get();
    // END REMOVE
    return FR_WRONG_RETURN;
  }

  return t;
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

//////////////////////////
// Builds laplace_exp_u term
//////////////////////////
template <class FieldType>
void build_laplacian_exp(
  Field<FieldType> & field,
  Field<FieldType> & laplace_exp_field,
  double dx)
{
  Site x(field.lattice());
  double dx2 = dx*dx;

  field.updateHalo();

  for(x.first(); x.test(); x.next())
  {
    // TODO: Check which definition works best

    // Option 1: laplace (exp(u))
    laplace_exp_field(x) = exp(field(x+0)) + exp(field(x-0)) + exp(field(x+1)) + exp(field(x-1)) + exp(field(x+2)) + exp(field(x-2)) - 6.*exp(field(x));
    // End Option 1

    // Option 2: div ( exp(field) * grad(field) )
    // double aplus, aminus;
    // laplace_exp_field(x) = 0.;
    // aplus = exp(field(x+0)) + exp(field(x));
    // aminus = exp(field(x)) + exp(field(x-0));
    // laplace_exp_field(x) += aplus * ( field(x+0) - field(x) ) - aminus * ( field(x) - field(x-0) );
    // aplus = exp(field(x+1)) + exp(field(x));
    // aminus = exp(field(x)) + exp(field(x-1));
    // laplace_exp_field(x) += aplus * ( field(x+1) - field(x) ) - aminus * ( field(x) - field(x-1) );
    // aplus = exp(field(x+2)) + exp(field(x));
    // aminus = exp(field(x)) + exp(field(x-2));
    // laplace_exp_field(x) += aplus * ( field(x+2) - field(x) ) - aminus * ( field(x) - field(x-2) );
    // laplace_exp_field(x) /= 2.;
    // End Option 2

    // Common for both options
    laplace_exp_field(x) /= dx2;
  }
  return;
}

/////////////////////////////////////////////////
// Builds laplacian from field
/////////////////////////////////////////////////
template <class FieldType>
void build_laplacian_scalaron(
  Field<FieldType> & scalaron,
  Field<FieldType> & laplacian,
  double dx,
  const metadata & sim)
{
  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    build_laplacian_exp(scalaron, laplacian, dx);
  }
  else
  {
    build_laplacian(scalaron, laplacian, dx);
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



#endif
