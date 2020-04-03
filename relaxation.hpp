//////////////////////////
// relaxation.hpp
//////////////////////////
//
// Implementation of relaxation method solver for f(R) gravity
//
// Authors: David Daverio (Cambridge University)
//          Lorenzo Reverberi (CEICO - Czech Academy of Sciences, Prague)
//
//////////////////////////
#include <unistd.h>
#ifndef RELAXATION_HEADER
#define RELAXATION_HEADER

//////////////////////////
// Euclidean Norm
//////////////////////////
template <class FieldType>
double euclidean_norm(
  MultiField<FieldType> & field,
  MultiGrid & engine,
  int level,
  long numpts3d)
{
  double error = 0,
         temp;

  Site x(field.lattice());
  if(engine.isPartLayer(level))
  {
    for(x.first(); x.test(); x.next())
    {
      temp = field(x);
      error += temp * temp;
    }

    if(level == 0)
    {
      parallel.sum(error);
    }
    else
    {
      parallel.layer_from_level(level).sum(error);
    }
  }

  error = sqrt(error)/numpts3d;

  return error;
}


//////////////////////////
// TODO: add description
//////////////////////////
template <class FieldType>
void restrict_to_level(
  MultiField<FieldType> * field,
  MultiGrid & engine,
  int target_level)
{
  for(int i=0; i<target_level; ++i)
  {
    engine.restrict(field, i);
  }
  return;
}

//////////////////////////
// Prepare the initial guess for the trace equation
// TODO: Details here
//////////////////////////
template <template<class> class FieldClass, class FieldType>
int prepare_initial_guess_trace_equation(
  FieldClass<FieldType> & deltaR,
  FieldClass<FieldType> & eightpiG_deltaT,
  FieldClass<FieldType> & xi,
  FieldClass<FieldType> & laplace_xi,
  FieldClass<FieldType> & rhs,
  double a,
  double dx,
  double Rbar,
  double fRbar,
  const metadata & sim)
{
  copy_field(eightpiG_deltaT, rhs, a*a/3.); // source term

  // TODO
  // METHOD 1 -- deltaR = - eightpiG_deltaT; xi, laplace_xi consistent;
  copy_field(eightpiG_deltaT, deltaR, -1.);
  convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim, "prepare_initial_guess_trace_equation");
  build_laplacian(xi, laplace_xi, dx);
  // END METHOD 1

  //// METHOD 2 -- all fields to zero
  // erase_field(xi);
  // erase_field(laplace_xi);
  // erase_field(deltaR);
  //// END METHOD 2

  return 0;
}

////////////////////////////////////////////////////////////
///////////////// ROUTINES AND TOOLS FOR xi ////////////////
////////////////////////////////////////////////////////////

//////////////////////////
// Residual Y
// of equation Y[xi] == 0
// TODO: Details here
//////////////////////////
double residual_xi(
  double phi,
  double xi,
  double xi_old,
  double laplace_xi,
  double deltaR,
  double rhs,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim,
  string message = ""
)
{
  double temp = deltaR;

  if(sim.relativistic_flag)
  {
    double hub_corr = sim.xi_Hubble ? 1. : 0.;

    temp = temp * (1. - fRbar - xi) + 2. * (f(Rbar + deltaR, sim, message+" residual_xi") - fbar) - Rbar * xi;
    temp = laplace_xi - a2_over_3 * temp - rhs;
    temp += hub_corr * (2. * phi * laplace_xi - two_Hubble_over_dtau * (xi - xi_old));
  }
  else
  {
    temp = laplace_xi - a2_over_3 * temp - rhs;
  }

  return temp;
}

//////////////////////////
// Calculation of (dY/dxi)
// from equation Y[xi] == 0
// Needed for next guess xi^{i+1} = xi^{i} - Y^{i}/(dY/dxi)^{i}
// TODO Correct the description above
//////////////////////////
double dresidual_dxi(
  double xi,
  double deltaR,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fRbar,
  const metadata & sim,
  string message = ""
)
{
  double R = Rbar + deltaR, temp = fRR(R, sim, message+" dresidual_dxi");

  if(!temp || temp > FR_WRONG)
  {
    cout << endl;
    cout << " dresidual_dxi" << endl;
    cout << " deltaR " << deltaR << endl;
    cout << " xi " << xi << endl;
    cout << " Rbar " << Rbar << endl;
    cout << " R " << R << endl;
    cout << " fRR " << temp << endl;
    cout << endl;
    parallel.abortForce();
  }

  temp = 1. / temp;

  if(sim.relativistic_flag)
  {
    // TODO Additional terms -- see if necessary
    temp = temp * (1. + fRbar + xi) - R;
    // End additional terms
  }

  temp = - a2_over_3 * temp;

  if(sim.xi_Hubble)
  {
    temp -= two_Hubble_over_dtau;
  }

  return temp;
}

//////////////////////////
// Computes new guess for xi
// Solves the equation Y[xi] == 0
// where Y[xi] := laplacian(xi) - a^2/3 * zeta(xi)
// TODO: Details here
//////////////////////////
double update_xi_single(
  double phi,
  double & xi,
  double xi_old,
  double laplace_xi,
  double deltaR,
  double rhs,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double coeff_laplacian,
  double overrelax,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double
    fr,
    denominator,
    resid,
    xi_prev = xi;

  if(std::isnan(deltaR) || std::isnan(xi) || Rbar + deltaR <= 0. || fabs(xi) >= 1.E20)
  {
    cout << endl;
    cout << " update_xi_single before" << endl;
    cout << " Rbar = " << Rbar << endl;
    cout << " R = " << Rbar + deltaR << endl;
    cout << " xi = " << xi << endl;
    cout << " fRbar = " << fRbar << endl;
    cout << " fR = " << xi + fRbar << endl;
    cout << " deltaR = " << deltaR << endl;
    deltaR = convert_xi_to_deltaR_single(deltaR, xi, Rbar, fRbar, sim);
    cout << " deltaR_new = " << deltaR << endl;
    cout << endl;
  }

  resid = residual_xi(phi, xi, xi_old, laplace_xi, deltaR, rhs, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
  denominator = dresidual_dxi(xi, deltaR, a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, "update_xi_and_deltaR 1");

  if(sim.relativistic_flag)
  {
    denominator += (1. + 2.*phi) * coeff_laplacian;
  }
  else
  {
    denominator += coeff_laplacian;
  }

  if(denominator)
  {
    if(sim.fR_model == FR_MODEL_HU_SAWICKI)
    {
      fr = xi + fRbar;
      if(fr)
      {
        resid = - overrelax * resid / fr / denominator;

        while(true)
        {
          xi = fr * exp(resid) - fRbar;
          if(!xi_allowed(xi, Rbar, fRbar, sim))
          {
            resid /= 2.;
          }
          else break;
        }
      }
      else
      {
        cout << " fR = 0 in update_xi_single." << endl;
      }
    }
    else
    {
      resid = - overrelax * resid / denominator;
      while(true)
      {
        xi = xi_prev + resid;
        if(!xi_allowed(xi, Rbar, fRbar, sim))
        {
          resid /= 2.;
        }
        else break;
      }
    }
  }

  return resid;
}

double update_xi_and_deltaR_single(
  double phi,
  double & xi,
  double xi_old,
  double laplace_xi,
  double & deltaR,
  double rhs,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double coeff_laplacian,
  double overrelax,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double temp, dxi;

  dxi = xi; // original value
  temp = update_xi_single(phi, xi, xi_old, laplace_xi, deltaR, rhs, a2_over_3, two_Hubble_over_dtau, coeff_laplacian, overrelax, Rbar, fbar, fRbar, sim);

  if(std::isnan(deltaR) || fabs(xi) >= 1.E20 || deltaR + Rbar <= 0.)
	{
    cout << endl;
    cout << " update_xi_and_deltaR_single 2" << endl;
		cout << " deltaR = " << deltaR << endl;
		cout << " Rbar = " << Rbar << endl;
		cout << " deltaR + Rbar = " << deltaR + Rbar << endl;
		cout << " xi = " << xi << endl;
    cout << endl;
	}

  dxi = (xi - dxi) / 10.;

  while(true)
  {
    temp = convert_xi_to_deltaR_single(deltaR, xi, Rbar, fRbar, sim);

    if(temp < FR_WRONG)
    {
      break;
    }
    else
    {
      cout << "W ";
      xi -= dxi;
    }
  }

  if(std::isnan(xi) || fabs(xi) > 1.E20)
  {
    cout << endl;
    cout << " update_xi_and_deltaR after" << endl;
    cout << " xi = " << xi << endl;
    cout << " dxi = " << laplace_xi << endl;
    cout << " Rbar = " << Rbar << endl;
    cout << " deltaR = " << deltaR << endl;
    cout << " R = " << Rbar + deltaR << endl;
    cout << endl;
    parallel.abortForce();
  }

  return deltaR;
}

////////////////////////////////////////////////////////////
// Relaxation -- Update values of xi and deltaR
template <template<class> class FieldClass, class FieldType>
void update_xi_and_deltaR(
  FieldClass<FieldType> & phi,
  FieldClass<FieldType> & xi,
  FieldClass<FieldType> & xi_old,
  FieldClass<FieldType> & laplace_xi,
  FieldClass<FieldType> & deltaR,
  FieldClass<FieldType> & eightpiG_deltaT,
  FieldClass<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim,
  string message = ""
)
{
  double temp, dxi, resid, denominator, fr, coeff_laplacian, overrelax, relat_flag;

  coeff_laplacian = -6./dx/dx;
  overrelax = sim.overrelaxation_factor;
  relat_flag = sim.relativistic_flag ? 1. : 0.;

  if(!sim.red_black)// Sequential sweep
  {
    Site x(xi.lattice());

    // Update xi
    if(sim.fR_model == FR_MODEL_HU_SAWICKI)
    {
      double n = sim.fR_params[2], c1 = sim.fR_params[3], c2 = sim.fR_params[1], xi_min;

      if(n == 1 || n == 1.)
      {
        xi_min = - c1 - fRbar;
      }
      else
      {
        xi_min = - (1.+n) * (1.+n) * pow((n-1.) / c2 / (1.+n), 1.-1./n) / n / 4. - fRbar;
      }

      for(x.first(); x.test(); x.next())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        fr = xi(x) + fRbar;
        if(fr)
        {
          resid = - overrelax * resid / fr / denominator;
          while(true)
          {
            xi(x) = fr * exp(resid) - fRbar;

            // Check if value of xi is allowed
            if(xi(x) <= xi_min || xi(x) + fRbar >= 0.) // reduce step if not allowed
            {
              resid /= 2.;
            }
            else break; // break if ok
          }
        }
        else
        {
          cout << " fR = 0 in update_xi_single." << endl;
        }
      }
    }
    else if(sim.fR_model == FR_MODEL_RN)
  	{
  		double xi_min = - fRbar;

      for(x.first(); x.test(); x.next())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;

        resid = - overrelax * resid / denominator;

        while(true)
        {
          if(xi(x) + resid <= xi_min)
          {
            resid /= 2.;
          }
          else
          {
            xi(x) += resid;
            break;
          }
        }
      }
  	}
    else if(sim.fR_model == FR_MODEL_DELTA)
  	{
  		double a = sim.fR_params[0], delta = sim.fR_params[1];
  		double xi_min = - a * (1. + delta) * pow(Rbar, delta);

      for(x.first(); x.test(); x.next())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        resid = - overrelax * resid / denominator;

        while(true)
        {
          if(xi(x) + resid <= xi_min)
          {
            resid /= 2.;
          }
          else
          {
            xi(x) += resid;
            break;
          }
        }
      }
  	}
  }
  else // Red-Black scheme
  {
    SiteRedBlack3d x(xi.lattice());

    if(sim.fR_model == FR_MODEL_HU_SAWICKI)
    {
      double n = sim.fR_params[2], c1 = sim.fR_params[3], c2 = sim.fR_params[1], xi_min;
      if(n == 1 || n == 1.)
      {
        xi_min = - c1 - fRbar;
      }
      else
      {
        xi_min = - (1.+n)*(1.+n) * pow((n-1.) / c2 / (1.+n), 1.-1./n) / n / 4.- fRbar;
      }

      // Red cells first
      for(x.firstRed(); x.testRed(); x.nextRed())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        fr = xi(x) + fRbar;
        if(fr)
        {
          resid = - overrelax * resid / fr / denominator;
          while(true)
          {
            xi(x) = fr * exp(resid) - fRbar;

            // Check if value of xi is allowed
            if(xi(x) <= xi_min || xi(x) + fRbar >= 0.) // reduce step if not allowed
            {
              resid /= 2.;
            }
            else break; // break if ok
          }
        }
        else
        {
          cout << " fR = 0 in update_xi_single." << endl;
        }
      }

      // convert_xi_to_deltaR(deltaR, xi, Rbar, fRbar, sim);
      build_laplacian(xi, laplace_xi, dx);

      // Then Black cells
      for(x.firstBlack(); x.testBlack(); x.nextBlack())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        fr = xi(x) + fRbar;
        if(fr)
        {
          resid = - overrelax * resid / fr / denominator;
          while(true)
          {
            xi(x) = fr * exp(resid) - fRbar;

            // Check if value of xi is allowed
            if(xi(x) <= xi_min || xi(x) + fRbar >= 0.) // reduce step if not allowed
            {
              resid /= 2.;
            }
            else break; // break if ok
          }
        }
        else
        {
          cout << " fR = 0 in update_xi_single." << endl;
        }
      }
    }
    else if(sim.fR_model == FR_MODEL_RN)
    {
      double xi_min = - fRbar;

      for(x.firstRed(); x.testRed(); x.nextRed())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        resid = - overrelax * resid / denominator;

        while(true)
        {
          if(xi(x) + resid <= xi_min)
          {
            resid /= 2.;
          }
          else
          {
            xi(x) += resid;
            break;
          }
        }
      }

      convert_xi_to_deltaR(deltaR, xi, Rbar, fRbar, sim);
      build_laplacian(xi, laplace_xi, dx);

      for(x.firstBlack(); x.testBlack(); x.nextBlack())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        resid = - overrelax * resid / denominator;

        while(true)
        {
          if(xi(x) + resid <= xi_min)
          {
            resid /= 2.;
          }
          else
          {
            xi(x) += resid;
            break;
          }
        }
      }
    }
    else if(sim.fR_model == FR_MODEL_DELTA)
    {
      double a = sim.fR_params[0], delta = sim.fR_params[1];
      double xi_min = - a * (1. + delta) * pow(Rbar, delta);

      for(x.firstRed(); x.testRed(); x.nextRed())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        resid = - overrelax * resid / denominator;

        while(true)
        {
          if(xi(x) + resid <= xi_min)
          {
            resid /= 2.;
          }
          else
          {
            xi(x) += resid;
            break;
          }
        }
      }

      convert_xi_to_deltaR(deltaR, xi, Rbar, fRbar, sim);
      build_laplacian(xi, laplace_xi, dx);

      for(x.firstBlack(); x.testBlack(); x.nextBlack())
      {
        resid = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator = dresidual_dxi(xi(x), deltaR(x), a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
        denominator += (1. + 2. * relat_flag * phi(x)) * coeff_laplacian;
        resid = - overrelax * resid / denominator;

        while(true)
        {
          if(xi(x) + resid <= xi_min)
          {
            resid /= 2.;
          }
          else
          {
            xi(x) += resid;
            break;
          }
        }
      }
    }
  }

  convert_xi_to_deltaR(deltaR, xi, Rbar, fRbar, sim, message+" update_xi_and_deltaR");
  build_laplacian(xi, laplace_xi, dx);

  return;
}

//////////////////////////////////////////////////
//////////// OTHER TOOLS AND ROUTINES ////////////
//////////////////////////////////////////////////

template <template<class> class FieldClass, class FieldType>
double compute_error(
  FieldClass<FieldType> & phi,
  FieldClass<FieldType> & xi,
  FieldClass<FieldType> & xi_old,
  FieldClass<FieldType> & laplace_xi,
  FieldClass<FieldType> & deltaR,
  FieldClass<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim,
  int level,
  string message = ""
)
{
  double temp, error;
  Site x(xi.lattice());

  error = 0.;

  for(x.first(); x.test(); x.next())
  {
    temp = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, message+" compute_error");
    error += temp*temp;
  }

  if(level == 0)
  {
    parallel.sum(error);
  }
  else
  {
    parallel.layer_from_level(level).sum(error);
  }

  error = sqrt(error)/numpts3d;

  return error;
}


//////////////////////////
// TODO Details
//////////////////////////
template <template<class> class FieldClass, class FieldType>
void build_residual(
  FieldClass<FieldType> & phi,
  FieldClass<FieldType> & xi,
  FieldClass<FieldType> & xi_old,
  FieldClass<FieldType> & laplace_xi,
  FieldClass<FieldType> & deltaR,
  FieldClass<FieldType> & rhs,
  FieldClass<FieldType> & destination_for_residual,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  Site x(xi.lattice());

  for(x.first(); x.test(); x.next())
  {
    destination_for_residual(x) = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
  }

  return;
}

///////////////////////////////////////////////////
///////////////////// SOLVERS /////////////////////
///////////////////////////////////////////////////

//////////////////////////
// Relaxation solver
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation_solver(
  MultiField<FieldType> * phi,
  MultiField<FieldType> * xi,
  MultiField<FieldType> * xi_old,
  MultiField<FieldType> * laplace_xi,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * residual,
  MultiField<FieldType> * err,
  MultiGrid & engine,
  double dx,
  double a,
  double two_Hubble_over_dtau,
  long numpts3d,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim
)
{
  double error;
  double a2_over_3 = a*a/3.;
  int max_level = sim.multigrid_n_grids - 1; // TODO Not really needed, just for clarity

  COUT << left << " z = " << setw(8) << 1./a - 1. << right << " -- " << flush;
  error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, 0, "relaxation_solver");
  COUT << "in.error = " << error << endl;

  if(sim.relaxation_method == METHOD_RELAX)
  {
    error = single_layer_solver(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], eightpiG_deltaT[0], rhs[0], dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, 0);
  }
  else if(sim.relaxation_method == METHOD_MULTIGRID)
  {
    error = multigrid_solver(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, err, engine, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
  }
  else if(sim.relaxation_method == METHOD_FMG)
  {
    error = FMG_solver(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, err, engine, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
  }

  // TODO REMOVE after debug
  if(error > FR_WRONG)
  {
    COUT << " Method " << sim.relaxation_method << " returned FR_WRONG_RETURN. Check what's going on. Press [Enter] to continue..." << endl;
    parallel.abortForce();
  }
  // END REMOVE

  COUT << "                fin.error = " << error << endl;
  return error;
}


//////////////////////////
// Single relaxation step
// TODO: Details here
//////////////////////////
template <template<class> class FieldClass, class FieldType>
double relaxation_step(
  FieldClass<FieldType> & phi,
  FieldClass<FieldType> & xi,
  FieldClass<FieldType> & xi_old,
  FieldClass<FieldType> & laplace_xi,
  FieldClass<FieldType> & deltaR,
  FieldClass<FieldType> & eightpiG_deltaT,
  FieldClass<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim,
  int level
)
{
  update_xi_and_deltaR(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, "relaxation_step");

  return compute_error(phi, xi, xi_old, laplace_xi, deltaR, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, level, "relaxation_step");
}



//////////////////////////
// Relaxation solver
// TODO: Details here
//////////////////////////
template <template<class> class FieldClass, class FieldType>
double single_layer_solver(
  FieldClass<FieldType> & phi,
  FieldClass<FieldType> & xi,
  FieldClass<FieldType> & xi_old,
  FieldClass<FieldType> & laplace_xi,
  FieldClass<FieldType> & deltaR,
  FieldClass<FieldType> & eightpiG_deltaT,
  FieldClass<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim,
  int level,
  string message = "",
  int smoothing = 0
)
{
  int count = 0;
  double error = 0.;
  double previous_error = compute_error(phi, xi, xi_old, laplace_xi, deltaR, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, level, message+" single_layer_solver 1");

  if(!level && previous_error < sim.relaxation_error)
  {
    return previous_error;
  }

  while(true)
  {
    relaxation_step(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, level);

    error = compute_error(phi, xi, xi_old, laplace_xi, deltaR, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, level, message+" single_layer_solver 2");

    if(error < sim.relaxation_error)
    {
      break;
    }

    if(error/previous_error > 1.)
    {
      scatter_field(deltaR, 0.5, eightpiG_deltaT, -0.5, deltaR);
      convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim, message+" single_layer_solver");
      build_laplacian(xi, laplace_xi, dx);
      error = compute_error(phi, xi, xi_old, laplace_xi, deltaR, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim, level, message+" single_layer_solver 3");
    }

    if(smoothing && count >= smoothing) break;

    previous_error = error;
    ++count;
  }

  return error;
}


//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double multigrid_solver(
  MultiField<FieldType> * phi,
  MultiField<FieldType> * xi,
  MultiField<FieldType> * xi_old,
  MultiField<FieldType> * laplace_xi,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * residual,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double deltax,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  int max_level = sim.multigrid_n_grids-1;
  long numpts3d[sim.multigrid_n_grids];
  double
    dx[sim.multigrid_n_grids],
    error = 0.,
    previous_error = 0.;

  numpts3d[0] = (long) sim.numpts * sim.numpts * sim.numpts;
  dx[0] = deltax;

  for(int i=1; i<sim.multigrid_n_grids; ++i)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8;
    dx[i] = dx[i-1] * 2.;
  }

  // Build eightpiG_deltaT, phi and xi_old on each level (they will not change)
  restrict_to_level(eightpiG_deltaT, engine, max_level);
  restrict_to_level(xi_old, engine, max_level);
  restrict_to_level(phi, engine, max_level);

  while(true)
  {
    previous_error = error;
    error = gamma_cycle(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, trunc, engine, dx, a2_over_3, two_Hubble_over_dtau, numpts3d, Rbar, fbar, fRbar, sim, 0, "multigrid_solver");

    if(error > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 1 trunc_error = FR_WRONG in multigrid_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    convert_xi_to_deltaR(deltaR[0], xi[0], Rbar, fRbar, sim);
    error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[0], sim, 0, "multigrid_solver");

    if(error < sim.relaxation_error)
    {
      break;
    }
    else
    {
      COUT << "             interm.error = " << error << endl;
    }
  }

  return error;
}


//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double FMG_solver(
  MultiField<FieldType> * phi,
  MultiField<FieldType> * xi,
  MultiField<FieldType> * xi_old,
  MultiField<FieldType> * laplace_xi,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * residual,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double deltax,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  int
    i,
    j,
    max_level = sim.multigrid_n_grids - 1; // Not exactly needed, but more redable this way
  double
   dx[sim.multigrid_n_grids],
   error = 0.,
   previous_error,
   check_fr_wrong;

  long numpts3d[sim.multigrid_n_grids];

  numpts3d[0] = (long) (sim.numpts * sim.numpts * sim.numpts);
  dx[0] = deltax;
  for(i=1; i<sim.multigrid_n_grids; ++i)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8.;
    dx[i] = 2. * dx[i-1];
  }

  // Build eightpiG_deltaT, phi and xi_old on each level (they will not change during the relaxation)
  restrict_to_level(eightpiG_deltaT, engine, max_level);
  restrict_to_level(xi_old, engine, max_level);
  restrict_to_level(phi, engine, max_level);

  while(true)
  {
    previous_error = error;

    // Build source term on coarsest and finest grid
    if(engine.isPartLayer(max_level))
    {
      copy_field(eightpiG_deltaT[max_level], rhs[max_level], a2_over_3);
    }

    copy_field(eightpiG_deltaT[0], rhs[0], a2_over_3);

    // build initial guess on coarsest grid -- xi and deltaR should be consistent
    if(sim.multigrid_restrict_mode == RESTRICT_XI)
    {
      restrict_to_level(xi, engine, max_level);
    }
    else
    {
      restrict_to_level(deltaR, engine, max_level);
    }

    // Compute initial error
    // TODO: remove after debugging
    build_laplacian(xi[0], laplace_xi[0], dx[0]);
    error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[0], sim, 0, "FMG_solver 1");

    // solve problem on coarsest grid
    if(engine.isPartLayer(max_level))
    {
      build_laplacian(xi[max_level], laplace_xi[max_level], dx[max_level]);

      if(sim.multigrid_restrict_mode == RESTRICT_XI)
      {
        check_fr_wrong = convert_xi_to_deltaR(deltaR[max_level], xi[max_level], Rbar, fRbar, sim);
      }
      else
      {
        check_fr_wrong = convert_deltaR_to_xi(deltaR[max_level], xi[max_level], Rbar, fRbar, sim);
      }

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 10 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }

      while(true)
      {
        update_xi_and_deltaR(phi[max_level], xi[max_level], xi_old[max_level], laplace_xi[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[max_level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim, "FMG_solver");

        error = compute_error(phi[max_level], xi[max_level], xi_old[max_level], laplace_xi[max_level], deltaR[max_level], rhs[max_level], dx[max_level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[max_level], sim, max_level, "FMG_solver 2");

        if(error < sim.relaxation_error)
        {
          break;
        }
      }
    }

    for(i=sim.multigrid_n_grids-2; i>=0; --i) // i is the starting level for the V or W mini-cycle
    {
      // prolong i+1 solution
      if(sim.multigrid_restrict_mode == RESTRICT_XI)
      {
        engine.prolong(xi, i+1);
      }
      else
      {
        engine.prolong(deltaR, i+1);
      }

      // build/copy source to the correct field (rhs)
      if(engine.isPartLayer(i))
      {
        copy_field(eightpiG_deltaT[i], rhs[i], a2_over_3);
      }

      for(j=0; j<sim.multigrid_n_cycles; ++j)
      {
        check_fr_wrong = gamma_cycle(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, trunc, engine, dx, a2_over_3, two_Hubble_over_dtau, numpts3d, Rbar, fbar, fRbar, sim, i);

        if(check_fr_wrong > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 12 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
          //END REMOVE
          return FR_WRONG_RETURN;
        }

        if(sim.multigrid_check_shape && i && j == sim.multigrid_n_cycles-1)
        {
          for(int k=0; k<i; ++k) COUT << "  ";
          COUT << i << endl;
          for(int k=0; k<i-1; ++k) COUT << "  ";
          COUT << " /" << endl;
        }
      }
    }

    if(sim.multigrid_restrict_mode == RESTRICT_XI)
    {
      check_fr_wrong = convert_xi_to_deltaR(deltaR[0], xi[0], Rbar, fRbar, sim);
    }
    else
    {
      check_fr_wrong = convert_deltaR_to_xi(deltaR[0], xi[0], Rbar, fRbar, sim);
    }

    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 13 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    build_laplacian(xi[0], laplace_xi[0], dx[0]);
    copy_field(eightpiG_deltaT[0], rhs[0], a2_over_3);

    error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[0], sim, 0, "FMG_solver 3");

    if(error < sim.relaxation_error)
    {
      break;
    }
  }

  return error;
}


//////////////////////////
// TODO Details here
//////////////////////////
template <class FieldType>
double gamma_cycle(
  MultiField<FieldType> * phi,
  MultiField<FieldType> * xi,
  MultiField<FieldType> * xi_old,
  MultiField<FieldType> * laplace_xi,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * residual,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double * dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  long * numpts3d,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim,
  int level,
  string message = ""
)
{
  double
    error = 0.,
    trunc_error = 0.,
    check_fr_wrong;
  int max_level = sim.multigrid_n_grids-1;

  if(engine.isPartLayer(level))
  {
    // TODO Maybe not necessary on finest grid? CHECK
    build_laplacian(xi[level], laplace_xi[level], dx[level]);

    error = single_layer_solver(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[level], sim, level, message+" gamma_cycle level "+to_string(level), sim.pre_smoothing);

    if(level == 0 && error < sim.relaxation_error)
    {
      return error;
    }

    // Writes residual on residual[level]
    build_residual(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], rhs[level], residual[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
  }

  if(sim.multigrid_check_shape)
  {
    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << level << endl;
    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << " \\" << endl;
  }

  // Restrict from level to level+1
  if(sim.multigrid_restrict_mode == RESTRICT_XI)
  {
    engine.restrict(xi, level);
  }
  else
  {
    engine.restrict(deltaR, level);
  }

  engine.restrict(residual, level);
  engine.restrict(rhs, level); // TODO Maybe not necessary? Perhaps we can do it just once like for eightpiG_deltaT

  if(engine.isPartLayer(level+1))
  {
    if(sim.multigrid_restrict_mode == RESTRICT_XI)
    {
      check_fr_wrong = convert_xi_to_deltaR(deltaR[level+1], xi[level+1], Rbar, fRbar, sim);
    }
    else
    {
      check_fr_wrong = convert_deltaR_to_xi(deltaR[level+1], xi[level+1], Rbar, fRbar, sim, message+" gamma cycle 1");
    }

    if(check_fr_wrong > FR_WRONG)
    {
      cout << parallel.rank() << " 2 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
      return FR_WRONG_RETURN;
    }

    build_laplacian(xi[level+1], laplace_xi[level+1], dx[level+1]);
    build_residual(phi[level+1], xi[level+1], xi_old[level+1], laplace_xi[level+1], deltaR[level+1], rhs[level+1], trunc[level+1], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
    subtract_fields(trunc[level+1], residual[level+1], trunc[level+1]);
    add_fields(rhs[level+1], trunc[level+1], rhs[level+1]);
  }

  // Coarsest grid
  if(level+1 == max_level)
  {
    if(engine.isPartLayer(level+1))
    {
      single_layer_solver(phi[level+1], xi[level+1], xi_old[level+1], laplace_xi[level+1], deltaR[level+1], eightpiG_deltaT[level+1], rhs[level+1], dx[level+1], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[level+1], sim, level+1, message+" gamma_cycle level "+to_string(level));
    }
  }
  else
  {
    for(int j=0; j<sim.multigrid_shape; ++j)
    {
      check_fr_wrong = gamma_cycle(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, trunc, engine, dx, a2_over_3, two_Hubble_over_dtau, numpts3d, Rbar, fbar, fRbar, sim, level+1, message+" gamma_cycle level "+to_string(level));

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 4 check_fr_wrong = FR_WRONG in gamma_cycle 4" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }
  }

  if(sim.multigrid_restrict_mode == RESTRICT_XI)
  {
    engine.restrict(xi, residual, level);
    if(engine.isPartLayer(level+1))
    {
      subtract_fields(xi[level+1], residual[level+1], residual[level+1]);
    }
  }
  else
  {
    engine.restrict(deltaR, residual, level);
    if(engine.isPartLayer(level+1))
    {
      subtract_fields(deltaR[level+1], residual[level+1], residual[level+1]);
    }
  }


  if(sim.multigrid_check_shape)
  {
    for(int k=0; k<level+1; ++k) COUT << "  ";
    COUT << level+1 << endl;

    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << " /" << endl;
    if(!level) COUT << "0" << endl;
  }

  // if(engine.isPartLayer(level+1)) check_field(trunc[level+1], "trunc["+to_string(level+1)+"]", numpts3d[level+1], sim, "BEFORE PROLONGATION", level+1);

  engine.prolong(residual, trunc, level+1);

  // if(engine.isPartLayer(level)) check_field(trunc[level], "trunc["+to_string(level)+"]", numpts3d[level], sim, "AFTER PROLONGATION", level);

  if(engine.isPartLayer(level))
  {
    if(sim.multigrid_restrict_mode == RESTRICT_XI)
    {
      add_fields(xi[level], trunc[level], xi[level]);
      check_fr_wrong = convert_xi_to_deltaR(deltaR[level], xi[level], Rbar, fRbar, sim, message+" gamma_cycle level "+to_string(level));
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 5 check_fr_wrong = FR_WRONG in gamma_cycle 5" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }
    else
    {
      add_fields(deltaR[level], trunc[level], deltaR[level]);
      check_fr_wrong = convert_deltaR_to_xi(deltaR[level], xi[level], Rbar, fRbar, sim, message+" gamma_cycle 2");
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 6 check_fr_wrong = FR_WRONG in gamma_cycle 6" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }

    build_laplacian(xi[level], laplace_xi[level], dx[level]);

    if(engine.isPartLayer(level))
    {
      error = single_layer_solver(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[level], sim, level, message+" gamma_cycle level "+to_string(level), sim.post_smoothing);

      if(level == 0 && error < sim.relaxation_error)
      {
        return error;
      }

      // Writes residual on residual[level]
      build_residual(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], rhs[level], residual[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
    }
  }

  trunc_error = euclidean_norm(trunc[level], engine, level, numpts3d[level]);
  return trunc_error;
}

#endif
