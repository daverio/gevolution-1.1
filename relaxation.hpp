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

#ifndef RELAXATION_HEADER
#define RELAXATION_HEADER

//////////////////////////
// Euclidean Norm
//////////////////////////
template <class FieldType>
double euclidean_norm(
  Field<FieldType> & field,
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
  }
  parallel.layer(engine.player(level)).sum(error);
  error = sqrt(error) / numpts3d;

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
template <class FieldType>
int prepare_initial_guess_trace_equation(
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & xi,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & rhs,
  double a,
  double dx,
  double Rbar,
  double fRbar,
  const metadata & sim)
{
  if(sim.fR_model == FR_MODEL_HU_SAWICKI)
  {
    erase_field(xi);
    erase_field(deltaR);
    erase_field(laplace_xi);
  }
  else if(sim.fR_model == FR_MODEL_RN || sim.fR_model == FR_MODEL_DELTA)
  {
    erase_field(xi);
    erase_field(deltaR);
    erase_field(laplace_xi);

    // copy_field(eightpiG_deltaT, deltaR, -1.);
    // convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim);
    // build_laplacian(xi, laplace_xi, dx);
  }

  copy_field(eightpiG_deltaT, rhs, a*a/3.);

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//                       ROUTINES AND TOOLS FOR xi
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// TODO Details?

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
  const metadata & sim)
{
  double temp = deltaR;

  if(sim.relativistic_flag || a2_over_3 * temp + rhs == 0.)
  {
    // TODO Additional terms -- see if necessary
    temp = temp * (1. - fRbar - xi) + 2. * (f(Rbar + deltaR, sim, 884) - fbar) - Rbar * xi;
    // End additional terms
  }

  temp = laplace_xi - a2_over_3 * temp - rhs;

  if(sim.xi_Hubble)
  {
    temp += 2. * phi * laplace_xi - two_Hubble_over_dtau * (xi - xi_old);
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
  const metadata & sim)
{
  double
    R = Rbar + deltaR,
    temp = fRR(R, sim, 738);

  if(!temp || temp > FR_WRONG)
  {
    cout << endl;
    cout << " dresidual_dxi" << endl;
    cout << " deltaR " << deltaR << endl;
    cout << " xi " << xi << endl;
    cout << " Rbar " << Rbar << endl;
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

  if(isnan(deltaR) || isnan(xi) || Rbar + deltaR <= 0. || fabs(xi) >= 1.E20)
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
  denominator = dresidual_dxi(xi, deltaR, a2_over_3, two_Hubble_over_dtau, Rbar, fRbar, sim);

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
            resid *= 0.5;
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
          resid *= 0.5;
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
  double xi_previous_timestep,
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
  temp = update_xi_single(phi, xi, xi_previous_timestep, laplace_xi, deltaR, rhs, a2_over_3, two_Hubble_over_dtau, coeff_laplacian, overrelax, Rbar, fbar, fRbar, sim);

  if(isnan(deltaR) || fabs(xi) >= 1.E20 || deltaR + Rbar <= 0.)
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

  if(isnan(xi) || fabs(xi) > 1.E20)
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


//////////////////////////
// Computes new guess for xi
// Solves the equation Y[xi] == 0
// where Y[xi] := laplacian(xi) - a^2/3 * zeta(xi)
// TODO: Details here
//////////////////////////
template <class FieldType>
double update_xi_and_deltaR(
  Field<FieldType> & phi,
  Field<FieldType> & xi,
  Field<FieldType> & xi_previous_timestep,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double temp, resid, coeff_laplacian, overrelax;

  coeff_laplacian = -6./dx/dx;
  overrelax = sim.overrelaxation_factor;

  if(sim.red_black)
  {
    SiteRedBlack3d x(xi.lattice());

    for(x.firstRed(); x.testRed(); x.nextRed())
    {
      update_xi_and_deltaR_single(phi(x), xi(x), xi_previous_timestep(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, coeff_laplacian, overrelax, Rbar, fbar, fRbar, sim);
    }

    build_laplacian(xi, laplace_xi, dx); // TODO: This might be optimised somehow -- only build_laplacian on the Blacks?

    for(x.firstBlack(); x.testBlack(); x.nextBlack())
    {
      update_xi_and_deltaR_single(phi(x), xi(x), xi_previous_timestep(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, coeff_laplacian, overrelax, Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    Site x(xi.lattice());

    for(x.first(); x.test(); x.next())
    {
      if(isnan(xi(x)) || fabs(xi(x)) > 1.E20)
      {
        cout << endl;
        cout << " update_xi_and_deltaR before" << endl;
        cout << " xi = " << xi(x) << endl;
        cout << " dxi = " << laplace_xi(x) << endl;
        cout << " Rbar = " << Rbar << endl;
        cout << " deltaR = " << deltaR(x) << endl;
        cout << " R = " << Rbar + deltaR(x) << endl;
        cout << endl;
        parallel.abortForce();
      }

      update_xi_and_deltaR_single(phi(x), xi(x), xi_previous_timestep(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, coeff_laplacian, overrelax, Rbar, fbar, fRbar, sim);

      if(isnan(xi(x)) || fabs(xi(x)) > 1.E20)
      {
        cout << endl;
        cout << " update_xi_and_deltaR after xi" << endl;
        cout << " xi = " << xi(x) << endl;
        cout << " dxi = " << laplace_xi(x) << endl;
        cout << " Rbar = " << Rbar << endl;
        cout << " deltaR = " << deltaR(x) << endl;
        cout << " R = " << Rbar + deltaR(x) << endl;
        cout << endl;
        parallel.abortForce();
      }
    }
  }

  build_laplacian(xi, laplace_xi, dx);

  // TODO: Why doesn't it work without this (apparently useless) call?
  temp = compute_error(phi, xi, xi_previous_timestep, laplace_xi, deltaR, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, 1, sim);

  return temp;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//                        OTHER TOOLS AND ROUTINES
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double compute_error(
  Field<FieldType> & phi,
  Field<FieldType> & xi,
  Field<FieldType> & xi_old,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim)
{
  double temp, error = 0.;

  Site x(xi.lattice());

  for(x.first(); x.test(); x.next())
  {
    temp = residual_xi(phi(x), xi(x), xi_old(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
    error += temp*temp;
  }

  parallel.sum(error);
  error = sqrt(error / numpts3d);

  return error;
}


//////////////////////////
// TODO Details
//////////////////////////
template <class FieldType>
void build_residual(
  Field<FieldType> & phi,
  Field<FieldType> & xi,
  Field<FieldType> & xi_old,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  Field<FieldType> & destination_for_residual,
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

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//                                SOLVERS
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////
// Relaxation solver
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation(
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
  const metadata & sim)
{
  double error;
  double a2_over_3 = a*a/3.;
  int max_level = sim.multigrid_n_grids - 1; // TODO Not really needed, just for clarity

  COUT << " z = " << 1./a - 1. << " -- " << flush;
  error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim);
  COUT << "in.error = " << error << " -- " << flush;

  if(sim.relaxation_method == METHOD_RELAX)
  {
    error = single_layer_solver(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], eightpiG_deltaT[0], rhs[0], dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim);
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

  COUT << "fin.error = " << error << endl;

  // build_residual(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], err[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
  // check_field(err[0], "residual", numpts3d, sim, "Final Residual");

  return error;
}


//////////////////////////
// Single relaxation step
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation_step(
  Field<FieldType> & phi,
  Field<FieldType> & xi,
  Field<FieldType> & xi_old,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim)
{
  double error;
  // Attempt new guess
  error = update_xi_and_deltaR(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);

  if(error > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << " relaxation_step Rbar = " << Rbar << endl;
    //END REMOVE
    return FR_WRONG_RETURN;
  }
  else
  {
    error = compute_error(phi, xi, xi_old, laplace_xi, deltaR, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim);
    return error;
  }
}



//////////////////////////
// Relaxation solver
// TODO: Details here
//////////////////////////
template <class FieldType>
double single_layer_solver(
  Field<FieldType> & phi,
  Field<FieldType> & xi,
  Field<FieldType> & xi_old,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim
)
{
  int count = 0;
  double error = 0., previous_error = 1.;

  if(sim.truncate_relaxation)
  {
    while(true)
    {
      error = relaxation_step(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim);

      if(error < sim.relaxation_error)
      {
        break;
      }

      if(fabs(error/previous_error - 1.) <= sim.relaxation_truncation_threshold)
      {
        ++count;
      }
      else
      {
        count = 0;
      }

      if(count >= sim.truncate_relaxation)
      {
        break;
      }

      previous_error = error;
    }
  }
  else
  {
    while(true)
    {
      error = relaxation_step(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d, sim);

      if(error < sim.relaxation_error)
      {
        break;
      }
    }
  }

  COUT << "(" << count << ") ";
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
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double DX,
  double a2_over_3,
  double two_Hubble_over_dtau,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  int
    max_level = sim.multigrid_n_grids-1;

  double
    dx[sim.multigrid_n_grids],
    error = 0.,
    trunc_error = 0.,
    previous_error = 0.;

  long numpts3d[sim.multigrid_n_grids];

  numpts3d[0] = (long) sim.numpts * sim.numpts * sim.numpts;
  dx[0] = DX;
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

    trunc_error = gamma_cycle(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, dx, a2_over_3, two_Hubble_over_dtau, numpts3d, Rbar, fbar, fRbar, sim, 0);

    if(trunc_error > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 1 trunc_error = FR_WRONG in multigrid_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    convert_xi_to_deltaR(deltaR[0], xi[0], Rbar, fRbar, sim);

    error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[0], sim);

    // TODO Check if this criterion makes sense -- test more
    if(3. * error / trunc_error <= 1. || error < sim.relaxation_error)
    {
      break;
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
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double DX,
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
  dx[0] = DX;
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
    error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[0], sim);

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
        check_fr_wrong = update_xi_and_deltaR(phi[max_level], xi[max_level], xi_old[max_level], laplace_xi[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[max_level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);

        if(check_fr_wrong > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 11 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
          //END REMOVE
          return FR_WRONG_RETURN;
        }

        error = compute_error(phi[max_level], xi[max_level], xi_old[max_level], laplace_xi[max_level], deltaR[max_level], rhs[max_level], dx[max_level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[max_level], sim);

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
        check_fr_wrong = gamma_cycle(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, dx, a2_over_3, two_Hubble_over_dtau, numpts3d, Rbar, fbar, fRbar, sim, i);

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

    error = compute_error(phi[0], xi[0], xi_old[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[0], sim);

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
  MultiField<FieldType> * temp,
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
  int level)
{
  double
    error = 0.,
    trunc_error = 0.,
    check_fr_wrong;
  int max_level = sim.multigrid_n_grids-1;


  if(engine.isPartLayer(level))
  {
    if(level) // Not necessary on finest grid? TODO CHECK
    {
      build_laplacian(xi[level], laplace_xi[level], dx[level]);
    }

    // TODO REMOVE
    // error = compute_error(xi[level], laplace_xi[level], deltaR[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, numpts3d[level], sim);
    // for(int jj=0; jj<level; ++jj) COUT << "     ";
    // COUT << " level " << level << " - i. error = " << error << endl;
    // END REMOVE

    for(int j=0; j<sim.pre_smoothing; ++j)
    {
      check_fr_wrong = update_xi_and_deltaR(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
        return FR_WRONG_RETURN;
      }
    }

    build_residual(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], rhs[level], temp[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
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

  engine.restrict(temp, level);
  engine.restrict(rhs, level); // TODO Maybe not necessary? Perhaps we can do it just once like for eightpiG_deltaT

  if(engine.isPartLayer(level+1))
  {
    if(sim.multigrid_restrict_mode == RESTRICT_XI)
    {
      check_fr_wrong = convert_xi_to_deltaR(deltaR[level+1], xi[level+1], Rbar, fRbar, sim);
    }
    else
    {
      check_fr_wrong = convert_deltaR_to_xi(deltaR[level+1], xi[level+1], Rbar, fRbar, sim);
    }

    if(check_fr_wrong > FR_WRONG)
    {
      cout << parallel.rank() << " 2 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
      return FR_WRONG_RETURN;
    }

    build_laplacian(xi[level+1], laplace_xi[level+1], dx[level+1]);
    build_residual(phi[level+1], xi[level+1], xi_old[level+1], laplace_xi[level+1], deltaR[level+1], rhs[level+1], trunc[level+1], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
    subtract_fields(trunc[level+1], temp[level+1], trunc[level+1]);
    add_fields(rhs[level+1], trunc[level+1], rhs[level+1]);
  }

  // Coarsest grid
  if(level+1 == max_level)
  {
    if(engine.isPartLayer(max_level))
    {
      error = single_layer_solver(phi[max_level], xi[max_level], xi_old[max_level], laplace_xi[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[max_level], sim);
    }
  }
  else
  {
    for(int j=0; j<sim.multigrid_shape; ++j)
    {
      check_fr_wrong = gamma_cycle(phi, xi, xi_old, laplace_xi, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, dx, a2_over_3, two_Hubble_over_dtau, numpts3d, Rbar, fbar, fRbar, sim, level+1);
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
    engine.restrict(xi, temp, level);
    if(engine.isPartLayer(level+1))
    {
      subtract_fields(xi[level+1], temp[level+1], temp[level+1]);
    }
  }
  else
  {
    engine.restrict(deltaR, temp, level);
    if(engine.isPartLayer(level+1))
    {
      subtract_fields(deltaR[level+1], temp[level+1], temp[level+1]);
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

  engine.prolong(temp, trunc, level+1);
  if(engine.isPartLayer(level))
  {
    if(sim.multigrid_restrict_mode == RESTRICT_XI)
    {
      add_fields(xi[level], trunc[level], xi[level]);
      check_fr_wrong = convert_xi_to_deltaR(deltaR[level], xi[level], Rbar, fRbar, sim);
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
      check_fr_wrong = convert_deltaR_to_xi(deltaR[level], xi[level], Rbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 6 check_fr_wrong = FR_WRONG in gamma_cycle 6" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }

    build_laplacian(xi[level], laplace_xi[level], dx[level]);

    // TODO In place of post-smoothing: requires that error be smaller than relaxation_error at each layer
    for(int j=0; j<sim.post_smoothing; ++j)
    {
      check_fr_wrong = update_xi_and_deltaR(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
        return FR_WRONG_RETURN;
      }
    }

    error = compute_error(phi[level], xi[level], xi_old[level], laplace_xi[level], deltaR[level], rhs[level], dx[level], a2_over_3, two_Hubble_over_dtau, Rbar, fbar, fRbar, numpts3d[level], sim);

    // TODO REMOVE
    // for(int jj=0; jj<level; ++jj) COUT << "     ";
    // COUT << " level " << level << " - f. error = " << error << endl;
    // END REMOVE
  }

  trunc_error = euclidean_norm(trunc[level], engine, level, numpts3d[level]);

  return trunc_error;
}

#endif
