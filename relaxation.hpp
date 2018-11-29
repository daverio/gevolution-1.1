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
// Initial guess for deltaR in relaxation solver
// TODO: add description
//////////////////////////
template <class FieldType>
void initial_guess_deltaR(
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  const metadata & sim) // TODO not necessary for these initial conditions, remove if not needed
{
  copy_field(eightpiG_deltaT, deltaR, -1.);
  // zero_field(deltaR);
  return;
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
// Prepare the initial conditions for the trace equation
// TODO: Details here
//////////////////////////
template <class FieldType>
int prepare_initial_conditions_trace_equation(
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & scalaron,
  Field<FieldType> & laplace_scalaron,
  Field<FieldType> & rhs,
  double a,
  double dx,
  double Rbar,
  double fRbar,
  const metadata & sim)
{
  // Build initial guess for deltaR -- common to all methods
  initial_guess_deltaR(deltaR, eightpiG_deltaT, sim);
  // Prepare source term = a^2 * eightpiG_deltaT / 3.
  copy_field(eightpiG_deltaT, rhs, a*a/3.);

  convert_deltaR_to_scalaron(deltaR, scalaron, Rbar, fRbar, sim); // deltaR and u are now consistent
  build_laplacian_scalaron(scalaron, laplace_scalaron, dx, sim);

  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//                       ROUTINES AND TOOLS FOR XI
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// TODO Details?

//////////////////////////
// Residual Y
// of equation Y[xi] == 0
// TODO: Details here
//////////////////////////
double residual_xi(
  double xi,
  double laplace_xi,
  double deltaR,
  double rhs,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double R = Rbar + deltaR,
         temp = deltaR;

  // TODO Additional terms -- see if necessary
  temp = temp * (1. - fRbar - xi) + 2. * (f(R, sim, 884) - fbar) - Rbar * xi;
  // End additional terms

  return laplace_xi - a2_over_3 * temp - rhs;
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
  double Rbar,
  double fRbar,
  const metadata & sim)
{
  double R = Rbar + deltaR,
         temp = 1. / fRR(R, sim, 738);

  // TODO Additional terms -- see if necessary
  temp = temp * (1. + fRbar + xi) - R;
  // End additional terms

  return - a2_over_3 * temp;
}

//////////////////////////
// Computes new guess for xi
// Solves the equation Y[xi] == 0
// where Y[xi] := laplacian(xi) - a^2/3 * zeta(xi)
// TODO: Details here
//////////////////////////
double update_xi_single(
  double & xi,
  double laplace_xi,
  double deltaR,
  double rhs,
  double a2_over_3,
  double coeff,
  double overrelax,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double fr,
         denomin,
         resid;

  fr = xi + fRbar;
  denomin = dresidual_dxi(xi, deltaR, a2_over_3, Rbar, fRbar, sim) + coeff;

  if(denomin && fr)
  {
    resid = residual_xi(xi, laplace_xi, deltaR, rhs, a2_over_3, Rbar, fbar, fRbar, sim);
    resid = - overrelax * resid / denomin;
    xi = fr * exp(resid / fr) - fRbar;
    return resid;
  }
  else if(fr)
  {
    cout << " denomin evaluated to 0 in update_xi_single -- xi, deltaR, fr = " << xi << ", " << deltaR << ", " << fr << endl;
    return 0.;
  }
  else
  {
    cout << " fr evaluated to 0 in update_xi_single -- xi, deltaR, denomin = " << xi << ", " << deltaR << ", " << denomin << endl;
    return 0.;
  }

}

//////////////////////////
// Computes new guess for xi
// Solves the equation Y[xi] == 0
// where Y[xi] := laplacian(xi) - a^2/3 * zeta(xi)
// TODO: Details here
//////////////////////////
template <class FieldType>
double update_xi_and_deltaR(
  Field<FieldType> & xi,
  Field<FieldType> & laplace_xi,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double temp,
         d_xi,
         coeff = -6./dx/dx,
         overrelax = sim.overrelaxation_coeff,
         resid;

  if(sim.red_black)
  {
    SiteRedBlack3d x(xi.lattice());

    for(x.firstRed(); x.testRed(); x.nextRed())
    {
      update_xi_single(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }

    build_laplacian(xi, laplace_xi, dx); // TODO: This might be optimised somehow -- only build_laplacian on the Blacks?

    for(x.firstBlack(); x.testBlack(); x.nextBlack())
    {
      update_xi_single(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    Site x(xi.lattice());
    for(x.first(); x.test(); x.next())
    {
      update_xi_single(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }
  }

  build_laplacian(xi, laplace_xi, dx);

  temp = convert_xi_to_deltaR(eightpiG_deltaT, deltaR, xi, Rbar, fRbar, sim);

  return temp;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//                       ROUTINES AND TOOLS FOR U
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// TODO Details?

//////////////////////////
// in-place calculation of residual Y
// TODO: Details here
//////////////////////////
double residual_u(
  double u,
  double laplace_exp_u,
  double deltaR,
  double rhs,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double R = Rbar + deltaR,
         temp = deltaR,
         xi = fRbar * (exp(u) - 1.);

  // TODO Additional terms -- see if necessary
  temp = temp * (1. - fRbar - xi) + 2. * (f(R, sim, 884) - fbar) - Rbar * xi;
  // End additional terms

  return fRbar * laplace_exp_u - a2_over_3 * temp - rhs;
}


//////////////////////////
// In-place calculation of (dY/du)
// TODO Correct this above!
//////////////////////////
double dresidual_du(
  double u,
  double deltaR,
  double a2_over_3,
  double Rbar,
  double fr,
  const metadata & sim)
{
  double R = Rbar + deltaR,
         temp = 1. / fRR(R, sim, 738);

  // TODO Additional terms -- see if necessary
  temp = temp * (1. + fr) - R;
  // End additional terms

  return - a2_over_3 * temp * fr;
}


//////////////////////////
// Computes new guess for u
// Solves the equation Y[u] == 0
// where Y[u] := laplacian( fRbar * exp(u) ) - a^2/3 * zeta(u)
// TODO: Details here
//////////////////////////
template <class FieldType>
double update_u_and_deltaR(
  Field<FieldType> & u,
  Field<FieldType> & laplace_exp_u,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double temp,
         d_xi,
         coeff = -6./dx/dx,
         overrelax = sim.overrelaxation_coeff,
         resid;

  if(sim.red_black)
  {
    SiteRedBlack3d x(u.lattice());

    for(x.firstRed(); x.testRed(); x.nextRed())
    {
      update_u_single(u(x), laplace_exp_u(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }

    build_laplacian_exp(u, laplace_exp_u, dx); // TODO: This might be optimised somehow -- only build_laplacian on the Blacks?

    for(x.firstBlack(); x.testBlack(); x.nextBlack())
    {
      update_u_single(u(x), laplace_exp_u(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    Site x(u.lattice());
    for(x.first(); x.test(); x.next())
    {
      update_u_single(u(x), laplace_exp_u(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }
  }

  build_laplacian_exp(u, laplace_exp_u, dx);

  temp = convert_u_to_deltaR(eightpiG_deltaT, deltaR, u, Rbar, fRbar, sim);

  return temp;
}


//////////////////////////
// Computes new guess for u
// Solves the equation Y[u] == 0
// where Y[u] := fRbar * laplacian(exp(u)) - a^2/3. * zeta(u)
// TODO: Details here
//////////////////////////
double update_u_single(
  double & u,
  double laplace_exp_u,
  double deltaR,
  double rhs,
  double a2_over_3,
  double coeff,
  double overrelax,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double fr,
         denomin,
         resid;

  fr = fRbar * exp(u);
  // TODO CHECK IF CORRECT
  denomin = dresidual_du(u, deltaR, a2_over_3, Rbar, fr, sim) + fr * coeff;


  if(denomin)
  {
    resid = residual_u(u, laplace_exp_u, deltaR, rhs, a2_over_3, Rbar, fbar, fRbar, sim);
    u -= overrelax * resid / denomin;
  }
  else
  {
    cout << " denomin evaluated to 0. -- u, deltaR, coeff = " << u << ", " << deltaR << ", " << coeff << endl;
    resid = 0.;
  }

  return resid;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//              TOOLS AND ROUTINES GENERIC FOR BOTH XI AND U
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// TODO Details?

//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double compute_error(
  Field<FieldType> & scalaron,
  Field<FieldType> & laplace_scalaron,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim)
{
  double temp,
         error = 0.;

  Site x(scalaron.lattice());

  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    for(x.first(); x.test(); x.next())
    {
      temp = residual_u(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, Rbar, fbar, fRbar, sim);
      error += temp*temp;
    }
  }
  else
  {
    for(x.first(); x.test(); x.next())
    {
      temp = residual_xi(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, Rbar, fbar, fRbar, sim);
      error += temp*temp;
    }
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
  Field<FieldType> & scalaron,
  Field<FieldType> & laplace_scalaron,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  Field<FieldType> & destination_for_residual,
  MultiGrid & engine,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{

  Site x(scalaron.lattice());

  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    for(x.first(); x.test(); x.next())
    {
      destination_for_residual(x) = residual_u(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    for(x.first(); x.test(); x.next())
    {
      destination_for_residual(x) = residual_xi(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, Rbar, fbar, fRbar, sim);
    }
  }

  return;
}

//////////////////////////
// Computes new guess for scalaron
// TODO: Details here
//////////////////////////
double update_scalaron_single(
  double & scalaron,
  double laplace_scalaron,
  double deltaR,
  double rhs,
  double a2_over_3,
  double coeff,
  double overrelax,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double t;

  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    t = update_u_single(scalaron, laplace_scalaron, deltaR, rhs, a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
  }
  else
  {
    t = update_xi_single(scalaron, laplace_scalaron, deltaR, rhs, a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
  }

  return t;
}


//////////////////////////
// Computes new guess for scalaron
// TODO: Details here
//////////////////////////
template <class FieldType>
double update_scalaron_and_deltaR(
  Field<FieldType> & scalaron,
  Field<FieldType> & laplace_scalaron,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double temp,
         coeff = -6./dx/dx,
         overrelax = sim.overrelaxation_coeff;

  if(sim.relaxation_variable == RELAXATION_VAR_U)
  {
    if(sim.red_black)
    {
      SiteRedBlack3d x(scalaron.lattice());

      for(x.firstRed(); x.testRed(); x.nextRed())
      {
        update_u_single(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
      }

      build_laplacian_exp(scalaron, laplace_scalaron, dx); // TODO: This might be optimised somehow -- only build_laplacian on the Blacks?

      for(x.firstBlack(); x.testBlack(); x.nextBlack())
      {
        update_u_single(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
      }
    }
    else
    {
      Site x(scalaron.lattice());
      for(x.first(); x.test(); x.next())
      {
        update_u_single(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
      }
    }

    build_laplacian_exp(scalaron, laplace_scalaron, dx);
    temp = convert_u_to_deltaR(eightpiG_deltaT, deltaR, scalaron, Rbar, fRbar, sim);
  }
  else
  {
    if(sim.red_black)
    {
      SiteRedBlack3d x(scalaron.lattice());

      for(x.firstRed(); x.testRed(); x.nextRed())
      {
        update_xi_single(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
      }

      build_laplacian(scalaron, laplace_scalaron, dx); // TODO: This might be optimised somehow -- only build_laplacian on the Blacks?

      for(x.firstBlack(); x.testBlack(); x.nextBlack())
      {
        update_xi_single(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
      }
    }
    else
    {
      Site x(scalaron.lattice());
      for(x.first(); x.test(); x.next())
      {
        update_xi_single(scalaron(x), laplace_scalaron(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
      }
    }

    build_laplacian(scalaron, laplace_scalaron, dx);
    temp = convert_xi_to_deltaR(eightpiG_deltaT, deltaR, scalaron, Rbar, fRbar, sim);
  }

return temp;
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
int relaxation(
  MultiField<FieldType> * scalaron,
  MultiField<FieldType> * laplace_scalaron,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * residual,
  MultiField<FieldType> * err,
  MultiGrid & engine,
  double a,
  double dx,
  long numpts3d,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double error;
  int max_level = sim.multigrid_n_grids - 1; // TODO Not really needed, just for clarity

  COUT << " z = " << 1./a - 1. << endl;

  error = compute_error(scalaron[0], laplace_scalaron[0], deltaR[0], rhs[0], dx, a*a/3., Rbar, fbar, fRbar, numpts3d, sim);

  if(error < sim.relaxation_error)
  {
    COUT << "  initial and final error = " << error << endl;
    return 0;
  }
  else
  {
    COUT << " initial error = " << error << endl;
  }

  if(sim.relaxation_method == METHOD_RELAX)
  {
    error = single_layer_solver(scalaron[0], laplace_scalaron[0], deltaR[0], eightpiG_deltaT[0], rhs[0], a*a/3., dx, numpts3d, Rbar, fbar, fRbar, sim);
  }
  else if(sim.relaxation_method == METHOD_MULTIGRID)
  {
    error = multigrid_solver(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, residual, err, engine, a*a/3., dx, Rbar, fbar, fRbar, sim);
  }
  else if(sim.relaxation_method == METHOD_FMG)
  {
    error = FMG_solver(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, residual, err, engine, a, dx, Rbar, fbar, fRbar, sim);
  }

  // TODO REMOVE after debug
  if(error > FR_WRONG)
  {
    COUT << " Method " << sim.relaxation_method << " returned FR_WRONG_RETURN. Check what's going on. Press [Enter] to continue..." << endl;
    cin.get();
  }
  // END REMOVE

  COUT << "   final error = " << error << endl;
  return 0;
}


//////////////////////////
// Single relaxation step
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation_step(
  Field<FieldType> & scalaron,
  Field<FieldType> & laplace_scalaron,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double dx,
  double a2_over_3,
  double Rbar,
  double fbar,
  double fRbar,
  long numpts3d,
  const metadata & sim)
{
  double error;
  // Attempt new guess
  error = update_scalaron_and_deltaR(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, Rbar, fbar, fRbar, sim);

  if(error > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << parallel.rank() << " 8 check_fr_wrong = FR_WRONG in relaxation_step" << endl;
    //END REMOVE
    return FR_WRONG_RETURN;
  }

  error = compute_error(scalaron, laplace_scalaron, deltaR, rhs, dx, a2_over_3, Rbar, fbar, fRbar, numpts3d, sim);

  return error;
}



//////////////////////////
// Relaxation solver
// TODO: Details here
//////////////////////////
template <class FieldType>
double single_layer_solver(
  Field<FieldType> & scalaron,
  Field<FieldType> & laplace_scalaron,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  double a2_over_3,
  double dx,
  long numpts3d,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim,
  double target_error = -1.)
{
  double error;

  if(target_error <= 0.)
  {
    target_error = sim.relaxation_error;
  }

  while(true)
  {
    error = relaxation_step(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, Rbar, fbar, fRbar, numpts3d, sim);

    if(error < target_error)
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
double multigrid_solver(
  MultiField<FieldType> * scalaron,
  MultiField<FieldType> * laplace_scalaron,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double a2_over_3,
  double DX,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  long numpts3d[sim.multigrid_n_grids];
  int max_level = sim.multigrid_n_grids-1; // TODO Not necessary, but good for clarity
  double dx[sim.multigrid_n_grids],
         error = 0.,
         trunc_error = 0.;

  numpts3d[0] = (long) sim.numpts * sim.numpts * sim.numpts;
  dx[0] = DX;
  for(int i=1; i<sim.multigrid_n_grids; ++i)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8;
    dx[i] = dx[i-1] * 2.;
  }

  // Build eightpiG_deltaT on each level (it will not change)
  restrict_to_level(eightpiG_deltaT, engine, max_level);

  while(true)
  {
    trunc_error = gamma_cycle(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a2_over_3, dx, numpts3d, Rbar, fbar, fRbar, sim, 0);

    if(trunc_error > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 1 trunc_error = FR_WRONG in multigrid_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    convert_scalaron_to_deltaR(eightpiG_deltaT[0], deltaR[0], scalaron[0], Rbar, fRbar, sim);

    error = compute_error(scalaron[0], laplace_scalaron[0], deltaR[0], rhs[0], dx[0], a2_over_3, Rbar, fbar, fRbar, numpts3d[0], sim);

    // TODO Check if this criterion makes sense -- test more
    if(3. * error / trunc_error <= 1. || error < sim.relaxation_error)
    {
      COUT << " (1/3) * trunc_error = " << trunc_error / 3. << endl;
      COUT << "               error = " << error << endl;
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
  MultiField<FieldType> * scalaron,
  MultiField<FieldType> * laplace_scalaron,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double a2_over_3,
  double DX,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  int i,
      j,
      max_level = sim.multigrid_n_grids - 1; // Not exactly needed, but more redable this way
  long numpts3d[sim.multigrid_n_grids];
  double dx[sim.multigrid_n_grids],
         error,
         check_fr_wrong;

  numpts3d[0] = (long) (sim.numpts * sim.numpts * sim.numpts);
  dx[0] = DX;
  for(i=1; i<sim.multigrid_n_grids; ++i)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8.;
    dx[i] = 2. * dx[i-1];
  }

  // Build eightpiG_deltaT on each level (it will not change during the relaxation)
  restrict_to_level(eightpiG_deltaT, engine, max_level);
  if(engine.isPartLayer(max_level))
  {
    copy_field(eightpiG_deltaT[max_level], rhs[max_level], a2_over_3);
  }
  copy_field(eightpiG_deltaT[0], rhs[0], a2_over_3);
  // build initial guess for u on coarsest grid
  restrict_to_level(scalaron, engine, max_level);

  // Compute initial error
  // TODO: remove after debugging
  build_laplacian_exp(scalaron[0], laplace_scalaron[0], dx[0]);
  error = compute_error(scalaron[0], laplace_scalaron[0], deltaR[0], rhs[0], dx[0], a2_over_3, Rbar, fbar, fRbar, numpts3d[0], sim);

  // solve problem on coarsest grid
  if(engine.isPartLayer(max_level))
  {
    build_laplacian_exp(scalaron[max_level], laplace_scalaron[max_level], dx[max_level]);

    check_fr_wrong = convert_scalaron_to_deltaR(eightpiG_deltaT[max_level], deltaR[max_level], scalaron[max_level], Rbar, fRbar, sim);
    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 10 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    while(true)
    {
      check_fr_wrong = update_scalaron_and_deltaR(scalaron[max_level], laplace_scalaron[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[max_level], a2_over_3, Rbar, fbar, fRbar, sim);

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 11 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }

      error = compute_error(scalaron[max_level], laplace_scalaron[max_level], deltaR[max_level], rhs[max_level], dx[max_level], a2_over_3, Rbar, fbar, fRbar, numpts3d[max_level], sim);

      if(error < sim.relaxation_error)
      {
        break;
      }
    }
  }

  for(i=sim.multigrid_n_grids-2; i>=0; i--) // i is the starting level for the V or W mini-cycle
  {
    // prolong i+1 solution
    engine.prolong(scalaron, i+1);
    // build/copy source to the correct field (rhs)
    if(engine.isPartLayer(i))
    {
      copy_field(eightpiG_deltaT[i], rhs[i], a2_over_3);
    }

    for(j=0; j<sim.multigrid_n_cycles; ++j)
    {
      check_fr_wrong = gamma_cycle(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a2_over_3, dx, numpts3d, Rbar, fbar, fRbar, sim, i);
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

  build_laplacian_exp(scalaron[0], laplace_scalaron[0], dx[0]);
  check_fr_wrong = convert_u_to_deltaR(eightpiG_deltaT[0], deltaR[0], scalaron[0], Rbar, fRbar, sim);
  if(check_fr_wrong > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << parallel.rank() << " 13 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
    //END REMOVE
    return FR_WRONG_RETURN;
  }

  copy_field(eightpiG_deltaT[0], rhs[0], a2_over_3);
  error = compute_error(scalaron[0], laplace_scalaron[0], deltaR[0], rhs[0], dx[0], a2_over_3, Rbar, fbar, fRbar, numpts3d[0], sim);

  return error; // TODO: Compute and divide by truncation error here! Fix!
}


//////////////////////////
// TODO Details here
//////////////////////////
template <class FieldType>
double gamma_cycle(
  MultiField<FieldType> * scalaron,
  MultiField<FieldType> * laplace_scalaron,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  double a2_over_3,
  double * dx,
  long * numpts3d,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim,
  int level)
{
  double error = 0.,
         trunc_error = 0.,
         check_fr_wrong;
  int max_level = sim.multigrid_n_grids-1;


  if(engine.isPartLayer(level))
  {
    if(level) // Not necessary on finest grid? TODO CHECK
    {
      build_laplacian_scalaron(scalaron[level], laplace_scalaron[level], dx[level], sim);
    }

    // TODO REMOVE
    // error = compute_error(scalaron[level], laplace_scalaron[level], deltaR[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, numpts3d[level], sim);
    // for(int jj=0; jj<level; ++jj) COUT << "     ";
    // COUT << " level " << level << " - i. error = " << error << endl;
    // END REMOVE

    for(int j=0; j<sim.pre_smoothing; ++j)
    {
      check_fr_wrong = update_scalaron_and_deltaR(scalaron[level], laplace_scalaron[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
        return FR_WRONG_RETURN;
      }
    }

    build_residual(scalaron[level], laplace_scalaron[level], deltaR[level], rhs[level], temp[level], engine, a2_over_3, Rbar, fbar, fRbar, sim);
  }

  if(sim.multigrid_check_shape)
  {
    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << level << endl;
    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << " \\" << endl;
  }

  // Restrict from level to level+1
  if(sim.multigrid_restrict_mode == RESTRICT_SCALARON)
  {
    engine.restrict(scalaron, level);
  }
  else
  {
    engine.restrict(deltaR, level);
  }

  engine.restrict(temp, level);
  engine.restrict(rhs, level); // TODO Maybe not necessary? Perhaps we can do it just once like for eightpiG_deltaT

  if(engine.isPartLayer(level+1))
  {
    if(sim.multigrid_restrict_mode == RESTRICT_SCALARON)
    {
      check_fr_wrong = convert_scalaron_to_deltaR(eightpiG_deltaT[level+1], deltaR[level+1], scalaron[level+1], Rbar, fRbar, sim);
    }
    else
    {
      check_fr_wrong = convert_deltaR_to_scalaron(deltaR[level+1], scalaron[level+1], Rbar, fRbar, sim);
    }

    if(check_fr_wrong > FR_WRONG)
    {
      cout << parallel.rank() << " 2 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
      return FR_WRONG_RETURN;
    }

    build_laplacian_scalaron(scalaron[level+1], laplace_scalaron[level+1], dx[level+1], sim);
    build_residual(scalaron[level+1], laplace_scalaron[level+1], deltaR[level+1], rhs[level+1], trunc[level+1], engine, a2_over_3, Rbar, fbar, fRbar, sim);
    subtract_fields(trunc[level+1], temp[level+1], trunc[level+1]);
    add_fields(rhs[level+1], trunc[level+1], rhs[level+1]);
  }

  // Coarsest grid
  if(level+1 == max_level)
  {
    if(engine.isPartLayer(max_level))
    {
      error = single_layer_solver(scalaron[max_level], laplace_scalaron[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], a2_over_3, dx[level], numpts3d[max_level], Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    for(int j=0; j<sim.multigrid_shape; ++j)
    {
      check_fr_wrong = gamma_cycle(scalaron, laplace_scalaron, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a2_over_3, dx, numpts3d, Rbar, fbar, fRbar, sim, level+1);
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 4 check_fr_wrong = FR_WRONG in gamma_cycle 4" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }
  }

  if(sim.multigrid_restrict_mode == RESTRICT_SCALARON)
  {
    engine.restrict(scalaron, temp, level);
    if(engine.isPartLayer(level+1))
    {
      subtract_fields(scalaron[level+1], temp[level+1], temp[level+1]);
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
    if(sim.multigrid_restrict_mode == RESTRICT_SCALARON)
    {
      add_fields(scalaron[level], trunc[level], scalaron[level]);
      check_fr_wrong = convert_scalaron_to_deltaR(eightpiG_deltaT[level], deltaR[level], scalaron[level], Rbar, fRbar, sim);
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
      check_fr_wrong = convert_deltaR_to_scalaron(deltaR[level], scalaron[level], Rbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 6 check_fr_wrong = FR_WRONG in gamma_cycle 6" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }

    build_laplacian_scalaron(scalaron[level], laplace_scalaron[level], dx[level], sim);

    // TODO In place of post-smoothing: requires that error be smaller than relaxation_error at each layer
    for(int j=0; j<sim.post_smoothing; ++j)
    {
      check_fr_wrong = update_scalaron_and_deltaR(scalaron[level], laplace_scalaron[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
        return FR_WRONG_RETURN;
      }
    }

    error = compute_error(scalaron[level], laplace_scalaron[level], deltaR[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, numpts3d[level], sim);

    // TODO REMOVE
    // for(int jj=0; jj<level; ++jj) COUT << "     ";
    // COUT << " level " << level << " - f. error = " << error << endl;
    // END REMOVE
  }

  trunc_error = euclidean_norm(trunc[level], engine, level, numpts3d[level]);

  return trunc_error;
}

#endif
