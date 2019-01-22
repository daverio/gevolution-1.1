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
  Field<FieldType> & xi,
  Field<FieldType> & laplace_xi,
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

  convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim); // deltaR and u are now consistent
  build_laplacian(xi, laplace_xi, dx);

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

  if(!sim.newtonian_fR)
  {
    // TODO Additional terms -- see if necessary
    temp = temp * (1. - fRbar - xi) + 2. * (f(R, sim, 884) - fbar) - Rbar * xi;
    // End additional terms
  }

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
         temp = fRR(R, sim, 738);

  if(!temp || temp > FR_WRONG) return FR_WRONG_RETURN;

  temp = 1. / temp;

  if(!sim.newtonian_fR)
  {
    // TODO Additional terms -- see if necessary
    temp = temp * (1. + fRbar + xi) - R;
    // End additional terms
  }

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


double update_xi_and_deltaR_single(
  double & xi,
  double laplace_xi,
  double & deltaR,
  double rhs,
  double a2_over_3,
  double coeff,
  double overrelax,
  double Rbar,
  double fbar,
  double fRbar,
  const metadata & sim)
{
  double temp, dxi;

  dxi = xi;
  update_xi_single(xi, laplace_xi, deltaR, rhs, a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
  dxi = (xi - dxi) / 10.;

  while(true)
  {
    temp = convert_xi_to_deltaR_single(deltaR, xi, Rbar, fRbar, sim);
    if(temp < FR_WRONG) break;
    xi -= dxi;
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
      update_xi_and_deltaR_single(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }

    build_laplacian(xi, laplace_xi, dx); // TODO: This might be optimised somehow -- only build_laplacian on the Blacks?

    for(x.firstBlack(); x.testBlack(); x.nextBlack())
    {
      update_xi_and_deltaR_single(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    Site x(xi.lattice());
    for(x.first(); x.test(); x.next())
    {
      update_xi_and_deltaR_single(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, coeff, overrelax, Rbar, fbar, fRbar, sim);
    }
  }

  build_laplacian(xi, laplace_xi, dx);

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
  Field<FieldType> & xi,
  Field<FieldType> & laplace_xi,
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

  Site x(xi.lattice());

  for(x.first(); x.test(); x.next())
  {
    temp = residual_xi(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, Rbar, fbar, fRbar, sim);
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
  Field<FieldType> & xi,
  Field<FieldType> & laplace_xi,
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

  Site x(xi.lattice());

  for(x.first(); x.test(); x.next())
  {
    destination_for_residual(x) = residual_xi(xi(x), laplace_xi(x), deltaR(x), rhs(x), a2_over_3, Rbar, fbar, fRbar, sim);
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
int relaxation(
  MultiField<FieldType> * xi,
  MultiField<FieldType> * laplace_xi,
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

  COUT << " z = " << 1./a - 1. << " -- " << flush;

  error = compute_error(xi[0], laplace_xi[0], deltaR[0], rhs[0], dx, a*a/3., Rbar, fbar, fRbar, numpts3d, sim);

  if(error < sim.relaxation_error)
  {
    COUT << "initial and final error = " << error << endl;
    return 0;
  }
  else
  {
    COUT << "initial error = " << error << " -- " << flush;
  }

  if(sim.relaxation_method == METHOD_RELAX)
  {
    error = single_layer_solver(xi[0], laplace_xi[0], deltaR[0], eightpiG_deltaT[0], rhs[0], a*a/3., dx, numpts3d, Rbar, fbar, fRbar, sim);
  }
  else if(sim.relaxation_method == METHOD_MULTIGRID)
  {
    error = multigrid_solver(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, err, engine, a*a/3., dx, Rbar, fbar, fRbar, sim);
  }
  else if(sim.relaxation_method == METHOD_FMG)
  {
    error = FMG_solver(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, residual, err, engine, a, dx, Rbar, fbar, fRbar, sim);
  }

  // TODO REMOVE after debug
  if(error > FR_WRONG)
  {
    COUT << " Method " << sim.relaxation_method << " returned FR_WRONG_RETURN. Check what's going on. Press [Enter] to continue..." << endl;
    cin.get();
  }
  // END REMOVE

  COUT << "final error = " << error << endl;
  return 1;
}


//////////////////////////
// Single relaxation step
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation_step(
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
  long numpts3d,
  const metadata & sim)
{
  double error;
  // Attempt new guess
  error = update_xi_and_deltaR(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, Rbar, fbar, fRbar, sim);

  if(error > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << parallel.rank() << " 8 check_fr_wrong = FR_WRONG in relaxation_step" << endl;
    //END REMOVE
    return FR_WRONG_RETURN;
  }

  error = compute_error(xi, laplace_xi, deltaR, rhs, dx, a2_over_3, Rbar, fbar, fRbar, numpts3d, sim);

  return error;
}



//////////////////////////
// Relaxation solver
// TODO: Details here
//////////////////////////
template <class FieldType>
double single_layer_solver(
  Field<FieldType> & xi,
  Field<FieldType> & laplace_xi,
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
    error = relaxation_step(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, dx, a2_over_3, Rbar, fbar, fRbar, numpts3d, sim);

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
  MultiField<FieldType> * xi,
  MultiField<FieldType> * laplace_xi,
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
    trunc_error = gamma_cycle(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a2_over_3, dx, numpts3d, Rbar, fbar, fRbar, sim, 0);

    if(trunc_error > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 1 trunc_error = FR_WRONG in multigrid_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    convert_xi_to_deltaR(deltaR[0], xi[0], Rbar, fRbar, sim);

    error = compute_error(xi[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, Rbar, fbar, fRbar, numpts3d[0], sim);

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
  MultiField<FieldType> * xi,
  MultiField<FieldType> * laplace_xi,
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
  restrict_to_level(xi, engine, max_level);

  // Compute initial error
  // TODO: remove after debugging
  build_laplacian(xi[0], laplace_xi[0], dx[0]);
  error = compute_error(xi[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, Rbar, fbar, fRbar, numpts3d[0], sim);

  // solve problem on coarsest grid
  if(engine.isPartLayer(max_level))
  {
    build_laplacian(xi[max_level], laplace_xi[max_level], dx[max_level]);

    check_fr_wrong = convert_xi_to_deltaR(deltaR[max_level], xi[max_level], Rbar, fRbar, sim);
    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 10 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    while(true)
    {
      check_fr_wrong = update_xi_and_deltaR(xi[max_level], laplace_xi[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[max_level], a2_over_3, Rbar, fbar, fRbar, sim);

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 11 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }

      error = compute_error(xi[max_level], laplace_xi[max_level], deltaR[max_level], rhs[max_level], dx[max_level], a2_over_3, Rbar, fbar, fRbar, numpts3d[max_level], sim);

      if(error < sim.relaxation_error)
      {
        break;
      }
    }
  }

  for(i=sim.multigrid_n_grids-2; i>=0; i--) // i is the starting level for the V or W mini-cycle
  {
    // prolong i+1 solution
    engine.prolong(xi, i+1);
    // build/copy source to the correct field (rhs)
    if(engine.isPartLayer(i))
    {
      copy_field(eightpiG_deltaT[i], rhs[i], a2_over_3);
    }

    for(j=0; j<sim.multigrid_n_cycles; ++j)
    {
      check_fr_wrong = gamma_cycle(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a2_over_3, dx, numpts3d, Rbar, fbar, fRbar, sim, i);
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

  build_laplacian(xi[0], laplace_xi[0], dx[0]);
  check_fr_wrong = convert_xi_to_deltaR(deltaR[0], xi[0], Rbar, fRbar, sim);
  if(check_fr_wrong > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << parallel.rank() << " 13 check_fr_wrong = FR_WRONG in FMG_solver" << endl;
    //END REMOVE
    return FR_WRONG_RETURN;
  }

  copy_field(eightpiG_deltaT[0], rhs[0], a2_over_3);
  error = compute_error(xi[0], laplace_xi[0], deltaR[0], rhs[0], dx[0], a2_over_3, Rbar, fbar, fRbar, numpts3d[0], sim);

  return error; // TODO: Compute and divide by truncation error here! Fix!
}


//////////////////////////
// TODO Details here
//////////////////////////
template <class FieldType>
double gamma_cycle(
  MultiField<FieldType> * xi,
  MultiField<FieldType> * laplace_xi,
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
      build_laplacian(xi[level], laplace_xi[level], dx[level]);
    }

    // TODO REMOVE
    // error = compute_error(xi[level], laplace_xi[level], deltaR[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, numpts3d[level], sim);
    // for(int jj=0; jj<level; ++jj) COUT << "     ";
    // COUT << " level " << level << " - i. error = " << error << endl;
    // END REMOVE

    for(int j=0; j<sim.pre_smoothing; ++j)
    {
      check_fr_wrong = update_xi_and_deltaR(xi[level], laplace_xi[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
        return FR_WRONG_RETURN;
      }
    }

    build_residual(xi[level], laplace_xi[level], deltaR[level], rhs[level], temp[level], engine, a2_over_3, Rbar, fbar, fRbar, sim);
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
    build_residual(xi[level+1], laplace_xi[level+1], deltaR[level+1], rhs[level+1], trunc[level+1], engine, a2_over_3, Rbar, fbar, fRbar, sim);
    subtract_fields(trunc[level+1], temp[level+1], trunc[level+1]);
    add_fields(rhs[level+1], trunc[level+1], rhs[level+1]);
  }

  // Coarsest grid
  if(level+1 == max_level)
  {
    if(engine.isPartLayer(max_level))
    {
      error = single_layer_solver(xi[max_level], laplace_xi[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], a2_over_3, dx[level], numpts3d[max_level], Rbar, fbar, fRbar, sim);
    }
  }
  else
  {
    for(int j=0; j<sim.multigrid_shape; ++j)
    {
      check_fr_wrong = gamma_cycle(xi, laplace_xi, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a2_over_3, dx, numpts3d, Rbar, fbar, fRbar, sim, level+1);
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
      check_fr_wrong = update_xi_and_deltaR(xi[level], laplace_xi[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
        return FR_WRONG_RETURN;
      }
    }

    error = compute_error(xi[level], laplace_xi[level], deltaR[level], rhs[level], dx[level], a2_over_3, Rbar, fbar, fRbar, numpts3d[level], sim);

    // TODO REMOVE
    // for(int jj=0; jj<level; ++jj) COUT << "     ";
    // COUT << " level " << level << " - f. error = " << error << endl;
    // END REMOVE
  }

  trunc_error = euclidean_norm(trunc[level], engine, level, numpts3d[level]);

  return trunc_error;
}

#endif
