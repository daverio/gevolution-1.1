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


template <class FieldType>
void initial_guess_deltaR(
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  const metadata & sim
)
{
  copy_field(eightpiG_deltaT, deltaR, -1.);
  // zero_field(deltaR);
  return;
}


//////////////////////////
// Builds diff_u term
//////////////////////////
template <class FieldType>
void build_diff_u(
  Field<FieldType> & u,
  Field<FieldType> & diff_u,
  const double dx,
  const double fRbar
)
{
  Site x(u.lattice());
  double dx2 = dx*dx;
  double aplus, aminus;

  for(x.first(); x.test(); x.next())
  {
    // TODO: Check which definition works best

    // Option 1: laplace (exp(u))
    diff_u(x) = exp(u(x+0)) + exp(u(x-0)) + exp(u(x+1)) + exp(u(x-1)) + exp(u(x+2)) + exp(u(x-2)) - 6.*exp(u(x));
    // End Option 1

    // Option 2: div ( exp(u) * grad(u) )
    // diff_u(x) = 0.;
    // aplus = exp(u(x+0)) + exp(u(x));
    // aminus = exp(u(x)) + exp(u(x-0));
    // diff_u(x) += aplus * ( u(x+0) - u(x) ) - aminus * ( u(x) - u(x-0) );
    // aplus = exp(u(x+1)) + exp(u(x));
    // aminus = exp(u(x)) + exp(u(x-1));
    // diff_u(x) += aplus * ( u(x+1) - u(x) ) - aminus * ( u(x) - u(x-1) );
    // aplus = exp(u(x+2)) + exp(u(x));
    // aminus = exp(u(x)) + exp(u(x-2));
    // diff_u(x) += aplus * ( u(x+2) - u(x) ) - aminus * ( u(x) - u(x-2) );
    // diff_u(x) /= 2.;
    // End Option 2

    // Common for both options
    diff_u(x) *= fRbar / dx2;
  }
  return;
}

//////////////////////////
// in-place calculation of residual Y
// of equation Y[xi] == 0
// TODO: Details here
//////////////////////////
double residual_xi(
  const double xi,
  const double laplace_xi,
  const double deltaR,
  const double rhs,
  const double coeff1,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim
)
{
  double R = Rbar + deltaR,
         temp;

  temp = deltaR;
  // TODO Additional terms -- see if necessary
  temp = temp * (1 - fRbar - xi) + 2. * (f(R, sim, 884) - fbar) - Rbar * xi;
  // End additional terms

  return laplace_xi - coeff1 * temp - rhs;
}

//////////////////////////
// in-place calculation of residual Y
// of equation Y[u] == 0
// TODO: Details here
//////////////////////////
double residual_u(
  const double u,
  const double diff_u,
  const double deltaR,
  const double rhs,
  const double coeff1,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim
)
{
  double R = Rbar + deltaR,
         temp,
         temp2;

  temp = deltaR;
  // TODO Additional terms -- see if necessary
  temp2 = exp(u);
  temp = temp * (1 - fRbar * temp2) + 2. * (f(R, sim, 884) - fbar) - Rbar * fRbar * (temp2 - 1.);
  // End additional terms

  return diff_u - coeff1 * temp - rhs;
}


//////////////////////////
// In-place calculation of d(residual)/d(xi) = dY/d(xi)
// from equation Y[xi] == 0
// Needed for next guess xi^{i+1} = xi^{i} - Y^{i}/(dY/d(xi))^{i}
//////////////////////////
inline double dresidual_dxi(
  const double xi,
  const double laplace_xi,
  const double deltaR,
  const double coeff1,
  const double Rbar,
  const double fRbar,
  const metadata & sim
)
{
  double temp,
         R = Rbar + deltaR;

  temp = 1. / fRR(R, sim, 738);
  // TODO Additional terms -- see if necessary
  temp = temp * (1. + fRbar + xi) - R;
  // End additional terms

  return - coeff1 * temp;
}


//////////////////////////
// In-place calculation of d(residual)/du = dY/du
// from equation Y[u] == 0
// Needed for next guess u^{i+1} = u^{i} - Y^{i}/(dY/du)^{i}
//////////////////////////
inline double dresidual_du(
  const double u,
  const double diff_u,
  const double deltaR,
  const double coeff1,
  const double Rbar,
  const double fRbar,
  const metadata & sim
)
{
  double temp,
         fr,
         R = Rbar + deltaR;

  fr = fRbar * exp(u);
  temp = fr / fRR(R, sim, 738);
  // TODO Additional terms -- see if necessary
  temp = temp * (1. + fr) - R * fr;
  // End additional terms

  return - coeff1 * temp;
}


//////////////////////////
// Euclidean Norm
//////////////////////////
template <class FieldType>
double euclidean_norm(
  Field<FieldType> & field,
  MultiGrid & engine,
  const int level,
  const long numpts3d
)
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
//
//////////////////////////
template <class FieldType>
double compute_error(
  Field<FieldType> & u,
  Field<FieldType> & diff_u,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  MultiGrid & engine,
  const double coeff1,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim,
  const int level,
  const long numpts3d
)
{
  double temp,
         error = 0.;
  Site x(u.lattice());

  if(sim.relaxation_error_method == RELAXATION_ERROR_METHOD_SUM) // Euclidean Norm
  {
    if(engine.isPartLayer(level))
    {
      for(x.first(); x.test(); x.next())
      {
        temp = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
        error += temp*temp;
      }
    }
    parallel.layer(engine.player(level)).sum(error);
    error = sqrt(error) / numpts3d;
  }
  else if(sim.relaxation_error_method == RELAXATION_ERROR_METHOD_MAX) // TODO: Probably wrong, just an attempt
  {
    for(x.first(); x.test(); x.next())
    {
      temp = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
      temp = fabs(temp);
      if(temp > error)
      {
        error = temp;
      }
      parallel.layer(engine.player(level)).max(error);
    }
  }

  return error;
}


//////////////////////////
// Error on finest layer -- no need for multigrid engine
//////////////////////////
template <class FieldType>
double compute_error(
  Field<FieldType> & u,
  Field<FieldType> & diff_u,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  const double coeff1,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim,
  const long numpts3d
)
{
  double temp,
         error = 0.;
  Site x(u.lattice());

  if(sim.relaxation_error_method == RELAXATION_ERROR_METHOD_SUM) // Euclidean Norm
  {
    for(x.first(); x.test(); x.next())
    {
      temp = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
      error += temp*temp;
    }
    parallel.sum(error);
    error = sqrt(error)/numpts3d;
  }
  else if(sim.relaxation_error_method == RELAXATION_ERROR_METHOD_MAX) // TODO: Probably wrong, just an attempt
  {
    for(x.first(); x.test(); x.next())
    {
      temp = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
      temp = fabs(temp);
      if(temp > error)
      {
        error = temp;
      }
      parallel.max(error);
    }
  }

  return error;
}


//////////////////////////
//
//////////////////////////
template <class FieldType>
void build_residual_u(
  Field<FieldType> & u,
  Field<FieldType> & diff_u,
  Field<FieldType> & deltaR,
  Field<FieldType> & rhs,
  Field<FieldType> & destination_for_residual,
  MultiGrid & engine,
  const double coeff1,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim,
  const int level
)
{

  double temp,
         error = 0.;
  Site x(u.lattice());

  if(engine.isPartLayer(level))
  {
    for(x.first(); x.test(); x.next())
    {
      destination_for_residual(x) = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
    }
  }
  return;
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
  Field<FieldType> & diff_u,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  const double dx,
  const double coeff1,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim
)
{
  double Y, dY, temp;

  if(sim.multigrid_red_black)
  {
    SiteRedBlack3d x(u.lattice());

    for(x.firstRed(); x.testRed(); x.nextRed())
    {
      Y = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
      dY = dresidual_du(u(x), diff_u(x), deltaR(x), coeff1, Rbar, fRbar, sim);
      if(dY)
      {
        temp = sim.relaxation_overrel_coeff * Y/dY;
        if(u(x) > temp)
        {
          u(x) -= temp;
        }
      }
      else
      {
        cout << " Warning: dY is zero!" << endl;
        cout << u(x) << ", " << diff_u(x) << ", " << deltaR(x) << ", " << Y << ", " << dY;
        cin.get();
      }
    }

    u.updateHalo();
    build_diff_u(u, diff_u, dx, fRbar); // TODO: This might be optimised somehow -- only build_diff_u on the Blacks

    for(x.firstBlack(); x.testBlack(); x.nextBlack())
    {
      Y = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
      dY = dresidual_du(u(x), diff_u(x), deltaR(x), coeff1, Rbar, fRbar, sim);
      if(dY)
      {
        temp = sim.relaxation_overrel_coeff * Y/dY;
        if(u(x) > temp)
        {
          u(x) -= temp;
        }
      }
      else
      {
        cout << " Warning: dY is zero!" << endl;
        cout << u(x) << ", " << diff_u(x) << ", " << deltaR(x) << ", " << Y << ", " << dY;
        cin.get();
      }
    }
  }
  else
  {
    Site x(u.lattice());
    for(x.first(); x.test(); x.next())
    {
      Y = residual_u(u(x), diff_u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
      dY = dresidual_du(u(x), diff_u(x), deltaR(x), coeff1, Rbar, fRbar, sim);
      if(dY)
      {
        temp = sim.relaxation_overrel_coeff * Y / dY;
        u(x) -= temp;
      }
      else
      {
        cout << " Warning: dY is zero!" << endl;
        cout << u(x) << ", " << diff_u(x) << ", " << deltaR(x) << ", " << Y << ", " << dY;
        cin.get();
      }
    }
  }

  u.updateHalo();
  build_diff_u(u, diff_u, dx, fRbar);
  temp = convert_u_to_deltaR(eightpiG_deltaT, deltaR, u, Rbar, fRbar, sim);
  return temp;

}

//////////////////////////
// Relaxation solver for u
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation_u(
  Field<FieldType> & u,
  Field<FieldType> & diff_u,
  Field<FieldType> & deltaR,
  Field<FieldType> & eightpiG_deltaT,
  Field<FieldType> & rhs,
  Field<FieldType> & u_temp,
  const double a,
  const double dx,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim,
  double err
)
{
  int i = 0, error_increased = 0;
  double coeff1 = a*a/3.,
         initial_error = err,
         check_fr_wrong,
         error = 0.;
  long numpts3d = (long) sim.numpts * sim.numpts * sim.numpts;

  copy_field(u, u_temp);

  while(true)
  {
    // Attempt new guess, writing on u_temp
    check_fr_wrong = update_u_and_deltaR(u_temp, diff_u, deltaR, eightpiG_deltaT, rhs, dx, coeff1, Rbar, fbar, fRbar, sim);

    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 8 check_fr_wrong = FR_WRONG in relaxation_u" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    error = compute_error(u_temp, diff_u, deltaR, rhs, coeff1, Rbar, fbar, fRbar, sim, numpts3d);

    // if(error < sim.relaxation_error / numpts3d || i > sim.multigrid_pre_smoothing) break;
    // TODO: Check if this condition is ok
    if(++i > sim.multigrid_pre_smoothing) break;
    else if(error_increased > 100)
    {
      COUT << " Tried to trim 100 times, accepting current error = " << error << endl;
      break;
    }
  }

  return error;
}

//////////////////////////
// TODO: add description
//////////////////////////
template <class FieldType>
void restrict_to_level(MultiField<FieldType> * field, MultiGrid & engine, const int target_level)
{
  for(int i=0; i<target_level; ++i)
  {
    engine.restrict(field, i);
  }
  return;
}


//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double multigrid_u(
  MultiField<FieldType> * u,
  MultiField<FieldType> * diff_u,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  const double a,
  const double DX,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim
)
{
  int i,
      max_level = sim.multigrid_n_grids-1,
      count = 0;
  long numpts3d[sim.multigrid_n_grids];
  double dx[sim.multigrid_n_grids],
         coeff1 = a*a/3.,
         initial_error = 0.,
         error = 0.,
         trunc_error = 0.,
         check_fr_wrong;
  numpts3d[0] = (long) sim.numpts * sim.numpts * sim.numpts;
  dx[0] = DX;
  for(i=1; i<sim.multigrid_n_grids; ++i)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8;
    dx[i] = dx[i-1]*2.;
  }

  // Build eightpiG_deltaT on each level (it will not change)
  restrict_to_level(eightpiG_deltaT, engine, max_level);
  copy_field(eightpiG_deltaT[0], rhs[0], coeff1);
  // Compute initial error
  // TODO: remove after debugging
  initial_error = compute_error(u[0], diff_u[0], deltaR[0], rhs[0], engine, coeff1, Rbar, fbar, fRbar, sim, 0, numpts3d[0]);

  COUT << " z = " << 1./a - 1. << endl;

  while(true)
  {
    trunc_error = gamma_cycle(u, diff_u, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim, 0);

    if(trunc_error > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 1 trunc_error = FR_WRONG in multigrid_u" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }
    u[0].updateHalo();
    build_diff_u(u[0], diff_u[0], dx[0], fRbar);

    check_fr_wrong = convert_u_to_deltaR(eightpiG_deltaT[0], deltaR[0], u[0], Rbar, fRbar, sim);
    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 9 check_fr_wrong = FR_WRONG in multigrid_u" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    error = compute_error(u[0], diff_u[0], deltaR[0], rhs[0], engine, coeff1, Rbar, fbar, fRbar, sim, 0, numpts3d[0]);
    ++count;

    if(3. * error / trunc_error <= 1.)
    {
      COUT << " Error small enough after " << count << " iterations of multigrid_u" << endl;
      break;
    }
  }

  COUT << " (1/3) * trunc_error = " << trunc_error / 3. << endl;
  COUT << "               error = " << error << endl;

  return 3. * error / trunc_error;
}



//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double multigrid_FMG(
  MultiField<FieldType> * u,
  MultiField<FieldType> * diff_u,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  const double a,
  const double DX,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim
)
{
  int i,
      j,
      max_level = sim.multigrid_n_grids - 1; // Not exactly needed, but more redable this way
  long numpts3d[sim.multigrid_n_grids];
  double dx[sim.multigrid_n_grids],
         error,
         check_fr_wrong,
         coeff1 = a*a/3.;

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
    copy_field(eightpiG_deltaT[max_level], rhs[max_level], coeff1);
  }
  copy_field(eightpiG_deltaT[0], rhs[0], coeff1);
  // build initial guess for u on coarsest grid
  restrict_to_level(u, engine, max_level);

  // Compute initial error
  // TODO: remove after debugging
  u[0].updateHalo();
  build_diff_u(u[0], diff_u[0], dx[0], fRbar);
  error = compute_error(u[0], diff_u[0], deltaR[0], rhs[0], engine, coeff1, Rbar, fbar, fRbar, sim, 0, numpts3d[0]);
  COUT << " z = " << 1./a - 1. << ", initial error = " << error << endl;

  // solve problem on coarsest grid
  if(engine.isPartLayer(max_level))
  {
    u[max_level].updateHalo();
    build_diff_u(u[max_level], diff_u[max_level], dx[max_level], fRbar);

    check_fr_wrong = convert_u_to_deltaR(eightpiG_deltaT[max_level], deltaR[max_level], u[max_level], Rbar, fRbar, sim);
    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      cout << parallel.rank() << " 10 check_fr_wrong = FR_WRONG in multigrid_FMG" << endl;
      //END REMOVE
      return FR_WRONG_RETURN;
    }

    while(true)
    {
      check_fr_wrong = update_u_and_deltaR(u[max_level], diff_u[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[max_level], coeff1, Rbar, fbar, fRbar, sim);

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 11 check_fr_wrong = FR_WRONG in multigrid_FMG" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }

      error = compute_error(u[max_level], diff_u[max_level], deltaR[max_level], rhs[max_level], engine, coeff1, Rbar, fbar, fRbar, sim, max_level, numpts3d[max_level]);

      if(error < sim.relaxation_error) break;
    }
  }

  for(i=sim.multigrid_n_grids-2; i>=0; i--) // i is the starting level for the V or W mini-cycle
  {
    // prolong i+1 solution
    engine.prolong(u, i+1);
    // build/copy source to the correct field (rhs)
    if(engine.isPartLayer(i))
    {
      copy_field(eightpiG_deltaT[i], rhs[i], coeff1);
    }

    for(j=0; j<sim.multigrid_n_cycles; ++j)
    {
      check_fr_wrong = gamma_cycle(u, diff_u, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim, i);
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 12 check_fr_wrong = FR_WRONG in multigrid_FMG" << endl;
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

  u[0].updateHalo();
  build_diff_u(u[0], diff_u[0], dx[0], fRbar);
  check_fr_wrong = convert_u_to_deltaR(eightpiG_deltaT[0], deltaR[0], u[0], Rbar, fRbar, sim);
  if(check_fr_wrong > FR_WRONG)
  {
    //TODO REMOVE after debugging
    cout << parallel.rank() << " 13 check_fr_wrong = FR_WRONG in multigrid_FMG" << endl;
    //END REMOVE
    return FR_WRONG_RETURN;
  }

  copy_field(eightpiG_deltaT[0], rhs[0], coeff1);
  error = compute_error(u[0], diff_u[0], deltaR[0], rhs[0], engine, coeff1, Rbar, fbar, fRbar, sim, 0, numpts3d[0]);
  COUT << " z = " << 1./a - 1. << ", error = " << error << endl;

  return error; // TODO: Compute and divide by truncation error here! Fix!
}


//////////////////////////
// // Solves: diff_u - coeff1 * deltaR == rhs
//////////////////////////
template <class FieldType>
double gamma_cycle(
  MultiField<FieldType> * u,
  MultiField<FieldType> * diff_u,
  MultiField<FieldType> * deltaR,
  MultiField<FieldType> * eightpiG_deltaT,
  MultiField<FieldType> * rhs,
  MultiField<FieldType> * temp,
  MultiField<FieldType> * trunc,
  MultiGrid & engine,
  const double a,
  double * dx,
  long * numpts3d,
  const double Rbar,
  const double fbar,
  const double fRbar,
  const metadata & sim,
  const int level
)
{
  double coeff1 = a*a/3.,
         error = 0.,
         trunc_error = 0.,
         check_fr_wrong;
  int max_level = sim.multigrid_n_grids-1, count = 0;

  if(engine.isPartLayer(level))
  {
    u[level].updateHalo();
    build_diff_u(u[level], diff_u[level], dx[level], fRbar);

    // TODO remove after debug
    error = compute_error(u[level], diff_u[level], deltaR[level], rhs[level], engine, coeff1, Rbar, fbar, fRbar, sim, level, numpts3d[level]);
    // for(int jj=0; jj<level; ++j) COUT << "     ";
    // COUT << " Level " << level << ",  in. error = " << error << endl;
    // End remove

    for(int j=0; j<sim.multigrid_pre_smoothing; ++j)
    {
      //TODO REMOVE after debugging
      // check_field(u[level], "u[level]", numpts3d[level], "debug 1");
      // check_field(diff_u[level], "diff_u[level]", numpts3d[level], "debug 1");
      // check_field(deltaR[level], "deltaR[level]", numpts3d[level], "debug 1");
      // check_field(rhs[level], "rhs[level]", numpts3d[level], "debug 1");
      //END REMOVE

      check_fr_wrong = update_u_and_deltaR(u[level], diff_u[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], coeff1, Rbar, fbar, fRbar, sim);

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        check_field(u[level], "u[level]", numpts3d[level], "debug 2");
        check_field(diff_u[level], "diff_u[level]", numpts3d[level], "debug 2");
        check_field(deltaR[level], "deltaR[level]", numpts3d[level], "debug 2");
        check_field(rhs[level], "rhs[level]", numpts3d[level], "debug 2");
        cout << parallel.rank() << " 1 check_fr_wrong = FR_WRONG in gamma_cycle" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }

    build_residual_u(u[level], diff_u[level], deltaR[level], rhs[level], temp[level], engine, coeff1, Rbar, fbar, fRbar, sim, level);

    // TODO remove after debug
    error = compute_error(u[level], diff_u[level], deltaR[level], rhs[level], engine, coeff1, Rbar, fbar, fRbar, sim, level, numpts3d[level]);
    // COUT << level << "     After pre-smoothing -- Error = " << error << endl;
    // End remove
  }

  if(sim.multigrid_check_shape)
  {
    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << level << endl;
    for(int k=0; k<level; ++k) COUT << "  ";
    COUT << " \\" << endl;
  }

  // Restrict from level to level+1
  if(sim.multigrid_restrict_mode == RESTRICT_U)
  {
    engine.restrict(u, level);
  }
  else
  {
    engine.restrict(deltaR, level);
  }

  engine.restrict(temp, level);
  engine.restrict(rhs, level);

  if(engine.isPartLayer(level+1))
  {
    if(sim.multigrid_restrict_mode == RESTRICT_U)
    {
      check_fr_wrong = convert_u_to_deltaR(eightpiG_deltaT[level+1], deltaR[level+1], u[level+1], Rbar, fRbar, sim);
    }
    else
    {
      check_fr_wrong = convert_deltaR_to_u(deltaR[level+1], u[level+1], Rbar, fRbar, sim);
    }

    if(check_fr_wrong > FR_WRONG)
    {
      //TODO REMOVE after debugging
      //END REMOVE
      cout << parallel.rank() << " 2 check_fr_wrong = FR_WRONG in gamma_cycle 2" << endl;
      return FR_WRONG_RETURN;
    }

    u[level+1].updateHalo();
    build_diff_u(u[level+1], diff_u[level+1], dx[level+1], fRbar);
    build_residual_u(u[level+1], diff_u[level+1], deltaR[level+1], rhs[level+1], trunc[level+1], engine, coeff1, Rbar, fbar, fRbar, sim, level+1);
    subtract_fields(trunc[level+1], temp[level+1], trunc[level+1]);
    add_fields(rhs[level+1], trunc[level+1], rhs[level+1]);
  }

  if(level+1 == max_level) // Coarsest grid
  {
    count = 0;

    if(engine.isPartLayer(max_level))
    {
      error = compute_error(u[max_level], diff_u[max_level], deltaR[max_level], rhs[max_level], engine, coeff1, Rbar, fbar, fRbar, sim, max_level, numpts3d[max_level]);
      while(true)
      {
        trunc_error = error; // Saving the original error

        check_fr_wrong = update_u_and_deltaR(u[max_level], diff_u[max_level], deltaR[max_level], eightpiG_deltaT[max_level], rhs[max_level], dx[max_level], coeff1, Rbar, fbar, fRbar, sim);

        if(check_fr_wrong > FR_WRONG)
        {
          //TODO REMOVE after debugging
          cout << parallel.rank() << " 3 check_fr_wrong = FR_WRONG in gamma_cycle 3" << endl;
          //END REMOVE
          return FR_WRONG_RETURN;
        }

        // //TODO REMOVE
        // if(count == 3)
        // {
        //   cout << " NEW GUESS " << count << endl;
        //   initial_guess_deltaR(deltaR[max_level], eightpiG_deltaT[max_level], sim);
        //   copy_field(eightpiG_deltaT[max_level], rhs[max_level], a*a/3.);
        //
        //   convert_deltaR_to_u(deltaR[max_level], u[max_level], Rbar, fRbar, sim); // deltaR and u are now consistent
        //   u[max_level].updateHalo();
        //   build_diff_u(u[max_level], diff_u[max_level], dx[max_level], fRbar);
        // }
        // // END REMOVE

        error = compute_error(u[max_level], diff_u[max_level], deltaR[max_level], rhs[max_level], engine, coeff1, Rbar, fbar, fRbar, sim, max_level, numpts3d[max_level]);

        // COUT << " Level " << level << ", error = " << error << endl;

        // TODO: Check if this convergence condition makes sense
        if(error < sim.relaxation_error/numpts3d[max_level])
        {
          break;
        }
      }
    }
  }
  else
  {
    for(int j=0; j<sim.multigrid_shape; ++j)
    {
      check_fr_wrong = gamma_cycle(u, diff_u, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim, level+1);
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 4 check_fr_wrong = FR_WRONG in gamma_cycle 4" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }
  }

  if(sim.multigrid_restrict_mode == RESTRICT_U)
  {
    engine.restrict(u, temp, level);

    if(engine.isPartLayer(level+1))
    {
      subtract_fields(u[level+1], temp[level+1], temp[level+1]);
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
    if(sim.multigrid_restrict_mode == RESTRICT_U)
    {
      add_fields(u[level], trunc[level], u[level]);
      check_fr_wrong = convert_u_to_deltaR(eightpiG_deltaT[level], deltaR[level], u[level], Rbar, fRbar, sim);
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
      check_fr_wrong = convert_deltaR_to_u(deltaR[level], u[level], Rbar, fRbar, sim);
      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 6 check_fr_wrong = FR_WRONG in gamma_cycle 6" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }

    u[level].updateHalo();
    build_diff_u(u[level], diff_u[level], dx[level], fRbar);

    for(int j=0; j<sim.multigrid_post_smoothing; ++j)
    {
      check_fr_wrong = update_u_and_deltaR(u[level], diff_u[level], deltaR[level], eightpiG_deltaT[level], rhs[level], dx[level], coeff1, Rbar, fbar, fRbar, sim);

      if(check_fr_wrong > FR_WRONG)
      {
        //TODO REMOVE after debugging
        cout << parallel.rank() << " 7 check_fr_wrong = FR_WRONG in gamma_cycle 7" << endl;
        //END REMOVE
        return FR_WRONG_RETURN;
      }
    }

    // TODO remove after debug
    // error = compute_error(u[level], diff_u[level], deltaR[level], rhs[level], engine, coeff1, Rbar, fbar, fRbar, sim, level, numpts3d[level]);
    // for(int jj=0; jj<level; ++jj) cout << "     ";
    // COUT << " Level " << level << ", fin. error = " << error << endl;
    // End remove
  }

  trunc_error = euclidean_norm(trunc[level], engine, level, numpts3d[level]);

  return trunc_error;
}



#endif
