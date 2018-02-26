//////////////////////////
// relaxation.hpp
//////////////////////////
//
// Implementation of relaxation method solver for F(R) gravity
//
// Authors: David Daverio (Cambridge University)
//          Lorenzo Reverberi (CEICO - Czech Academy of Sciences, Prague
//                             University of Cape Town)
//
//////////////////////////

#ifndef RELAXATION_HEADER
#define RELAXATION_HEADER

//////////////////////////
// Builds diff_u term
//////////////////////////
template <class FieldType>
void build_diff_u(Field<FieldType> & u,
                  Field<FieldType> & diff_u,
                  const double dx,
                  const double fRbar)
{
  Site x(u.lattice());
  double dx2 = dx*dx;

  for(x.first(); x.test(); x.next())
  {
    diff_u(x) = exp(u(x+0)) + exp(u(x-0)) + exp(u(x+1)) + exp(u(x-1)) + exp(u(x+2)) + exp(u(x-2)) - 6.*exp(u(x));
    diff_u(x) *= fRbar / dx2;
  }
  return;
}

//////////////////////////
// Computes the error (locally)
//////////////////////////
inline double error_u(const double err, const double val)
{
  if(val)
  {
    return fabs(err/val);
  }
  else return 1.;
}


//////////////////////////
// in-place calculation of residual Y
// of equation Y[u] == 0
// TODO: Details here
//////////////////////////
inline double residual(const double diff_u,
                       const double u,
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
  // Additional terms
  temp = temp * (1 - fR(R, sim, 885)) + 2. * (f(R, sim, 884) - fbar) - Rbar * fRbar * (exp(u) - 1.);
  return diff_u - coeff1 * temp - rhs;
}

//////////////////////////
// in-place calculation of d(residual)/du = dY/du
// from equation Y[u] == 0
// Needed for next guess u^{i+1} = u^{i} - Y^{i}/(dY/du)^{i}
//////////////////////////
inline double dresidual_du(const double diff_u,
                           const double u,
                           const double deltaR,
                           const double coeff1,
                           const double Rbar,
                           const double fRbar,
                           const metadata & sim
                           )
{
  double temp,
         R = Rbar + deltaR;
  temp = fR(R, sim, 737);
  // Additional terms
  temp *= (1. + temp)/fRR(R, sim, 738) - deltaR;
  temp -= Rbar * fRbar * exp(u);

  return diff_u - coeff1 * temp;
}



//////////////////////////
//
//////////////////////////
template <class FieldType>
double compute_residual_u(Field<FieldType> & diff_u,
                          Field<FieldType> & u,
                          Field<FieldType> & deltaR,
                          Field<FieldType> & rhs,
                          MultiGrid & engine,
                          const double coeff1,
                          const double Rbar,
                          const double fbar,
                          const double fRbar,
                          const metadata & sim,
                          const int level)
{
  double temp,
         error = 0.;
  Site x(diff_u.lattice());

  for(x.first(); x.test(); x.next())
  {
    temp = residual(diff_u(x), u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
    temp = error_u(temp, coeff1 * Rbar);
    if(temp > error)
    {
      error = temp;
    }
  }
  parallel.layer(engine.player(level)).max(error);
  return error;
}


//////////////////////////
// Computes new guess for u
// Solves the equation Y[u] == 0
// where Y[u] := laplacian( fRbar * exp(u) ) - a^2/3 * zeta(u)
// TODO: Details here
//////////////////////////

template <class FieldType>
void update_u(Field<FieldType> & u,
              Field<FieldType> & diff_u,
              Field<FieldType> & deltaR,
              Field<FieldType> & rhs,
              const double coeff1,
              const double Rbar,
              const double fbar,
              const double fRbar,
              const metadata & sim)
{
  Site x(u.lattice());
  double Y, dY, temp, temp2;
  for(x.first(); x.test(); x.next())
  {
    Y = residual(diff_u(x), u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
    dY = dresidual_du(diff_u(x), u(x), deltaR(x), coeff1, Rbar, fRbar, sim);
    temp = dY ? Y/dY : 0.;
    temp = Y/dY;
    u(x) -= temp;
  }
  return;
}


template <class FieldType>
void update_u_red_black(Field<FieldType> & u,
                        Field<FieldType> & deltaR,
                        Field<FieldType> & rhs,
                        const double coeff1,
                        const double Rbar,
                        const double fbar,
                        const double fRbar,
                        const metadata & sim)
{
  SiteRedBlack3d x(u.lattice());
  double Y,
         dY,
         temp,
         diff_u;
  for(x.first(); x.test(); x.next())
  {
    diff_u = exp(u(x+0)) + exp(u(x-0)) + exp(u(x+1)) + exp(u(x-1)) + exp(u(x+2)) + exp(u(x-2)) - 6. * exp(u(x));
    diff_u *= fRbar * (double) sim.numpts * (double) sim.numpts;
    Y = residual(diff_u, u(x), deltaR(x), rhs(x), coeff1, Rbar, fbar, fRbar, sim);
    dY = dresidual_du(diff_u, u(x), deltaR(x), coeff1, Rbar, fRbar, sim);
    temp = dY ? Y/dY : 0.;
    temp = Y/dY;
    u(x) -= temp;
  }
  return;
}





//////////////////////////
// Relaxation solver for u
// TODO: Details here
//////////////////////////
template <class FieldType>
double relaxation_u(Field<FieldType> & u,
                    Field<FieldType> & diff_u,
                    Field<FieldType> & deltaR,
                    Field<FieldType> & eightpiG_deltaT,
                    Field<FieldType> & zeta,
                    Field<FieldType> & residual,
                    const double a,
                    const double dx,
                    const double Rbar,
                    const double fRbar,
                    const metadata & sim)
{
  // TODO: Rewrite!!!
}

//////////////////////////
// TODO: add description
//////////////////////////
template <class FieldType>
void restrict_to_level(MultiField<FieldType> * field,
                       MultiGrid & engine,
                       const int target_level)
{
  for(int i=0; i<target_level; i++)
  {
    engine.restrict(field, i);
  }
  return;
}


//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double multigrid_u(MultiField<FieldType> * u,
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
                   const metadata & sim)
{
  int i,
      mg_gamma = sim.multigrid_shape,
      n_grids = engine.nl(); // Min between chosen # of grids and maximum allowed by engine
  long numpts3d[n_grids];
  double dx[n_grids],
         coeff1 = a*a/3.;

  numpts3d[0] = (long) sim.numpts * sim.numpts * sim.numpts;
  dx[0] = DX;
  for(i=1; i<n_grids; i++)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8;
    dx[i] = dx[i-1]*2.;
  }

  // Build eightpiG_deltaT on each level (it will not change)
  restrict_to_level(eightpiG_deltaT, engine, n_grids-1);
  copy_field(eightpiG_deltaT[0], rhs[0], coeff1);

  gamma_cycle(u, diff_u, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim, mg_gamma, n_grids, 0);

  u[0].updateHalo();
  build_diff_u(u[0], diff_u[0], dx[0], fRbar);
  convert_u_to_deltaR(eightpiG_deltaT[0], deltaR[0], u[0], Rbar, fRbar, sim);
  copy_field(eightpiG_deltaT[0], rhs[0], coeff1);
  //
  double error = compute_residual_u(diff_u[0], u[0], deltaR[0], rhs[0], engine, coeff1, Rbar, fbar, fRbar, sim, 0);
  COUT << " z = " << 1./a - 1. << ", error = " << error << endl;

  return error;
}



//////////////////////////
// TODO: Details here
//////////////////////////
template <class FieldType>
double multigrid_FMG(MultiField<FieldType> * u,
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
                     const metadata & sim)
{
  int i,
      j,
      mg_gamma = sim.multigrid_shape,
      n_grids = engine.nl(), // Min between chosen # of grids and maximum allowed by engine
      max_level = n_grids - 1;
  long numpts3d[n_grids];
  double dx[n_grids],
         error,
         coeff1 = a*a/3.;

  numpts3d[0] = (long) sim.numpts * sim.numpts * sim.numpts;
  dx[0] = DX;
  for(i=1; i<n_grids; i++)
  {
    numpts3d[i] = (long) numpts3d[i-1] / 8;
    dx[i] = dx[i-1]*2.;
  }

  // // Build eightpiG_deltaT on each level (it will not change)
  restrict_to_level(eightpiG_deltaT, engine, max_level);
  copy_field(eightpiG_deltaT[max_level], rhs[max_level], coeff1);
  // build initial guess for u on coarsest grid
  restrict_to_level(u, engine, max_level);
  // solve problem on coarsest grid
  if(engine.isPartLayer(max_level))
  {
    u[max_level].updateHalo();
    build_diff_u(u[max_level], diff_u[max_level], dx[max_level], fRbar);
    convert_u_to_deltaR(eightpiG_deltaT[max_level], deltaR[max_level], u[max_level], Rbar, fRbar, sim);
    while(true)
    {
      if(sim.multigrid_red_black)
      {
        update_u_red_black(u[max_level], deltaR[max_level], rhs[max_level], coeff1, Rbar, fbar, fRbar, sim);
      }
      else
      {
        update_u(u[max_level], diff_u[max_level], deltaR[max_level], rhs[max_level], coeff1, Rbar, fbar, fRbar, sim);
      }
      u[max_level].updateHalo();
      build_diff_u(u[max_level], diff_u[max_level], dx[max_level], fRbar);
      convert_u_to_deltaR(eightpiG_deltaT[max_level], deltaR[max_level], u[max_level], Rbar, fRbar, sim);
      error = compute_residual_u(diff_u[max_level], u[max_level], deltaR[max_level], rhs[max_level], engine, coeff1, Rbar, fbar, fRbar, sim, max_level);

      if(error < sim.relaxation_error)
      {
        break;
      }
    }
  }

  for(i=n_grids-2; i>=0; i--) // i is the starting level for the V or W mini-cycle
  {
    // prolong i+1 solution
    engine.prolong(u, i+1);
    // build/copy source to the correct field (rhs)
    if(engine.isPartLayer(i))
    {
      copy_field(eightpiG_deltaT[i], rhs[i], coeff1);
    }

    for(j=0; j<sim.multigrid_n_cycles; j++)
    {
      gamma_cycle(u, diff_u, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim, mg_gamma, n_grids, i);
      if(sim.multigrid_check_shape && i && j == sim.multigrid_n_cycles-1)
      {
        for(int k=0; k<i; k++) COUT << "  ";
        COUT << i << endl;
        // if()
        for(int k=0; k<i-1; k++) COUT << "  ";
        COUT << " /" << endl;
      }
    }
  }

  u[0].updateHalo();
  build_diff_u(u[0], diff_u[0], dx[0], fRbar);
  convert_u_to_deltaR(eightpiG_deltaT[0], deltaR[0], u[0], Rbar, fRbar, sim);
  copy_field(eightpiG_deltaT[0], rhs[0], coeff1);
  error = compute_residual_u(diff_u[0], u[0], deltaR[0], rhs[0], engine, coeff1, Rbar, fbar, fRbar, sim, 0);
  COUT << " z = " << 1./a - 1. << ", error = " << error << endl;

  return error;
}


//////////////////////////
// // Solves: diff_u - coeff1 * deltaR == rhs
//////////////////////////
template <class FieldType>
double gamma_cycle(MultiField<FieldType> * u,
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
                   const int mg_gamma,
                   const int n_grids,
                   const int level)
{
  double coeff1 = a*a/3.,
         error;
  int max_level = n_grids-1;

  if(engine.isPartLayer(level))
  {
    u[level].updateHalo();
    build_diff_u(u[level], diff_u[level], dx[level], fRbar);
    for(int j=0; j<sim.multigrid_pre_smoothing; j++)
    {
      if(sim.multigrid_red_black)
      {
        update_u_red_black(u[level], deltaR[level], rhs[level], coeff1, Rbar, fbar, fRbar, sim);
      }
      else
      {
        update_u(u[level], diff_u[level], deltaR[level], rhs[level], coeff1, Rbar, fbar, fRbar, sim);
      }
      convert_u_to_deltaR(eightpiG_deltaT[level], deltaR[level], u[level], Rbar, fRbar, sim);
      u[level].updateHalo();
      build_diff_u(u[level], diff_u[level], dx[level], fRbar);
    }
    add_fields(diff_u[level], 1., deltaR[level], -coeff1, temp[level]);
  }

  // Restrict from level to level+1
  if(sim.multigrid_check_shape)
  {
    for(int k=0; k<level; k++) COUT << "  ";
    COUT << level << endl;
    for(int k=0; k<level; k++) COUT << "  ";
    COUT << " \\" << endl;
  }

  engine.restrict(temp, level);
  engine.restrict(u, level);
  engine.restrict(rhs, level);

  if(engine.isPartLayer(level+1))
  {
    u[level+1].updateHalo();
    build_diff_u(u[level+1], diff_u[level+1], dx[level+1], fRbar);
    convert_u_to_deltaR(eightpiG_deltaT[level+1], deltaR[level+1], u[level+1], Rbar, fRbar, sim);
    add_fields(diff_u[level+1], 1., deltaR[level+1], -coeff1, trunc[level+1]);
    subtract_fields(trunc[level+1], temp[level+1], trunc[level+1]);
    add_fields(rhs[level+1], trunc[level+1], rhs[level+1]);
  }

  if(level+1 == n_grids-1) // Coarsest grid
  {
    if(engine.isPartLayer(max_level))
    {
      u[max_level].updateHalo();
      build_diff_u(u[max_level], diff_u[max_level], dx[max_level], fRbar);
      convert_u_to_deltaR(eightpiG_deltaT[max_level], deltaR[max_level], u[max_level], Rbar, fRbar, sim);
      while(true)
      {
        if(sim.multigrid_red_black)
        {
          update_u_red_black(u[max_level], deltaR[max_level], rhs[max_level], coeff1, Rbar, fbar, fRbar, sim);
        }
        else
        {
          update_u(u[max_level], diff_u[max_level], deltaR[max_level], rhs[max_level], coeff1, Rbar, fbar, fRbar, sim);
        }
        u[max_level].updateHalo();
        build_diff_u(u[max_level], diff_u[max_level], dx[max_level], fRbar);
        convert_u_to_deltaR(eightpiG_deltaT[max_level], deltaR[max_level], u[max_level], Rbar, fRbar, sim);
        error = compute_residual_u(diff_u[max_level], u[max_level], deltaR[max_level], rhs[max_level], engine, coeff1, Rbar, fbar, fRbar, sim, max_level);

        if(error < sim.relaxation_error)
        {
          break;
        }
      }
    }
  }
  else
  {
    for(int j=0; j<mg_gamma; j++)
    {
      gamma_cycle(u, diff_u, deltaR, eightpiG_deltaT, rhs, temp, trunc, engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim, mg_gamma, n_grids, level+1);
    }
  }

  engine.restrict(u, temp, level);
  if(engine.isPartLayer(level+1))
  {
    subtract_fields(u[level+1], temp[level+1], temp[level+1]);
  }

  if(sim.multigrid_check_shape)
  {
    for(int k=0; k<level+1; k++) COUT << "  ";
    COUT << level+1 << endl;
    for(int k=0; k<level; k++) COUT << "  ";
    COUT << " /" << endl;
    if(!level) COUT << "0" << endl;
  }

  engine.prolong(temp, trunc, level+1);
  if(engine.isPartLayer(level))
  {
    add_fields(u[level], trunc[level], u[level]);
    u[level].updateHalo();
    build_diff_u(u[level], diff_u[level], dx[level], fRbar);
    convert_u_to_deltaR(eightpiG_deltaT[level], deltaR[level], u[level], Rbar, fRbar, sim);
    for(int j=0; j<sim.multigrid_post_smoothing; j++)
    {
      if(sim.multigrid_red_black)
      {
        update_u_red_black(u[level], deltaR[level], rhs[level], coeff1, Rbar, fbar, fRbar, sim);
      }
      else
      {
        update_u(u[level], diff_u[level], deltaR[level], rhs[level], coeff1, Rbar, fbar, fRbar, sim);
      }
      u[level].updateHalo();
      build_diff_u(u[level], diff_u[level], dx[level], fRbar);
      convert_u_to_deltaR(eightpiG_deltaT[level], deltaR[level], u[level], Rbar, fRbar, sim);
    }
  }
}




#endif
