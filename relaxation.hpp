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


// Builds diff_u term
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
    diff_u(x) = (u(x+0) - u(x-0))*(u(x+0) - u(x-0)) +
    (u(x+1) - u(x-1))*(u(x+1) - u(x-1)) +
    (u(x+2) - u(x-2))*(u(x+2) - u(x-2));
    diff_u(x) /= 4.;
    diff_u(x) += u(x+0) + u(x-0) + u(x+1) + u(x-1) + u(x+2) + u(x-2) - 6.*u(x);
    diff_u(x) /= dx2;
    diff_u(x) *= fRbar * exp(u(x));
  }
  return;
}


// Compute residual_u
template <class FieldType>
double compute_residual_u(Field<FieldType> & u,
                          Field<FieldType> & diff_u,
                          Field<FieldType> & deltaR,
                          Field<FieldType> & eightpiG_deltaT,
                          const double coeff1)
{
  double Y, dY, temp, error = 0.;
  Site x(u.lattice());

  for(x.first(); x.test(); x.next())
  {
    Y = diff_u(x);
    dY = coeff1 * (deltaR(x) + eightpiG_deltaT(x));
    temp = Y - dY;
    temp /= max(fabs(Y), fabs(dY));
    error += temp*temp;
  }
  parallel.sum(error);
  return error;
}


// Computes new guess for u
template <class FieldType>
void update_u(Field<FieldType> & u,
              Field<FieldType> & diff_u,
              Field<FieldType> & deltaR,
              Field<FieldType> & eightpiG_deltaT,
              const double coeff1,
              const double Rbar,
              const metadata & sim,
              int is_it_multifield = 0)
{
  Site x(u.lattice());
  double Y, dY;
  for(x.first(); x.test(); x.next())
  {
    Y = diff_u(x) - coeff1 * (deltaR(x) + eightpiG_deltaT(x));
    dY = diff_u(x) - coeff1 * fR(Rbar + deltaR(x), sim, 885) / fRR(Rbar + deltaR(x), sim, 886);
    u(x) -= dY ? Y/dY : 0.;
  }
  if(!is_it_multifield)
  {
    u.updateHalo();
  }
  return;
}


// Same as relax_xi, but using the redefined field u explicitly
template <class FieldType>
double relax_u(Field<FieldType> & u,
               Field<FieldType> & diff_u,
               Field<FieldType> & deltaR,
               Field<FieldType> & eightpiG_deltaT,
               Field<FieldType> & zeta,
               const double a,
               const double dx,
               const double Rbar,
               const double fRbar,
               const metadata & sim)
{
  Site x(u.lattice());
  double error, Y, dY, temp, coeff1 = a*a/3., dx2 = dx*dx;
  int count = 0, n=1;

  // Build initial guesses for deltaR and u
  add_fields(eightpiG_deltaT, -1., zeta, 1., deltaR);
  convert_deltaR_to_u(u, deltaR, Rbar, fRbar, sim);

  while(true)
  {
    // COUT << "count = " << count << endl;
    // check_field(deltaR, "deltaR", 64*64*64);
    // check_field(zeta, "zeta", 64*64*64);
    // check_field(u, "u", 64*64*64);
    // check_field(diff_u, "diff_u", 64*64*64);
    // cin.get();

    // Compute new guess for u and deltaR
    update_u(u, diff_u, deltaR, eightpiG_deltaT, coeff1, Rbar, sim);
    convert_u_to_deltaR(eightpiG_deltaT, deltaR, u, Rbar, fRbar, sim);
    add_fields(deltaR, eightpiG_deltaT, zeta);
    // Build diff_u = fRbar * exp(u) * ( laplace u + (grad u)**2 ) -- corresponds to laplace(xi)
    build_diff_u(u, diff_u, dx, fRbar);
    // compute residual at each point
    error = compute_residual_u(u, diff_u, deltaR, eightpiG_deltaT, coeff1);
    error /= pow(sim.numpts, 3);

    COUT << "error, relax_error = " << error << ", " << sim.fR_relax_error << endl;

    if(error <= sim.fR_relax_error)
    {
      break;
    }
    if(count % (n*10) == 0)
    {
      COUT << ".";
      n++;
    }
    count++;
  }

  COUT << " Error = " << error << " after " << count << " cycles" << endl;

  // Re-convert u -> xi
  convert_u_to_xi(u, u, fRbar);

  return error;
}


// TODO: add description
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



// MultiGrid Relaxation for u
template <class FieldType>
double MultiGrid_relax_u(MultiField<FieldType> * u,
                         MultiField<FieldType> * diff_u,
                         MultiField<FieldType> * deltaR,
                         MultiField<FieldType> * eightpiG_deltaT,
                         MultiField<FieldType> * residual,
                         MultiField<FieldType> * err,
                         MultiGrid & engine,
                         Field<FieldType> & zeta,
                         const double a,
                         const double dx,
                         const double Rbar,
                         const double fRbar,
                         const metadata & sim)
{
  double error, Y, dY, temp, coeff1 = a*a/3., dx2 = dx*dx, old_dx, new_dx;
  int old_level, new_level;
  double w = 0.;
  int num_steps = 8;
  int level[] = {2, 1, 2, 1, 0, 1, 2, 1, 0};
  int max_level = 0;
  for(int g=0; g<num_steps; g++)
  {
    if(level[g] > max_level)
    {
      max_level = level[g];
    }
  }

  // Build initial guesses for deltaR and u
  add_fields(eightpiG_deltaT[0], -1., zeta, 1., deltaR[0]);
  convert_deltaR_to_u(u[0], deltaR[0], Rbar, fRbar, sim);

  // Move to desired initial grid (restrict if necessary)
  restrict_to_level(eightpiG_deltaT, engine, max_level);

  error = 2 * sim.fR_relax_error;
  while(error > sim.fR_relax_error)
  {
    restrict_to_level(deltaR, engine, max_level);
    restrict_to_level(u, engine, max_level);

    check_field(deltaR[0], "deltaR[0]", 64*64*64);
    check_field(deltaR[1], "deltaR[1]", 64*64*64/8);
    check_field(deltaR[2], "deltaR[2]", 64*64*64/64);
    check_field(u[0], "u[0]", 64*64*64);
    check_field(u[1], "u[1]", 64*64*64/8);
    check_field(u[2], "u[2]", 64*64*64/64);

    for(int j=0; j<num_steps; j++)
    {
      old_level = level[j];
      new_level = level[j+1];
      // Rescale dx for coarser/finer grids
      old_dx = dx * pow(2, old_level);
      new_dx = dx * pow(2, new_level);


      // Restricting:
      if(old_level - new_level == -1)
      {
        if(parallel.layer(engine.player(old_level)).isPartLayer()) // Makes sure we are on the correct layer
        {// Build diff_u and compute residual (on old level)
          u[old_level].updateHalo();
          build_diff_u(u[old_level], diff_u[old_level], old_dx, fRbar);
          compute_residual_u(u[old_level], diff_u[old_level], deltaR[old_level], eightpiG_deltaT[old_level], coeff1);
        }
        // Restrict guess and residual
        engine.restrict(u, old_level);
        engine.restrict(residual, old_level);
        // Relaxation steps (if required) -- See relax_u
        if(parallel.layer(engine.player(new_level)).isPartLayer())
        {// Copy the old guess on {error} -- we will need it later
          copy_field(u[new_level], err[new_level]);
          for(int h = 0; h < sim.relax_steps; h++)
          {
            u[new_level].updateHalo();
            build_diff_u(u[new_level], diff_u[new_level], new_dx, fRbar);
            update_u(u[new_level], diff_u[new_level], deltaR[new_level], eightpiG_deltaT[new_level], coeff1, Rbar, sim, 1);
            convert_u_to_deltaR(eightpiG_deltaT[new_level], deltaR[new_level], u[new_level], Rbar, fRbar, sim);
            // Updates the error = (new guess) - (old guess)
            add_fields(u[new_level], 1., err[new_level], -1., err[new_level]);
            error = compute_residual_u(u[new_level], diff_u[new_level], deltaR[new_level], eightpiG_deltaT[new_level], coeff1);
            if(error <= sim.fR_relax_error)
            {
              break;
            }
          }
        }
      }
      // Prolonging
      else if(old_level - new_level == 1)
      {
        // Prolong the error
        engine.prolong(err, old_level);
        // build new guess from this updated error
        if(parallel.layer(engine.player(new_level)).isPartLayer())
        {
          add_fields(u[new_level], err[new_level], u[new_level]);
          convert_u_to_deltaR(eightpiG_deltaT[new_level], deltaR[new_level], u[new_level], Rbar, fRbar, sim);
        }
      }
      else
      {
        COUT << "/!\\ There's something wrong with the MultiGrid steps!!! /!\\" << endl;
      }
    }

    if(parallel.layer(engine.player(0)).isPartLayer()) // Final level must be zero!
    {
      u[0].updateHalo();
      build_diff_u(u[0], diff_u[0], dx, fRbar);
      error = compute_residual_u(u[0], diff_u[0], deltaR[0], eightpiG_deltaT[0], coeff1);
    }

    COUT << " Error, relax_error = " << error << ", " << sim.fR_relax_error << endl;
    cin.get();

  }

  return error;
}


// Find a solution to the quasi-static trace equation with an iterative FFT scheme
template <class FieldType>
double relax_u_FFT(Field<FieldType> & u,
                   Field<FieldType> & deltaR,
                   Field<FieldType> & zeta,
                   Field<FieldType> & eightpiG_deltaT,
                   Field<FieldType> & source,
                   Field<Cplx> & scalarFT,
                   PlanFFT<Cplx> * plan_source,
                   PlanFFT<Cplx> * plan_u,
                   const double a,
                   const double dx,
                   const double Rbar,
                   const double fRbar,
                   const metadata & sim)
{
  double error, coeff1 = a*a/3./fRbar, dx2 = dx*dx, temp;
  int count = 0;
  Site x(u.lattice());
  error = 2. * sim.fR_relax_error;

  // convert xi to u (xi is passed to function in place of u)
  convert_xi_to_u(u, u, fRbar);

  while(error > sim.fR_relax_error)
  {
    error = 0.;

    // Builds source term
    for(x.first(); x.test(); x.next())
    {
      source(x) = (u(x+0) - u(x-0))*(u(x+0) - u(x-0)) +
                  (u(x+1) - u(x-1))*(u(x+1) - u(x-1)) +
                  (u(x+2) - u(x-2))*(u(x+2) - u(x-2));
      source(x) /= - 4. * dx2;
      source(x) += coeff1 * zeta(x);
      source(x) /= exp(u(x));
    }


    plan_source->execute(FFT_FORWARD);
    solveModifiedPoissonFT(scalarFT, scalarFT, 1.); // Solves laplace equation for u
    plan_u->execute(FFT_BACKWARD);
    u.updateHalo();

    convert_u_to_deltaR(eightpiG_deltaT, deltaR, u, Rbar, fRbar, sim);
    add_fields(deltaR, eightpiG_deltaT, zeta);

    // Check equation
    for(x.first(); x.test(); x.next())
    {
      temp = (u(x+0) - u(x-0))*(u(x+0) - u(x-0)) +
             (u(x+1) - u(x-1))*(u(x+1) - u(x-1)) +
             (u(x+2) - u(x-2))*(u(x+2) - u(x-2));
      temp /= 4.;
      temp += u(x+0) + u(x-0) + u(x+1) + u(x-1) + u(x+2) + u(x-2) - 6.*u(x);
      temp /= dx2;
      temp -= coeff1 * zeta(x);
      temp *= temp;
      if(temp > error)
      {
        error = temp;
      }
    }
    count++;
    parallel.max(error);
  }

  // Convert back to xi
  convert_u_to_xi(u, u, fRbar);

  return error;
}






#endif
