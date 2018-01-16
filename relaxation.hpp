//////////////////////////
// relaxation.hpp
//////////////////////////
//
// Implementation of relaxation method solver for F(R) gravity
//
// Authors: David Daverio (Cambridge University)
//          Lorenzo Reverberi (University of Cape Town; CEICO - Czech Academy of Sciences, Prague)
//
//////////////////////////

#ifndef RELAXATION_HEADER
#define RELAXATION_HEADER


// Relaxation method to solve the quasi-static trace equation for xi
template <class FieldType>
double relax_xi(Field<FieldType> & xi,
                Field<FieldType> & laplace_xi,
                Field<FieldType> & deltaR,
                Field<FieldType> & zeta,
                Field<FieldType> & eightpiG_deltaT,
                const double a,
                const double dx,
                const double Rbar,
                const double fRbar,
                const metadata & sim)
{
  Site x(xi.lattice());
  double error, Y, dY, dxi, temp, coeff1 = a*a/3., dx2 = dx*dx;
  double maxY, maxdY, maxdxi, maxexpon;
  double signY, signdY, signdxi, signexpon;
  int count = 0, n = 1;

  // Build initial guesses for xi and laplace_xi
  for(x.first(); x.test(); x.next())
  {
    deltaR(x) = - eightpiG_deltaT(x) + zeta(x);
    xi(x) = fR(Rbar + deltaR(x), sim, 652) - fRbar;
  }
  deltaR.updateHalo();
  xi.updateHalo();

  build_laplacian(xi, laplace_xi, dx);

  do
  {
    error = 0.;
    maxY = maxdY = maxdxi = 0.;
    // Update xi with next guess
    check_field(xi, "xi", 1, "first");
    check_field(deltaR, "deltaR", 1, "first");
    check_field(zeta, "zeta", 1, "first");
    check_field(laplace_xi, "laplace_xi", 1, "first");
    for(x.first(); x.test(); x.next())
    {
      Y = laplace_xi(x) - coeff1 * zeta(x);
      dY = laplace_xi(x) - coeff1 * fR(deltaR(x), sim, 589) / fRR(deltaR(x), sim, 590);
      dxi = exp(-Y/dY) * (xi(x) + fRbar) - fRbar - xi(x);
      if(fabs(Y) > fabs(maxY))
      {
        maxY = Y;
        signY = maxY > 0 ? 1 : -1;
      }
      if(fabs(dY) > fabs(maxdY))
      {
        maxdY = dY;
        signdY = maxdY > 0 ? 1 : -1;
      }
      if(fabs(dxi) > fabs(maxdxi))
      {
        maxdxi = dxi;
        signdxi = maxdxi > 0 ? 1 : -1;
        maxexpon = exp(-Y/dY);
      }

      xi(x) += dxi;
    }
    if(signY > 0) parallel.max(maxY);
    else parallel.min(maxY);
    if(signdY > 0) parallel.max(maxdY);
    else parallel.min(maxY);
    if(signdxi > 0) parallel.max(maxdxi);
    else parallel.min(maxY);
    parallel.max(maxexpon);
    COUT << " maxY = " << maxY << "  maxdY = " << maxdY << "  maxdxi = " << maxdxi << "  maxexpon = " <<  maxexpon << endl;

    xi.updateHalo();
    check_field(xi, "xi", 1, "second");

    // Update zeta and deltaR given the new xi
    computeDRzeta(eightpiG_deltaT, deltaR, zeta, xi, Rbar, fRbar, sim);
    check_field(deltaR, "deltaR", 1, "second");
    check_field(zeta, "zeta", 1, "second");

    // Compute laplace_xi and total error
    for(x.first(); x.test(); x.next())
    {
      laplace_xi(x) = xi(x+0) + xi(x-0) + xi(x+1) + xi(x-1) + xi(x+2) + xi(x-2) - 6.*xi(x);
      laplace_xi(x) /= dx2;
      if(laplace_xi(x) && zeta(x))
      {
        temp = coeff1 * zeta(x) / laplace_xi(x) - 1.;
        temp *= temp;
        if(temp >= error)
        {
          error = temp;
        }
      }
      else
      {
        error = 1.;
      }
    }
    laplace_xi.updateHalo();

    parallel.max(error);

    count ++;
    if(count >= n*100)
    {
      COUT << " Count reached " << n*100 << "  error = " << error << endl;
      n++;
    }
  } while(error > sim.fR_relax_error);

  return error;
}


// // Relaxation method to solve the quasi-static trace equation for xi
// template <class FieldType>
// double relax_xi_Fourier(Field<FieldType> & xi,
//                         Field<FieldType> & laplace_xi,
//                         Field<FieldType> & deltaR,
//                         Field<FieldType> & zeta,
//                         Field<FieldType> & eightpiG_deltaT,
//                         PlanFFT<Cplx> * plan_xi,
//                         const double a,
//                         const double dx,
//                         const double Rbar,
//                         const double fRbar,
//                         const metadata & sim)
// {
//   Site x(xi.lattice());
//   double error, Y, dY, dxi, temp, coeff1 = a*a/3., dx2 = dx*dx;
//   int count = 0, n = 1;
//
//   // Build initial guesses for xi and laplace_xi
//   for(x.first(); x.test(); x.next())
//   {
//     xi(x) = fR(Rbar - eightpiG_deltaT(x) + zeta(x), sim, 653);
//   }
//   xi.updateHalo();
//
//   build_laplacian(xi, laplace_xi, dx);
//
//   do
//   {
//     error = 0.;
//
//     // Update xi with next guess
//     for(x.first(); x.test(); x.next())
//     {
//       temp = Rbar - eightpiG_deltaT(x) + zeta(x);
//       Y = laplace_xi(x) - coeff1 * zeta(x);
//       dY = laplace_xi(x) - coeff1 * fR(temp, sim, 589) / fRR(temp, sim, 590);
//       dxi = exp(-Y/dY) * (xi(x) + fRbar) - fRbar - xi(x);
//       xi(x) += dxi;
//     }
//     xi.updateHalo();
//
//     // Update zeta and deltaR given the new xi
//     computeDRzeta(eightpiG_deltaT, deltaR, zeta, xi, Rbar, fRbar, sim);
//
//     // Compute laplace_xi and total error
//     for(x.first(); x.test(); x.next())
//     {
//       laplace_xi(x) = xi(x+0) + xi(x-0) + xi(x+1) + xi(x-1) + xi(x+2) + xi(x-2) - 6.*xi(x);
//       laplace_xi(x) /= dx2;
//       if(laplace_xi(x) && zeta(x))
//       {
//         temp = coeff1 * zeta(x) / laplace_xi(x) - 1.;
//         temp *= temp;
//         if(temp >= error)
//         {
//           error = temp;
//         }
//       }
//       else
//       {
//         error = 1.;
//       }
//     }
//
//     parallel.max(error);
//     laplace_xi.updateHalo();
//
//     count ++;
//     if(count >= n*100)
//     {
//       COUT << " Count reached " << n*100 << "  error = " << error << endl;
//       n++;
//     }
//   } while(error > sim.fR_relax_error);
//
//
// }
//
//
//













#endif
