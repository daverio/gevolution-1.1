//////////////////////////
// gevolution.hpp
//////////////////////////
//
// Geneva algorithms for evolution of metric perturbations
// and relativistic free-streaming particles (gevolution)
//
// 1. Suite of Fourier-based methods for the computation of the
//    relativistic scalar (Phi, Phi-Psi) and vector modes [see J. Adamek,
//    R. Durrer, and M. Kunz, Class. Quant. Grav. 31, 234006 (2014)]
//
// 2. Collection of "update position" and "update velocity/momentum" methods
//    [see J. Adamek, D. Daverio, R. Durrer, and M. Kunz, JCAP 1607, 053 (2016)]
//
// 3. Collection of projection methods for the construction of the
//    stress-energy-tensor
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: December 2016
//
//////////////////////////

#ifndef GEVOLUTION_HEADER
#define GEVOLUTION_HEADER

#ifndef Cplx
#define Cplx Imag
#endif

using namespace std;
using namespace LATfield2;


//////////////////////////
// prepareFTsource (1)
//////////////////////////
// Description:
//   construction of real-space source tensor for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration
//   Tij        reference to symmetric tensor field containing the space-space
//              components of the stress-energy tensor (rescaled by a^3)
//   Sij        reference to allocated symmetric tensor field which will contain
//              the source tensor (may be identical to Tji)
//   coeff      scaling coefficient for Tij ("8 pi G dx^2 / a")
//
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource(Field<FieldType> & phi, Field<FieldType> & Tij, Field<FieldType> & Sij, const double coeff)
{
	Site x(phi.lattice());

	for (x.first(); x.test(); x.next())
	{
		// 0-0-component:
		Sij(x, 0, 0) = coeff * Tij(x, 0, 0);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 0, 0) -= 4. * phi(x) * (phi(x-0) + phi(x+0) - 2. * phi(x));
		Sij(x, 0, 0) -= 0.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0));
#else
		Sij(x, 0, 0) += 0.5 * (phi(x+0) - phi(x-0)) * (phi(x+0) - phi(x-0));
#endif
#endif

		// 1-1-component:
		Sij(x, 1, 1) = coeff * Tij(x, 1, 1);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 1, 1) -= 4. * phi(x) * (phi(x-1) + phi(x+1) - 2. * phi(x));
		Sij(x, 1, 1) -= 0.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1));
#else
		Sij(x, 1, 1) += 0.5 * (phi(x+1) - phi(x-1)) * (phi(x+1) - phi(x-1));
#endif
#endif

		// 2-2-component:
		Sij(x, 2, 2) = coeff * Tij(x, 2, 2);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		Sij(x, 2, 2) -= 4. * phi(x) * (phi(x-2) + phi(x+2) - 2. * phi(x));
		Sij(x, 2, 2) -= 0.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2));
#else
		Sij(x, 2, 2) += 0.5 * (phi(x+2) - phi(x-2)) * (phi(x+2) - phi(x-2));
#endif
#endif

		// 0-1-component:
		Sij(x, 0, 1) = coeff * Tij(x, 0, 1);
#ifdef PHINONLINEAR
		Sij(x, 0, 1) += phi(x+0) * phi(x+1) - phi(x) * phi(x+0+1);
#ifdef ORIGINALMETRIC
		Sij(x, 0, 1) -= 1.5 * phi(x) * phi(x);
		Sij(x, 0, 1) += 1.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 1) += 1.5 * phi(x+1) * phi(x+1);
		Sij(x, 0, 1) -= 1.5 * phi(x+0+1) * phi(x+0+1);
#else
		Sij(x, 0, 1) += 0.5 * phi(x) * phi(x);
		Sij(x, 0, 1) -= 0.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 1) -= 0.5 * phi(x+1) * phi(x+1);
		Sij(x, 0, 1) += 0.5 * phi(x+0+1) * phi(x+0+1);
#endif
#endif

		// 0-2-component:
		Sij(x, 0, 2) = coeff * Tij(x, 0, 2);
#ifdef PHINONLINEAR
		Sij(x, 0, 2) += phi(x+0) * phi(x+2) - phi(x) * phi(x+0+2);
#ifdef ORIGINALMETRIC
		Sij(x, 0, 2) -= 1.5 * phi(x) * phi(x);
		Sij(x, 0, 2) += 1.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 2) += 1.5 * phi(x+2) * phi(x+2);
		Sij(x, 0, 2) -= 1.5 * phi(x+0+2) * phi(x+0+2);
#else
		Sij(x, 0, 2) += 0.5 * phi(x) * phi(x);
		Sij(x, 0, 2) -= 0.5 * phi(x+0) * phi(x+0);
		Sij(x, 0, 2) -= 0.5 * phi(x+2) * phi(x+2);
		Sij(x, 0, 2) += 0.5 * phi(x+0+2) * phi(x+0+2);
#endif
#endif

		// 1-2-component:
		Sij(x, 1, 2) = coeff * Tij(x, 1, 2);
#ifdef PHINONLINEAR
		Sij(x, 1, 2) += phi(x+1) * phi(x+2) - phi(x) * phi(x+1+2);
#ifdef ORIGINALMETRIC
		Sij(x, 1, 2) -= 1.5 * phi(x) * phi(x);
		Sij(x, 1, 2) += 1.5 * phi(x+1) * phi(x+1);
		Sij(x, 1, 2) += 1.5 * phi(x+2) * phi(x+2);
		Sij(x, 1, 2) -= 1.5 * phi(x+1+2) * phi(x+1+2);
#else
		Sij(x, 1, 2) += 0.5 * phi(x) * phi(x);
		Sij(x, 1, 2) -= 0.5 * phi(x+1) * phi(x+1);
		Sij(x, 1, 2) -= 0.5 * phi(x+2) * phi(x+2);
		Sij(x, 1, 2) += 0.5 * phi(x+1+2) * phi(x+1+2);
#endif
#endif
	}
}

template <class FieldType>
void prepareFTsource_S00(Field<FieldType> & source,
												 Field<FieldType> & phi,
												 Field<FieldType> & chi,
												 Field<FieldType> & xi,
											 	 Field<FieldType> & deltaR,
												 Field<FieldType> & T00,
												 double bgmodel,
												 double dx2,
												 double dtau,
												 double Hubble,
												 double a2,
												 double eightpiG,
												 double Rbar,
												 double Fbar,
												 double FRbar,
												 double * paramF,
												 const int Ftype)
{
	Site x(phi.lattice());
	double laplace;
	double grad[3];

	for (x.first(); x.test(); x.next())
	{
		source(x) = eightpiG * a2 * (T00(x) - bgmodel) + 6.0*Hubble*Hubble*chi(x) + 3.0*Hubble * (2.0*phi(x)-xi(x))/dtau;
		//add laplace phi
		laplace=0.0;
		for(int i =0;i<3;i++)laplace += phi(x+i)+phi(x-i);
		laplace -= 6.0*phi(x);
		source(x) += (8.0*phi(x) + xi(x) + FRbar)*laplace/dx2;
		//add gradient^2 phi
		for(int i =0;i<3;i++)
		{
			grad[i] = phi(x+i)-phi(x-i);
			grad[i] *= grad[i];
		}
		//f(R) terms
		source(x) -= 0.5 * a2 * (Rbar*xi(x) + Fbar - F(Rbar+ deltaR(x),paramF,Ftype));
		source(x) += (grad[0] + grad[1] + grad[2]) * 0.75/dx2;
	}
}

template <class FieldType>
void prepareFTsource_S0i(Field<FieldType> & S0i,
												 Field<FieldType> & phi,
												 Field<FieldType> & Bi,
												 double eightpiGa2,
											 	 double a2)
{
	Site x(phi.lattice());

	double phicjk[3][3];

	for (x.first(); x.test(); x.next())
	{
			S0i(x) *= eightpiGa2;
			//for(int i=0;i<3;i++)phicjk[i][i] =
	}
}

template <class FieldType>
void projectFTsource_S0i(Field<FieldType> & S0i,
												 Field<FieldType> & output)
{

	const int linesize = S0i.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	Real * sinc;
	rKSite k(S0i.lattice());
	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));
	double temp = (long)linesize * (long)linesize * (long)linesize;
	Cplx im = Cplx(0.0,1.0);


	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		output(k) = Cplx(0.,0.);
		k.next();
	}

	for (; k.test(); k.next())
	{
		output(k) = S0i(k,0)*kshift[k.coord(0)] + S0i(k,1)*kshift[k.coord(1)] + S0i(k,2)*kshift[k.coord(2)];
		output(k) /= gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
		output(k) /= temp;
		output(k) *= im;
	}

	free(gridk2);
	free(kshift);
}

template <class FieldType>
void stepXi(Field<FieldType> & xi,
							Field<FieldType> & source,
							Field<FieldType> & phidot,
							Field<FieldType> & phi,
							Field<FieldType> & chi,
							double H,
							double dtau)
{
	Site x(xi.lattice());
	FieldType * pointer;
	for(x.first();x.test();x.next())
	{
		xi(x) = (source(x) + (phidot(x)-xi(x)+2.0*phi(x))/dtau )/H + 2.0 *chi(x);

		phidot(x) = 0.25*(xi(x)-phidot(x));
		phi(x) = (phidot(x) - phi(x))/dtau;
		xi(x) = xi(x) - 2.0* phidot(x);
	}
	pointer = phidot.data_;
	phidot.data_ = phi.data_;
	phi.data_ = pointer;
}

template <class FieldType>
void computeTtrace(Field<FieldType> & T,
							Field<FieldType> & T00,
							Field<FieldType> & Tij)
{
	Site x(T.lattice());
	for(x.first();x.test();x.next())
	{
		T(x) = T00(x) + Tij(x,0,0) + Tij(x,1,1) + Tij(x,2,2);
	}
}

template <class FieldType>
void addXi(Field<FieldType> & chi,
							Field<FieldType> & xi)
{
	Site x(chi.lattice());
	for(x.first();x.test();x.next())
	{
		chi(x)+=xi(x);
	}
}

template <class FieldType>
void computeDRzeta(Field<FieldType> & deltaR,
									 Field<FieldType> & zeta,
									 Field<FieldType> & xi,
									 double Rbar,
								 	 double FRbar,
								 	 double eightpiG,
								 	 double Tbar,
								 	 double * params,
									 const int fofR_type)
{
	Site x(deltaR.lattice());
	double temp,temp2;
	for(x.first();x.test();x.next())
	{
		temp = eightpiG * (deltaR(x)-Tbar);
		temp2 = FRR(Rbar - temp, params,fofR_type);
		if(temp2!=0)
		{
			zeta(x) = (xi(x) + FRbar + FR(Rbar - temp, params,fofR_type))/temp2;
			deltaR(x) = zeta(x) - temp;
		}
		else
		{
			COUT<<"FRR is equal to 0, aborting"<<endl;
			exit(2);
		}

	}
}

//////////////////////////
// prepareFTsource (2)
//////////////////////////
// Description:
//   construction of real-space source field for Fourier-based solvers
//
// Arguments:
//   phi        reference to field configuration (first Bardeen potential)
//   chi        reference to field configuration (difference between Bardeen potentials, phi-psi)
//   source     reference to fully dressed source field (rescaled by a^3)
//   bgmodel    background model of the source (rescaled by a^3) to be subtracted
//   result     reference to allocated field which will contain the result (may be identical to source)
//   coeff      diffusion coefficient ("3 H_conformal dx^2 / dtau")
//   coeff2     scaling coefficient for the source ("4 pi G dx^2 / a")
//   coeff3     scaling coefficient for the psi-term ("3 H_conformal^2 dx^2")
//
// Returns:
//
//////////////////////////

template <class FieldType>
void prepareFTsource(Field<FieldType> & phi, Field<FieldType> & chi, Field<FieldType> & source, const FieldType bgmodel, Field<FieldType> & result, const double coeff, const double coeff2, const double coeff3)
{
	Site x(phi.lattice());

	for (x.first(); x.test(); x.next())
	{
		result(x) = coeff2 * (source(x) - bgmodel);
#ifdef PHINONLINEAR
#ifdef ORIGINALMETRIC
		result(x) *= 1. - 4. * phi(x);
		result(x) -= 0.375 * (phi(x-0) - phi(x+0)) * (phi(x-0) - phi(x+0));
		result(x) -= 0.375 * (phi(x-1) - phi(x+1)) * (phi(x-1) - phi(x+1));
		result(x) -= 0.375 * (phi(x-2) - phi(x+2)) * (phi(x-2) - phi(x+2));
#else
		result(x) *= 1. - 2. * phi(x);
		result(x) += 0.125 * (phi(x-0) - phi(x+0)) * (phi(x-0) - phi(x+0));
		result(x) += 0.125 * (phi(x-1) - phi(x+1)) * (phi(x-1) - phi(x+1));
		result(x) += 0.125 * (phi(x-2) - phi(x+2)) * (phi(x-2) - phi(x+2));
#endif
#endif
		result(x) += (coeff3 - coeff) * phi(x) - coeff3 * chi(x);
	}
}


#ifdef FFT3D
//////////////////////////
// projectFTscalar
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the trace-free
//   longitudinal (scalar) component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   chiFT      reference to allocated field which will contain the Fourier
//              image of the trace-free longitudinal (scalar) component
//
// Returns:
//
//////////////////////////

void projectFTscalar(Field<Cplx> & SijFT, Field<Cplx> & chiFT, const int add = 0)
{
	const int linesize = chiFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(chiFT.lattice());

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		chiFT(k) = Cplx(0.,0.);
		k.next();
	}

	if (add)
	{
		for (; k.test(); k.next())
		{
			chiFT(k) += ((gridk2[k.coord(1)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(0)]) * SijFT(k, 0, 0) +
						(gridk2[k.coord(0)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(1)]) * SijFT(k, 1, 1) +
						(gridk2[k.coord(0)] + gridk2[k.coord(1)] - 2. * gridk2[k.coord(2)]) * SijFT(k, 2, 2) -
						6. * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1) -
						6. * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2) -
						6. * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) /
						(2. * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * linesize);
		}
	}
	else
	{
		for (; k.test(); k.next())
		{
			chiFT(k) = ((gridk2[k.coord(1)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(0)]) * SijFT(k, 0, 0) +
						(gridk2[k.coord(0)] + gridk2[k.coord(2)] - 2. * gridk2[k.coord(1)]) * SijFT(k, 1, 1) +
						(gridk2[k.coord(0)] + gridk2[k.coord(1)] - 2. * gridk2[k.coord(2)]) * SijFT(k, 2, 2) -
						6. * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1) -
						6. * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2) -
						6. * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2)) /
						(2. * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)]) * linesize);
		}
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// evolveFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a tensor field on the spin-1 component
//   used as a source for the evolution of the vector perturbation
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   BiFT       reference to the Fourier image of the vector perturbation
//   a2dtau     conformal time step times scale factor squared (a^2 * dtau)
//
// Returns:
//
//////////////////////////

void evolveFTvector(Field<Cplx> & SijFT, Field<Cplx> & BiFT, const Real a2dtau)
{
	const int linesize = BiFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(BiFT.lattice());
	Real k4;

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		BiFT(k, 0) = Cplx(0.,0.);
		BiFT(k, 1) = Cplx(0.,0.);
		BiFT(k, 2) = Cplx(0.,0.);
		k.next();
	}

	for (; k.test(); k.next())
	{
		k4 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
		k4 *= k4;

		BiFT(k, 0) += Cplx(0.,-2.*a2dtau/k4) * (kshift[k.coord(0)].conj() * ((gridk2[k.coord(1)] + gridk2[k.coord(2)]) * SijFT(k, 0, 0)
				- gridk2[k.coord(1)] * SijFT(k, 1, 1) - gridk2[k.coord(2)] * SijFT(k, 2, 2) - 2. * kshift[k.coord(1)] * kshift[k.coord(2)] * SijFT(k, 1, 2))
				+ (gridk2[k.coord(1)] + gridk2[k.coord(2)] - gridk2[k.coord(0)]) * (kshift[k.coord(1)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 0, 2)));
		BiFT(k, 1) += Cplx(0.,-2.*a2dtau/k4) * (kshift[k.coord(1)].conj() * ((gridk2[k.coord(0)] + gridk2[k.coord(2)]) * SijFT(k, 1, 1)
				- gridk2[k.coord(0)] * SijFT(k, 0, 0) - gridk2[k.coord(2)] * SijFT(k, 2, 2) - 2. * kshift[k.coord(0)] * kshift[k.coord(2)] * SijFT(k, 0, 2))
				+ (gridk2[k.coord(0)] + gridk2[k.coord(2)] - gridk2[k.coord(1)]) * (kshift[k.coord(0)] * SijFT(k, 0, 1) + kshift[k.coord(2)] * SijFT(k, 1, 2)));
		BiFT(k, 2) += Cplx(0.,-2.*a2dtau/k4) * (kshift[k.coord(2)].conj() * ((gridk2[k.coord(0)] + gridk2[k.coord(1)]) * SijFT(k, 2, 2)
				- gridk2[k.coord(0)] * SijFT(k, 0, 0) - gridk2[k.coord(1)] * SijFT(k, 1, 1) - 2. * kshift[k.coord(0)] * kshift[k.coord(1)] * SijFT(k, 0, 1))
				+ (gridk2[k.coord(0)] + gridk2[k.coord(1)] - gridk2[k.coord(2)]) * (kshift[k.coord(0)] * SijFT(k, 0, 2) + kshift[k.coord(1)] * SijFT(k, 1, 2)));
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// projectFTvector
//////////////////////////
// Description:
//   projects the Fourier image of a vector field on the transverse component
//   and solves the constraint equation for the vector perturbation
//
// Arguments:
//   SiFT       reference to the Fourier image of the input vector field
//   BiFT       reference to the Fourier image of the vector perturbation (can be identical to input)
//   coeff      rescaling coefficient (default 1)
//   modif      modification k^2 -> k^2 + modif (default 0)
//
// Returns:
//
//////////////////////////

void projectFTvector(Field<Cplx> & SiFT, Field<Cplx> & BiFT, const Real coeff = 1., const Real modif = 0.)
{
	const int linesize = BiFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(BiFT.lattice());
	Real k2;
	Cplx tmp(0., 0.);

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		BiFT(k, 0) = Cplx(0.,0.);
		BiFT(k, 1) = Cplx(0.,0.);
		BiFT(k, 2) = Cplx(0.,0.);
		k.next();
	}

	for (; k.test(); k.next())
	{
		k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];

		tmp = (kshift[k.coord(0)] * SiFT(k, 0) + kshift[k.coord(1)] * SiFT(k, 1) + kshift[k.coord(2)] * SiFT(k, 2)) / k2;

		BiFT(k, 0) = (SiFT(k, 0) - kshift[k.coord(0)].conj() * tmp) * 4. * coeff / (k2 + modif);
		BiFT(k, 1) = (SiFT(k, 1) - kshift[k.coord(1)].conj() * tmp) * 4. * coeff / (k2 + modif);
		BiFT(k, 2) = (SiFT(k, 2) - kshift[k.coord(2)].conj() * tmp) * 4. * coeff / (k2 + modif);
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// projectFTtensor
//////////////////////////
// Description:
//   projection of the Fourier image of a tensor field on the transverse
//   trace-free tensor component
//
// Arguments:
//   SijFT      reference to the Fourier image of the input tensor field
//   hijFT      reference to allocated field which will contain the Fourier
//              image of the transverse trace-free tensor component
//
// Returns:
//
//////////////////////////

void projectFTtensor(Field<Cplx> & SijFT, Field<Cplx> & hijFT)
{
	const int linesize = hijFT.lattice().size(1);
	int i;
	Real * gridk2;
	Cplx * kshift;
	rKSite k(hijFT.lattice());
	Cplx SxxFT, SxyFT, SxzFT, SyyFT, SyzFT, SzzFT;
	Real k2, k6;

	gridk2 = (Real *) malloc(linesize * sizeof(Real));
	kshift = (Cplx *) malloc(linesize * sizeof(Cplx));

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		kshift[i] = gridk2[i] * Cplx(cos(M_PI * (Real) i / (Real) linesize), -sin(M_PI * (Real) i / (Real) linesize));
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		for (i = 0; i < hijFT.components(); i++)
			hijFT(k, i) = Cplx(0.,0.);

		k.next();
	}

	for (; k.test(); k.next())
	{
		SxxFT = SijFT(k, 0, 0);
		SxyFT = SijFT(k, 0, 1);
		SxzFT = SijFT(k, 0, 2);
		SyyFT = SijFT(k, 1, 1);
		SyzFT = SijFT(k, 1, 2);
		SzzFT = SijFT(k, 2, 2);

		k2 = gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)];
		k6 = k2 * k2 * k2 * linesize;

		hijFT(k, 0, 0) = ((gridk2[k.coord(0)] - k2) * ((gridk2[k.coord(0)] - k2) * SxxFT + 2. * kshift[k.coord(0)] * (kshift[k.coord(1)] * SxyFT + kshift[k.coord(2)] * SxzFT))
				+ ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SyyFT
				+ ((gridk2[k.coord(0)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SzzFT
				+ 2. * (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)] * kshift[k.coord(2)] * SyzFT) / k6;

		hijFT(k, 0, 1) = (2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(1)] - k2) * SxyFT + (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(1)].conj() * SzzFT
				+ (gridk2[k.coord(0)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(0)].conj() * SxxFT + 2. * kshift[k.coord(2)] * SxzFT)
				+ (gridk2[k.coord(1)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(1)].conj() * SyyFT + 2. * kshift[k.coord(2)] * SyzFT)) / k6;

		hijFT(k, 0, 2) = (2. * (gridk2[k.coord(0)] - k2) * (gridk2[k.coord(2)] - k2) * SxzFT + (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)].conj() * kshift[k.coord(2)].conj() * SyyFT
				+ (gridk2[k.coord(0)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(0)].conj() * SxxFT + 2. * kshift[k.coord(1)] * SxyFT)
				+ (gridk2[k.coord(2)] - k2) * kshift[k.coord(0)].conj() * (kshift[k.coord(2)].conj() * SzzFT + 2. * kshift[k.coord(1)] * SyzFT)) / k6;

		hijFT(k, 1, 1) = ((gridk2[k.coord(1)] - k2) * ((gridk2[k.coord(1)] - k2) * SyyFT + 2. * kshift[k.coord(1)] * (kshift[k.coord(0)] * SxyFT + kshift[k.coord(2)] * SyzFT))
				+ ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SxxFT
				+ ((gridk2[k.coord(1)] + k2) * (gridk2[k.coord(2)] + k2) - 2. * k2 * k2) * SzzFT
				+ 2. * (gridk2[k.coord(1)] + k2) * kshift[k.coord(0)] * kshift[k.coord(2)] * SxzFT) / k6;

		hijFT(k, 1, 2) = (2. * (gridk2[k.coord(1)] - k2) * (gridk2[k.coord(2)] - k2) * SyzFT + (gridk2[k.coord(0)] + k2) * kshift[k.coord(1)].conj() * kshift[k.coord(2)].conj() * SxxFT
				+ (gridk2[k.coord(1)] - k2) * kshift[k.coord(2)].conj() * (kshift[k.coord(1)].conj() * SyyFT + 2. * kshift[k.coord(0)] * SxyFT)
				+ (gridk2[k.coord(2)] - k2) * kshift[k.coord(1)].conj() * (kshift[k.coord(2)].conj() * SzzFT + 2. * kshift[k.coord(0)] * SxzFT)) / k6;

		hijFT(k, 2, 2) = ((gridk2[k.coord(2)] - k2) * ((gridk2[k.coord(2)] - k2) * SzzFT + 2. * kshift[k.coord(2)] * (kshift[k.coord(0)] * SxzFT + kshift[k.coord(1)] * SyzFT))
				+ ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(0)] + k2) - 2. * k2 * k2) * SxxFT
				+ ((gridk2[k.coord(2)] + k2) * (gridk2[k.coord(1)] + k2) - 2. * k2 * k2) * SyyFT
				+ 2. * (gridk2[k.coord(2)] + k2) * kshift[k.coord(0)] * kshift[k.coord(1)] * SxyFT) / k6;
	}

	free(gridk2);
	free(kshift);
}


//////////////////////////
// solveModifiedPoissonFT
//////////////////////////
// Description:
//   Modified Poisson solver using the standard Fourier method
//
// Arguments:
//   sourceFT   reference to the Fourier image of the source field
//   potFT      reference to the Fourier image of the potential
//   coeff      coefficient applied to the source ("4 pi G / a")
//   modif      modification k^2 -> k^2 + modif (default 0 gives standard Poisson equation)
//
// Returns:
//
//////////////////////////

void solveModifiedPoissonFT(Field<Cplx> & sourceFT, Field<Cplx> & potFT, Real coeff, const Real modif = 0.)
{
	const int linesize = potFT.lattice().size(1);
	int i;
	Real * gridk2;
	Real * sinc;
	rKSite k(potFT.lattice());

	gridk2 = (Real *) malloc(linesize * sizeof(Real));

	coeff /= -((long) linesize * (long) linesize * (long) linesize);

	for (i = 0; i < linesize; i++)
	{
		gridk2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
		gridk2[i] *= gridk2[i];
	}

	k.first();
	if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
	{
		if (modif == 0.)
			potFT(k) = Cplx(0.,0.);
		else
			potFT(k) = sourceFT(k) * coeff / modif;
		k.next();
	}

	for (; k.test(); k.next())
	{
		potFT(k) = sourceFT(k) * coeff / (gridk2[k.coord(0)] + gridk2[k.coord(1)] + gridk2[k.coord(2)] + modif);
	}

	free(gridk2);
}
#endif


//////////////////////////
// update_q
//////////////////////////
// Description:
//   Update momentum method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   Sites      array of Sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

Real update_q(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * Sites, int nfield, double * params, double * outputs, int noutputs)
{
#define phi (*fields[0])
#define chi (*fields[1])
#define Bi (*fields[2])
#define xphi (Sites[0])
#define xchi (Sites[1])
#define xB (Sites[2])

	Real gradphi[3]={0,0,0};
	Real pgradB[3]={0,0,0};
	Real v2 = (*part).vel[0] * (*part).vel[0] + (*part).vel[1] * (*part).vel[1] + (*part).vel[2] * (*part).vel[2];
	Real e2 = v2 + params[0] * params[0];

	gradphi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * (phi(xphi+0) - phi(xphi));
	gradphi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * (phi(xphi+1) - phi(xphi));
	gradphi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * (phi(xphi+2) - phi(xphi));
	gradphi[0] += ref_dist[1] * (1.-ref_dist[2]) * (phi(xphi+1+0) - phi(xphi+1));
	gradphi[1] += ref_dist[0] * (1.-ref_dist[2]) * (phi(xphi+1+0) - phi(xphi+0));
	gradphi[2] += ref_dist[0] * (1.-ref_dist[1]) * (phi(xphi+2+0) - phi(xphi+0));
	gradphi[0] += (1.-ref_dist[1]) * ref_dist[2] * (phi(xphi+2+0) - phi(xphi+2));
	gradphi[1] += (1.-ref_dist[0]) * ref_dist[2] * (phi(xphi+2+1) - phi(xphi+2));
	gradphi[2] += (1.-ref_dist[0]) * ref_dist[1] * (phi(xphi+2+1) - phi(xphi+1));
	gradphi[0] += ref_dist[1] * ref_dist[2] * (phi(xphi+2+1+0) - phi(xphi+2+1));
	gradphi[1] += ref_dist[0] * ref_dist[2] * (phi(xphi+2+1+0) - phi(xphi+2+0));
	gradphi[2] += ref_dist[0] * ref_dist[1] * (phi(xphi+2+1+0) - phi(xphi+1+0));

	gradphi[0] *= (v2 + e2) / e2;
	gradphi[1] *= (v2 + e2) / e2;
	gradphi[2] *= (v2 + e2) / e2;

	if (nfield >= 2 && fields[1] != NULL)
	{
		gradphi[0] -= (1.-ref_dist[1]) * (1.-ref_dist[2]) * (chi(xchi+0) - chi(xchi));
		gradphi[1] -= (1.-ref_dist[0]) * (1.-ref_dist[2]) * (chi(xchi+1) - chi(xchi));
		gradphi[2] -= (1.-ref_dist[0]) * (1.-ref_dist[1]) * (chi(xchi+2) - chi(xchi));
		gradphi[0] -= ref_dist[1] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+1));
		gradphi[1] -= ref_dist[0] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+0));
		gradphi[2] -= ref_dist[0] * (1.-ref_dist[1]) * (chi(xchi+2+0) - chi(xchi+0));
		gradphi[0] -= (1.-ref_dist[1]) * ref_dist[2] * (chi(xchi+2+0) - chi(xchi+2));
		gradphi[1] -= (1.-ref_dist[0]) * ref_dist[2] * (chi(xchi+2+1) - chi(xchi+2));
		gradphi[2] -= (1.-ref_dist[0]) * ref_dist[1] * (chi(xchi+2+1) - chi(xchi+1));
		gradphi[0] -= ref_dist[1] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+1));
		gradphi[1] -= ref_dist[0] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+0));
		gradphi[2] -= ref_dist[0] * ref_dist[1] * (chi(xchi+2+1+0) - chi(xchi+1+0));
	}

	e2 = sqrt(e2);

	if (nfield >= 3 && fields[2] != NULL)
	{
		pgradB[0] = ((1.-ref_dist[2]) * (Bi(xB+0,1) - Bi(xB,1)) + ref_dist[2] * (Bi(xB+2+0,1) - Bi(xB+2,1))) * (*part).vel[1];
		pgradB[0] += ((1.-ref_dist[1]) * (Bi(xB+0,2) - Bi(xB,2)) + ref_dist[1] * (Bi(xB+1+0,2) - Bi(xB+1,2))) * (*part).vel[2];
		pgradB[0] += (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((ref_dist[0]-1.) * Bi(xB-0,0) + (1.-2.*ref_dist[0]) * Bi(xB,0) + ref_dist[0] * Bi(xB+0,0)) * (*part).vel[0];
		pgradB[0] += ref_dist[1] * (1.-ref_dist[2]) * ((ref_dist[0]-1.) * Bi(xB+1-0,0) + (1.-2.*ref_dist[0]) * Bi(xB+1,0) + ref_dist[0] * Bi(xB+1+0,0)) * (*part).vel[0];
		pgradB[0] += (1.-ref_dist[1]) * ref_dist[2] * ((ref_dist[0]-1.) * Bi(xB+2-0,0) + (1.-2.*ref_dist[0]) * Bi(xB+2,0) + ref_dist[0] * Bi(xB+2+0,0)) * (*part).vel[0];
		pgradB[0] += ref_dist[1] * ref_dist[2] * ((ref_dist[0]-1.) * Bi(xB+2+1-0,0) + (1.-2.*ref_dist[0]) * Bi(xB+2+1,0) + ref_dist[0] * Bi(xB+2+1+0,0)) * (*part).vel[0];

		pgradB[1] = ((1.-ref_dist[0]) * (Bi(xB+1,2) - Bi(xB,2)) + ref_dist[0] * (Bi(xB+1+0,2) - Bi(xB+0,2))) * (*part).vel[2];
		pgradB[1] += ((1.-ref_dist[2]) * (Bi(xB+1,0) - Bi(xB,0)) + ref_dist[2] * (Bi(xB+1+2,0) - Bi(xB+2,0))) * (*part).vel[0];
		pgradB[1] += (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((ref_dist[1]-1.) * Bi(xB-1,1) + (1.-2.*ref_dist[1]) * Bi(xB,1) + ref_dist[1] * Bi(xB+1,1)) * (*part).vel[1];
		pgradB[1] += ref_dist[0] * (1.-ref_dist[2]) * ((ref_dist[1]-1.) * Bi(xB+0-1,1) + (1.-2.*ref_dist[1]) * Bi(xB+0,1) + ref_dist[1] * Bi(xB+0+1,1)) * (*part).vel[1];
		pgradB[1] += (1.-ref_dist[0]) * ref_dist[2] * ((ref_dist[1]-1.) * Bi(xB+2-1,1) + (1.-2.*ref_dist[1]) * Bi(xB+2,1) + ref_dist[1] * Bi(xB+2+1,1)) * (*part).vel[1];
		pgradB[1] += ref_dist[0] * ref_dist[2] * ((ref_dist[1]-1.) * Bi(xB+2+0-1,1) + (1.-2.*ref_dist[1]) * Bi(xB+2+0,1) + ref_dist[1] * Bi(xB+2+0+1,1)) * (*part).vel[1];

		pgradB[2] = ((1.-ref_dist[1]) * (Bi(xB+2,0) - Bi(xB,0)) + ref_dist[1] * (Bi(xB+2+1,0) - Bi(xB+1,0))) * (*part).vel[0];
		pgradB[2] += ((1.-ref_dist[0]) * (Bi(xB+2,1) - Bi(xB,1)) + ref_dist[0] * (Bi(xB+2+0,1) - Bi(xB+0,1))) * (*part).vel[1];
		pgradB[2] += (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((ref_dist[2]-1.) * Bi(xB-2,2) + (1.-2.*ref_dist[2]) * Bi(xB,2) + ref_dist[2] * Bi(xB+2,2)) * (*part).vel[2];
		pgradB[2] += ref_dist[0] * (1.-ref_dist[1]) * ((ref_dist[2]-1.) * Bi(xB+0-2,2) + (1.-2.*ref_dist[2]) * Bi(xB+0,2) + ref_dist[2] * Bi(xB+0+2,2)) * (*part).vel[2];
		pgradB[2] += (1.-ref_dist[0]) * ref_dist[1] * ((ref_dist[2]-1.) * Bi(xB+1-2,2) + (1.-2.*ref_dist[2]) * Bi(xB+1,2) + ref_dist[2] * Bi(xB+2+1,2)) * (*part).vel[2];
		pgradB[2] += ref_dist[0] * ref_dist[1] * ((ref_dist[2]-1.) * Bi(xB+1+0-2,2) + (1.-2.*ref_dist[2]) * Bi(xB+1+0,2) + ref_dist[2] * Bi(xB+1+0+2,2)) * (*part).vel[2];

		gradphi[0] += pgradB[0] / params[1] / e2;
		gradphi[1] += pgradB[1] / params[1] / e2;
		gradphi[2] += pgradB[2] / params[1] / e2;
	}

	v2 = 0.;
	for (int i = 0; i < 3; i++)
	{
		(*part).vel[i] -= dtau * e2 * gradphi[i] / dx;
		v2 += (*part).vel[i] * (*part).vel[i];
	}

	return v2 / params[0] / params[0];

#undef phi
#undef chi
#undef Bi
#undef xphi
#undef xchi
#undef xB
}


//////////////////////////
// update_q_Newton
//////////////////////////
// Description:
//   Update momentum method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = psi
//              fields[1] = chi
//   Sites      array of Sites on the respective lattices
//   nfield     number of fields (should be 1)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns: squared velocity of particle after update
//
//////////////////////////

Real update_q_Newton(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * Sites, int nfield, double * params, double * outputs, int noutputs)
{
#define psi (*fields[0])
#define xpsi (Sites[0])
#define chi (*fields[1])
#define xchi (Sites[1])

	Real gradpsi[3]={0,0,0};

	gradpsi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * (psi(xpsi+0) - psi(xpsi));
	gradpsi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * (psi(xpsi+1) - psi(xpsi));
	gradpsi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * (psi(xpsi+2) - psi(xpsi));
	gradpsi[0] += ref_dist[1] * (1.-ref_dist[2]) * (psi(xpsi+1+0) - psi(xpsi+1));
	gradpsi[1] += ref_dist[0] * (1.-ref_dist[2]) * (psi(xpsi+1+0) - psi(xpsi+0));
	gradpsi[2] += ref_dist[0] * (1.-ref_dist[1]) * (psi(xpsi+2+0) - psi(xpsi+0));
	gradpsi[0] += (1.-ref_dist[1]) * ref_dist[2] * (psi(xpsi+2+0) - psi(xpsi+2));
	gradpsi[1] += (1.-ref_dist[0]) * ref_dist[2] * (psi(xpsi+2+1) - psi(xpsi+2));
	gradpsi[2] += (1.-ref_dist[0]) * ref_dist[1] * (psi(xpsi+2+1) - psi(xpsi+1));
	gradpsi[0] += ref_dist[1] * ref_dist[2] * (psi(xpsi+2+1+0) - psi(xpsi+2+1));
	gradpsi[1] += ref_dist[0] * ref_dist[2] * (psi(xpsi+2+1+0) - psi(xpsi+2+0));
	gradpsi[2] += ref_dist[0] * ref_dist[1] * (psi(xpsi+2+1+0) - psi(xpsi+1+0));

	if (nfield >= 2 && fields[1] != NULL)
	{
		gradpsi[0] -= (1.-ref_dist[1]) * (1.-ref_dist[2]) * (chi(xchi+0) - chi(xchi));
		gradpsi[1] -= (1.-ref_dist[0]) * (1.-ref_dist[2]) * (chi(xchi+1) - chi(xchi));
		gradpsi[2] -= (1.-ref_dist[0]) * (1.-ref_dist[1]) * (chi(xchi+2) - chi(xchi));
		gradpsi[0] -= ref_dist[1] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+1));
		gradpsi[1] -= ref_dist[0] * (1.-ref_dist[2]) * (chi(xchi+1+0) - chi(xchi+0));
		gradpsi[2] -= ref_dist[0] * (1.-ref_dist[1]) * (chi(xchi+2+0) - chi(xchi+0));
		gradpsi[0] -= (1.-ref_dist[1]) * ref_dist[2] * (chi(xchi+2+0) - chi(xchi+2));
		gradpsi[1] -= (1.-ref_dist[0]) * ref_dist[2] * (chi(xchi+2+1) - chi(xchi+2));
		gradpsi[2] -= (1.-ref_dist[0]) * ref_dist[1] * (chi(xchi+2+1) - chi(xchi+1));
		gradpsi[0] -= ref_dist[1] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+1));
		gradpsi[1] -= ref_dist[0] * ref_dist[2] * (chi(xchi+2+1+0) - chi(xchi+2+0));
		gradpsi[2] -= ref_dist[0] * ref_dist[1] * (chi(xchi+2+1+0) - chi(xchi+1+0));
	}

	Real v2 = 0.;
	for (int i = 0; i < 3; i++)
	{
		(*part).vel[i] -= dtau * params[0] * gradpsi[i] / dx;
		v2 += (*part).vel[i] * (*part).vel[i];
	}

	return v2 / params[0] / params[0];

#undef psi
#undef xpsi
#undef chi
#undef xchi
}


//////////////////////////
// update_pos
//////////////////////////
// Description:
//   Update position method (arbitrary momentum)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that as q^2 << m^2 a^2 the meaning of vel[3]
//   is ~ v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation
//              fields[0] = phi
//              fields[1] = chi
//              fields[2] = Bi
//   Sites      array of Sites on the respective lattices
//   nfield     number of fields
//   params     array of additional parameters
//              params[0] = a
//              params[1] = scaling coefficient for Bi
//   outputs    array of reduction variables
//   noutputs   number of reduction variables
//
// Returns:
//
//////////////////////////

void update_pos(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * Sites, int nfield, double * params, double * outputs, int noutputs)
{
	Real v[3];
	Real v2 = (*part).vel[0] * (*part).vel[0] + (*part).vel[1] * (*part).vel[1] + (*part).vel[2] * (*part).vel[2];
	Real e2 = v2 + params[0] * params[0];
	Real phi = 0;
	Real chi = 0;

	if (nfield >= 1)
	{
		phi = (*fields[0])(Sites[0]) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		phi += (*fields[0])(Sites[0]+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		phi += (*fields[0])(Sites[0]+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
		phi += (*fields[0])(Sites[0]+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
		phi += (*fields[0])(Sites[0]+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
		phi += (*fields[0])(Sites[0]+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
		phi += (*fields[0])(Sites[0]+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
		phi += (*fields[0])(Sites[0]+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
	}

	if (nfield >= 2)
	{
		chi = (*fields[1])(Sites[1]) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		chi += (*fields[1])(Sites[1]+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		chi += (*fields[1])(Sites[1]+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
		chi += (*fields[1])(Sites[1]+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
		chi += (*fields[1])(Sites[1]+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
		chi += (*fields[1])(Sites[1]+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
		chi += (*fields[1])(Sites[1]+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
		chi += (*fields[1])(Sites[1]+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
	}

	v2 = (1. + (3. - v2 / e2) * phi - chi) / sqrt(e2);

	v[0] = (*part).vel[0] * v2;
	v[1] = (*part).vel[1] * v2;
	v[2] = (*part).vel[2] * v2;

	if (nfield >= 3)
	{
		Real b[3];

		b[0] = (*fields[2])(Sites[2], 0) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
		b[1] = (*fields[2])(Sites[2], 1) * (1.-ref_dist[0]) * (1.-ref_dist[2]);
		b[2] = (*fields[2])(Sites[2], 2) * (1.-ref_dist[0]) * (1.-ref_dist[1]);
		b[1] += (*fields[2])(Sites[2]+0, 1) * ref_dist[0] * (1.-ref_dist[2]);
		b[2] += (*fields[2])(Sites[2]+0, 2) * ref_dist[0] * (1.-ref_dist[1]);
		b[0] += (*fields[2])(Sites[2]+1, 0) * ref_dist[1] * (1.-ref_dist[2]);
		b[2] += (*fields[2])(Sites[2]+1, 2) * (1.-ref_dist[0]) * ref_dist[1];
		b[0] += (*fields[2])(Sites[2]+2, 0) * (1.-ref_dist[1]) * ref_dist[2];
		b[1] += (*fields[2])(Sites[2]+2, 1) * (1.-ref_dist[0]) * ref_dist[2];
		b[1] += (*fields[2])(Sites[2]+2+0, 1) * ref_dist[0] * ref_dist[2];
		b[0] += (*fields[2])(Sites[2]+2+1, 0) * ref_dist[1] * ref_dist[2];
		b[2] += (*fields[2])(Sites[2]+1+0, 2) * ref_dist[0] * ref_dist[1];

		for (int i = 0; i < 3; i++) (*part).pos[i] += dtau * (v[i] + b[i] / params[1]);
	}
	else
	{
		for (int i = 0; i < 3 ; i++) (*part).pos[i] += dtau * v[i];
	}
}


//////////////////////////
// update_pos_Newton
//////////////////////////
// Description:
//   Update position method (Newtonian version)
//   Note that vel[3] in the particle structure is used to store q[3] in units
//   of the particle mass, such that the meaning of vel[3] is v*a.
//
// Arguments:
//   dtau       time step
//   dx         lattice unit (unused)
//   part       pointer to particle structure
//   ref_dist   distance vector to reference point (unused)
//   partInfo   global particle properties (unused)
//   fields     array of pointers to fields appearing in geodesic equation (unused)
//   Sites      array of Sites on the respective lattices (unused)
//   nfield     number of fields (unused)
//   params     array of additional parameters
//              params[0] = a
//   outputs    array of reduction variables (unused)
//   noutputs   number of reduction variables (unused)
//
// Returns:
//
//////////////////////////

void update_pos_Newton(double dtau, double dx, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * Sites, int nfield, double * params, double * outputs, int noutputs)
{
	for (int i = 0; i < 3; i++) (*part).pos[i] += dtau * (*part).vel[i] / params[0];
}


//////////////////////////
// projection_T00_project
//////////////////////////
// Description:
//   Particle-mesh projection for T00, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   T00        pointer to target field
//   a          scale factor at projection (needed in order to convert
//              canonical momenta to energies)
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_T00_project(Particles<part, part_info, part_dataType> * pcls, Field<Real> * T00, double a = 1., Field<Real> * phi = NULL, double coeff = 1.)
{
	if (T00->lattice().halo() == 0)
	{
		cout<< "projection_T00_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xField(T00->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridUp[3];
	Real weightScalarGridDown[3];
	Real dx = pcls->res();

	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());
	mass /= a;

	Real e = a, f = 0.;
	Real * q;
	size_t offset_q = offsetof(part,vel);

	Real localCube[8]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111
	Real localCubePhi[8];

	for (int i = 0; i < 8; i++) localCubePhi[i] = 0.0;

	for (xPart.first(), xField.first(); xPart.test(); xPart.next(), xField.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for(int i = 0; i < 3; i++) referPos[i] = xPart.coord(i)*dx;
			for(int i = 0; i < 8; i++) localCube[i] = 0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xField);
				localCubePhi[1] = (*phi)(xField+2);
				localCubePhi[2] = (*phi)(xField+1);
				localCubePhi[3] = (*phi)(xField+1+2);
				localCubePhi[4] = (*phi)(xField+0);
				localCubePhi[5] = (*phi)(xField+0+2);
				localCubePhi[6] = (*phi)(xField+0+1);
				localCubePhi[7] = (*phi)(xField+0+1+2);
			}

			for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i = 0; i < 3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				if (phi != NULL)
				{
					q = (Real*)((char*)&(*it)+offset_q);

					f = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
					e = sqrt(f + a * a);
					f = 3. * e + f / e;
				}

				//000
				localCube[0] += weightScalarGridDown[0]*weightScalarGridDown[1]*weightScalarGridDown[2]*(e+f*localCubePhi[0]);
				//001
				localCube[1] += weightScalarGridDown[0]*weightScalarGridDown[1]*weightScalarGridUp[2]*(e+f*localCubePhi[1]);
				//010
				localCube[2] += weightScalarGridDown[0]*weightScalarGridUp[1]*weightScalarGridDown[2]*(e+f*localCubePhi[2]);
				//011
				localCube[3] += weightScalarGridDown[0]*weightScalarGridUp[1]*weightScalarGridUp[2]*(e+f*localCubePhi[3]);
				//100
				localCube[4] += weightScalarGridUp[0]*weightScalarGridDown[1]*weightScalarGridDown[2]*(e+f*localCubePhi[4]);
				//101
				localCube[5] += weightScalarGridUp[0]*weightScalarGridDown[1]*weightScalarGridUp[2]*(e+f*localCubePhi[5]);
				//110
				localCube[6] += weightScalarGridUp[0]*weightScalarGridUp[1]*weightScalarGridDown[2]*(e+f*localCubePhi[6]);
				//111
				localCube[7] += weightScalarGridUp[0]*weightScalarGridUp[1]*weightScalarGridUp[2]*(e+f*localCubePhi[7]);
			}

			(*T00)(xField)       += localCube[0] * mass;
			(*T00)(xField+2)     += localCube[1] * mass;
			(*T00)(xField+1)     += localCube[2] * mass;
			(*T00)(xField+1+2)   += localCube[3] * mass;
			(*T00)(xField+0)     += localCube[4] * mass;
			(*T00)(xField+0+2)   += localCube[5] * mass;
			(*T00)(xField+0+1)   += localCube[6] * mass;
			(*T00)(xField+0+1+2) += localCube[7] * mass;
		}
	}
}

#define projection_T00_comm scalarProjectionCIC_comm


//////////////////////////
// projection_T0i_project
//////////////////////////
// Description:
//   Particle-mesh projection for T0i, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   T0i        pointer to target field
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_T0i_project(Particles<part,part_info,part_dataType> * pcls, Field<Real> * T0i, Field<Real> * phi = NULL, double coeff = 1.)
{
	if (T0i->lattice().halo() == 0)
	{
		cout<< "projection_T0i_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xT0i(T0i->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridDown[3];
	Real weightScalarGridUp[3];
	Real dx = pcls->res();

	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());

    Real w;
	Real * q;
	size_t offset_q = offsetof(part,vel);

	Real  qi[12];
	Real  localCubePhi[8];

	for (int i = 0; i < 8; i++) localCubePhi[i] = 0;

	for(xPart.first(), xT0i.first(); xPart.test(); xPart.next(), xT0i.next())
	{
		if(pcls->field()(xPart).size!=0)
        {
        	for(int i=0; i<3; i++)
        		referPos[i] = xPart.coord(i)*dx;

            for(int i = 0; i < 12; i++) qi[i]=0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xT0i);
				localCubePhi[1] = (*phi)(xT0i+2);
				localCubePhi[2] = (*phi)(xT0i+1);
				localCubePhi[3] = (*phi)(xT0i+1+2);
				localCubePhi[4] = (*phi)(xT0i+0);
				localCubePhi[5] = (*phi)(xT0i+0+2);
				localCubePhi[6] = (*phi)(xT0i+0+1);
				localCubePhi[7] = (*phi)(xT0i+0+1+2);
			}

			for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i = 0; i < 3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				q = (Real*)((char*)&(*it)+offset_q);

				w = mass * q[0];

				qi[0] +=  w * weightScalarGridDown[1] * weightScalarGridDown[2];
				qi[1] +=  w * weightScalarGridUp[1]   * weightScalarGridDown[2];
				qi[2] +=  w * weightScalarGridDown[1] * weightScalarGridUp[2];
				qi[3] +=  w * weightScalarGridUp[1]   * weightScalarGridUp[2];

				w = mass * q[1];

				qi[4] +=  w * weightScalarGridDown[0] * weightScalarGridDown[2];
				qi[5] +=  w * weightScalarGridUp[0]   * weightScalarGridDown[2];
				qi[6] +=  w * weightScalarGridDown[0] * weightScalarGridUp[2];
				qi[7] +=  w * weightScalarGridUp[0]   * weightScalarGridUp[2];

                w = mass * q[2];

				qi[8] +=  w * weightScalarGridDown[0] * weightScalarGridDown[1];
				qi[9] +=  w * weightScalarGridUp[0]   * weightScalarGridDown[1];
				qi[10]+=  w * weightScalarGridDown[0] * weightScalarGridUp[1];
				qi[11]+=  w * weightScalarGridUp[0]   * weightScalarGridUp[1];
			}

			(*T0i)(xT0i,0) += qi[0] * (1. + localCubePhi[0] + localCubePhi[4]);
			(*T0i)(xT0i,1) += qi[4] * (1. + localCubePhi[0] + localCubePhi[2]);
			(*T0i)(xT0i,2) += qi[8] * (1. + localCubePhi[0] + localCubePhi[1]);

            (*T0i)(xT0i+0,1) += qi[5] * (1. + localCubePhi[4] + localCubePhi[6]);
            (*T0i)(xT0i+0,2) += qi[9] * (1. + localCubePhi[4] + localCubePhi[5]);

            (*T0i)(xT0i+1,0) += qi[1] * (1. + localCubePhi[2] + localCubePhi[6]);
            (*T0i)(xT0i+1,2) += qi[10] * (1. + localCubePhi[2] + localCubePhi[3]);

            (*T0i)(xT0i+2,0) += qi[2] * (1. + localCubePhi[1] + localCubePhi[5]);
            (*T0i)(xT0i+2,1) += qi[6] * (1. + localCubePhi[1] + localCubePhi[3]);

            (*T0i)(xT0i+1+2,0) += qi[3] * (1. + localCubePhi[3] + localCubePhi[7]);
            (*T0i)(xT0i+0+2,1) += qi[7] * (1. + localCubePhi[5] + localCubePhi[7]);
            (*T0i)(xT0i+0+1,2) += qi[11] * (1. + localCubePhi[6] + localCubePhi[7]);
		}
	}
}

#define projection_T0i_comm vectorProjectionCICNGP_comm


//////////////////////////
// projection_Tij_project
//////////////////////////
// Description:
//   Particle-mesh projection for Tij, including geometric corrections
//
// Arguments:
//   pcls       pointer to particle handler
//   Tij        pointer to target field
//   a          scale factor at projection (needed in order to convert
//              canonical momenta to energies)
//   phi        pointer to Bardeen potential which characterizes the
//              geometric corrections (volume distortion); can be set to
//              NULL which will result in no corrections applied
//   coeff      coefficient applied to the projection operation (default 1)
//
// Returns:
//
//////////////////////////

template<typename part, typename part_info, typename part_dataType>
void projection_Tij_project(Particles<part, part_info, part_dataType> * pcls, Field<Real> * Tij, double a = 1., Field<Real> * phi = NULL, double coeff = 1.)
{
	if (Tij->lattice().halo() == 0)
	{
		cout<< "projection_Tij_project: target field needs halo > 0" << endl;
		exit(-1);
	}

	Site xPart(pcls->lattice());
	Site xTij(Tij->lattice());

	typename std::list<part>::iterator it;

	Real referPos[3];
	Real weightScalarGridDown[3];
	Real weightScalarGridUp[3];
	Real dx = pcls->res();

	double mass = coeff / (dx*dx*dx);
	mass *= *(double*)((char*)pcls->parts_info() + pcls->mass_offset());
	mass /= a;

	Real e, f, w;
	Real * q;
	size_t offset_q = offsetof(part,vel);

	Real  tij[6];           // local cube
	Real  tii[24];          // local cube
	Real  localCubePhi[8];

	for (int i = 0; i < 8; i++) localCubePhi[i] = 0;

	for (xPart.first(), xTij.first(); xPart.test(); xPart.next(), xTij.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for (int i = 0; i < 3; i++)
				referPos[i] = (double)xPart.coord(i)*dx;

			for (int i = 0; i < 6; i++)  tij[i]=0.0;
			for (int i = 0; i < 24; i++) tii[i]=0.0;

			if (phi != NULL)
			{
				localCubePhi[0] = (*phi)(xTij);
				localCubePhi[1] = (*phi)(xTij+2);
				localCubePhi[2] = (*phi)(xTij+1);
				localCubePhi[3] = (*phi)(xTij+1+2);
				localCubePhi[4] = (*phi)(xTij+0);
				localCubePhi[5] = (*phi)(xTij+0+2);
				localCubePhi[6] = (*phi)(xTij+0+1);
				localCubePhi[7] = (*phi)(xTij+0+1+2);
			}

			for (it = (pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				for (int i = 0; i < 3; i++)
				{
					weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / dx;
					weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
				}

				q = (Real*)((char*)&(*it)+offset_q);
				f = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
				e = sqrt(f + a * a);
				f = 4. + a * a / (f + a * a);

				// diagonal components
				for (int i = 0; i < 3; i++)
				{
					w = mass * q[i] * q[i] / e;
					//000
					tii[0+i*8] += w * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridDown[2] * (1. + f * localCubePhi[0]);
					//001
					tii[1+i*8] += w * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridUp[2]   * (1. + f * localCubePhi[1]);
					//010
					tii[2+i*8] += w * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridDown[2] * (1. + f * localCubePhi[2]);
					//011
					tii[3+i*8] += w * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridUp[2]   * (1. + f * localCubePhi[3]);
					//100
					tii[4+i*8] += w * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridDown[2] * (1. + f * localCubePhi[4]);
					//101
					tii[5+i*8] += w * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridUp[2]   * (1. + f * localCubePhi[5]);
					//110
					tii[6+i*8] += w * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridDown[2] * (1. + f * localCubePhi[6]);
					//111
					tii[7+i*8] += w * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridUp[2]   * (1. + f * localCubePhi[7]);
				}

				w = mass * q[0] * q[1] / e;
				tij[0] +=  w * weightScalarGridDown[2] * (1. + f * 0.25 * (localCubePhi[0] + localCubePhi[2] + localCubePhi[4] + localCubePhi[6]));
				tij[1] +=  w * weightScalarGridUp[2] * (1. + f * 0.25 * (localCubePhi[1] + localCubePhi[3] + localCubePhi[5] + localCubePhi[7]));

				w = mass * q[0] * q[2] / e;
				tij[2] +=  w * weightScalarGridDown[1] * (1. + f * 0.25 * (localCubePhi[0] + localCubePhi[1] + localCubePhi[4] + localCubePhi[5]));
				tij[3] +=  w * weightScalarGridUp[1] * (1. + f * 0.25 * (localCubePhi[2] + localCubePhi[3] + localCubePhi[6] + localCubePhi[7]));

				w = mass * q[1] * q[2] / e;
				tij[4] +=  w * weightScalarGridDown[0] * (1. + f * 0.25 * (localCubePhi[0] + localCubePhi[1] + localCubePhi[2] + localCubePhi[3]));
				tij[5] +=  w * weightScalarGridUp[0] * (1. + f * 0.25 * (localCubePhi[4] + localCubePhi[5] + localCubePhi[6] + localCubePhi[7]));

			}


			for (int i = 0; i < 3; i++) (*Tij)(xTij,i,i) += tii[8*i];
			(*Tij)(xTij,0,1) += tij[0];
			(*Tij)(xTij,0,2) += tij[2];
			(*Tij)(xTij,1,2) += tij[4];

			for (int i = 0; i < 3; i++) (*Tij)(xTij+0,i,i) += tii[4+8*i];
			(*Tij)(xTij+0,1,2) += tij[5];

			for (int i = 0; i < 3; i++) (*Tij)(xTij+1,i,i) += tii[2+8*i];
			(*Tij)(xTij+1,0,2) += tij[3];

			for (int i = 0; i < 3; i++) (*Tij)(xTij+2,i,i) += tii[1+8*i];
			(*Tij)(xTij+2,0,1) += tij[1];

			for (int i = 0; i < 3; i++) (*Tij)(xTij+0+1,i,i) += tii[6+8*i];
			for (int i = 0; i < 3; i++) (*Tij)(xTij+0+2,i,i) += tii[5+8*i];
			for (int i = 0; i < 3; i++) (*Tij)(xTij+1+2,i,i) += tii[3+8*i];
			for (int i = 0; i < 3; i++) (*Tij)(xTij+0+1+2,i,i) += tii[7+8*i];
		}
	}
}

#ifndef projection_Tij_comm
#define projection_Tij_comm symtensorProjectionCICNGP_comm
#endif

#endif
