//////////////////////////
// output.hpp
//////////////////////////
//
// Output of snapshots and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: December 2016
//
//////////////////////////

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;



//////////////////////////
// writeSnapshots
//////////////////////////
// Description:
//   output of snapshots
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   hdr            Gadget-2 header structure
//   a              scale factor
//   snapcount      snapshot index
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
//
//////////////////////////

void writeSnapshots(metadata & sim, cosmology & cosmo, const double fourpiG, gadget2_header & hdr, const double a, const int snapcount, string h5filename, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * deltaR, Field<Real> * eightpiG_deltaT, Field<Real> * xi, Field<Real> * laplace_xi, Field<Real> * zeta, Field<Real> * phi_ddot, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_deltaR, PlanFFT<Cplx> * plan_eightpiG_deltaT, PlanFFT<Cplx> * plan_xi, PlanFFT<Cplx> * plan_laplace_xi, PlanFFT<Cplx> * plan_zeta, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[3*PARAM_MAX_LENGTH+40];
	char buffer[64];
	int i;
	Site x(phi->lattice());
	Real divB, curlB, divh, traceh, normh;

	sprintf(filename, "%03d", snapcount);

#ifdef EXTERNAL_IO
	while (ioserver.openOstream()== OSTREAM_FAIL);

	if(sim.out_snapshot & MASK_PCLS)
	{
		pcls_cdm->saveHDF5_server_open(h5filename + filename + "_cdm");
		if(sim.baryon_flag)
			pcls_b->saveHDF5_server_open(h5filename + filename + "_b");
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5_server_open(h5filename + filename + buffer);
		}
	}

	if(sim.out_snapshot & MASK_T00)
		source->saveHDF5_server_open(h5filename + filename + "_T00");

	if(sim.out_snapshot & MASK_B)
		Bi->saveHDF5_server_open(h5filename + filename + "_B");

	if(sim.out_snapshot & MASK_PHI)
		phi->saveHDF5_server_open(h5filename + filename + "_phi");

	if(sim.out_snapshot & MASK_XI)
		xi->saveHDF5_server_open(h5filename + filename + "_xi");

	if(sim.out_snapshot & MASK_LAPLACE_XI)
			laplace_xi->saveHDF5_server_open(h5filename + filename + "_laplace_xi");

	if(sim.out_snapshot & MASK_PHI_DDOT)
			phi_ddot->saveHDF5_server_open(h5filename + filename + "_phi_ddot");

	if(sim.out_snapshot & MASK_ZETA)
		zeta->saveHDF5_server_open(h5filename + filename + "_zeta");

	if(sim.out_snapshot & MASK_DELTAR)
		deltaR->saveHDF5_server_open(h5filename + filename + "_deltaR");

	if(sim.out_snapshot & MASK_DELTAT)
		deltaR->saveHDF5_server_open(h5filename + filename + "_eightpiG_deltaT");

	if(sim.out_snapshot & MASK_CHI)
		chi->saveHDF5_server_open(h5filename + filename + "_chi");

	if(sim.out_snapshot & MASK_HIJ)
		Sij->saveHDF5_server_open(h5filename + filename + "_hij");

#ifdef CHECK_B
	if(sim.out_snapshot & MASK_B)
		Bi_check->saveHDF5_server_open(h5filename + filename + "_B_check");
#endif
#endif

	if(sim.out_snapshot & MASK_RBARE || sim.out_snapshot & MASK_POT)
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if(sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		scalarProjectionCIC_comm(source);
	}

	if(sim.out_snapshot & MASK_RBARE)
	{
		if(sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_rhoN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_rhoN.h5");
	}

	if(sim.out_snapshot & MASK_POT)
	{
		plan_source->execute(FFT_FORWARD);
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
		plan_source->execute(FFT_BACKWARD);
		if(sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_psiN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_psiN.h5");
	}

	if(sim.out_snapshot & MASK_T00)
	{
		projection_init(source);
		if(sim.relativistic_flag > 0)
		{
			projection_T00_project(pcls_cdm, source, a, phi);
			if(sim.baryon_flag)
				projection_T00_project(pcls_b, source, a, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T00_project(pcls_ncdm+i, source, a, phi);
		}
		else
		{
			scalarProjectionCIC_project(pcls_cdm, source);
			if(sim.baryon_flag)
				scalarProjectionCIC_project(pcls_b, source);
			for (i = 0; i < cosmo.num_ncdm; i++)
				scalarProjectionCIC_project(pcls_ncdm+i, source);
		}
		projection_T00_comm(source);
#ifdef EXTERNAL_IO
		source->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_T00.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_T00.h5");
#endif
	}

	if(sim.out_snapshot & MASK_B)
	{
		if(sim.relativistic_flag == 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
		}
		for (x.first(); x.test(); x.next())
		{
			(*Bi)(x,0) /= a * a * sim.numpts;
			(*Bi)(x,1) /= a * a * sim.numpts;
			(*Bi)(x,2) /= a * a * sim.numpts;
		}
		Bi->updateHalo();

		computeVectorDiagnostics(*Bi, divB, curlB);
		COUT << " B diagnostics: max |divB| = " << divB << ", max |curlB| = " << curlB << endl;

#ifdef EXTERNAL_IO
		Bi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_B.h5", sim.downgrade_factor);
		else
			Bi->saveHDF5(h5filename + filename + "_B.h5");
#endif

		if(sim.relativistic_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}

	if(sim.out_snapshot & MASK_PHI)
	{
#ifdef EXTERNAL_IO
		phi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			phi->saveHDF5_coarseGrain3D(h5filename + filename + "_phi.h5", sim.downgrade_factor);
		}
		else
		{
			phi->saveHDF5(h5filename + filename + "_phi.h5");
		}
#endif
	}

	if(sim.out_snapshot & MASK_XI)
	{
#ifdef EXTERNAL_IO
		xi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			xi->saveHDF5_coarseGrain3D(h5filename + filename + "_xi.h5", sim.downgrade_factor);
		}
		else
		{
			xi->saveHDF5(h5filename + filename + "_xi.h5");
		}
#endif
	}

	if(sim.out_snapshot & MASK_LAPLACE_XI)
	{
#ifdef EXTERNAL_IO
		laplace_xi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			laplace_xi->saveHDF5_coarseGrain3D(h5filename + filename + "_laplace_xi.h5", sim.downgrade_factor);
		}
		else
		{
			laplace_xi->saveHDF5(h5filename + filename + "_laplace_xi.h5");
		}
#endif
	}

	if(sim.out_snapshot & MASK_ZETA)
	{
#ifdef EXTERNAL_IO
		zeta->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			zeta->saveHDF5_coarseGrain3D(h5filename + filename + "_zeta.h5", sim.downgrade_factor);
		}
		else
		{
			zeta->saveHDF5(h5filename + filename + "_zeta.h5");
		}
#endif
	}

	if(sim.out_snapshot & MASK_PHI_DDOT)
	{
		#ifdef EXTERNAL_IO
		phi_ddot->saveHDF5_server_write(NUMBER_OF_IO_FILES);
		#else
		if(sim.downgrade_factor > 1)
		{
			phi_ddot->saveHDF5_coarseGrain3D(h5filename + filename + "_phi_ddot.h5", sim.downgrade_factor);
		}
		else
		{
			phi_ddot->saveHDF5(h5filename + filename + "_phi_ddot.h5");
		}
		#endif
	}

	if(sim.out_snapshot & MASK_DELTAR)
	{
#ifdef EXTERNAL_IO
		deltaR->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			deltaR->saveHDF5_coarseGrain3D(h5filename + filename + "_deltaR.h5", sim.downgrade_factor);
		}
		else
		{
			deltaR->saveHDF5(h5filename + filename + "_deltaR.h5");
		}
#endif
	}

		if(sim.out_snapshot & MASK_DELTAT)
		{
	#ifdef EXTERNAL_IO
			eightpiG_deltaT->saveHDF5_server_write(NUMBER_OF_IO_FILES);
	#else
			if(sim.downgrade_factor > 1)
			{
				eightpiG_deltaT->saveHDF5_coarseGrain3D(h5filename + filename + "_eightpiG_deltaT.h5", sim.downgrade_factor);
			}
			else
			{
				eightpiG_deltaT->saveHDF5(h5filename + filename + "_eightpiG_deltaT.h5");
			}
	#endif
		}

	if(sim.out_snapshot & MASK_CHI)
	{
#ifdef EXTERNAL_IO
		chi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			chi->saveHDF5_coarseGrain3D(h5filename + filename + "_chi.h5", sim.downgrade_factor);
		}
		else
		{
			chi->saveHDF5(h5filename + filename + "_chi.h5");
		}
#endif
	}

	if(sim.out_snapshot & MASK_HIJ)
	{
		projectFTtensor(*SijFT, *SijFT);
		plan_Sij->execute(FFT_BACKWARD);
		Sij->updateHalo();

		computeTensorDiagnostics(*Sij, divh, traceh, normh);
		COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
		Sij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij.h5", sim.downgrade_factor);
		}
		else
		{
			Sij->saveHDF5(h5filename + filename + "_hij.h5");
		}
#endif
	}

	if(sim.out_snapshot & MASK_TIJ)
	{
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if(sim.baryon_flag)
		{
			projection_Tij_project(pcls_b, Sij, a, phi);
		}
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		}

		projection_Tij_comm(Sij);

		if(sim.downgrade_factor > 1)
		{
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_Tij.h5", sim.downgrade_factor);
		}
		else
		{
			Sij->saveHDF5(h5filename + filename + "_Tij.h5");
		}
	}

	if(sim.out_snapshot & MASK_P)
	{
		projection_init(Bi);
		projection_T0i_project(pcls_cdm, Bi, phi);
		if(sim.baryon_flag)
		{
			projection_T0i_project(pcls_b, Bi, phi);
		}
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			projection_T0i_project(pcls_ncdm+i, Bi, phi);
		}

		projection_T0i_comm(Bi);

		if(sim.downgrade_factor > 1)
		{
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_p.h5", sim.downgrade_factor);
		}
		else
		{
			Bi->saveHDF5(h5filename + filename + "_p.h5");
		}

		if(sim.relativistic_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}

#ifdef CHECK_B
	if(sim.out_snapshot & MASK_B)
	{
		if(sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if(sim.baryon_flag)
			{
				projection_T0i_project(pcls_b, Bi_check, phi);
			}
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			}

			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		plan_Bi_check->execute(FFT_BACKWARD);

		for (x.first(); x.test(); x.next())
		{
			(*Bi_check)(x,0) /= a * a * sim.numpts;
			(*Bi_check)(x,1) /= a * a * sim.numpts;
			(*Bi_check)(x,2) /= a * a * sim.numpts;
		}
#ifdef EXTERNAL_IO
		Bi_check->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if(sim.downgrade_factor > 1)
		{
			Bi_check->saveHDF5_coarseGrain3D(h5filename + filename + "_B_check.h5", sim.downgrade_factor);
		}
		else
		{
			Bi_check->saveHDF5(h5filename + filename + "_B_check.h5");
		}
#endif
	}
#endif

	if(sim.out_snapshot & MASK_GADGET)
	{
		hdr.time = a;
		hdr.redshift = (1./a) - 1.;

		hdr.npart[1] = (unsigned int) (sim.numpcl[0] / sim.tracer_factor[0]);
		hdr.npartTotal[1] = hdr.npart[1];
		if(sim.baryon_flag)
		{
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
		}
		else
		{
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
		}
		pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.tracer_factor[0]);

		if(sim.baryon_flag)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[1] / sim.tracer_factor[1]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for(i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (unsigned int) (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
		}
	}

	if(sim.out_snapshot & MASK_PCLS)
	{
#ifdef EXTERNAL_IO
		pcls_cdm->saveHDF5_server_write();
		if(sim.baryon_flag)
			pcls_b->saveHDF5_server_write();
		for (i = 0; i < cosmo.num_ncdm; i++)
			pcls_ncdm[i].saveHDF5_server_write();
#else
		pcls_cdm->saveHDF5(h5filename + filename + "_cdm", 1);
		if(sim.baryon_flag)
			pcls_b->saveHDF5(h5filename + filename + "_b", 1);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5(h5filename + filename + buffer, 1);
		}
#endif
	}

#ifdef EXTERNAL_IO
	ioserver.closeOstream();
#endif
}


//////////////////////////
// writeSpectra
//////////////////////////
// Description:
//   output of spectra
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   pkcount        spectrum output index
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
//
//////////////////////////

void writeSpectra(
	metadata & sim,
	cosmology & cosmo,
	const double fourpiG,
	const double a,
	const int pkcount,
	const int cycle,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm,
	Field<Real> * phi,
	Field<Real> * deltaR,
	Field<Real> * eightpiG_deltaT,
	Field<Real> * xi,
	Field<Real> * laplace_xi,
	Field<Real> * zeta,
	Field<Real> * phi_ddot,
	Field<Real> * chi,
	Field<Real> * Bi,
	Field<Real> * source,
	Field<Real> * Sij,
	Field<Cplx> * scalarFT,
	Field<Cplx> * BiFT,
	Field<Cplx> * SijFT,
	PlanFFT<Cplx> * plan_phi,
	PlanFFT<Cplx> * plan_deltaR,
	PlanFFT<Cplx> * plan_eightpiG_deltaT,
	PlanFFT<Cplx> * plan_xi,
	PlanFFT<Cplx> * plan_laplace_xi,
	PlanFFT<Cplx> * plan_zeta,
	PlanFFT<Cplx> * plan_phi_ddot,
	PlanFFT<Cplx> * plan_phi_effective,
	PlanFFT<Cplx> * plan_chi,
	PlanFFT<Cplx> * plan_Bi,
	PlanFFT<Cplx> * plan_source,
	PlanFFT<Cplx> * plan_Sij,
	Field<Real> * Bi_check = NULL,
	Field<Cplx> * BiFT_check = NULL,
	PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[3*PARAM_MAX_LENGTH+40];
	char buffer[64];
	int i, j;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
  long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
	Cplx tempk;
	double Omega_ncdm;

	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;

	kbin = (Real *) malloc(sim.numbins * sizeof(Real));
	power = (Real *) malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	occupation = (int *) malloc(sim.numbins * sizeof(int));

	if(sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || sim.out_pk & MASK_POT || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.relativistic_flag == 0))
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if(sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		scalarProjectionCIC_comm(source);
		plan_source->execute(FFT_FORWARD);

		if(sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.relativistic_flag == 0))
		{
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
		}

		if(sim.out_pk & MASK_RBARE)
		{
			sprintf(filename, "%s%s_%s%03d_rhoN.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of rho_N", a, cycle);
		}

		if(sim.out_pk & MASK_DBARE)
		{
			sprintf(filename, "%s%s_%s%03d_deltaN.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta_N", a, cycle);
		}

		if(sim.out_pk & MASK_T00 && sim.relativistic_flag == 0)
		{
			sprintf(filename, "%s%s_%s%03d_T00.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, cycle);
		}

		if(sim.out_pk & MASK_DELTA && sim.relativistic_flag == 0)
		{
			sprintf(filename, "%s%s_%s%03d_delta.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta", a, cycle);
		}

		if(sim.out_pk & MASK_POT)
		{
			solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%s_%s%03d_psiN.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of psi_N", a, cycle);
		}

		if((cosmo.num_ncdm > 0 || sim.baryon_flag) && (sim.out_pk & MASK_DBARE || (sim.out_pk & MASK_DELTA && sim.relativistic_flag == 0)))
		{
			projection_init(source);
			scalarProjectionCIC_project(pcls_cdm, source);
			scalarProjectionCIC_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s_%s%03d_cdm.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta_N for cdm", a, cycle);
			if(sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if(sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				scalarProjectionCIC_project(pcls_b, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s_%s%03d_b.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta_N for baryons", a, cycle);
				if(sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s_%s%03d_cdmxb.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_cdm * cosmo.Omega_b, filename, "cross power spectrum of delta_N for cdm x baryons", a, cycle);
				}
			}
			Omega_ncdm = 0.;
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				scalarProjectionCIC_project(pcls_ncdm+i, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s_%s%03d_ncdm%d.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount, i);
				sprintf(buffer, "power spectrum of delta_N for ncdm %d", i);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[i], filename, buffer, a);
				Omega_ncdm += cosmo.Omega_ncdm[i];
				// store k-space information for cross-spectra using SijFT as temporary array
				if(cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}
			}
			if(cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s_%s%03d_ncdm.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * Omega_ncdm * Omega_ncdm, filename, "power spectrum of delta_N for total ncdm", a, cycle);
			}
			if(cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if(sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							sprintf(filename, "%s%s_%s%03d_ncdm%dx%d.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount, i, j);
							sprintf(buffer, "cross power spectrum of delta_N for ncdm %d x %d", i, j);
							writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[j], filename, buffer, a);
						}
					}
				}
			}
		}
	}

	if(sim.out_pk & MASK_PHI)
	{
		plan_phi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_phi.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a, cycle);
	}

	if(sim.out_pk & MASK_XI)
	{
		plan_xi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_xi.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of xi", a, cycle);
	}

	if(sim.out_pk & MASK_LAPLACE_XI)
	{
		plan_laplace_xi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_laplace_xi.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of laplace_xi", a, cycle);
	}

	if(sim.out_pk & MASK_PHI_DDOT)
	{
		plan_phi_ddot->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_phi_ddot.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi_ddot", a, cycle);
	}

	if(sim.out_pk & MASK_PHI_EFFECTIVE)
	{
		plan_phi_effective->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_phi_effective.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi_effective", a, cycle);
	}

	if(sim.out_pk & MASK_ZETA)
	{
		plan_zeta->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_zeta.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of zeta", a, cycle);
	}

	if(sim.out_pk & MASK_DELTAR)
	{
		plan_deltaR->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_deltaR.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of deltaR", a, cycle);
	}

	if(sim.out_pk & MASK_DELTAT)
	{
		plan_eightpiG_deltaT->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_eightpiG_deltaT.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of eightpiG_deltaT", a, cycle);
	}

	if(sim.out_pk & MASK_CHI)
	{
		plan_chi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_chi.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of chi", a, cycle);
	}

	if(sim.out_pk & MASK_HIJ)
	{
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if(sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		projection_Tij_comm(Sij);

		prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / (double) sim.numpts / (double) sim.numpts / a);
		plan_Sij->execute(FFT_FORWARD);
		projectFTtensor(*SijFT, *SijFT);

		extractPowerSpectrum(*SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_hij.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij", a, cycle);
	}

	if((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.relativistic_flag > 0)
	{
		projection_init(source);
		projection_T00_project(pcls_cdm, source, a, phi);
		if(sim.baryon_flag)
			projection_T00_project(pcls_b, source, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_T00_project(pcls_ncdm+i, source, a, phi);
		projection_T00_comm(source);

		plan_source->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);

		if(sim.out_pk & MASK_T00)
		{
			sprintf(filename, "%s%s_%s%03d_T00.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, cycle);
		}

		if(sim.out_pk & MASK_DELTA)
		{
			sprintf(filename, "%s%s_%s%03d_delta.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)), filename, "power spectrum of delta", a, cycle);
		}

		if(cosmo.num_ncdm > 0 || sim.baryon_flag)
		{
			projection_init(source);
			projection_T00_project(pcls_cdm, source, a, phi);
			projection_T00_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			if(sim.out_pk & MASK_T00)
			{
				sprintf(filename, "%s%s_%s%03d_T00cdm.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for cdm", a, cycle);
			}
			if(sim.out_pk & MASK_DELTA)
			{
				sprintf(filename, "%s%s_%s%03d_deltacdm.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta for cdm", a, cycle);
			}
			if(sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if(sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				projection_T00_project(pcls_b, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if(sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s_%s%03d_T00b.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for baryons", a, cycle);
				}
				if(sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s_%s%03d_deltab.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta for baryons", a, cycle);
				}
				if(sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s_%s%03d_deltacdmxb.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_cdm, filename, "cross power spectrum of delta for cdm x baryons", a, cycle);
				}
			}
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				projection_T00_project(pcls_ncdm+i, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if(sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s_%s%03d_T00ncdm%d.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of T00 for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a);
				}
				if(sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s_%s%03d_deltancdm%d.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of delta for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, i), filename, buffer, a);
				}
				// store k-space information for cross-spectra using SijFT as temporary array
				if(cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}
			}
			if(cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if(sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s_%s%03d_T00ncdm.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for total ncdm", a, cycle);
				}
				if(sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s_%s%03d_deltancdm.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo) * bg_ncdm(a, cosmo), filename, "power spectrum of delta for total ncdm", a, cycle);
				}
			}
			if(cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if(sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							if(sim.out_pk & MASK_T00)
							{
								sprintf(filename, "%s%s_%s%03d_T00ncdm%dx%d.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of T00 for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a);
							}
							if(sim.out_pk & MASK_DELTA)
							{
								sprintf(filename, "%s%s_%s%03d_deltancdm%dx%d.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of delta for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, j), filename, buffer, a);
							}
						}
					}
				}
			}
		}
	}

	if(sim.out_pk & MASK_B)
	{
		extractPowerSpectrum(*BiFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_B.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, cycle);

#ifdef CHECK_B
		if(sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if(sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		extractPowerSpectrum(*BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s_%s%03d_B_check.dat", sim.output_path, sim.basename_generic, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, cycle);
#endif
	}

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}








/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Print background info.
// Simplifies things in f(R) when we need more points in the background file
// than the quasi_static timesteps
// TODO More info here?
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
bool print_background(
	const int cycle,
	std::ofstream &bgoutfile,
	const string bgfilename,
	const double tau,
	const double a,
	const double Hubble,
	const double Rbar
)
{
	if(parallel.rank() == 0)
	{
		bgoutfile.open(bgfilename, std::ofstream::app);
		if(!bgoutfile.is_open())
		{
			cout << " error opening file for background output!" << endl;
			return 1;
		}

		int wid = 6;
		bgoutfile << scientific << setprecision(wid);
		wid += 9;

		bgoutfile
		<< setw(9) << cycle
		<< setw(wid) << tau
		<< setw(wid) << a
		<< setw(wid) << Hubble
		<< setw(wid) << Rbar
		<< endl;

		bgoutfile.close();
	}

	return 0;
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
bool print_background(
	const int cycle,
	std::ofstream &bgoutfile,
	const string bgfilename,
	const double tau,
	const double a,
	const double Hubble,
	const double Rbar,
	const double phi_hom,
	const double T00_hom_rescaled_a3
)
{
	if(parallel.rank() == 0)
	{
		bgoutfile.open(bgfilename, std::ofstream::app);
		if(!bgoutfile.is_open())
		{
			cout << " error opening file for background output!" << endl;
			return 1;
		}

		int wid = 6;
		bgoutfile << scientific << setprecision(wid);
		wid += 9;

		bgoutfile
		<< setw(9) << cycle
		<< setw(wid) << tau
		<< setw(wid) << a
		<< setw(wid) << Hubble
		<< setw(wid) << Rbar
		<< setw(wid) << phi_hom
		<< setw(wid) << -T00_hom_rescaled_a3
		<< endl;

		bgoutfile.close();
	}

	return 0;
}


int output_background_data(
	double tau,
	double a,
	double Hubble,
	double dtau_old,
	double Rbar,
	double dot_Rbar,
	cosmology cosmo,
	double T00_hom,
	double T00_hom_rescaled_a3,
	double fbar,
	double fRbar,
	double fRRbar,
	int modified_gravity_flag)
{
	COUT << "              tau = " << tau << endl;
	COUT << "                z = " << (1./a) - 1. << endl;
	COUT << "                H = " << Hubble << endl;
	COUT << "         dtau_old = " << dtau_old << endl;
	COUT << "             Rbar = " << Rbar << endl;
	COUT << "         dot_Rbar = " << dot_Rbar << endl;
	COUT << " background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
	COUT << " avg/rescaled T00 = " << T00_hom*a*a*a << " / " << T00_hom_rescaled_a3 << endl;
	if(modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
	{
		COUT << "             fbar = " << fbar << endl;
		COUT << "            fRbar = " << fRbar << endl;
		COUT << "           fRRbar = " << fRRbar << endl;
	}
	COUT << endl;

	return 1;
}


///black hole output:
/// function output lines (for scalars)
///

template <typename FType>
void output_lines(Field<FType> & field, double resolution, string filename)
{
	fstream file;
	Site x(field.lattice());

	double line[field.lattice().size(0)];
	double line_y[field.lattice().size(1)];

	double line_dist[field.lattice().size(0)];
	double line_y_dist[field.lattice().size(1)];

	int point_coord[3];
	double part_coord[3];

	int local_desc[2];
	int xofs,zofs;

	int i,j,k;

	bool proc_have = false;




	int md_local_desc[parallel.size()*2];


	for(i=0;i<3;i++)
	{
		point_coord[i] = field.lattice().size(i) / 2;
		part_coord[i] = (double)point_coord[i] + 0.5*resolution;
	}

	xofs = point_coord[0] -  point_coord[1];
	zofs = point_coord[2] -  point_coord[1];

	//creating line x...

	for(i=0;i<field.lattice().size(0);i++)
	{
		line[i] = 0;
		line_dist[i] = 0;
	}

	if(x.setCoord(point_coord) )
	{
		for(i=0;i<field.lattice().size(0);i++)
		{
				x.setCoord(i,point_coord[1],point_coord[2]);
				line[i] = ( field(x) + field(x+1) + field(x+2) + field(x+1+2) )/4.0;
				line_dist[i] =  i*resolution - part_coord[0];
		}

		//this proc write in the file the first line


    file.open( filename.c_str() , std::fstream::out | std::fstream::trunc);
    if( !file.is_open() )
    {
      cout << "cannot open file: " << filename << ", exiting" << endl;
      exit(145);
    }
		file<<line_dist[0];
		for(i=1;i<field.lattice().size(0);i++)file<<","<<line_dist[i];
		file<<endl;
		file<<line[0];
		for(i=1;i<field.lattice().size(0);i++)file<<","<<line[i];
		file<<endl;
		file.close();
	}

	//creating line xy (here index is the y_coord)

	local_desc[0]=-1;
	local_desc[1]=0;
	for(i=0;i<parallel.size();i++)
	{
		md_local_desc[2*i] = -1;
		md_local_desc[2*i+1] = 0;
	}
	for(j=0;j<field.lattice().size(1);j++)
	{
		line_y[j] = 0;
		line_y_dist[j] = 0;
	}


	for(j=0;j<field.lattice().size(1);j++)
	{
		i = j + xofs;
		if(i<0)i += field.lattice().size(0);
		if(i>=field.lattice().size(0)) i -= field.lattice().size(0);

		if(x.setCoord(i,j,point_coord[2]) )
		{
			line_y[j] = (field(x) + field(x+2))/2.0;
			local_desc[1]++;
			if(local_desc[0]==-1)local_desc[0] = j;

			line_y_dist[j] = sqrt( (part_coord[0]- i*resolution)*(part_coord[0]- i*resolution) +
													(part_coord[1]- j*resolution)*(part_coord[1]- j*resolution) );
			if(j<point_coord[1]) line_y_dist[j] *= -1;
		}
	}

	//MPI gather all data..
	MPI_Gather(local_desc, 2*sizeof(int), MPI_BYTE,
			   		 md_local_desc, 2*sizeof(int), MPI_BYTE, 0, parallel.lat_world_comm());

	//send data to root:
	if(local_desc[0]!=-1 && parallel.rank()!=0)
	{
		parallel.send(&line_y[local_desc[0]],local_desc[1],0);
		parallel.send(&line_y_dist[local_desc[0]],local_desc[1],0);
	}

	if(parallel.rank()==0)
	{
		for(i=1;i<parallel.size();i++)
		{
			if(md_local_desc[2*i]!=-1)
			{
				parallel.receive(&line_y[md_local_desc[2*i]],md_local_desc[2*i+1],i);
				parallel.receive(&line_y_dist[md_local_desc[2*i]],md_local_desc[2*i+1],i);
			}
		}

		file.open( filename.c_str() , std::fstream::out | std::fstream::app);
    if( !file.is_open() )
    {
      cout << "cannot open file: " << filename << ", exiting" << endl;
      exit(145);
    }
		file<<line_dist[0];
		for(i=1;i<field.lattice().size(0);i++)file<<","<<line_y_dist[i];
		file<<endl;
		file<<line[0];
		for(i=1;i<field.lattice().size(0);i++)file<<","<<line_y[i];
		file<<endl;
		file.close();

	}


	///XYZ direction
	local_desc[0]=-1;
	local_desc[1]=0;
	for(i=0;i<parallel.size();i++)
	{
		md_local_desc[2*i] = -1;
		md_local_desc[2*i+1] = 0;
	}
	for(j=0;j<field.lattice().size(1);j++)
	{
		line_y[j] = 0;
		line_y_dist[j] = 0;
	}

	for(j=0;j<field.lattice().size(1);j++)
	{
		i = j + xofs;
		if(i<0)i += field.lattice().size(0);
		if(i>=field.lattice().size(0)) i -= field.lattice().size(0);
		k = j + zofs;
		if(k<0)k += field.lattice().size(2);
		if(k>=field.lattice().size(2)) k -= field.lattice().size(2);

		if(x.setCoord(i,j,k) )
		{
			line_y[j] = field(x);
			local_desc[1]++;
			if(local_desc[0]==-1)local_desc[0] = j;

			line_y_dist[j] = sqrt( (part_coord[0]- i*resolution)*(part_coord[0]- i*resolution) +
													(part_coord[1]- j*resolution)*(part_coord[1]- j*resolution) +
												  (part_coord[2]- k*resolution)*(part_coord[2]- k*resolution) );
			if(j<point_coord[1]) line_y_dist[j] *= -1;

		}
	}

	MPI_Gather(local_desc, 2*sizeof(int), MPI_BYTE,
			   		 md_local_desc, 2*sizeof(int), MPI_BYTE, 0, parallel.lat_world_comm());

	//send data to root:
	if(local_desc[0]!=-1 && parallel.rank()!=0)
	{
		parallel.send(&line_y[local_desc[0]],local_desc[1],0);
		parallel.send(&line_y_dist[local_desc[0]],local_desc[1],0);
	}

	if(parallel.rank()==0)
	{
		for(i=1;i<parallel.size();i++)
		{
			if(md_local_desc[2*i]!=-1)
			{
				parallel.receive(&line_y[md_local_desc[2*i]],md_local_desc[2*i+1],i);
				parallel.receive(&line_y_dist[md_local_desc[2*i]],md_local_desc[2*i+1],i);
			}
		}

		file.open( filename.c_str() , std::fstream::out | std::fstream::app);
    if( !file.is_open() )
    {
      cout << "cannot open file: " << filename << ", exiting" << endl;
      exit(145);
    }
		file<<line_dist[0];
		for(i=1;i<field.lattice().size(0);i++)file<<","<<line_y_dist[i];
		file<<endl;
		file<<line[0];
		for(i=1;i<field.lattice().size(0);i++)file<<","<<line_y[i];
		file<<endl;
		file.close();

	}
}



#endif
