//////////////////////////
// Copyright (c) 2015-2016 Julian Adamek (Université de Genève)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////

//////////////////////////
// main.cpp
//////////////////////////
//
// main control sequence of Geneva N-body code with evolution of metric perturbations (gevolution)
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: December 2016
//
//////////////////////////


#include <iomanip>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX			// due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "class_tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"

// f(R) headers
#include "fR_tools.hpp"
#include "background_fR.hpp"

#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "tools.hpp"
#include "output.hpp"
#include "hibernation.hpp"
#include "relaxation.hpp"

using namespace std;
using namespace LATfield2;

int main(int argc, char **argv)
{
	Site ref_site;

	#ifdef BENCHMARK
	//benchmarking variables

	double ref_time, ref2_time, cycle_start_time;


	double initialization_time;
	double run_time;
	double cycle_time = 0;
	double projection_time = 0;
	double snapshot_output_time = 0;
	double spectra_output_time = 0;
	double gravity_solver_time = 0;
	double fft_time = 0;
	int fft_count = 0;
	double update_q_time = 0;
	int update_q_count = 0;
	double moveParts_time = 0;
	int moveParts_count = 0;

	#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;

	int i, j, g, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numsteps_bg, numspecies;
	int temr; // TODO REMOVE AFTER DEBUGGING
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, dx, tau, a, fourpiG, tau_Lambda, tmp, start_time;
	double dtau_print, tau_print;
	double error; // TODO Remove after debugging
	double maxvel[MAX_PCL_SPECIES];
	char filename[2*PARAM_MAX_LENGTH+24];
	string bgfilename;
	int wid;
	string h5filename;
	char * settingsfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	gadget2_header hdr;
	Real T00_hom;
	Real T00_hom_rescaled_a3;
	Real Tii_hom;
	Real phi_hom;
	Real max_fRR;
	std::ofstream bgoutfile;

	// f(R) background description
	double Hubble;
	double Rbar;
	double dot_Rbar;
	double fbar;
	double fRbar;
	double dot_fRbar;
	double fRRbar;
	double coeffm2;

	bool do_I_check;

	#ifndef H5_DEBUG
	H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
	#endif

	for(i=1; i<argc; ++i)
	{
		if(argv[i][0] != '-')
		{
			continue;
		}
		switch(argv[i][1])
		{
			case 's':
			settingsfile = argv[++i]; //settings file name
			break;
			case 'n':
			n = atoi(argv[++i]); //size of the dim 1 of the processor grid
			break;
			case 'm':
			m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
			break;
			case 'i':
			// #ifndef EXTERNAL_IO
			// 				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server" << endl;
			// 				exit(-1000);
			// #endif
			io_size = atoi(argv[++i]);
			break;
			case 'g':
			// #ifndef EXTERNAL_IO
			// 	cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server" << endl;
			// 	exit(-1000);
			// #endif
			io_group_size = atoi(argv[++i]);
		}
	}

	string settingsfile_bin(settingsfile);
	settingsfile_bin += ".bin";

#ifndef EXTERNAL_IO
	parallel.initialize(n, m);
#else
	parallel.initialize(n, m, io_size, io_group_size);
	if(parallel.isIO())
	{
		ioserver.start();
	}
	else
#endif
	{
		COUT
		<< "   __  ___                 _        _    _            " << endl
		<< "  / _|| _ \\ ___ __ __ ___ | | _  _ | |_ (_) ___  _ _  " << endl
		<< " |  _||   // -_)\\ V // _ \\| || || ||  _|| |/ _ \\| ' \\ " << endl
		<< " |_|  |_|_\\\\___| \\_/ \\___/|_| \\_,_| \\__||_|\\___/|_||_| version 1.0 running on " << n*m << " cores." << endl;

		COUT << COLORTEXT_RESET << endl;

		if(settingsfile == NULL)
		{
			COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
			parallel.abortForce();
		}

		COUT << " initializing..." << endl;
		start_time = MPI_Wtime();

		numparam = loadParameterFile(settingsfile, params);
		usedparams = parseMetadata(params, numparam, sim, cosmo, ic);
		COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;
		sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
		saveParameterFile(filename, params, numparam);
		free(params);

		fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
		a = 1.;
		dx = 1.0 / (double) sim.numpts;
		//TODO check particleHorizon
		tau = particleHorizon(a, fourpiG, cosmo);

		bgfilename.reserve(2*PARAM_MAX_LENGTH + 24);
		bgfilename.assign(sim.output_path);
		bgfilename += sim.basename_generic;
		bgfilename += "_background.dat";

		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			if(ic.generator != ICGEN_READ_FROM_DISK) fR_details(cosmo, &sim, fourpiG);

			Rbar = R_initial_fR(a, fourpiG, cosmo);
			fbar = f(Rbar, sim, 100);
			fRbar = fR(Rbar, sim, 101);
		}
		else
		{
			Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
			fbar = fRbar = fRRbar = 0.;
		}

		//=========================================== Full evolution -- Background + Perturbations ===========================================//
		//========= TODO: A lot of stuff unnecessary (background, for instance), remove all of it
		//
		// #ifdef HAVE_CLASS
		// 		background class_background;
		// 		perturbs class_perturbs;
		// 		spectra class_spectra;
		// #endif

		h5filename.reserve(2*PARAM_MAX_LENGTH);
		h5filename.assign(sim.output_path);
		h5filename += sim.basename_generic;
		h5filename += "_";
		h5filename += sim.basename_snapshot;

		box[0] = sim.numpts;
		box[1] = sim.numpts;
		box[2] = sim.numpts;

		Lattice lat(3,box,2);
		Lattice latFT;
		latFT.initializeRealFFT(lat,0);

		MultiGrid mg_engine;
		mg_engine.initialize(&lat, sim.multigrid_n_grids, 4);

		if(sim.multigrid_or_not)
		{
			sim.multigrid_n_grids = mg_engine.nl();
			COUT << " Initialized multigrid with " << sim.multigrid_n_grids << " layers." << endl;
			COUT << " Dimensions of coarsest grid: " << sim.numpts/pow(2,sim.multigrid_n_grids-1) << " x " << sim.numpts/pow(2,sim.multigrid_n_grids-1) << " x " << sim.numpts/pow(2,sim.multigrid_n_grids-1) << endl;
		}

		Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
		Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_b;
		Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES-2];
		Field<Real> * update_cdm_fields[3];
		Field<Real> * update_b_fields[3];
		Field<Real> * update_ncdm_fields[3];
		double f_params[5];

		Field<Real> phi;
		Field<Real> chi;
		Field<Real> Bi;
		Field<Real> Bi_old;
		Field<Real> source;
		Field<Real> Sij;
		Field<Cplx> scalarFT;
		Field<Cplx> SijFT;
		Field<Cplx> BiFT;
		Field<Cplx> Bi_oldFT;

		//additional f(R) fields
		Field<Real> xi;
		Field<Real> xi_old;
		Field<Real> zeta;
		Field<Real> deltaR;
		Field<Real> dot_deltaR;
		Field<Real> deltaR_prev;
		Field<Real> eightpiG_deltaT;
		Field<Real> phi_dot;
		Field<Real> chi_dot;
		Field<Real> phi_ddot;
		Field<Real> phi_effective;
		Field<Real> xi_dot;
		Field<Real> laplace_xi;
		Field<Real> rhs;
		Field<Real> residual;
		Field<Real> err;

		// MultiFields
		MultiField<Real> * mg_deltaR;
		MultiField<Real> * mg_eightpiG_deltaT;
		MultiField<Real> * mg_phi;
		MultiField<Real> * mg_xi;
		MultiField<Real> * mg_xi_old;
		MultiField<Real> * mg_laplace_xi;
		MultiField<Real> * mg_rhs;
		MultiField<Real> * mg_residual;
		MultiField<Real> * mg_err;

		PlanFFT<Cplx> plan_source;
		PlanFFT<Cplx> plan_phi;
		PlanFFT<Cplx> plan_phi_dot;
		PlanFFT<Cplx> plan_phi_ddot;
		PlanFFT<Cplx> plan_phi_effective;
		PlanFFT<Cplx> plan_chi;
		PlanFFT<Cplx> plan_zeta;
		PlanFFT<Cplx> plan_deltaR; // TODO: needed only to output power spectra
		PlanFFT<Cplx> plan_eightpiG_deltaT;
		PlanFFT<Cplx> plan_xi;
		PlanFFT<Cplx> plan_laplace_xi;
		PlanFFT<Cplx> plan_Bi;
		PlanFFT<Cplx> plan_Bi_old;
		PlanFFT<Cplx> plan_Sij;

		// #ifdef CHECK_B
		// 		PlanFFT<Cplx> plan_Bi_check;
		// #endif

		// TODO: For the moment, initialize and alloc everything. After debugging, only stuff that is needed
		scalarFT.initialize(latFT,1);
		BiFT.initialize(latFT,3);
		SijFT.initialize(latFT,3,3,symmetric);

		phi.initialize(lat,1);
		plan_phi.initialize(&phi, &scalarFT);
		phi_dot.initialize(lat,1);
		plan_phi_dot.initialize(&phi_dot, &scalarFT);
		phi_ddot.initialize(lat,1);
		plan_phi_ddot.initialize(&phi_ddot, &scalarFT);
		phi_effective.initialize(lat,1);
		plan_phi_effective.initialize(&phi_effective, &scalarFT);

		chi.initialize(lat,1);
		plan_chi.initialize(&chi, &scalarFT);
		chi_dot.initialize(lat,1);
		chi_dot.alloc();

		xi.initialize(lat,1);
		plan_xi.initialize(&xi, &scalarFT);
		xi_old.initialize(lat,1);
		xi_old.alloc();
		xi_dot.initialize(lat,1);
		xi_dot.alloc();
		laplace_xi.initialize(lat,1);
		plan_laplace_xi.initialize(&laplace_xi, &scalarFT);

		deltaR.initialize(lat,1);
		plan_deltaR.initialize(&deltaR, &scalarFT);
		deltaR_prev.initialize(lat,1);
		deltaR_prev.alloc();
		dot_deltaR.initialize(lat, 1);
		dot_deltaR.alloc();
		zeta.initialize(lat,1);
		plan_zeta.initialize(&zeta, &scalarFT);

		source.initialize(lat,1);
		plan_source.initialize(&source, &scalarFT);
		eightpiG_deltaT.initialize(lat,1);
		plan_eightpiG_deltaT.initialize(&eightpiG_deltaT, &scalarFT);

		residual.initialize(lat, 1);
		residual.alloc();
		err.initialize(lat, 1);
		err.alloc();
		rhs.initialize(lat, 1);
		rhs.alloc();

		Bi.initialize(lat,3);
		plan_Bi.initialize(&Bi, &BiFT);
		Bi_old.initialize(lat,3);
		plan_Bi_old.initialize(&Bi_old, &BiFT);

		Sij.initialize(lat,3,3,symmetric);
		plan_Sij.initialize(&Sij, &SijFT);

		// Initialize MultiGrid fields
		mg_engine.initialize_Field(&phi, mg_phi);
		mg_engine.initialize_Field(&xi, mg_xi);
		mg_engine.initialize_Field(&xi_old, mg_xi_old);
		mg_engine.initialize_Field(&laplace_xi, mg_laplace_xi);
		mg_engine.initialize_Field(&deltaR, mg_deltaR);
		mg_engine.initialize_Field(&eightpiG_deltaT, mg_eightpiG_deltaT);
		mg_engine.initialize_Field(&residual, mg_residual);
		mg_engine.initialize_Field(&err, mg_err);
		mg_engine.initialize_Field(&rhs, mg_rhs);

		update_cdm_fields[0] = &phi;
		update_cdm_fields[1] = &chi;
		update_cdm_fields[2] = &Bi;

		update_b_fields[0] = &phi;
		update_b_fields[1] = &chi;
		update_b_fields[2] = &Bi;

		update_ncdm_fields[0] = &phi;
		update_ncdm_fields[1] = &chi;
		update_ncdm_fields[2] = &Bi;

		Site x(lat);
		rKSite kFT(latFT);

		numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
		COUT << " numpts = " << sim.numpts << "   linesize = " << phi.lattice().size(1) << "   fourpiG = " << fourpiG << endl;

		for(i=0; i<3; ++i) // particles may never move farther than to the adjacent domain
		{
			if(lat.sizeLocal(i)-1 < sim.movelimit)
			{
				sim.movelimit = lat.sizeLocal(i)-1;
			}
		}
		parallel.min(sim.movelimit);

		tau_Lambda = -1.0;
		numsteps_bg = 1;

		// Generate Initial Conditions
		if(ic.generator == ICGEN_BASIC)
		{
			generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij); // generates ICs on the fly
		}
		else if(ic.generator == ICGEN_BASIC_BH)
		{
			generateIC_basic_BH(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
		}
// 		else if(ic.generator == ICGEN_READ_FROM_DISK)
// 		{
// 			if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
// 			{
// 				readIC_fR(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar, dot_Rbar, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &xi, &xi_old, &zeta, &deltaR, &deltaR_prev, &dot_deltaR, &eightpiG_deltaT, &phi_dot, &phi_ddot, &xi_dot, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount, settingsfile_bin);
//
// 				fbar = f(Rbar, sim, 400);
// 				fRbar = fR(Rbar, sim, 401);
// 				fRRbar = fRR(Rbar, sim, 402);
//
// 				dot_fRbar = 0.;
// 				dot_Rbar = 0.;
// 			}
// 			else
// 			{
// 				readIC_GR(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
// 				Hubble = Hconf(a, fourpiG, cosmo);
// 				Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
// 			}
// 		}
// #ifdef ICGEN_PREVOLUTION
// 		else if(ic.generator == ICGEN_PREVOLUTION)
// 		{
// 			generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
// 		}
// #endif
// #ifdef ICGEN_FALCONIC
// 		else if(ic.generator == ICGEN_FALCONIC)
// 		{
// 			maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
// 		}
// #endif
		else
		{
			COUT << " error: IC generator not implemented!" << endl;
			parallel.abortForce();
		}

		if(sim.baryon_flag > 1)
		{
			COUT << " error: baryon_flag > 1 after IC generation, something went wrong in IC generator!" << endl;
			parallel.abortForce();
		}

		numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
		parallel.max<double>(maxvel, numspecies);

		// if(sim.relativistic_flag > 0)
		// {
		// 	for(i=0; i<numspecies; ++i)
		// 	{
		// 		maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		// 	}
		// }

		for(i=0; i<6; ++i)
		{
			hdr.npart[i] = 0;
			hdr.npartTotal[i] = 0;
			hdr.mass[i] = 0.;
		}
		hdr.num_files = 1;
		hdr.Omega0 = cosmo.Omega_m;
		hdr.OmegaLambda = cosmo.Omega_Lambda;
		hdr.HubbleParam = cosmo.h;
		hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
		hdr.flag_sfr = 0;
		hdr.flag_cooling = 0;
		hdr.flag_feedback = 0;
		for(i=0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; ++i)
		{
			hdr.fill[i] = 0;
		}

#ifdef BENCHMARK
		initialization_time = MPI_Wtime() - start_time;
		parallel.sum(initialization_time);
		COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
		COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

// #ifdef HAVE_CLASS
// 		if(sim.radiation_flag > 0)
// 		{
// 			initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra);
// 			if(sim.relativistic_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
// 			{
// 				prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
// 				plan_source.execute(FFT_BACKWARD);
// 				for(x.first(); x.test(); x.next())
// 				{
// 					chi(x) += source(x);
// 				}
// 				chi.updateHalo();
// 			}
// 		}
// #endif

		//=========================================== Main loop ===========================================//

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif

		// construct stress-energy tensor -- Writes -a^3 T00 in {source} (note the minus sign!)
		projection_init(&source);

#ifdef HAVE_CLASS
		if(sim.radiation_flag > 0)
		{
			projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
		}
#endif

		// Project T00
		if(sim.relativistic_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if(sim.baryon_flag)
			{
				projection_T00_project(&pcls_b, &source, a, &phi);
			}

			for(i=0; i<cosmo.num_ncdm; ++i)
			{
				if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
				{
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
				}
				else if(sim.radiation_flag == 0)
				{
					tmp = bg_ncdm(a, cosmo, i);
					for(x.first(); x.test(); x.next())
					{
						source(x) += tmp;
					}
				}
			}
			projection_T00_comm(&source);
		}
		else // Newtonian gravity and f(R) Newtonian
		{//n-body gauge
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if(sim.baryon_flag)
			{
				scalarProjectionCIC_project(&pcls_b, &source);
			}

			for(i=0; i<cosmo.num_ncdm; ++i)
			{
				if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
				{
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
				}
			}
			projection_T00_comm(&source);
		}

		// // Project T0i -- writes a^4 T0i in {Bi}
		if(sim.vector_flag == VECTOR_ELLIPTIC)
		{
			copy_vector_field(Bi, Bi_old);

			projection_init(&Bi);
			projection_T0i_project(&pcls_cdm, &Bi, &phi);
			if(sim.baryon_flag)
			{
				projection_T0i_project(&pcls_b, &Bi, &phi);
			}
			for(i=0; i<cosmo.num_ncdm; ++i)
			{
				if(a >= 1. / (sim.z_switch_Bncdm[i] + 1.))
				{
					projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
				}
			}
			projection_T0i_comm(&Bi);
		}

		// Project Tij -- Write a^3 Tij in {Sij}
		projection_init(&Sij);
		projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
		if(sim.baryon_flag)
		{
			projection_Tij_project(&pcls_b, &Sij, a, &phi);
		}
		if(a >= 1. / (sim.z_switch_linearchi + 1.))
		{
			for(i=0; i<cosmo.num_ncdm; ++i)
			{
				projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
			}
		}
		projection_Tij_comm(&Sij);

#ifdef BENCHMARK
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif

		// Computes 8piG*deltaT = 8piG( (T00 + T11 + T22 + T33) - (-rho_backg + 3P_backg) ), to be stored in {eightpiG_deltaT}
		// At the moment: {source} = -a^3*T00, {Sij} = a^3*Tij
		// NB: in Newtonian f(R), consider only delta_rho instead of full T_mn
		compute_eightpiG_deltaT(eightpiG_deltaT, source, Sij, a, -cosmo.Omega_m, fourpiG, sim);

		if(sim.modified_gravity_flag != MODIFIED_GRAVITY_FLAG_FR || sim.gr_curvature) // in GR/Newton, deltaR tracks exactly eightpiG_deltaT
		{
			copy_field(eightpiG_deltaT, deltaR, -1.);
		}

		build_homogeneous_terms(source, phi, Sij, T00_hom, Tii_hom, phi_hom, T00_hom_rescaled_a3, a, numpts3d);

		//==========================================================================================================//
		//========================================= EVOLVE deltaR, zeta, xi ========================================//
		//==========================================================================================================//
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			if(sim.fR_model == FR_MODEL_R2) // First let's work out the R + R^2 case
			{
				coeffm2 = a*a/6./sim.fR_params[0];
				// Writes source term in {zeta}, then to Fourier space, {zeta} -> {scalarFT}
				copy_field(eightpiG_deltaT, zeta, coeffm2);
				plan_zeta.execute(FFT_FORWARD);
				// Solve wave-like equation transformed into Poisson-like then back to real space, {scalarFT} -> {zeta}
				solveModifiedPoissonFT(scalarFT, scalarFT, 1., coeffm2); // TODO CHECK
				plan_zeta.execute(FFT_BACKWARD);
				// deltaR --> deltaR_prev, Writes zeta --> deltaR, deltaR + eightpiG_deltaT --> zeta
				flip_fields(deltaR, deltaR_prev, zeta);
				add_fields(deltaR, eightpiG_deltaT, zeta);
				convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim); // Needs value of deltaR
			}
			else // Generic f(R) model (not R + R^2), quasi-static or Newtonian
			{
				// Initial conditions for the relaxation solver
				prepare_initial_guess_trace_equation(deltaR, eightpiG_deltaT, xi, laplace_xi, rhs, a, dx, Rbar, fRbar, sim);
				// Relaxation solver
				relaxation(mg_phi, mg_xi, mg_xi_old, mg_laplace_xi, mg_deltaR, mg_eightpiG_deltaT, mg_rhs, mg_residual, mg_err, mg_engine, dx, a, 0., numpts3d, Rbar, fbar, fRbar, sim);
				add_fields(deltaR, eightpiG_deltaT, zeta); // Build zeta from deltaR and eightpiG_deltaT
			}
			// TODO: Avoid updating the halo twice (see later)
			// Build laplacian of xi
			build_laplacian(xi, laplace_xi, dx);
		}

		//==========================================================================================================//
		//=============================================== EVOLVE phi ===============================================//
		//==========================================================================================================//
		// GR and relativistic f(R)
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			// Write source term for the 00 equation in {source}
			prepareFTsource_S00_fR_BH<Real>(source, phi, chi, xi, deltaR, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), dx*dx, fourpiG, Rbar, fbar, fRbar, sim);
		}
		else if(sim.relativistic_flag) // GR
		{
			prepareFTsource_BH<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, fourpiG * dx * dx);  // prepare nonlinear source for phi update
		}

#ifdef BENCHMARK
		ref2_time = MPI_Wtime();
#endif

		plan_source.execute(FFT_FORWARD);  // go to k-space, {source} --> FFT_FORWARD --> {scalarFT}

#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif

		// phi update (k-space)
		solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx));

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		// go back to position space: {scalarFT} --> {phi_dot}
		plan_phi.execute(FFT_BACKWARD);

#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif

		//==========================================================================================================//
		//=============================================== EVOLVE chi ===============================================//
		//==========================================================================================================//
		copy_field(chi, chi_dot); // Before evolving chi

		// prepare nonlinear source for additional equations
		prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx);

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif

		plan_Sij.execute(FFT_FORWARD);  // go to k-space

#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count += 6;
#endif

		// #ifdef HAVE_CLASS
		// if(sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
		// {
		// 	prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
		// 	projectFTscalar(SijFT, scalarFT, 1);
		// }
		// else
		// #endif
		{
			projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)
		}

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		plan_chi.execute(FFT_BACKWARD);	 // go back to position space

#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif

		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR) // TODO: check if it should be done only after the 0th time step or not
		{
			add_fields(chi, xi, chi);
		}

		chi.updateHalo();  // communicate halo values


		//==========================================================================================================//
		//=============================================== EVOLVE Bi ================================================//
		//==========================================================================================================//

		if(sim.vector_flag == VECTOR_ELLIPTIC) // solve B using elliptic constraint; TODO: check for f(R) -- should be the same
		{

#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			plan_Bi.execute(FFT_FORWARD);

#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
			projectFTvector(BiFT, BiFT, fourpiG * dx * dx);

			// #ifdef CHECK_B
			// 				evolveFTvector(SijFT, BiFT_check, 0.);
			// #endif
		}
		else // evolve B using vector projection
		{
			evolveFTvector(SijFT, BiFT, 0.);
		}

		if(sim.relativistic_flag)
		{
#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			plan_Bi.execute(FFT_BACKWARD);  // go back to position space

#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count += 3;
#endif
			Bi.updateHalo();  // communicate halo values
		}

#ifdef BENCHMARK
		gravity_solver_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		check_all_fields(phi, xi, laplace_xi, chi, deltaR, eightpiG_deltaT, zeta, phi_effective, phi_ddot, Bi, numpts3d, sim);
		add_fields(xi, fRbar, xi);
		check_field(xi, "fR", numpts3d, sim);
		add_fields(xi, -fRbar, xi);

		// snapshot output
		COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

		writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &phi_ddot, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);

#ifdef BENCHMARK
		snapshot_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// power spectra
		writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &phi_ddot, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_phi_ddot, &plan_phi_effective, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);

#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif

		COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef BENCHMARK
		run_time = MPI_Wtime() - start_time;

// #ifdef HAVE_CLASS
// 		if(sim.radiation_flag > 0)
// 		{
// 			freeCLASSstructures(class_background, class_perturbs, class_spectra);
// 		}
// #endif

		parallel.sum(run_time);
		parallel.sum(cycle_time);
		parallel.sum(projection_time);
		parallel.sum(snapshot_output_time);
		parallel.sum(spectra_output_time);
		parallel.sum(gravity_solver_time);
		parallel.sum(fft_time);
		parallel.sum(update_q_time);
		parallel.sum(moveParts_time);

		COUT
		<< endl << "BENCHMARK" << endl
		<< "total execution time  : "<< hourMinSec(run_time) << endl
		<< "total number of cycles: "<< cycle << endl
		<< "time consumption breakdown:" << endl
		<< "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time << "%." << endl
		<< "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time << "%." << endl
		<< "----------- main loop: components -----------" << endl
		<< "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time << "%." << endl
		<< "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time << "%." << endl
		<< "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time << "%." << endl
		<< "update momenta (count: "<< update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time << "%." << endl
		<< "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time << "%." << endl
		<< "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time << "%." << endl
		<< "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time << "%." << endl;
#endif

#ifdef EXTERNAL_IO
		ioserver.stop();
#endif
	}

	return 0;
}
