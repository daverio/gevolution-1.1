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
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, dx, tau, a, fourpiG, tau_Lambda, tmp, start_time;
	double dtau_print, tau_print;
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
	double fRRbar;
	double coeffm2;

	bool do_I_check;

	#ifndef H5_DEBUG
		H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
	#endif

	for(i=1; i<argc; i++)
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
			#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server" << endl;
				exit(-1000);
			#endif
				io_size = atoi(argv[++i]);
				break;
			case 'g':
			#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server" << endl;
				exit(-1000);
			#endif
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
		a = 1. / (1. + sim.z_in);
		dx = 1.0 / (double) sim.numpts;
		//TODO check particleHorizon
		tau = particleHorizon(a, fourpiG, cosmo);

		if(sim.mg_flag == FLAG_FR) fR_details(cosmo, &sim, fourpiG);

		bgfilename.reserve(2*PARAM_MAX_LENGTH + 24);
		bgfilename.assign(sim.output_path);
		bgfilename += sim.basename_generic;
		bgfilename += "_background.dat";

		if(sim.mg_flag == FLAG_FR)
		{
			Rbar = R_initial_fR(a, fourpiG, cosmo);
			fbar = f(Rbar, sim, 100);
			fRbar = fR(Rbar, sim, 101);
			fRRbar = fRR(Rbar, sim, 102);
			Hubble = H_initial_fR(a, Hconf(a, fourpiG, cosmo), Rbar, fbar, fRbar,	6. * fourpiG * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * fRRbar / a / a / a); // TODO: Check Omega_m term
			dot_Rbar = dot_R_initial_fR(a, Hubble, fourpiG, cosmo, sim);
			if(sim.Cf * dx < sim.steplimit / Hubble)
			{
				dtau = sim.Cf * dx;
			}
			else
			{
				dtau = sim.steplimit / Hubble;
			}

			// If perturbations_are_GR or quasi-static, dtau for fields dictated by usual steplimit (background is still split in smaller timesteps)
			if(!sim.perturbations_are_GR && !sim.quasi_static)
			{
				dtau_osci = sim.fR_epsilon_fields * sqrt(3.0 * fabs(fRRbar))/a;
				if(dtau > dtau_osci)
				{
					dtau = dtau_osci;
				}
			}
		}
		else
		{
			Hubble = Hconf(a, fourpiG, cosmo);
			Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
			if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
			{
				dtau = sim.Cf * dx;
			}
			else
			{
				dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
			}
		}

		dtau_old = dtau_old_2 = dtau_print = 0.;
		tau_print = tau;

		if(sim.mg_flag == FLAG_FR && sim.BACKGROUND_NUMPTS > 0)
		{
			dtau_print = tau * sqrt(sim.z_in + 1.) / sqrt(sim.z_fin + 1.) / sim.BACKGROUND_NUMPTS;
		}

		if(parallel.rank() == 0)
		{
			bgoutfile.open(bgfilename);
			if(!bgoutfile.is_open())
			{
				cout << " error opening file for background output!" << endl;
				return -1;
			}

			wid = 6;
			bgoutfile << scientific << setprecision(wid);
			wid += 9;

			bgoutfile
			<< "# background statistics\n#"
			<< setw(8) 	 << "cycle"
			<< setw(wid) << "tau/boxsize"
			<< setw(wid) << "a"
			<< setw(wid) << "conformal H"
			<< setw(wid) << "R";

			if(!sim.background_only)
			{
				bgoutfile
				<< setw(wid) << "phi(k=0)"
				<< setw(wid) << "T00(k=0)";
			}

			bgoutfile << endl;
			bgoutfile.close();
		}

		//=========================================== Background Only ===========================================//
		if(sim.background_only)
		{
			while(a < 1./(1. + sim.z_fin))
			{
				if(sim.mg_flag == FLAG_FR)
				{
					dtau_old = dtau;
					fbar = f(Rbar, sim, 300);
					fRbar = fR(Rbar, sim, 301);
					fRRbar = fRR(Rbar, sim, 302);
					dtau = sim.fR_epsilon_bg * sqrt(3. * fabs(fRRbar)) / a;
					dtau = rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau, sim);
					tau += dtau; // TODO: CHECK
				}
				else
				{
					dtau = sim.steplimit / Hubble;
					rungekutta4bg(a, fourpiG, cosmo, dtau);  // GR -- evolve background
					Hubble = Hconf(a, fourpiG, cosmo);
					Rbar = Rbar_GR(a, fourpiG, cosmo);
					dot_Rbar = dot_Rbar_GR(a, Hubble, fourpiG, cosmo);
					tau += dtau;
				}

				// Write background data on file
				if(sim.mg_flag != FLAG_FR || tau >= tau_print)
				{
					if(parallel.rank() == 0)
					{
						cout
						<< " cycle " << cycle
						<< "\tz = " << 1./a - 1.
						<< "\ttau = " << tau
						<< "\tdtau = " << dtau
						<< endl;

						if(print_background(cycle, bgoutfile, bgfilename, tau, a, Hubble, Rbar))
						{
							return -1;
						}
					}

					tau_print += dtau_print;
				}

				cycle++;
			}

			// Write last step for f(R)
			if(parallel.rank() == 0)
			{
				bgoutfile.open(bgfilename, std::ofstream::app);
				if(!bgoutfile.is_open())
				{
					cout << " error opening file for background output!" << endl;
					return -1;
				}

				if(sim.mg_flag == FLAG_FR)
				{
					cout
					<< " FINAL STEP:\n"
					<< " cycle " << cycle
					<< "\tz = " << 1./a - 1.
					<< "\ttau = " << tau
					<< "\tdtau = " << dtau
					<< endl;

					bgoutfile
					<< setw(9) << cycle
					<< setw(wid) << tau
					<< setw(wid) << a
					<< setw(wid) << Hubble
					<< setw(wid) << Rbar
					<< endl;
				}

				bgoutfile.close();
				cout << "\n Background evolution complete.\n";
			}

			return 0;
		}


		//=========================================== Full evolution -- Background + Perturbations ===========================================//

#ifdef HAVE_CLASS
		background class_background;
		perturbs class_perturbs;
		spectra class_spectra;
#endif

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
		Field<Real> source;
		Field<Real> Sij;
		Field<Cplx> scalarFT;
		Field<Cplx> SijFT;
		Field<Cplx> BiFT;

		//additional f(R) fields
		Field<Real> xi;
		Field<Real> xi_prev;
		Field<Real> zeta;
		Field<Real> deltaR;
		Field<Real> dot_deltaR;
		Field<Real> deltaR_prev;
		Field<Real> eightpiG_deltaT;
		Field<Real> phidot;
		Field<Real> lensing;
		Field<Real> xidot;
		Field<Real> laplace_xi;
		Field<Real> rhs;
		Field<Real> residual;
		Field<Real> err;
		Field<Real> debug_field; // TODO: Just for debugging
		Field<Real> debug_field2; // TODO: Just for debugging
		Field<Real> debug_field3; // TODO: Just for debugging

		// MultiFields
		MultiField<Real> * mg_deltaR;
		MultiField<Real> * mg_eightpiG_deltaT;
		MultiField<Real> * mg_xi;
		MultiField<Real> * mg_laplace_xi;
		MultiField<Real> * mg_rhs;
		MultiField<Real> * mg_residual;
		MultiField<Real> * mg_err;

		PlanFFT<Cplx> plan_source;
		PlanFFT<Cplx> plan_phi;
		PlanFFT<Cplx> plan_phidot;
		PlanFFT<Cplx> plan_lensing;
		PlanFFT<Cplx> plan_chi;
		PlanFFT<Cplx> plan_zeta;
		PlanFFT<Cplx> plan_deltaR; // TODO: needed only to output power spectra
		PlanFFT<Cplx> plan_deltaR_prev;
		PlanFFT<Cplx> plan_eightpiG_deltaT;
		PlanFFT<Cplx> plan_xi;
		PlanFFT<Cplx> plan_laplace_xi;
		PlanFFT<Cplx> plan_Bi;
		PlanFFT<Cplx> plan_Sij;

		PlanFFT<Cplx> plan_debug_field; // TODO: Just for debugging
		PlanFFT<Cplx> plan_debug_field2; // TODO: Just for debugging
		PlanFFT<Cplx> plan_debug_field3; // TODO: Just for debugging

#ifdef CHECK_B
		PlanFFT<Cplx> plan_Bi_check;
#endif

		// TODO: For the moment, initialize and alloc everything. After debugging, only stuff that's needed
		phi.initialize(lat,1);
		chi.initialize(lat,1);
		source.initialize(lat,1);
		scalarFT.initialize(latFT,1);
		xi.initialize(lat,1);
		xi_prev.initialize(lat,1);
		xi_prev.alloc();
		xidot.initialize(lat,1);
		xidot.alloc();
		laplace_xi.initialize(lat,1);
		zeta.initialize(lat,1);
		deltaR.initialize(lat,1);
		dot_deltaR.initialize(lat, 1);
		dot_deltaR.alloc();
		residual.initialize(lat, 1);
		residual.alloc();
		err.initialize(lat, 1);
		err.alloc();
		rhs.initialize(lat, 1);
		rhs.alloc();
		debug_field.initialize(lat, 1);
		debug_field2.initialize(lat, 1);
		debug_field3.initialize(lat, 1);
		deltaR_prev.initialize(lat,1);
		eightpiG_deltaT.initialize(lat,1);
		phidot.initialize(lat,1);
		lensing.initialize(lat,1);
		plan_source.initialize(&source, &scalarFT);
		plan_phi.initialize(&phi, &scalarFT);
		plan_phidot.initialize(&phidot, &scalarFT);
		plan_lensing.initialize(&lensing, &scalarFT);
		plan_chi.initialize(&chi, &scalarFT);
		plan_xi.initialize(&xi, &scalarFT);
		plan_laplace_xi.initialize(&laplace_xi, &scalarFT);
		plan_deltaR.initialize(&deltaR, &scalarFT);
		plan_deltaR_prev.initialize(&deltaR_prev, &scalarFT);
		plan_eightpiG_deltaT.initialize(&eightpiG_deltaT, &scalarFT);
		plan_zeta.initialize(&zeta, &scalarFT);
		plan_debug_field.initialize(&debug_field, &scalarFT);
		plan_debug_field2.initialize(&debug_field2, &scalarFT);
		plan_debug_field3.initialize(&debug_field3, &scalarFT);
		Bi.initialize(lat,3);
		BiFT.initialize(latFT,3);
		plan_Bi.initialize(&Bi, &BiFT);
		Sij.initialize(lat,3,3,symmetric);
		SijFT.initialize(latFT,3,3,symmetric);
		plan_Sij.initialize(&Sij, &SijFT);

		// Initialize MultiGrid fields
		mg_engine.initialize_Field(&xi, mg_xi);
		mg_engine.initialize_Field(&laplace_xi, mg_laplace_xi);
		mg_engine.initialize_Field(&deltaR, mg_deltaR);
		mg_engine.initialize_Field(&eightpiG_deltaT, mg_eightpiG_deltaT);
		mg_engine.initialize_Field(&residual, mg_residual);
		mg_engine.initialize_Field(&err, mg_err);
		mg_engine.initialize_Field(&rhs, mg_rhs);

#ifdef CHECK_B
		Field<Real> Bi_check;
		Field<Cplx> BiFT_check;
		Bi_check.initialize(lat,3);
		BiFT_check.initialize(latFT,3);
		plan_Bi_check.initialize(&Bi_check, &BiFT_check);
#endif

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
		COUT << " numpts = " << sim.numpts << "   linesize = " << phi.lattice().size(1) << endl;

		for(i=0; i<3; i++) // particles may never move farther than to the adjacent domain
		{
			if(lat.sizeLocal(i)-1 < sim.movelimit)
			{
				sim.movelimit = lat.sizeLocal(i)-1;
			}
		}
		parallel.min(sim.movelimit);

		tau_Lambda = -1.0;
		numsteps_bg = 1;

		if(sim.mg_flag == FLAG_FR && !sim.lcdm_background) // set timesteps for f(R) background evolution
		{
			dtau_bg = sim.fR_epsilon_bg * sqrt(3. * fRRbar) / a; // This would be the "ideal" timestep to evolve the f(R) background
			if(dtau >= 2. * dtau_bg)// numsteps_bg must be even -- we split the background evolution in two steps of dtau/2 each
			{
				numsteps_bg = (int) (dtau / dtau_bg / 2.);
				numsteps_bg = (int) (2 * numsteps_bg);
				dtau_bg = (double) (dtau / numsteps_bg);
			}
			else dtau_bg = dtau;
		}
		else // keep same dtau
		{
			dtau_bg = dtau;
		}

		// Generate Initial Conditions
		if(ic.generator == ICGEN_BASIC)
		{
			generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij); // generates ICs on the fly
		}
		else if(ic.generator == ICGEN_READ_FROM_DISK)
		{
			if(sim.mg_flag == FLAG_FR)
			{
				readIC_fR(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar, dot_Rbar, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &xi, &xi_prev, &zeta, &deltaR, &deltaR_prev, &dot_deltaR, &eightpiG_deltaT, &phidot, &xidot, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount, settingsfile_bin);

				fbar = f(Rbar, sim, 400);
				fRbar = fR(Rbar, sim, 401);
				fRRbar = fRR(Rbar, sim, 402);
			}
			else
			{
				readIC_GR(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
				Hubble = Hconf(a, fourpiG, cosmo);
				Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
			}
		}
#ifdef ICGEN_PREVOLUTION
		else if(ic.generator == ICGEN_PREVOLUTION)
		{
			generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
		}
#endif
#ifdef ICGEN_FALCONIC
		else if(ic.generator == ICGEN_FALCONIC)
		{
			maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
		}
#endif
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

		if(sim.gr_flag > 0)
		{
			for(i=0; i<numspecies; i++)
			{
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
			}
		}

#ifdef CHECK_B
		if(sim.vector_flag == VECTOR_ELLIPTIC)
		{
			for(kFT.first(); kFT.test(); kFT.next())
			{
				BiFT_check(kFT, 0) = BiFT(kFT, 0);
				BiFT_check(kFT, 1) = BiFT(kFT, 1);
				BiFT_check(kFT, 2) = BiFT(kFT, 2);
			}
		}
#endif

		for(i=0; i<6; i++)
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
		for(i=0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; i++)
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

#ifdef HAVE_CLASS
		if(sim.radiation_flag > 0)
		{
			initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra);
			if(sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
			{
				prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
				plan_source.execute(FFT_BACKWARD);
				for(x.first(); x.test(); x.next())
				{
					chi(x) += source(x);
				}
				chi.updateHalo();
			}
		}
#endif

		//=========================================== Main loop ===========================================//
		while(true)
		{
			// TODO: REMOVE AFTER DEBUGGING
			do_I_check = (bool) (cycle % sim.CYCLE_INFO_INTERVAL == 0 && sim.check_fields);

			if(cycle % sim.CYCLE_INFO_INTERVAL == 0) // output some info
			{
				COUT << endl << "===================================  CYCLE " << cycle << "  ";
				if(cycle < 10)
				{
					COUT << "==================================" << endl;
				}
				else if(cycle < 100)
				{
					COUT << "=================================" << endl;
				}
				else if(cycle < 1000)
				{
					COUT << "================================" << endl;
				}
				else if(cycle < 10000)
				{
					COUT << "===============================" << endl;
				}
			}

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
			if(sim.gr_flag > 0)
			{
				projection_T00_project(&pcls_cdm, &source, a, &phi);
				if(sim.baryon_flag)
				{
					projection_T00_project(&pcls_b, &source, a, &phi);
				}

				for(i=0; i<cosmo.num_ncdm; i++)
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
			else
			{//n-body gauge
				scalarProjectionCIC_project(&pcls_cdm, &source);
				if(sim.baryon_flag)
				{
					scalarProjectionCIC_project(&pcls_b, &source);
				}

				for(i=0; i<cosmo.num_ncdm; i++)
				{
					if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					{
						scalarProjectionCIC_project(pcls_ncdm+i, &source);
					}
				}
				projection_T00_comm(&source);
			}

			// Project T0i -- actually writes a^4 T0i in {Bi}
			if(sim.vector_flag == VECTOR_ELLIPTIC)
			{
				projection_init(&Bi);
				projection_T0i_project(&pcls_cdm, &Bi, &phi);
				if(sim.baryon_flag)
				{
					projection_T0i_project(&pcls_b, &Bi, &phi);
				}
				for(i=0; i<cosmo.num_ncdm; i++)
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
				for(i=0; i<cosmo.num_ncdm; i++)
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
			computeTtrace(eightpiG_deltaT, source, Sij, a*a*a, -cosmo.Omega_m/a/a/a, fourpiG); // TODO: Should this include the explicit Lambda term?
			if(sim.mg_flag != FLAG_FR) // in GR/Newton, deltaR tracks exactly eightpiG_deltaT
			{
				copy_field(eightpiG_deltaT, deltaR, -1.);
			}

			// EVOLUTION FOR GR or f(R) -- NEWTONIAN GRAVITY LATER
			if(sim.gr_flag > 0)
			{
				T00_hom = 0.;
				Tii_hom = 0.;
				phi_hom = 0.;
				for(x.first(); x.test(); x.next())
				{
					T00_hom += -source(x); // source = - a^3 * T00
					Tii_hom += Sij(x,0,0) + Sij(x,1,1) + Sij(x,2,2); // Sij = a^3 Tij
					phi_hom += phi(x);
				}
				parallel.sum<Real>(T00_hom);
				parallel.sum<Real>(Tii_hom);
				parallel.sum<Real>(phi_hom);
				T00_hom /= (Real) numpts3d;
				Tii_hom /= (Real) numpts3d;
				phi_hom /= (Real) numpts3d;
				T00_hom_rescaled_a3 = T00_hom / (1. + 3. * phi_hom);
				T00_hom /= a*a*a; // T00_hom contained a^3 * <T00>, now contains the correct <T00>
				Tii_hom /= a*a*a; // same for Tii

				if(cycle % sim.CYCLE_INFO_INTERVAL == 0) // output some info
				{
					COUT
					<< "              tau = " << tau << endl
					<< "                z = " << (1./a) - 1. << endl
					<< "                H = " << Hubble << endl
					<< "         dtau_old = " << dtau_old << endl
					<< "             Rbar = " << Rbar << endl
					<< "         dot_Rbar = " << dot_Rbar << endl
					<< "             fbar = " << fbar << endl
					<< "            fRbar = " << fRbar << endl
					<< "           fRRbar = " << fRRbar << endl
					<< " background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl
					<< " avg/rescaled T00 = " << T00_hom*a*a*a << " / " << T00_hom_rescaled_a3 << endl
					<< endl;
				}

				//=========================================== EVOLVE deltaR, zeta, xi ===========================================//

				if(sim.mg_flag == FLAG_FR)
				{
					if(sim.perturbations_are_GR) // If perturbations_are_GR, assign deltaR = eightpiG_deltaT, zeta = 0, xi = 0 at all times
					{
						copy_field(eightpiG_deltaT, deltaR, -1.);
					}
					else if(sim.fR_type == FR_TYPE_R2) // First let's work out the R + R^2 case
					{
						coeffm2 = a * a / 6. / sim.fR_params[0];
						if(sim.quasi_static)
						{
							// Writes source term in {zeta}, then to Fourier space, {zeta} -> {scalarFT}
							copy_field(eightpiG_deltaT, zeta, coeffm2);
							plan_zeta.execute(FFT_FORWARD);
							// Solve wave-like equation transformed into Poisson-like then back to real space, {scalarFT} -> {zeta}
							solveModifiedPoissonFT(scalarFT, scalarFT, 1., coeffm2); // TODO CHECK
							plan_zeta.execute(FFT_BACKWARD);
							// deltaR --> deltaR_prev, Writes zeta --> deltaR, deltaR + eightpiG_deltaT --> zeta
							flip_fields(deltaR, deltaR_prev, zeta);
							add_fields(deltaR, eightpiG_deltaT, zeta);
							// Writes old and new xi
							copy_field(xi, xi_prev);
							convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim); // Needs value of deltaR
						}
						else // TODO: LEAPFROG
						{
							if(cycle <= 1) // Initial conditions for non-quasi-static are the same as quasi-static!
							{
								eightpiG_deltaT.updateHalo();
								prepareFTsource_leapfrog_R2(eightpiG_deltaT, zeta, dx);
								plan_zeta.execute(FFT_FORWARD);
								// Solve wave-like equation transformed into Poisson-like, then back to real space, {scalarFT} -> {zeta}
								solveModifiedPoissonFT(scalarFT, scalarFT, 1., coeffm2); // TODO CHECK
								plan_zeta.execute(FFT_BACKWARD);
								flip_fields(deltaR, deltaR_prev);
								subtract_fields(zeta, eightpiG_deltaT, deltaR);
							}
							else // Else builds deltaR_i from deltaR_{i-1} and dot_deltaR_{i-1/2}, for i >= 2
							{
								add_fields(deltaR, 1., dot_deltaR, dtau_old, deltaR);
							}
							// Done at each step:
							copy_field(xi, xi_prev);// TODO: is this needed?
							convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim); // Needs value of deltaR
							// Builds dot_deltaR_{1/2} from deltaR_1 and deltaR_0 (only first step!)
							if(cycle == 1)
							{
								add_fields(deltaR, 1./dtau_old, deltaR_prev, -1./dtau_old, dot_deltaR);
							}
							// Actual leapfrog step (not for initial cycle)
							if(cycle)
							{
								deltaR.updateHalo();
								leapfrog_dotR(deltaR, eightpiG_deltaT, dot_deltaR, dot_deltaR, Hubble, coeffm2, dtau_old, dtau, dx * dx);
							}
						}
					}
					else if(sim.quasi_static) // Generic f(R) model (not R + R^2), quasi-static
					{
						// TODO: Do we need something different for the first timestep?

						copy_field(xi, xi_prev); // Copy old value of xi on xi_prev

						if(sim.relaxation_variable == RELAXATION_VAR_U)
						{
							convert_xi_to_u(xi, xi, fRbar); // Convert back xi -> u
						}

						// Initial conditions for the relaxation solver
						prepare_initial_conditions_trace_equation(deltaR, eightpiG_deltaT, xi, laplace_xi, rhs, a, dx, Rbar, fRbar, sim);
						// Relaxation solver
						relaxation(mg_xi, mg_laplace_xi, mg_deltaR, mg_eightpiG_deltaT, mg_rhs, mg_residual, mg_err, mg_engine, a, dx, numpts3d, Rbar, fbar, fRbar, sim);

						add_fields(deltaR, eightpiG_deltaT, zeta); // Build zeta from deltaR and eightpiG_deltaT

						if(sim.relaxation_variable == RELAXATION_VAR_U)
						{
							convert_u_to_xi(xi, xi, fRbar); // Convert back u -> xi
						}

						sim.multigrid_check_shape = 0; // Just check the first time (at most)
					}
					else // Generic f(R) model (not R + R^2), all time derivatives
					{
						if(cycle < 2) // Initial conditions are deltaR = - eightpiG_deltaT
						{
							copy_field(deltaR, deltaR_prev);
							copy_field(xi, xi_prev);
							copy_field(eightpiG_deltaT, deltaR, -1.);
							convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim);
							if(cycle == 1) // Builds xidot_{1/2} from xi_1 and xi_0
							{
								add_fields(xi, 1./dtau_old, xi_prev, -1./dtau_old, xidot);
							}
						}
						else // Else builds xi_i from xi_{i-1} and xidot_{i-1/2}, for i >= 2
						{
							add_fields(xi, 1., xidot, dtau_old, xi);
						}

						// Compute deltaR from a given xi
						convert_xi_to_deltaR(eightpiG_deltaT, deltaR, xi, Rbar, fRbar, sim);
						// Build zeta from deltaR and eightpiG_deltaT
						add_fields(deltaR, eightpiG_deltaT, zeta);

						if(cycle) // Update xidot
						{
							xi.updateHalo(); // TODO: Avoid updating the halo twice (see later)
							leapfrog_dotxi(xi, zeta, xidot, xidot, Hubble, a*a, dtau_old, dtau, dx*dx);
						}
					}

					// TODO: Avoid updating the halo twice (see later)
					// Build laplacian of xi
					build_laplacian(xi, laplace_xi, dx);
				}

				//=========================================== EVOLVE phi ===========================================//
				if(cycle)
				{
					if(sim.mg_flag == FLAG_FR && !sim.perturbations_are_GR)
					{
						// Write source term for the 00 equation in {source}
						prepareFTsource_S00<Real>(source, phi, chi, xi, xi_prev, deltaR, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), dx*dx, dtau_old, Hubble, a, fourpiG / a, Rbar, fbar, fRbar, sim);
					}
					else
					{
						prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update
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
					if(sim.mg_flag == FLAG_FR && !sim.perturbations_are_GR)
					{
						solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), sim.quasi_static ? 3. * Hubble / dtau_old : 3. * Hubble / dtau_old);
					}
					else  // GR
					{
						solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);
					}

#ifdef BENCHMARK
					ref2_time= MPI_Wtime();
#endif
					// go back to position space: {scalarFT} --> {phidot}
					plan_phidot.execute(FFT_BACKWARD);
					flip_fields(phi, phidot); // Now {phi} contains phi_new, {phidot} contains phi_old
					add_fields(phi, 1./dtau_old, phidot, -1./dtau_old, phidot); // Computes phidot = (phi_new - phi_old)/dtau

#ifdef BENCHMARK
					fft_time += MPI_Wtime() - ref2_time;
					fft_count++;
#endif
				}
				// Finished evolving deltaR, xi, phi in GR or f(R)
			}
			else // Newtonian Evolution
			{
#ifdef BENCHMARK
				ref2_time = MPI_Wtime();
#endif
				plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space

#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif
				solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif

				plan_phidot.execute(FFT_BACKWARD);	 // go back to position space
				flip_fields(phi, phidot); // Now {phi} contains phi_new, {phidot} contains phi_old
				add_fields(phi, 1./dtau_old, phidot, -1./dtau_old, phidot); // Computes phidot = (phi_new - phi_old)/dtau

#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif
			}

			phi.updateHalo();  // communicate halo values
			phidot.updateHalo();

			if(sim.mg_flag == FLAG_FR && sim.perturbations_are_GR)
			{
				convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim);
				build_laplacian(xi, laplace_xi, dx);
			}


			// Background timesteps in f(R)
			if(sim.mg_flag == FLAG_FR && !sim.lcdm_background) // set timesteps for f(R) background evolution
			{
				dtau_bg = sim.fR_epsilon_bg * sqrt(3. * fRRbar) / a; // This would be the "ideal" timestep to evolve the f(R) background
				if(dtau >= 2. * dtau_bg)// numsteps_bg must be even -- we split the background evolution in two steps of dtau/2 each
				{
					numsteps_bg = (int) (dtau / dtau_bg / 2.);
					numsteps_bg = (int) (2 * numsteps_bg);
					dtau_bg = (double) (dtau / numsteps_bg);
				}
				else dtau_bg = dtau;
			}
			else // keep same dtau
			{
				dtau_bg = dtau;
			}

			// record background data in Newton, GR or LCDM -- see later for f(R)
			if(kFT.setCoord(0, 0, 0))
			{
				if(sim.mg_flag != FLAG_FR || sim.lcdm_background || numsteps_bg == 1)
				{
					if(print_background(cycle, bgoutfile, bgfilename, tau, a, Hubble, Rbar, phi_hom, -T00_hom_rescaled_a3)) return -1;
				}
			}

			// prepare nonlinear source for additional equations
			prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);

#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif

			plan_Sij.execute(FFT_FORWARD);  // go to k-space

#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count += 6;
#endif

#ifdef HAVE_CLASS
			if(sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
			{
				prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
				projectFTscalar(SijFT, scalarFT, 1);
			}
			else
#endif
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

			if(sim.mg_flag == FLAG_FR && cycle && !sim.perturbations_are_GR) // TODO: check if it should be done only after the 0th time step or not
			{
				add_fields(chi, xi, chi);
			}

			chi.updateHalo();  // communicate halo values

			if(do_I_check) check_all_fields(phi, xi, chi, deltaR, eightpiG_deltaT, zeta, Bi, numpts3d, sim);

			// Build lensing potential, proportional to phidot + psidot = 2*phidot - chidot
			// TODO Remove after debugging
			subtract_fields(chi, lensing, lensing);
			add_fields(phidot, 2., lensing, -1./dtau_old, lensing);
			// End Remove

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

#ifdef CHECK_B
				evolveFTvector(SijFT, BiFT_check, a * a * dtau_old);
#endif
			}
			else // evolve B using vector projection
			{
				evolveFTvector(SijFT, BiFT, a * a * dtau_old);
			}

			if(sim.gr_flag > 0)
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

			// snapshot output
			if(snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
			{
				COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;


#ifdef CHECK_B
				writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
				writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
				snapcount++;
			}

#ifdef BENCHMARK
			snapshot_output_time += MPI_Wtime() - ref_time;
			ref_time = MPI_Wtime();
#endif

			// power spectra
			if(pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
			{
				COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
				//TODO REMOVE
				copy_field(zeta, debug_field);
				add_fields(phi, 1., xi, -0.5, zeta);
				// END REMOVE

#ifdef CHECK_B
				writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_phidot, &plan_lensing, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
				writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_phidot, &plan_lensing, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
				pkcount++;
				//TODO REMOVE
				copy_field(debug_field, zeta);
				// END REMOVE
			}

#ifdef BENCHMARK
			spectra_output_time += MPI_Wtime() - ref_time;
#endif

			if(pkcount >= sim.num_pk && snapcount >= sim.num_snapshot) break; // simulation complete

			// For when the full background evolution only starts at sim.z_switch_fR_background
			// TODO: Remove this?
			if(1./a - 1. < sim.z_switch_fR_background && sim.mg_flag == FLAG_FR && sim.lcdm_background)
			{
				COUT << endl;
				COUT << " =================================== ";
				COUT << " Beginning full f(R) background evolution at z = " << 1./a - 1.;
				COUT << " ===================================";
				COUT << endl << endl;

				if(sim.fR_type == FR_TYPE_HU_SAWICKI)
				{
					cosmo.Omega_Lambda = 0.;
				}
				sim.lcdm_background = 0;
			}

			// compute number of step subdivisions for particle updates
			numsteps = 1;
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				if(dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				{
					numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
				}
				else
				{
					numsteps_ncdm[i] = 1;
				}

				if(numsteps < numsteps_ncdm[i])
				{
					numsteps = numsteps_ncdm[i];
				}
			}

			if(numsteps > 1 && numsteps % 2 > 0) // if numsteps > 1, make it an even number
			{
				numsteps++;
			}

			for(i=0; i<cosmo.num_ncdm; i++)
			{
				if(numsteps / numsteps_ncdm[i] <= 1)
				{
					numsteps_ncdm[i] = numsteps;
				}
				else if(numsteps_ncdm[i] > 1)
				{
					numsteps_ncdm[i] = numsteps / 2;
				}
			}

			if(cycle % sim.CYCLE_INFO_INTERVAL == 0) // Some outputs
			{
				COUT << "           max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx << ")\n";
				if(sim.baryon_flag)
				{
					COUT << "           baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx << ")\n";
				}

				if(sim.mg_flag == FLAG_FR)
				{
					COUT << "           Hubble * dtau = " <<  Hubble * dtau << endl
					<< "                      dx = " << dx << endl;

					if(!sim.lcdm_background)
					{
						COUT << " numsteps for background = " << numsteps_bg;
					}
				}
				else
				{
					COUT << "           time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau;
				}

				for(i=0; i<cosmo.num_ncdm; i++)
				{
					if(i == 0)
					{
						COUT << endl << " time step subdivision for ncdm species: ";
					}
					COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i+1+sim.baryon_flag] << ")";
					if(i<cosmo.num_ncdm-1)
					{
						COUT << ", ";
					}
				}
				COUT << endl;
			}

			for(j=0; j<numsteps; j++) // particle update
			{
#ifdef BENCHMARK
				ref2_time = MPI_Wtime();
#endif
				f_params[0] = a;
				f_params[1] = a * a * sim.numpts;
				if(j == 0)
				{
					if(sim.gr_flag > 0)
					{
						maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
						if(sim.baryon_flag)
						{
							maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
						}
					}
					else
					{
						maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
						if(sim.baryon_flag)
						{
							maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
						}
					}

#ifdef BENCHMARK
					update_q_count++;
#endif
				}

				for(i=0; i<cosmo.num_ncdm; i++)
				{
					if(j % (numsteps / numsteps_ncdm[i]) == 0)
					{
						if(sim.gr_flag > 0)
						{
							maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
						}
						else
						{
							maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
						}

#ifdef BENCHMARK
						update_q_count++;
#endif
					}
				}

#ifdef BENCHMARK
				update_q_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif

				for(i=0; i<cosmo.num_ncdm; i++) // move particles
				{
					if(numsteps > 1 && ((numsteps_ncdm[i] == 1 && j == numsteps / 2) || (numsteps_ncdm[i] == numsteps / 2 && j % 2 > 0)))
					{
						if(sim.gr_flag > 0)
						{
							pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
						}
						else
						{
							pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
						}
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
					}
				}

				// Background Evolution (first half of a time-step)
				if(numsteps == 1)
				{
					if(sim.mg_flag != FLAG_FR || sim.lcdm_background) // GR or Newtonian evolution
					{
						rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);
						Hubble = Hconf(a, fourpiG, cosmo);
						Rbar = Rbar_GR(a, fourpiG, cosmo);
						dot_Rbar = dot_Rbar_GR(a, Hubble, fourpiG, cosmo);
					}
					else // f(R) gravity, non-LCDM background
					{
						if(numsteps_bg == 1)
						{
							rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, 0.5 * dtau, sim);
						}
						else
						{
							for(g=0; g<numsteps_bg/2; g++)
							{
								rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau_bg, sim);
								if(tau + (g + 1) * dtau_bg > tau_print)
								{
									tau_print += dtau_print;
									if(print_background(cycle, bgoutfile, bgfilename, tau, a, Hubble, Rbar, phi_hom, -T00_hom_rescaled_a3)) return -1;
								}
							}
						}
					}
				}

				f_params[0] = a;
				f_params[1] = a * a * sim.numpts;
				if(numsteps == 1 || j == numsteps / 2) // Move particles
				{
					if(sim.gr_flag > 0)
					{
						pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
						if(sim.baryon_flag)
						{
							pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
						}
					}
					else
					{
						pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
						if(sim.baryon_flag)
						{
							pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
						}
					}

#ifdef BENCHMARK
					moveParts_count++;
					moveParts_time += MPI_Wtime() - ref2_time;
					ref2_time = MPI_Wtime();
#endif
				}

				// Background evolution (first half of a time-step)
				if(numsteps != 1) //  TODO: Doesn't this do the same as if numsteps == 1?
				{
					if(sim.mg_flag != FLAG_FR || sim.lcdm_background) // GR or Newtonian evolution
					{
						rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);
						Hubble = Hconf(a, fourpiG, cosmo);
						Rbar = Rbar_GR(a, fourpiG, cosmo);
						dot_Rbar = dot_Rbar_GR(a, Hubble, fourpiG, cosmo);
					}
					else // f(R) gravity, non-LCDM background
					{
						if(numsteps_bg == 1)
						{
							rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, 0.5 * dtau, sim);
						}
						else
						{
							for(g=0; g<numsteps_bg/2; g++)
							{
								rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau_bg, sim);
								if(tau + (g + 1) * dtau_bg > tau_print)
								{
									tau_print += dtau_print;
									if(print_background(cycle, bgoutfile, bgfilename, tau, a, Hubble, Rbar, phi_hom, -T00_hom_rescaled_a3)) return -1;
								}
							}
						}
					}
				}

				f_params[0] = a;
				f_params[1] = a * a * sim.numpts;
				for(i=0; i<cosmo.num_ncdm; i++)
				{
					if(numsteps_ncdm[i] == numsteps) // Move particles
					{
						if(sim.gr_flag > 0)
						{
							pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
						}
						else
						{
							pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
						}
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
					}
				}

				// Background evolution (second half of a time-step)
				if(sim.mg_flag != FLAG_FR || sim.lcdm_background) // GR or Newtonian evolution
				{
					rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);
					Hubble = Hconf(a, fourpiG, cosmo);
					Rbar = Rbar_GR(a, fourpiG, cosmo);
					dot_Rbar = dot_Rbar_GR(a, Hubble, fourpiG, cosmo);
				}
				else // f(R) gravity, non-LCDM background
				{
					if(numsteps_bg == 1)
					{
						rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, 0.5 * dtau, sim);
					}
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau_bg, sim);
							if(tau + (numsteps_bg/2 + g + 1) * dtau_bg > tau_print)
							{
								tau_print += dtau_print;
								if(print_background(cycle, bgoutfile, bgfilename, tau, a, Hubble, Rbar, phi_hom, -T00_hom_rescaled_a3)) return -1;
							}
						}
					}
				}
			}

			// Build background quantities after evolving background
			if(sim.mg_flag == FLAG_FR)
			{
				fbar = f(Rbar, sim, 700);
				fRbar = fR(Rbar, sim, 701);
				fRRbar = fRR(Rbar, sim, 702);
			}

			parallel.max<double>(maxvel, numspecies);
			if(sim.gr_flag > 0)
			{
				for(i=0; i<numspecies; i++)
				{
					maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
				}
			}

			tau += dtau;

			if(tau_Lambda < 0. && (cosmo.Omega_m / a / a / a) < cosmo.Omega_Lambda) // matter-dark energy equality
			{
				tau_Lambda = tau;
				COUT << "matter-dark energy equality at z=" << ((1./a) - 1.) << endl;
			}

			if(sim.mg_flag == FLAG_FR) // TODO: set new dtau in f(R) gravity -- must be done before hibernations!
			{
				dtau_old_2 = dtau_old;
				dtau_old = dtau;

				if(sim.Cf * dx < sim.steplimit / Hubble)
				{
					dtau = sim.Cf * dx;
				}
				else
				{
					dtau = sim.steplimit / Hubble;
				}
				if(!sim.perturbations_are_GR && !sim.quasi_static)
				{
					max_fRR = compute_max_fRR(deltaR, Rbar, sim);
					dtau_osci = sim.fR_epsilon_fields * sqrt(3. * max_fRR) / a;
					if(dtau > dtau_osci) dtau = dtau_osci;
				}
			}
			else
			{
				dtau_old = dtau;
				if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo)) dtau = sim.Cf * dx;
				else dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
			}

			if(sim.wallclocklimit > 0. && a <= 1.)   // check for wallclock time limit
			{
				tmp = MPI_Wtime() - start_time;
				parallel.max(tmp);
				if(tmp > sim.wallclocklimit)   // hibernate
				{
					COUT << COLORTEXT_YELLOW << " reaching hibernation wallclock limit, hibernating..." << COLORTEXT_RESET << endl;
					COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
					if(sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
					{
						plan_Bi.execute(FFT_BACKWARD);
					}

#ifdef CHECK_B
					if(sim.vector_flag == VECTOR_ELLIPTIC)
					{
						plan_Bi_check.execute(FFT_BACKWARD);
						if(sim.mg_flag == FLAG_FR)
						{
							hibernate_fR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, eightpiG_deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar, dot_Rbar, cycle);
						}
						else
						{
							hibernate_GR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle);
						}
					}
					else
					{
#endif
						if(sim.mg_flag == FLAG_FR)
						{
							hibernate_fR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, eightpiG_deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar,	dot_Rbar, cycle);
						}
						else
						{
							hibernate_GR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle);
						}
#ifdef CHECK_B
					}
#endif
					return 0; //  TODO: Maybe something more needed?
				}
			}

			if(restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
			{
				COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
				if(sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
				{
					plan_Bi.execute(FFT_BACKWARD);
				}
#ifdef CHECK_B
				if(sim.vector_flag == VECTOR_ELLIPTIC)
				{
					plan_Bi_check.execute(FFT_BACKWARD);
					if(sim.mg_flag == FLAG_FR)
					{
						hibernate_fR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, eightpiG_deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar, dot_Rbar, cycle, restartcount);
					}
					else
					{
						hibernate_GR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle, restartcount);
					}
				}
				else
				{
#endif
					if(sim.mg_flag == FLAG_FR)
					{
						hibernate_fR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, eightpiG_deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar, dot_Rbar, cycle, restartcount);
					}
					else
					{
						hibernate_GR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle, restartcount);
					}
#ifdef CHECK_B
				}
#endif
				restartcount++;
				return 0; //  TODO: Maybe something more needed?
			}

			cycle++;

#ifdef BENCHMARK
			cycle_time += MPI_Wtime() - cycle_start_time;
#endif
		}

		COUT
		<< endl
		<< "===================================  CYCLE " << cycle << "  ";
		if(cycle < 10)
		{
			COUT << "==================================" << endl;
		}
		else if(cycle < 100)
		{
			COUT << "=================================" << endl;
		}
		else if(cycle < 1000)
		{
			COUT << "================================" << endl;
		}
		else if(cycle < 10000)
		{
			COUT << "===============================" << endl;
		}

		// COUT << "              tau = " << tau << endl;
		COUT << "                z = " << (1./a) - 1. << endl;
		COUT << "                H = " << Hubble << endl;
		// COUT << "         dtau_old = " << dtau_old << endl;
		COUT << "             Rbar = " << Rbar << endl;
		// COUT << "         dot_Rbar = " << dot_Rbar << endl;
		COUT << "             fbar = " << fbar << endl;
		COUT << "            fRbar = " << fRbar << endl;
		COUT << "           fRRbar = " << fRRbar << endl;
		// COUT << " background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
		COUT << " avg/rescaled T00 = " << T00_hom*a*a*a << " / " << T00_hom_rescaled_a3 << endl;
		COUT << endl;

		check_all_fields(phi, xi, chi, deltaR, eightpiG_deltaT, zeta, Bi, numpts3d, sim);
		bgoutfile.close();
		COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef BENCHMARK
		run_time = MPI_Wtime() - start_time;

#ifdef HAVE_CLASS
		if(sim.radiation_flag > 0)
		{
			freeCLASSstructures(class_background, class_perturbs, class_spectra);
		}
#endif

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
