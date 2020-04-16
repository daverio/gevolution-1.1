//////////////////////////
// Copyright (c) 2015-2020 Julian Adamek
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

#include <iomanip>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX			// due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "class_tools.hpp"
#include "tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "fR_tools.hpp"
#include "background_fR.hpp"
#include "relaxation.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#include "ic_read_fR.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "output.hpp"
#include "hibernation.hpp"
#include "hibernation_fR.hpp"
#ifdef VELOCITY
#include "velocity.hpp"
#endif

using namespace std;
using namespace LATfield2;

int main(int argc, char **argv)
{
#ifdef BENCHMARK
	//benchmarking variables
	double ref_time, ref2_time, cycle_start_time;
	double initialization_time;
	double run_time;
	double cycle_time=0;
	double projection_time = 0;
	double snapshot_output_time = 0;
	double spectra_output_time = 0;
	double lightcone_output_time = 0;
	double gravity_solver_time = 0;
	double fft_time = 0;
	int fft_count = 0;
	double update_q_time = 0;
	int update_q_count = 0;
	double moveParts_time = 0;
	int  moveParts_count =0;
#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;
	int i, j, cycle, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numspecies, done_hij;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dx, tau, a, fourpiG, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	string bgfilename;
	string h5filename;
	char * settingsfile = NULL;
	char * precisionfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	gadget2_header hdr;
	Real T00_hom;
	Real phi_hom;

	///////////////////// additional fRevolution variables  /////////////////////
	double Hubble;
	double Rbar;
	double dot_Rbar;
	double fbar;
	double fRbar;
	double dot_fRbar;
	double fRRbar;
	double dtau_osci, dtau_bg;
	int numsteps_bg, wid, g;
	int numsteps_bg_ncdm[MAX_PCL_SPECIES-2];
	bool do_I_check;
	double tmpa, tmpHubble, tmpRbar, tmpdot_Rbar;
	Real T00_hom_rescaled_a3;
	Real Tii_hom;
	double Trace_hom;
	double max_fRR;
	std::ofstream bgoutfile;

#ifndef H5_DEBUG
	H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif

	for(i=1 ; i<argc ; i++ )
	{
		if( argv[i][0] != '-' )
			continue;
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
			case 'p':
#ifndef HAVE_CLASS
				cout << "HAVE_CLASS needs to be set at compilation to use CLASS precision files" << endl;
				exit(-100);
#endif
				precisionfile = argv[++i];
				break;
			case 'i':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_size =  atoi(argv[++i]);
				break;
			case 'g':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_group_size = atoi(argv[++i]);
		}
	}

#ifndef EXTERNAL_IO
	parallel.initialize(n,m);
#else
	if(!io_size || !io_group_size)
	{
		cout << "invalid number of I/O tasks and group sizes for I/O server (-DEXTERNAL_IO)" << endl;
		exit(-1000);
	}
	parallel.initialize(n,m,io_size,io_group_size);
	if(parallel.isIO()) ioserver.start();
	else
	{
#endif

	COUT << "   __  ___                 _        _    _            " << endl;
	COUT << "  / _|| _ \\ ___ __ __ ___ | | _  _ | |_ (_) ___  _ _  " << endl;
	COUT << " |  _||   // -_)\\ V // _ \\| || || ||  _|| |/ _ \\| ' \\ " << endl;
	COUT << " |_|  |_|_\\\\___| \\_/ \\___/|_| \\_,_| \\__||_|\\___/|_||_| version 1.1 running on " << n*m << " cores." << endl;
	COUT << COLORTEXT_RESET << endl;

	if(settingsfile == NULL)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		parallel.abortForce();
	}

	COUT << " initialising..." << endl;
	start_time = MPI_Wtime();
	numparam = loadParameterFile(settingsfile, params);
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);
	COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;
	sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
	saveParameterFile(filename, params, numparam);
	free(params);

#ifdef HAVE_CLASS
	background class_background;
	perturbs class_perturbs;
	spectra class_spectra;

	if(precisionfile != NULL)
	{
		numparam = loadParameterFile(precisionfile, params);
	}
	else
#endif
	{
		numparam = 0;
	}

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
		COUT << " Initialized multigrid with " << sim.multigrid_n_grids << " layers" << endl;
		COUT << " Dimensions of coarsest grid: " << sim.numpts/pow(2,sim.multigrid_n_grids-1) << " x " << sim.numpts/pow(2,sim.multigrid_n_grids-1) << " x " << sim.numpts/pow(2,sim.multigrid_n_grids-1) << endl;
	}

	Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> pcls_cdm;
	Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> pcls_b;
	Particles_gevolution<part_simple, part_simple_info, part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES-2];
	Field<Real> * update_cdm_fields[3];
	Field<Real> * update_b_fields[3];
	Field<Real> * update_ncdm_fields[3];
	double f_params[5];
	set<long> IDbacklog[MAX_PCL_SPECIES];

	///////////////////////// Standard gevolution fields /////////////////////////
	Field<Real> phi;
	Field<Real> chi;
	Field<Real> Bi;
	Field<Real> source;
	Field<Real> Sij;
	Field<Cplx> scalarFT;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;

  // FFT
	PlanFFT<Cplx> plan_phi;
	PlanFFT<Cplx> plan_chi;
	PlanFFT<Cplx> plan_source;
	PlanFFT<Cplx> plan_Sij;
	PlanFFT<Cplx> plan_Bi;

	// Initialising
	phi.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	plan_phi.initialize(&phi, &scalarFT);
	chi.initialize(lat,1);
	plan_chi.initialize(&chi, &scalarFT);
	source.initialize(lat,1);
	plan_source.initialize(&source, &scalarFT);
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	plan_Sij.initialize(&Sij, &SijFT);
	Bi.initialize(lat,3);
	BiFT.initialize(latFT,3);
	plan_Bi.initialize(&Bi, &BiFT);

	/////////////////////// additional fRevolution fields ///////////////////////
	Field<Real> deltaR;
	Field<Real> eightpiG_deltaT;
	Field<Real> zeta;
	Field<Real> phi_effective;
	Field<Real> phi_dot;
	Field<Real> phi_ddot;
	Field<Real> chi_dot;
	Field<Real> xi;
	Field<Real> xi_old;
	Field<Real> xi_dot;
	Field<Real> laplace_xi;
	Field<Real> rhs;
	Field<Real> residual;
	Field<Real> err;
	Field<Real> Bi_old;
	Field<Cplx> Bi_oldFT;

	// FFT
	PlanFFT<Cplx> plan_deltaR;
	PlanFFT<Cplx> plan_eightpiG_deltaT;
	PlanFFT<Cplx> plan_zeta;
	PlanFFT<Cplx> plan_phi_effective;
	PlanFFT<Cplx> plan_phi_dot;
	PlanFFT<Cplx> plan_phi_ddot;
	PlanFFT<Cplx> plan_xi;
	PlanFFT<Cplx> plan_laplace_xi;
	PlanFFT<Cplx> plan_Bi_old;

	// Initialising
	deltaR.initialize(lat,1);
	plan_deltaR.initialize(&deltaR, &scalarFT);
	eightpiG_deltaT.initialize(lat,1);
	plan_eightpiG_deltaT.initialize(&eightpiG_deltaT, &scalarFT);
	zeta.initialize(lat,1);
	plan_zeta.initialize(&zeta, &scalarFT);
	phi_effective.initialize(lat,1);
	plan_phi_effective.initialize(&phi_effective, &scalarFT);
	phi_dot.initialize(lat,1);
	plan_phi_dot.initialize(&phi_dot, &scalarFT);
	phi_ddot.initialize(lat,1);
	plan_phi_ddot.initialize(&phi_ddot, &scalarFT);
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
	rhs.initialize(lat, 1);
	rhs.alloc();
	residual.initialize(lat, 1);
	residual.alloc();
	err.initialize(lat, 1);
	err.alloc();
	Bi_old.initialize(lat,3);
	plan_Bi_old.initialize(&Bi_old, &BiFT);

	//////////////////////////////// MultiFields ////////////////////////////////
	MultiField<Real> * mg_deltaR;
	MultiField<Real> * mg_eightpiG_deltaT;
	MultiField<Real> * mg_phi;
	MultiField<Real> * mg_xi;
	MultiField<Real> * mg_xi_old;
	MultiField<Real> * mg_laplace_xi;
	MultiField<Real> * mg_rhs;
	MultiField<Real> * mg_residual;
	MultiField<Real> * mg_err;

	// Initialize
	mg_engine.initialize_Field(&phi, mg_phi);
	mg_engine.initialize_Field(&xi, mg_xi);
	mg_engine.initialize_Field(&xi_old, mg_xi_old);
	mg_engine.initialize_Field(&laplace_xi, mg_laplace_xi);
	mg_engine.initialize_Field(&deltaR, mg_deltaR);
	mg_engine.initialize_Field(&eightpiG_deltaT, mg_eightpiG_deltaT);
	mg_engine.initialize_Field(&residual, mg_residual);
	mg_engine.initialize_Field(&err, mg_err);
	mg_engine.initialize_Field(&rhs, mg_rhs);

#ifdef VELOCITY
	Field<Real> vi;
	Field<Cplx> viFT;
	vi.initialize(lat,3);
	viFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_vi(&vi, &viFT);
	double a_old;
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

	dx = 1.0 / (double) sim.numpts;
	numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

	for(i=0; i<3; i++) // particles may never move farther than to the adjacent domain
	{
		if(lat.sizeLocal(i)-1 < sim.movelimit)
		{
			sim.movelimit = lat.sizeLocal(i)-1;
		}
	}
	parallel.min(sim.movelimit);

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
	a = 1. / (1. + sim.z_in);
	tau = particleHorizon(a, fourpiG, cosmo);
	numsteps_bg = 1;

	///////////// Setup background quantities and initial timesteps /////////////
	if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
	{
		if(ic.generator != ICGEN_READ_FROM_DISK)
		{
			fR_details(cosmo, &sim, fourpiG);
			Rbar = R_initial_fR(a, fourpiG, cosmo);
			fbar = f(Rbar, sim, "initial");
			fRbar = fR(Rbar, sim, "initial");
			fRRbar = fRR(Rbar, sim, "initial");
			dot_fRbar = 0.;
			T00_hom = T00_hom_rescaled_a3 = 0.;
			Hubble = H_initial_fR(a, Hconf(a, fourpiG, cosmo), Rbar, fbar, fRbar,	6. * fourpiG * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * fRRbar / a / a / a);
			dot_Rbar = dot_R_initial_fR(a, Hubble, fourpiG, cosmo, sim);
		}
	}
	else
	{
		Hubble = Hconf(a, fourpiG, cosmo);
		Rbar = 2. * fourpiG * ( cosmo.Omega_m / a / a / a + 4. * cosmo.Omega_Lambda );
		dot_Rbar = 0.;
		fbar = fRbar = fRRbar = 0.;
	}

	if(sim.Cf * dx < sim.steplimit / Hubble)
	{
		dtau = sim.Cf * dx;
	}
	else
	{
		dtau = sim.steplimit / Hubble;
	}

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);
	h5filename += sim.basename_generic;
	h5filename += "_";
	h5filename += sim.basename_snapshot;

	bgfilename.reserve(2*PARAM_MAX_LENGTH + 24);
	bgfilename.assign(sim.output_path);
	bgfilename += sim.basename_generic;
	bgfilename += "_background.dat";

	///////////////////// Start writing the background file /////////////////////
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

		bgoutfile << "# background statistics\n#";
		bgoutfile << setw(8) 	 << "cycle";
		bgoutfile << setw(wid) << "tau/boxsize";
		bgoutfile << setw(wid) << "a";
		bgoutfile << setw(wid) << "conformal H";
		bgoutfile << setw(wid) << "R";;
		bgoutfile << setw(wid) << "phi(k=0)";
		bgoutfile << setw(wid) << "T00(k=0)";;
		bgoutfile << endl;
		bgoutfile.close();
	}

	//////////////////////// Generate initial conditions ////////////////////////
	if(ic.generator == ICGEN_BASIC)
	{
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam); // generates ICs on the fly
	}
	else if(ic.generator == ICGEN_READ_FROM_DISK)
	{
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			readIC_fR(sim, ic, cosmo, fourpiG, a, Hubble, Rbar, dot_Rbar, fbar, fRbar, fRRbar, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &xi, &xi_old, &deltaR, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
		}
		else
		{
			readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
		}
	}
#ifdef ICGEN_PREVOLUTION
	else if(ic.generator == ICGEN_PREVOLUTION)
	{
		generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam);
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

	if(sim.relativistic_flag > 0)
	{
		for(i=0; i<numspecies; i++)
			maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
	}

#ifdef VELOCITY
	a_old = a;
	projection_init(&vi);
#endif

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
	// for(i=0; i < 256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8; ++i)
	// {
	// 	hdr.fill[i] = '0';
	// }

#ifdef BENCHMARK
	initialization_time = MPI_Wtime() - start_time;
	parallel.sum(initialization_time);
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

#ifdef HAVE_CLASS
	if(sim.radiation_flag > 0 || sim.fluid_flag > 0)
	{
		initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra, params, numparam);
		if(sim.relativistic_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
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
	if(numparam > 0) free(params);
#endif

	///////////////////////////////// Main Loop /////////////////////////////////
	while(true)
	{
#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif

		do_I_check = (bool) (cycle % sim.CYCLE_INFO_INTERVAL == 0 && sim.check_fields);

		if(do_I_check) // output some info
		{
			COUT << endl << "================================ CYCLE " << cycle << " =====================";
			for(i=0; i<10-log10(cycle+1); ++i)
			{
				COUT << "=";
			}
			COUT << endl;
		}

		// construct stress-energy tensor -- Writes -a^3 T00 in {source} (note the minus sign!)
		projection_init(&source);

#ifdef HAVE_CLASS
		if(sim.radiation_flag > 0 || sim.fluid_flag > 0)
		{
			projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
		}
#endif

		if(sim.relativistic_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if(sim.baryon_flag)
			{
				projection_T00_project(&pcls_b, &source, a, &phi);
			}
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
				{
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
				}
				else if(sim.radiation_flag == 0 || (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] == 0))
				{
					tmp = bg_ncdm(a, cosmo, i);
					for(x.first(); x.test(); x.next())
					{
						source(x) += tmp;
					}
				}
			}
		}
		else // Newtonian gravity and f(R) Newtonian
		{
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if(sim.baryon_flag)
			{
				scalarProjectionCIC_project(&pcls_b, &source);
			}
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
				{
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
				}
			}
		}

		projection_T00_comm(&source);

#ifdef VELOCITY
		if((sim.out_pk & MASK_VEL) || (sim.out_snapshot & MASK_VEL))
		{
			projection_init(&Bi);
			projection_Ti0_project(&pcls_cdm, &Bi, &phi, &chi);
			vertexProjectionCIC_comm(&Bi);
			compute_vi_rescaled(cosmo, &vi, &source, &Bi, a, a_old);
			a_old = a;
		}
#endif

		// Project T0i -- writes a^4 T0i in {Bi}
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
				if(a >= 1. / (sim.z_switch_Bncdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
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
				if(sim.numpcl[1+sim.baryon_flag+i] > 0)
				{
					projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
				}
			}
		}
		projection_Tij_comm(&Sij);

#ifdef BENCHMARK
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif

		build_homogeneous_terms(source, phi, Sij, T00_hom, Tii_hom, Trace_hom, phi_hom, T00_hom_rescaled_a3, a, numpts3d);

		if(cycle % sim.CYCLE_INFO_INTERVAL == 0)
		{
			output_background_data(tau, a, Hubble, dtau_old, Rbar, dot_Rbar, cosmo, T00_hom, T00_hom_rescaled_a3, fbar, fRbar, fRRbar, sim.modified_gravity_flag);
		}

		compute_eightpiG_deltaT(eightpiG_deltaT, source, Sij, a, T00_hom, Trace_hom, fourpiG, sim);
		if(sim.modified_gravity_flag != MODIFIED_GRAVITY_FLAG_FR) // in GR/Newton, deltaR tracks exactly eightpiG_deltaT
		{
			copy_field(eightpiG_deltaT, deltaR, -1.);
		}

		if(do_I_check)
		{
			check_all_fields(phi, xi, laplace_xi, chi, deltaR, eightpiG_deltaT, zeta, phi_effective, phi_ddot, Bi, numpts3d, sim);
		}

		///////////////////////// EVOLVE deltaR, zeta, xi /////////////////////////
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			if(sim.fR_model == FR_MODEL_R2) // R + R^2 can be solved with FFT
			{
				// Writes source term in {zeta}, then to Fourier space, {zeta} -> {scalarFT}
				copy_field(eightpiG_deltaT, zeta, a*a/6./sim.fR_params[0]);
				plan_zeta.execute(FFT_FORWARD);
				solveModifiedPoissonFT(scalarFT, scalarFT, 1., a*a/6./sim.fR_params[0]);
				plan_deltaR.execute(FFT_BACKWARD);
				add_fields(deltaR, eightpiG_deltaT, zeta);
				// Writes old and new xi
				copy_field(xi, xi_old);
				convert_deltaR_to_xi(deltaR, xi, Rbar, fRbar, sim); // Needs value of deltaR
			}
			else // Generic f(R) model (not R + R^2), quasi-static
			{
				copy_field(xi, xi_old);
				prepare_initial_guess_trace_equation(deltaR, eightpiG_deltaT, xi, laplace_xi, rhs, a, dx, Rbar, fRbar, sim);

				if(!cycle)
				{
					copy_field(xi, xi_old);
				}

				relaxation_solver(mg_phi, mg_xi, mg_xi_old, mg_laplace_xi, mg_deltaR, mg_eightpiG_deltaT, mg_rhs, mg_residual, mg_err, mg_engine, dx, a, dtau_old ? 2. * Hubble / dtau_old : 0., numpts3d, Rbar, fbar, fRbar, sim);
				add_fields(deltaR, eightpiG_deltaT, zeta); // Build zeta from deltaR and eightpiG_deltaT
				sim.multigrid_check_shape = 0; // Just check the first time (at most)
			}
			// TODO: Avoid updating the halo twice (see later)
			build_laplacian(xi, laplace_xi, dx);
		}

		//////////////////////////////// Evolve phi ////////////////////////////////
		if(sim.relativistic_flag) // GR and relativistic f(R)
		{
			if(cycle)
			{
				if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
				{
					prepareFTsource_S00_fR_rel<Real>(source, phi, chi, xi, xi_old, deltaR, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), dx*dx, dtau_old, Hubble, a, fourpiG / a, Rbar, fbar, fRbar, sim);
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
				if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR) // Relativistic f(R)
				{
					solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hubble / dtau_old);
				}
				else // GR perturbations
				{
					solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);
				}

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif

				plan_phi_dot.execute(FFT_BACKWARD);	 // go back to position space

#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif

				flip_fields(phi, phi_dot); // Now {phi} contains phi_new, {phi_dot} contains 	phi_old
				add_fields(phi, 1./dtau_old, phi_dot, -1./dtau_old, phi_dot); // Computes phi_dot = (phi_new - phi_old)/dtau
				add_fields(phi, 1., xi, -0.5, phi_effective);
				phi_dot.updateHalo();
				phi.updateHalo(); // communicate halo values
			}
		}
		else // Newton or Newtonian f(R)
		{
			if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
			{
				prepareFTsource_S00_fR_Newtonian<Real>(source, phi, chi, xi, xi_old, deltaR, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), dx*dx, dtau_old, Hubble, a, fourpiG / a, Rbar, fbar, fRbar, sim);
			}
#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
			if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
			{
				solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx));  // Newton: phi update (k-space)
			}
			else
			{
				solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)
			}

#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif

			plan_phi_dot.execute(FFT_BACKWARD);	 // go back to position space

#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif

			flip_fields(phi, phi_dot); // Now {phi} contains phi_new, {phi_dot} contains 	phi_old
			add_fields(phi, 1./dtau_old, phi_dot, -1./dtau_old, phi_dot); // Computes phi_dot = (phi_new - phi_old)/dtau
			add_fields(phi, 1., xi, -0.5, phi_effective);
			phi_dot.updateHalo();
			phi.updateHalo(); // communicate halo values
		}

		// prepare nonlinear source for additional equations
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			prepareFTsource_Sij_fR<Real>(phi, chi, xi, Sij, Sij, 2. * fourpiG * dx * dx / a);
		}
		else
		{
			prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);
		}

		//////////////////////////////// Evolve chi ////////////////////////////////
		copy_field(chi, chi_dot); // Before evolving chi, keep copy of old value

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif

		plan_Sij.execute(FFT_FORWARD); // go to k-space

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

		chi.updateHalo(); // communicate halo values

		if(dtau_old)
		{
			add_fields(chi, 1./dtau_old, chi_dot, -1./dtau_old, chi_dot);
			chi_dot.updateHalo();
		}

		//////////////////////////////// Evolve Bi ////////////////////////////////
		if(sim.vector_flag == VECTOR_ELLIPTIC) // solve B using elliptic constraint
		{
			if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR && sim.relativistic_flag)
			{
				build_laplacian(phi, laplace_xi, dx);
				add_fR_source_S0i(Bi_old, phi, laplace_xi, chi, xi, xi_dot, Bi, dot_fRbar, dx, fourpiG);
				build_laplacian(xi, laplace_xi, dx);
			}

#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			plan_Bi.execute(FFT_FORWARD);

#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
			projectFTvector(BiFT, BiFT, fourpiG * dx * dx);
		}
		else // evolve B using vector projection
		{
			// TODO: Verify and fix for f(R) gravity!
			evolveFTvector(SijFT, BiFT, a * a * dtau_old);
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

		// Builds phi_ddot - to be compared with e.g. deltaR for a test of the quasi-static approximation
		build_phi_ddot(phi, chi, phi_dot, chi_dot, phi_ddot, deltaR, dx, a, Hubble, Rbar);

		if(do_I_check)
		{
			check_all_fields(phi, xi, laplace_xi, chi, deltaR, eightpiG_deltaT, zeta, phi_effective, phi_ddot, Bi, numpts3d, sim);
		}

		// record background data
		if(kFT.setCoord(0, 0, 0))
		{
			if(print_background(cycle, bgoutfile, bgfilename, tau, a, Hubble, Rbar, phi_hom, -T00_hom_rescaled_a3))
			{
				return -1;
			}
		}

		// Background timesteps in f(R)
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR && !sim.lcdm_background)
		{
			// This would be the "ideal" timestep to evolve the f(R) background
			dtau_bg = sim.fR_epsilon_bg * sqrt(3.0 * fRRbar) / a;

			if(dtau >= dtau_bg)
			{
				numsteps_bg = (int) 2 * ceil(dtau/2./dtau_bg) ; // must be even
				dtau = dtau_bg * numsteps_bg;
			}
			else
			{
				numsteps_bg = 2;
			}
		}
		else
		{
			numsteps_bg = 2;
		}

		///////////////////////////// lightcone output /////////////////////////////
		if(sim.num_lightcone > 0)
		{
			writeLightcones(sim, cosmo, fourpiG, a, tau, dtau, dtau_old, maxvel[0], cycle, (string) sim.output_path + sim.basename_lightcone, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &Sij, &BiFT, &SijFT, &plan_Bi, &plan_Sij, done_hij, IDbacklog);
		}
		else
		{
			done_hij = 0;
		}

#ifdef BENCHMARK
		lightcone_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		///////////////////////////// snapshot output /////////////////////////////
		if(snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

			writeSnapshots(sim, cosmo, fourpiG, a, dtau_old, done_hij, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &phi_ddot, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef VELOCITY
				, &vi
#endif
			);

			snapcount++;
		}

#ifdef BENCHMARK
		snapshot_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		////////////////////////////// power spectra //////////////////////////////
		if(pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

			writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle,
#ifdef HAVE_CLASS
				class_background, class_perturbs, class_spectra, ic,
#endif
				&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &phi_ddot, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_phi_ddot, &plan_phi_effective, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef VELOCITY
				, &vi, &viFT, &plan_vi
#endif
			);

			pkcount++;
		}

#ifdef EXACT_OUTPUT_REDSHIFTS
		tmpa = a;
		tmpHubble = Hubble;
		tmpRbar = Rbar;
		tmpdot_Rbar = dot_Rbar;

		rungekutta_background(tmpa, tmpHubble, tmpRbar, tmpdot_Rbar, fourpiG, 0.5 * dtau, Trace_hom, cosmo, sim, dtau / numsteps_bg, numsteps_bg / 2);
		rungekutta_background(tmpa, tmpHubble, tmpRbar, tmpdot_Rbar, fourpiG, 0.5 * dtau, Trace_hom, cosmo, sim, dtau / numsteps_bg, numsteps_bg / 2);

		if(pkcount < sim.num_pk && 1. / tmpa < sim.z_pk[pkcount] + 1.)
		{
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle,
#ifdef HAVE_CLASS
				class_background, class_perturbs, class_spectra, ic,
#endif
				&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &eightpiG_deltaT, &xi, &laplace_xi, &zeta, &phi_ddot, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_eightpiG_deltaT, &plan_xi, &plan_laplace_xi, &plan_zeta, &plan_phi_ddot, &plan_phi_effective, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef VELOCITY
				, &vi, &viFT, &plan_vi
#endif
			);
		}
#endif // EXACT_OUTPUT_REDSHIFTS

#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif

		if(pkcount >= sim.num_pk && snapcount >= sim.num_snapshot)
		{
			for(i=0; i<sim.num_lightcone; i++)
			{
				if(sim.lightcone[i].z + 1. < 1. / a)
				{
					i = sim.num_lightcone + 1;
				}
			}

			if(i == sim.num_lightcone)
			{
				break; // simulation complete
			}
		}

		// For when the full background evolution only starts at sim.z_switch_fR_background
		if(1./a - 1. < sim.z_switch_fR_background && sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR && sim.lcdm_background)
		{
			COUT << endl;
			COUT << " ===============================================================" << endl;
			COUT << " Beginning full f(R) background evolution at z = " << 1./a - 1. << endl;
			COUT << " ===============================================================" << endl;
			COUT << endl;

			if(sim.fR_model == FR_MODEL_HU_SAWICKI)
			{
				cosmo.Omega_Lambda = 0.;
			}
			sim.lcdm_background = 0;
		}

		// Rescale phi -> psi in Newtonian f(R) gravity -- Needed to evolve particles correctly
		if(!sim.relativistic_flag && sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			subtract_fields(phi, chi, phi);
			phi.updateHalo();
		}

		// compute number of step subdivisions for ncdm particle updates
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
		}

		if(cycle % sim.CYCLE_INFO_INTERVAL == 0)
		{
			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
			if(sim.baryon_flag)
			{
				COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
			}
			COUT << "), time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau;
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

#ifdef BENCHMARK
		ref2_time = MPI_Wtime();
#endif

		for(i=0; i<cosmo.num_ncdm; i++) // non-cold DM particle update
		{
			if(sim.numpcl[1+sim.baryon_flag+i] == 0) continue;

			tmpa = a;
			tmpHubble = Hubble;
			tmpRbar = Rbar;
			tmpdot_Rbar = dot_Rbar;

			if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR && numsteps_bg > numsteps_ncdm[i])
			{
				numsteps_bg_ncdm[i] = (int) ceil(numsteps_bg / numsteps_ncdm[i]); // Extra division for the timestep, if ncdm particles need more resolution than the f(R) background
			}
			else
			{
				numsteps_bg_ncdm[i] = 1;
			}

			for(j=0; j<numsteps_ncdm[i]; j++)
			{
				f_params[0] = tmpa;
				f_params[1] = tmpa * tmpa * sim.numpts;
				if(sim.relativistic_flag > 0)
				{
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				}
				else
				{
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
				}

#ifdef BENCHMARK
				update_q_count++;
				update_q_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif

				rungekutta_background(tmpa, tmpHubble, tmpRbar, tmpdot_Rbar, fourpiG, 0.5 * dtau / numsteps_ncdm[i], Trace_hom, cosmo, sim, dtau / (numsteps_bg * numsteps_bg_ncdm[i]), numsteps_bg * numsteps_ncdm[i] / 2);

				f_params[0] = tmpa;
				f_params[1] = tmpa * tmpa * sim.numpts;

				if(sim.relativistic_flag > 0)
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

				rungekutta_background(tmpa, tmpHubble, tmpRbar, tmpdot_Rbar, fourpiG, 0.5 * dtau / numsteps_ncdm[i], Trace_hom, cosmo, sim, dtau / (numsteps_bg * numsteps_bg_ncdm[i]), numsteps_bg * numsteps_ncdm[i] / 2);

			}
		}

		// cdm and baryon particle update
		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
		if(sim.relativistic_flag > 0)
		{
			maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
			if(sim.baryon_flag)
			{
				maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
			}
		}
		else
		{
			maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
			if(sim.baryon_flag)
			{
				maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
			}
		}

#ifdef BENCHMARK
		update_q_count++;
		update_q_time += MPI_Wtime() - ref2_time;
		ref2_time = MPI_Wtime();
#endif

		rungekutta_background(a, Hubble, Rbar, dot_Rbar, fourpiG, 0.5 * dtau, Trace_hom - 4.*cosmo.Omega_Lambda, cosmo, sim, dtau / numsteps_bg, numsteps_bg / 2);

		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;

		if(sim.relativistic_flag > 0)
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
#endif

		rungekutta_background(a, Hubble, Rbar, dot_Rbar, fourpiG, 0.5 * dtau, Trace_hom - 4.*cosmo.Omega_Lambda, cosmo, sim, dtau / numsteps_bg, numsteps_bg / 2);

		// Build background quantities after evolving background
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			fbar = f(Rbar, sim, "building background cycle "+to_string(cycle));
			dot_fRbar = - fRbar / dtau;
			fRbar = fR(Rbar, sim, "building background cycle "+to_string(cycle));
			dot_fRbar += fRbar / dtau;
			fRRbar = fRR(Rbar, sim, "building background cycle "+to_string(cycle));
		}

		parallel.max<double>(maxvel, numspecies);
		if(sim.relativistic_flag > 0)
		{
			for(i=0; i<numspecies; i++)
			{
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
			}
		}
		// done particle update

		tau += dtau;
		dtau_old = dtau;
		if(sim.Cf * dx < sim.steplimit / Hubble)
		{
			dtau = sim.Cf * dx;
		}
		else
		{
			dtau = sim.steplimit / Hubble;
		}

		if(sim.wallclocklimit > 0.)   // check for wallclock time limit
		{
			tmp = MPI_Wtime() - start_time;
			parallel.max(tmp);
			if(tmp > sim.wallclocklimit)   // hibernate
			{
				COUT << COLORTEXT_YELLOW << " reaching hibernation wallclock limit, 2ting..." << COLORTEXT_RESET << endl;
				COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
				if(sim.vector_flag == VECTOR_PARABOLIC && sim.relativistic_flag == 0)
				{
					plan_Bi.execute(FFT_BACKWARD);
				}

				if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
				{
					hibernate_fR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_old, deltaR, a, Hubble, Rbar, dot_Rbar, fbar, fRbar, fRRbar, tau, dtau, dtau_old, cycle);
				}
				else
				{
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle);
				}

				break;
			}
		}

		if(restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

			if(sim.vector_flag == VECTOR_PARABOLIC && sim.relativistic_flag == 0)
			{
				plan_Bi.execute(FFT_BACKWARD);
			}

			if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
			{
				hibernate_fR(sim, ic, cosmo, hdr, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_old, deltaR, a, Hubble, Rbar, dot_Rbar, fbar, fRbar, fRRbar, tau, dtau, dtau_old, cycle, restartcount);
			}
			else
			{
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle, restartcount);
			}

			restartcount++;
		}

		cycle++;

#ifdef BENCHMARK
		cycle_time += MPI_Wtime() - cycle_start_time;
#endif
	}

	COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef BENCHMARK
		ref_time = MPI_Wtime();
#endif

#ifdef HAVE_CLASS
	if(sim.radiation_flag > 0 || sim.fluid_flag > 0)
	{
		freeCLASSstructures(class_background, class_perturbs, class_spectra);
	}
#endif

#ifdef BENCHMARK
	lightcone_output_time += MPI_Wtime() - ref_time;
	run_time = MPI_Wtime() - start_time;

	parallel.sum(run_time);
	parallel.sum(cycle_time);
	parallel.sum(projection_time);
	parallel.sum(snapshot_output_time);
	parallel.sum(spectra_output_time);
	parallel.sum(lightcone_output_time);
	parallel.sum(gravity_solver_time);
	parallel.sum(fft_time);
	parallel.sum(update_q_time);
	parallel.sum(moveParts_time);

	COUT << endl << "BENCHMARK" << endl;
	COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
	COUT << "total number of cycles: "<< cycle << endl;
	COUT << "time consumption breakdown:" << endl;
	COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
	COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;

	COUT << "----------- main loop: components -----------"<<endl;
	COUT << "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
	COUT << "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
	COUT << "lightcone outputs          : "<< hourMinSec(lightcone_output_time) << " ; " << 100. * lightcone_output_time/cycle_time <<"%."<<endl;
	COUT << "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
	COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time <<"%."<<endl;
	COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time <<"%."<<endl;
	COUT << "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
	COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
#endif

#ifdef EXTERNAL_IO
		ioserver.stop();
	}
#endif

	return 0;
}
