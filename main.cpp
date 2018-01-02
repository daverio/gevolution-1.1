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
#include <stdlib.h>
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
#include "fofR_tools.hpp"
#include "background_fofR.hpp"

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
	double cycle_time=0;
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
	double dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, dx, tau, tau_temp, a, fourpiG, tau_Lambda, tmp, start_time;
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
	Real Trace_hom;
	Real phi_hom;
	Real max_FRR;

#ifndef H5_DEBUG
	H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif

	for(i=1 ; i < argc ; i++ ){
		if( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
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


	string settingsfile_bin(settingsfile);
	settingsfile_bin += ".bin";

#ifndef EXTERNAL_IO
	parallel.initialize(n,m);
#else
	parallel.initialize(n,m,io_size,io_group_size);
	if(parallel.isIO()) ioserver.start();
	else
	{
#endif

	COUT << COLORTEXT_WHITE << endl;
	COUT << "  __  ___                 _        _    _  " << endl
	 		 << " / _|| _ \\ ___ __ __ ___ | | _  _ | |_ (_) ___  _ _ " << endl
			 << "|  _||   // -_)\\ V // _ \\| || || ||  _|| |/ _ \\| ' \\ " << endl
			 << "|_|  |_|_\\\\___| \\_/ \\___/|_| \\_,_| \\__||_|\\___/|_||_|   version 1.0 running on " << n*m << " cores." << endl;

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

	bgfilename.reserve(2*PARAM_MAX_LENGTH + 24);
	bgfilename.assign(sim.output_path);
	bgfilename += sim.basename_generic;
	bgfilename += "_background.dat";

	std::ofstream bgoutfile;

	box[0] = sim.numpts;
	box[1] = sim.numpts;
	box[2] = sim.numpts;

	Lattice lat(3,box,2);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);

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

	//additional F(R) fields
	Field<Real> xi;
	Field<Real> xi_prev;
	Field<Real> zeta;
	Field<Real> deltaR;
	Field<Real> deltaR_plus_deltaT_minus_zeta;
	Field<Real> dot_deltaR;
	Field<Real> deltaR_prev;
	Field<Real> deltaT;
	Field<Real> phidot;
	Field<Real> xidot;
	Field<Real> laplace_xi;

  PlanFFT<Cplx> plan_source;
  PlanFFT<Cplx> plan_phi;
  PlanFFT<Cplx> plan_phidot;
  PlanFFT<Cplx> plan_chi;
  PlanFFT<Cplx> plan_zeta;
  PlanFFT<Cplx> plan_deltaR; // TODO: needed only to output power spectra
  PlanFFT<Cplx> plan_deltaR_prev;
  PlanFFT<Cplx> plan_deltaT;
  PlanFFT<Cplx> plan_xi;

  PlanFFT<Cplx> plan_laplace_xi;

	PlanFFT<Cplx> plan_Bi;
  PlanFFT<Cplx> plan_Sij;

#ifdef CHECK_B
  PlanFFT<Cplx> plan_Bi_check;
#endif


// TODO: For the moment, initialize and alloc everything. After debugging, only stuff that's needed
phi.initialize(lat,1);
phi.alloc();
chi.initialize(lat,1);
source.initialize(lat,1);
scalarFT.initialize(latFT,1);
xi.initialize(lat,1);
xi_prev.initialize(lat,1);
xi_prev.alloc();
xidot.initialize(lat,1);
xidot.alloc();

laplace_xi.initialize(lat,1);
deltaR_plus_deltaT_minus_zeta.initialize(lat,1);
deltaR_plus_deltaT_minus_zeta.alloc();

zeta.initialize(lat,1);
zeta.initialize(lat,1);
deltaR.initialize(lat,1);
dot_deltaR.initialize(lat, 1);
dot_deltaR.alloc();
deltaR_prev.initialize(lat,1);
deltaT.initialize(lat,1);
phidot.initialize(lat,1);
plan_source.initialize(&source, &scalarFT);
plan_phi.initialize(&phi, &scalarFT);
plan_phidot.initialize(&phidot, &scalarFT);
plan_chi.initialize(&chi, &scalarFT);
plan_xi.initialize(&xi, &scalarFT);

plan_laplace_xi.initialize(&laplace_xi, &scalarFT);

plan_deltaR.initialize(&deltaR, &scalarFT);
plan_deltaR_prev.initialize(&deltaR_prev, &scalarFT);
plan_deltaT.initialize(&deltaT, &scalarFT);
plan_zeta.initialize(&zeta, &scalarFT);
Bi.initialize(lat,3);
BiFT.initialize(latFT,3);
plan_Bi.initialize(&Bi, &BiFT);
Sij.initialize(lat,3,3,symmetric);
SijFT.initialize(latFT,3,3,symmetric);
plan_Sij.initialize(&Sij, &SijFT);

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

	//f(R) background description and gsl spline
	double Hubble;
	double Rbar;
	double dot_Rbar;
	double Fbar;
	double FRbar;
	double FRRbar;
	double t1, t2, t3; // TODO: Temp variables, remove after debugging
	double coeffm2;

	gsl_spline * a_spline;
	gsl_spline * H_spline;
	gsl_spline * Rbar_spline;
	gsl_interp_accel * gsl_inpl_acc = gsl_interp_accel_alloc ();
	const gsl_interp_type *gsl_inpl_type = gsl_interp_linear;

	Site x(lat);
	rKSite kFT(latFT);

	dx = 1.0 / (double) sim.numpts;
	numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

	COUT << "numpts = " << sim.numpts << "   linesize = " << phi.lattice().size(1) << endl;

	for(i = 0; i < 3; i++) // particles may never move farther than to the adjacent domain
	{
		if(lat.sizeLocal(i)-1 < sim.movelimit)
			sim.movelimit = lat.sizeLocal(i)-1;
	}
	parallel.min(sim.movelimit);

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
	a = 1. / (1. + sim.z_in);
	//TODO check particleHorizon
	tau = particleHorizon(a, fourpiG, cosmo);
	COUT << " initial tau/boxsize = " << tau << "\n";
	COUT << " fourpiG = " << fourpiG << "\n";

	tau_Lambda = -1.0;

	if(sim.mg_flag == FOFR && ic.generator != ICGEN_READ_FROM_DISK)
	{
		fR_details(cosmo, &sim, fourpiG);

		if(!sim.lcdm_background)
		{
			Rbar = R_initial_fR(a, fourpiG, cosmo);
			Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 100);
			FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 101);
			FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 102);
			Hubble = H_initial_fR(a, Hconf(a, fourpiG, cosmo), Rbar, Fbar, FRbar,	6. * fourpiG * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * FRRbar / a / a / a); // TODO: Check Omega_m term

			if(sim.background_trace)
			{
				dot_Rbar = dot_R_initial_fR(a, Hubble, fourpiG, cosmo, sim.fofR_params, sim.fofR_type);
			}

			if(sim.Cf * dx < sim.steplimit / Hubble)
			{
				dtau = sim.Cf * dx;
			}
			else
			{
				dtau = sim.steplimit / Hubble;
			}
		}
		else
		{
			Hubble = Hconf(a, fourpiG, cosmo);
			Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
			Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 200);
			FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 201);
			FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 202);
			if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo)) dtau = sim.Cf * dx;
		  else dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
		}
		if(!sim.back_to_GR && !sim.quasi_static) // If back_to_GR or quasi-static, dtau for fields dictated by usual steplimit (background is still split in smaller timesteps)
		{
			dtau_osci = sim.fofR_epsilon_fields * sqrt(3.0 * fabs(FRRbar))/a;
			if(dtau > dtau_osci) dtau = dtau_osci;
		}
	}
	else if(sim.mg_flag != FOFR)
	{
		Hubble = Hconf(a, fourpiG, cosmo);
		Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
		if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo)) dtau = sim.Cf * dx;
	  else dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
	}

	dtau_old = 0.;
	dtau_old_2 = 0.;

	//////////////////////////////////////////////////////// Background Only
	if(sim.background_only)
	{
		if(sim.mg_flag == FOFR)
		{
			dtau_print = tau * sqrt(sim.bg_initial_redshift + 1.) / sqrt(sim.bg_final_redshift + 1.) / sim.BACKGROUND_NUMPTS;
			tau_print = tau;
		}

		if(parallel.rank() == 0)
		{
			bgoutfile.open(bgfilename);
			if(bgoutfile == NULL)
			{
				cout << " error opening file for background output!" << endl;
				return -1;
			}

			wid = 6;
			bgoutfile << scientific << setprecision(wid);
			wid += 9;

			bgoutfile << "# background statistics\n#" << setw(8)
				    		<< "cycle" << setw(wid)
								<< "tau/boxsize" << setw(wid)
								<< "a" << setw(wid)
								<< "conformal H" << setw(wid)
								<< "R";
			if(sim.mg_flag == FOFR)
			{
				bgoutfile << setw(wid) << "dR/R";
			}
			bgoutfile << endl;
		}

		while(a < 1./(1. + sim.bg_final_redshift))
		{
			if(sim.mg_flag == FOFR && !sim.lcdm_background)
			{
				dtau_old = dtau;
				dtau = sim.fofR_epsilon_bg * sqrt(3. * fabs(FRRbar)) / a;
				if(sim.background_trace)
				{
					dtau = rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau, sim.fofR_params, sim.fofR_type);
				}
				else
				{
					dtau = rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, dtau, sim.fofR_params, sim.fofR_type);
				}
				Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 300);
				FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 301);
				FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 302);
				tau += dtau; // TODO: CHECK
			}
			else
			{
				if(sim.Cf * dx < sim.steplimit / Hubble)
				{
					dtau = sim.Cf * dx;
				}
				else
				{
					dtau = sim.steplimit / Hubble;
				}

				rungekutta4bg(a, fourpiG, cosmo, dtau);  // GR -- evolve background
				Hubble = Hconf(a, fourpiG, cosmo);
				Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4. * cosmo.Omega_Lambda);
				tau += dtau;
			}

			// Write background data on file
			if(parallel.rank() == 0)
			{
				if(sim.mg_flag != FOFR || tau >= tau_print)
				{
					cout << " cycle " << cycle << "\tz = " << 1./a - 1. << "\ttau = " << tau << "\tdtau = " << dtau << endl;
					bgoutfile << setw(9) << cycle << setw(wid) << tau	<< setw(wid) << a	<< setw(wid) << Hubble << setw(wid) << Rbar;
					if(sim.mg_flag == FOFR)
					{
						tau_print += dtau_print;
						t1 = Rbar/R_GR(a, fourpiG, cosmo) - 1.;
						bgoutfile << setw(wid) << fabs(t1);
					}
					bgoutfile << endl;
				}
			}

			cycle++;
		}

		// Write last step for f(R)
		if(parallel.rank() == 0)
		{
			if(sim.mg_flag == FOFR)
			{
				cout << " FINAL STEP:\n";
				cout << " cycle " << cycle << "\tz = " << 1./a - 1. << "\ttau = " << tau << "\tdtau = " << dtau << endl;
				bgoutfile << setw(9) << cycle	<< setw(wid) << tau	<< setw(wid) << a	<< setw(wid) << Hubble << setw(wid) << Rbar	<< "\n";
			}

			bgoutfile.close();
			cout << "\n Background evolution complete.\n";
		}

		return 0;
	}


		//////////////////////////////////////////////////////// Full evolution -- Background + Perturbations

		// Generate Initial Conditions
	if(ic.generator == ICGEN_BASIC)
	{
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij); // generates ICs on the fly
	}
	else if(ic.generator == ICGEN_READ_FROM_DISK)
	{
		if(sim.mg_flag == FOFR)
		{
				readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, dtau_old_2, dtau_osci, dtau_bg, Hubble, Rbar, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &xi, &xi_prev, &zeta, &deltaR, &deltaR_prev, &dot_deltaR, &deltaT, &phidot, &xidot, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount, settingsfile_bin);

				Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 400);
				FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 401);
				FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 402);
		}
		else
		{
			readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
			Hubble = Hconf(a, fourpiG, cosmo);
			Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
		}
	}
#ifdef ICGEN_PREVOLUTION
	else if(ic.generator == ICGEN_PREVOLUTION)
		generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
#ifdef ICGEN_FALCONIC
	else if(ic.generator == ICGEN_FALCONIC)
		maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
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
		for(i = 0; i < numspecies; i++)
			maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
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

	for(i = 0; i < 6; i++)
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
	for(i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; i++)
		hdr.fill[i] = 0;


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
				chi(x) += source(x);
			chi.updateHalo();
		}
	}
#endif

	//////////////////////////////////////////////////////// Main loop
	while(true)
	{
		if(1./a - 1. < sim.z_check) // TODO: Remove after debugging
		{
			sim.check_fields = 1;
		}

		if(1./a - 1. < sim.z_switch_fR_background && sim.mg_flag == FOFR && sim.lcdm_background)
		{
			COUT << "\n                              ----- Beginning full f(R) background evolution at z = " << 1./a - 1. << " -----\n";
			cosmo.Omega_Lambda = 0.;
			sim.lcdm_background = 0;
		}

		if(cycle % sim.CYCLE_INFO_INTERVAL == 0 && sim.check_fields)
		{
			check_field(phi, "phi", numpts3d);
			check_field(deltaR, "deltaR", numpts3d);
			check_field(deltaT, "deltaT", numpts3d);
			if(sim.mg_flag == FOFR)
			{
				check_field(zeta, "zeta", numpts3d);
				check_field(deltaR_plus_deltaT_minus_zeta, "deltaR_plus_deltaT_minus_zeta", numpts3d);
			}
		}

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor -- Write -a^3 T00 in {source} (note the minus sign!)
		projection_init(&source);

#ifdef HAVE_CLASS
		if(sim.radiation_flag > 0)
			projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
#endif

		if(sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if(sim.baryon_flag)
				projection_T00_project(&pcls_b, &source, a, &phi);
				for(i = 0; i < cosmo.num_ncdm; i++)
				{
					if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
					else if(sim.radiation_flag == 0)
					{
						tmp = bg_ncdm(a, cosmo, i);
						for(x.first(); x.test(); x.next())
						source(x) += tmp;
					}
				}
				projection_T00_comm(&source);
		}
		else
		{//n-body gauge
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if(sim.baryon_flag)
				scalarProjectionCIC_project(&pcls_b, &source);
			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
			}
			projection_T00_comm(&source);
		}

		if(sim.vector_flag == VECTOR_ELLIPTIC)
		{
			projection_init(&Bi);
			projection_T0i_project(&pcls_cdm, &Bi, &phi);
			if (sim.baryon_flag)
			projection_T0i_project(&pcls_b, &Bi, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.))
				projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
			}
			projection_T0i_comm(&Bi);
		}

		projection_init(&Sij); // Write a^3 Tij in {Sij}
		projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
		if(sim.baryon_flag)
			projection_Tij_project(&pcls_b, &Sij, a, &phi);
		if(a >= 1. / (sim.z_switch_linearchi + 1.))
		{
			for(i = 0; i < cosmo.num_ncdm; i++)
				projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
		}
		projection_Tij_comm(&Sij);

#ifdef BENCHMARK
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif

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

			// Computes 8piG*deltaT = 8piG( (T00 + T11 + T22 + T33) - (-rho_backg + 3P_backg) ), to be stored in {deltaT}
			// At the moment: {source} = -a^3*T00, {Sij} = a^3*Tij
			computeTtrace(deltaT, source, Sij, a*a*a, cosmo.Omega_m/a/a/a, fourpiG); // TODO: Should this include the explicit Lambda term?

			T00_hom /= a*a*a; // T00_hom contained a^3 * <T00>, now contains the correct <T00>
			Tii_hom /= a*a*a; // same for Tii
			Trace_hom = T00_hom + Tii_hom; // same for Trace

			if(sim.mg_flag != FOFR)
			{
				copy_field(deltaT, deltaR, -1.);
			}

			if(cycle % sim.CYCLE_INFO_INTERVAL == 0)
			{
				COUT << "\n cycle " << cycle << ", background information:\n"
				     << "           z = " << (1./a) - 1. << "\n"
						 << "           background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << "\n"
						 << "           average/rescaled T00 = " << T00_hom*a*a*a << " / " << T00_hom_rescaled_a3 << "\n"
						 << "           phi_hom = " << phi_hom << "\n";
			}

			if(sim.mg_flag == FOFR)
			{
				if(sim.back_to_GR) // If back_to_GR, assign deltaR = deltaT, zeta = 0, xi = 0 at all times
				{
					flip_fields(deltaR, deltaR_prev, zeta);
					copy_field(deltaT, deltaR, -1.);
					if(cycle) add_fields(deltaR, deltaR_prev, dot_deltaR, 1./dtau_old, -1./dtau_old); // Computes dot_deltaR at time (t-1/2)
					flip_fields(xi, xi_prev);
					copy_field(xi, xi, 0.); // Sets xi = 0.
					copy_field(zeta, zeta, 0.); // Sets zeta = 0.

					// TODO: Are these necessary?
					deltaR.updateHalo();
					deltaR_prev.updateHalo();
					zeta.updateHalo();
					dot_deltaR.updateHalo();
				}
				else if(sim.fofR_type == FOFR_TYPE_R2)
				{
					coeffm2 = a * a / 6. / sim.fofR_params[0];
					if(sim.quasi_static)
					{
						copy_field(deltaT, zeta, coeffm2); // Writes source term in {zeta}
						plan_zeta.execute(FFT_FORWARD); // Fourier space, {zeta} -> {scalarFT}
						// Solve wave-like equation transformed into Poisson-like
						solveModifiedPoissonFT(scalarFT, scalarFT, 1., coeffm2); // TODO CHECK
						plan_zeta.execute(FFT_BACKWARD); // Back to real space, {scalarFT} -> {zeta}
						flip_fields(deltaR, deltaR_prev, zeta);
						deltaR.updateHalo();
						deltaR_prev.updateHalo(); // TODO: Probably unnecessary
						add_fields(deltaR, deltaT, zeta);
						copy_field(xi, xi_prev);
						xi_set(xi, deltaR, Rbar, FRbar, sim.fofR_params, sim.fofR_type); // Needs value of deltaR
					}
					else // TODO: LEAPFROG
					{
						if(cycle <= 1) // Initial conditions for non-quasi-static are the same as quasi-static!
						{
							// copy_field(deltaT, zeta, coeffm2); // Writes source term in {zeta}
							Site xx(zeta.lattice());
							for(xx.first(); xx.test(); xx.next())
							{
								zeta(xx) = deltaT(xx+0) + deltaT(xx-0) + deltaT(xx+1) + deltaT(xx-1) + deltaT(xx+2) + deltaT(xx-2) - 6. * deltaT(xx);
								zeta(xx) /= dx*dx;
								// zeta(xx) = 0.;
							}
							plan_zeta.execute(FFT_FORWARD); // Fourier space, {zeta} -> {scalarFT}
							// Solve wave-like equation transformed into Poisson-like
							solveModifiedPoissonFT(scalarFT, scalarFT, 1., coeffm2); // TODO CHECK
							plan_zeta.execute(FFT_BACKWARD); // Back to real space, {scalarFT} -> {zeta}
							flip_fields(deltaR, deltaR_prev);
							add_fields(zeta, deltaT, deltaR, 1., -1.);
							zeta.updateHalo();
							deltaR.updateHalo();
						}
						else // Else builds deltaR_i from deltaR_{i-1} and dot_deltaR_{i-1/2}, for i >= 2
						{
							add_fields(deltaR, dot_deltaR, deltaR, 1., dtau_old);
							add_fields(deltaR, deltaT, zeta);
							deltaR.updateHalo();
							zeta.updateHalo();
						}
						// Done at each step:
						copy_field(xi, xi_prev);// TODO: is this needed?
						xi_set(xi, deltaR, Rbar, FRbar, sim.fofR_params, sim.fofR_type); // Needs value of deltaR
						if(cycle == 1) // Builds dot_deltaR_{1/2} from deltaR_1 and deltaR_0
						{
							add_fields(deltaR, deltaR_prev, dot_deltaR, 1./dtau_old, -1./dtau_old);
						}
						if(cycle)
						{
							lp_update_dotR(deltaR, deltaT, dot_deltaR, dot_deltaR, Hubble, coeffm2, dtau_old, dtau, dx * dx);
						}
					}
				}
				else // Generic f(R) model (not R + R^2)
				{
					if(sim.quasi_static)
					{
						COUT << "Only f(R) = R + R^2 is implemented in quasi-static so far, sorry! :)\n";
						exit(0);
					}
					else
					{
						if(cycle <= 1) // Initial conditions are deltaR = - deltaT
						{
							copy_field(deltaR, deltaR_prev);
							copy_field(xi, xi_prev);
							copy_field(deltaT, deltaR, -1.);
							xi_set(xi, deltaR, Rbar, FRbar, sim.fofR_params, sim.fofR_type);
							deltaR.updateHalo();
							deltaR_prev.updateHalo(); // TODO: Probably unnecessary
							if(cycle == 1) // Builds xidot_{1/2} from xi_1 and xi_0
							{
								add_fields(xi, xi_prev, xidot, 1./dtau_old, -1./dtau_old);
							}
						}
						else // Else builds xi_i from xi_{i-1} and xidot_{i-1/2}, for i >= 2
						{
							add_fields(xi, xidot, xi, 1., dtau_old);
							xi.updateHalo();
						}

						computeDRzeta(deltaT, deltaR, zeta, xi, Rbar, FRbar, sim.fofR_params, sim.fofR_type, sim.fofR_target_precision);

						if(cycle)
						{
							lp_update_dotxi(xi, zeta, xidot, xidot, Hubble, a*a, dtau_old, dtau, dx*dx);
						}
					}
				}
			}

			if(sim.mg_flag == FOFR) //TODO: Remove after debugging
			{
				add_fields(deltaR, deltaT, deltaR_plus_deltaT_minus_zeta);
				add_fields(deltaR_plus_deltaT_minus_zeta, zeta, deltaR_plus_deltaT_minus_zeta, 1., -1.);
			}

			if(cycle) // Only after the 0th cycle
			{
				if(sim.mg_flag == FOFR)
				{
					// Write source term for the 00 equation in {source}
					prepareFTsource_S00<Real>(source, phi, chi, xi, xi_prev, deltaR, source, -(cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo))/a/a/a, dx*dx, dtau_old, Hubble, a, fourpiG / a, Rbar, Fbar, FRbar, sim.fofR_params, sim.fofR_type, sim.back_to_GR, sim.quasi_static);
					// prepareFTsource_S00<Real>(source, phi, chi, xi, xi_prev, deltaR, source, T00_hom - bg_ncdm(a, cosmo), dx*dx, dtau_old, Hubble, a, fourpiG / a, Rbar, Fbar, FRbar, sim.fofR_params, sim.fofR_type, sim.back_to_GR, sim.quasi_static);
				}
				else if(sim.gr_flag)
				{
					prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update
				}

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				plan_source.execute(FFT_FORWARD);  // go to k-space, {source} --> FFT_FORWARD --> {scalarFT}
#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif
				// phi update (k-space)
				if(sim.mg_flag == FOFR)
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
				add_fields(phi, phidot, phidot, 1./dtau_old, -1./dtau_old); // Computes phidot = (phi_new - phi_old)/dtau
				// Update halos for phi and phidot
				phi.updateHalo();
				phidot.updateHalo();

#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif

			}
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
			plan_phi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
			phi.updateHalo();  // communicate halo values
			phidot.updateHalo();
		}

		// record some background data
		if(kFT.setCoord(0, 0, 0))
		{
			if(!cycle)
			{
				bgoutfile.open(bgfilename);
				if(bgoutfile == NULL)
				{
					cout << " error opening file for background output!" << endl;
					return -1;
				}

				wid = 6;
				bgoutfile << scientific << setprecision(wid);
				wid += 9;

				bgoutfile << "# background statistics\n#"
									<< setw(8) 	 << "cycle"
									<< setw(wid) << "tau/boxsize"
									<< setw(wid) << "a"
									<< setw(wid) << "conformal H"
									<< setw(wid) << "R"
							 		<< setw(wid) << "phi(k=0)"
									<< setw(wid) << "T00(k=0)"
									<< "\n";
			}
			else
			{
				bgoutfile.open(bgfilename, std::ofstream::app);
				if(bgoutfile == NULL)
				{
					cout << " error opening file for background output!" << endl;
					return -1;
				}
			}

			bgoutfile << setw(9) << cycle << setw(wid) << tau	<< setw(wid) << a	<< setw(wid) << Hubble << setw(wid) << Rbar	<< setw(wid) << phi_hom << setw(wid) << -T00_hom_rescaled_a3 << "\n";

			bgoutfile.close();
		}

		// done recording background data
		prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations

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
		projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		plan_chi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif

		if(sim.mg_flag == FOFR && cycle) // TODO: check if it should be done only after the 0th time step
		{
			add_fields(chi, xi, chi);
		}
		chi.updateHalo();  // communicate halo values

		if(sim.vector_flag == VECTOR_ELLIPTIC) // TODO: check for f(R)!
		{
			#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
			#endif
			plan_Bi.execute(FFT_FORWARD);
			#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
			#endif
			projectFTvector(BiFT, BiFT, fourpiG * dx * dx); // solve B using elliptic constraint (k-space)
			#ifdef CHECK_B
			evolveFTvector(SijFT, BiFT_check, a * a * dtau_old);
			#endif
		}
		else
		evolveFTvector(SijFT, BiFT, a * a * dtau_old);  // evolve B using vector projection (k-space)

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
			writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &deltaT, &xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_deltaT, &plan_xi, &plan_zeta, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &deltaT, &xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_deltaT, &plan_xi, &plan_zeta, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
			snapcount++;
		}

#ifdef BENCHMARK
		snapshot_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// power spectra
		if(pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		// if(cycle < 6)
		{

			for(x.first(); x.test(); x.next())
			{
				laplace_xi(x) = xi(x+0) + xi(x-0) + xi(x+1) + xi(x-1) + xi(x+2) + xi(x-2) - 6.*xi(x);
				laplace_xi(x) /= dx*dx;
			}

			COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &deltaT, &laplace_xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_deltaT, &plan_laplace_xi, &plan_zeta, &plan_phidot, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, cycle, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &deltaR, &deltaT, &laplace_xi, &zeta, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_deltaR, &plan_deltaT, &plan_laplace_xi, &plan_zeta, &plan_phidot, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
			pkcount++;
		}

#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif
if(pkcount >= sim.num_pk && snapcount >= sim.num_snapshot) break; // simulation complete

		// compute number of step subdivisions for particle updates
		numsteps = 1;
		for(i = 0; i < cosmo.num_ncdm; i++)
		{
			if(dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;

			if(numsteps < numsteps_ncdm[i]) numsteps = numsteps_ncdm[i];
		}
		if(numsteps > 1 && numsteps % 2 > 0) numsteps++;   // if >1, make it an even number

		for(i = 0; i < cosmo.num_ncdm; i++)
		{
			if(numsteps / numsteps_ncdm[i] <= 1) numsteps_ncdm[i] = numsteps;
			else if(numsteps_ncdm[i] > 1) numsteps_ncdm[i] = numsteps / 2;
		}

		numsteps_bg = 1;

		if(sim.mg_flag == FOFR && !sim.lcdm_background)
		{
			dtau_bg = sim.fofR_epsilon_bg * sqrt(3. * FRRbar) / a; // This would be the "ideal" timestep to evolve the f(R) background
			if(dtau >= 2. * dtau_bg)// numsteps_bg must be even -- we split the background evolution in two steps of dtau/2 each
			{
				numsteps_bg = (int) (dtau / dtau_bg / 2.);
				numsteps_bg = (int) (2 * numsteps_bg);
				dtau_bg = (double) (dtau / numsteps_bg);
			}
			else dtau_bg = dtau;
		}
		else
		{
			dtau_bg = dtau;
		}

		if(cycle % sim.CYCLE_INFO_INTERVAL == 0)
		{
			COUT << "           max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx << ")\n";
			if(sim.baryon_flag)
			{
				COUT << "           baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx << ")\n";
			}
			if(sim.mg_flag == FOFR)
			{
				COUT << "           Hubble * dtau = " <<  Hubble * dtau << ", dx = " << dx << "\n";
				if(!sim.lcdm_background) COUT << "           numsteps_bg = " << numsteps_bg;
			}
			else
			{
				COUT << "           time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau;
			}

			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(i == 0)
				{
					COUT << endl << " time step subdivision for ncdm species: ";
				}
				COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i+1+sim.baryon_flag] << ")";
				if(i < cosmo.num_ncdm-1)
				{
					COUT << ", ";
				}
			}
			COUT << endl;
		}

		for(j = 0; j < numsteps; j++) // particle update
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
						maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				}
				else
				{
					maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
					if(sim.baryon_flag)
						maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
				}

#ifdef BENCHMARK
				update_q_count++;
#endif
			}

			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(j % (numsteps / numsteps_ncdm[i]) == 0)
				{
					if(sim.gr_flag > 0)
						maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

#ifdef BENCHMARK
					update_q_count++;
#endif
				}
			}

#ifdef BENCHMARK
			update_q_time += MPI_Wtime() - ref2_time;
			ref2_time = MPI_Wtime();
#endif

			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(numsteps > 1 && ((numsteps_ncdm[i] == 1 && j == numsteps / 2) || (numsteps_ncdm[i] == numsteps / 2 && j % 2 > 0)))
				{
					if(sim.gr_flag > 0)
						pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
				}
			}

			if(numsteps == 1)
			{
				if(sim.mg_flag != FOFR || sim.lcdm_background) // GR or Newtonian evolution
				{
					rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);  // GR -- evolve background by half a time step
					Hubble = Hconf(a, fourpiG, cosmo);
					Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4. * cosmo.Omega_Lambda);
				}
				else if(sim.background_trace) // f(R) gravity, non-LCDM background, trace equation
				{
					if(numsteps_bg == 1)
					{
						rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, 0.5 * dtau, sim.fofR_params, sim.fofR_type);
					}
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
						}
					}
				}
				else // f(R) gravity, non-LCDM background, Friedmann equation
				{
					if(numsteps_bg == 1)
					{
						rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, 0.5 * dtau, sim.fofR_params, sim.fofR_type);
					}
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
						}
					}
				}

				if(sim.mg_flag == FOFR)
				{
					Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 500);
					FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 501);
					FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 502);
				}
			}
			// if numsteps == 1 ends here

			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			if(numsteps == 1 || j == numsteps / 2)
			{
				if(sim.gr_flag > 0)
				{
					pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
					if(sim.baryon_flag)
						pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
				}
				else
				{
					pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
					if(sim.baryon_flag)
						pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
				}

#ifdef BENCHMARK
				moveParts_count++;
				moveParts_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif
			}

			if(numsteps != 1) //  TODO: Doesn't this do the same as if numsteps == 1?
			{
				if(sim.mg_flag != FOFR || sim.lcdm_background) // GR or Newtonian evolution
				{
					rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);  // GR -- evolve background by half a time step
					Hubble = Hconf(a, fourpiG, cosmo);
					Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4. * cosmo.Omega_Lambda);
				}
				else if(sim.background_trace) // f(R) gravity, non-LCDM background, trace equation
				{
					if(numsteps_bg == 1)
					{
						rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, 0.5 * dtau, sim.fofR_params, sim.fofR_type);
					}
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
						}
					}
				}
				else // f(R) gravity, non-LCDM background, Friedmann equation
				{
					if(numsteps_bg == 1)
					{
						rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, 0.5 * dtau, sim.fofR_params, sim.fofR_type);
					}
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
						}
					}
				}

				if(sim.mg_flag == FOFR)
				{
					Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 600);
					FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 601);
					FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 602);
				}
			} // if (numsteps != 1) ends here

			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(numsteps_ncdm[i] == numsteps)
				{
					if(sim.gr_flag > 0)
															pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
				}
			}

			if(sim.mg_flag != FOFR || sim.lcdm_background) // GR or Newtonian evolution
			{
				rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);  // GR -- evolve background by half a time step
				Hubble = Hconf(a, fourpiG, cosmo);
				Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4. * cosmo.Omega_Lambda);
			}
			else if(sim.background_trace) // f(R) gravity, non-LCDM background, trace equation
			{
				if(numsteps_bg == 1)
				{
					rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, 0.5 * dtau, sim.fofR_params, sim.fofR_type);
				}
				else
				{
					for(g=0; g<numsteps_bg/2; g++)
					{
						rungekutta_fR_trace(a, Hubble, Rbar, dot_Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
					}
				}
			}
			else // f(R) gravity, non-LCDM background, Friedmann equation
			{
				if(numsteps_bg == 1)
				{
					rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, 0.5 * dtau, sim.fofR_params, sim.fofR_type);
				}
				else
				{
					for(g=0; g<numsteps_bg/2; g++)
					{
						rungekutta_fR_45(a, Hubble, Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
					}
				}
			}

			if(sim.mg_flag == FOFR)
			{
				Fbar = F(Rbar, sim.fofR_params, sim.fofR_type, 700);
				FRbar = FR(Rbar, sim.fofR_params, sim.fofR_type, 701);
				FRRbar = FRR(Rbar, sim.fofR_params, sim.fofR_type, 702);
			}
		}   // particle update done
			// for(j = 0; j < numsteps; j++) ends here

		parallel.max<double>(maxvel, numspecies);

		if(sim.gr_flag > 0)
		{
			for(i = 0; i < numspecies; i++)
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		}

		tau += dtau;

		if(sim.mg_flag == FOFR)
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
			if(!sim.back_to_GR && !sim.quasi_static)
			{
				max_FRR = compute_max_FRR(deltaR, Rbar, sim.fofR_params, sim.fofR_type);
				dtau_osci = sim.fofR_epsilon_fields * sqrt(3. * max_FRR) / a;
				if(dtau > dtau_osci) dtau = dtau_osci;
			}
		}

		if(tau_Lambda < 0. && (cosmo.Omega_m / a / a / a) < cosmo.Omega_Lambda)
		{
			tau_Lambda = tau;
			COUT << "matter-dark energy equality at z=" << ((1./a) - 1.) << endl;
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
					plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
				if(sim.vector_flag == VECTOR_ELLIPTIC)
				{
					plan_Bi_check.execute(FFT_BACKWARD);
					if(sim.mg_flag == FOFR)
						hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, Hubble, Rbar, cycle);
					else
						hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle);
				}
				else
				{
#endif
					if(sim.mg_flag == FOFR)
						hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, deltaT, phidot, xidot,	a, tau, dtau, dtau_old, dtau_old_2, Hubble, Rbar,	cycle);
					else
						hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle);
#ifdef CHECK_B
				}
#endif
				break;
			}
		}

		if(restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			if(sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
				plan_Bi.execute(FFT_BACKWARD);

#ifdef CHECK_B
			if(sim.vector_flag == VECTOR_ELLIPTIC)
			{
				plan_Bi_check.execute(FFT_BACKWARD);
				if(sim.mg_flag == FOFR)
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, Hubble, Rbar, cycle, restartcount);
				else
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau, cycle, restartcount);
			}
			else
			{
#endif
				if(sim.mg_flag == FOFR)
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, xi, xi_prev, zeta, deltaR, deltaR_prev, dot_deltaR, deltaT, phidot, xidot, a, tau, dtau, dtau_old, dtau_old_2, Hubble, Rbar, cycle, restartcount);
				else
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau, cycle, restartcount);
#ifdef CHECK_B
			}
#endif
			restartcount++;
		}

		if(sim.mg_flag != FOFR)
		{
			dtau_old = dtau;
			if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo)) dtau = sim.Cf * dx;
			else dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
		}

		cycle++;

#ifdef BENCHMARK
		cycle_time += MPI_Wtime() - cycle_start_time;
#endif
	}

	bgoutfile.close();

	COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef BENCHMARK
	run_time = MPI_Wtime() - start_time;

#ifdef HAVE_CLASS
	if(sim.radiation_flag > 0)
		freeCLASSstructures(class_background, class_perturbs, class_spectra);
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

	COUT << endl << "BENCHMARK" << endl;
	COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
	COUT << "total number of cycles: "<< cycle << endl;
	COUT << "time consumption breakdown:" << endl;
	COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
	COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;

	COUT << "----------- main loop: components -----------"<<endl;
	COUT << "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
	COUT << "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
	COUT << "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
	COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time <<"%."<<endl;
	COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time <<"%."<<endl;
	COUT << "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
	COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time << "%." << endl;
#endif

#ifdef EXTERNAL_IO
		ioserver.stop();
	}
#endif

	return 0;
}
