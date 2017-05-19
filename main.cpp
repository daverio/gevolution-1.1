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
	int  moveParts_count = 0;

#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;

	int i, j, g, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numsteps_bg, numspecies;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dtau_new, dtau_osci, dtau_bg, dx, tau, tau_temp, a, fourpiG, tau_Lambda, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	string h5filename;
	char * settingsfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	gadget2_header hdr;
	Real T00hom;
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
	h5filename += '_';
	h5filename += sim.basename_snapshot;

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
	Field<Real> T00_old;
	Field<Real> T0i_old;
	Field<Real> T00_new;
	Field<Real> T0i_new;
	Field<Real> Si;
	Field<Real> Sij;
	Field<Cplx> scalarFT;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;
	//additional F(R) fields

	Field<Real> xi;
	Field<Real> xi_prev;
	Field<Real> xi_next;
	Field<Real> zeta;
	Field<Real> deltaR;
	Field<Real> deltaT;
	Field<Real> phidot;
	Field<Real> xidot;
	Field<Cplx> SiFT;

	double temp_val;// TODO Remove after debugging

PlanFFT<Cplx> plan_source;
PlanFFT<Cplx> plan_phi;
PlanFFT<Cplx> plan_phidot;
PlanFFT<Cplx> plan_chi;
PlanFFT<Cplx> plan_xi_next;
PlanFFT<Cplx> plan_Bi;
PlanFFT<Cplx> plan_Si;
PlanFFT<Cplx> plan_Sij;
#ifdef CHECK_B
PlanFFT<Cplx> plan_Bi_check;
#endif

if(sim.mg_flag == FOFR)
{
	phi.initialize(lat,1);
	phi.alloc();
	chi.initialize(lat,1);
	source.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	xi.initialize(lat,1);
	xi.alloc();
	xi_prev.initialize(lat,1);
	xi_prev.alloc();
	xi_next.initialize(lat,1);
	xi_next.alloc();
	zeta.initialize(lat,1);
	zeta.alloc();
	deltaR.initialize(lat,1);
	deltaR.alloc();
	deltaT.initialize(lat,1);
	deltaT.alloc();
	phidot.initialize(lat,1);
	plan_source.initialize(&source, &scalarFT);
	plan_phi.initialize(&phi, &scalarFT);
	plan_phidot.initialize(&phidot, &scalarFT);
	plan_chi.initialize(&chi, &scalarFT);
	plan_xi_next.initialize(&xi_next, &scalarFT);
	Bi.initialize(lat,3);
	Si.initialize(lat,3);
	BiFT.initialize(latFT,3);
	SiFT.initialize(latFT,3);
	plan_Bi.initialize(&Bi, &BiFT);
	plan_Si.initialize(&Si, &SiFT);
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	plan_Sij.initialize(&Sij, &SijFT);
}
else
{
	source.initialize(lat,1);
	phi.initialize(lat,1);
	phidot.initialize(lat,1);
	chi.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	plan_source.initialize(&source, &scalarFT);
	plan_phidot.initialize(&phidot, &scalarFT);
	plan_phi.initialize(&phi, &scalarFT);
	plan_chi.initialize(&chi, &scalarFT);
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	plan_Sij.initialize(&Sij, &SijFT);
	Bi.initialize(lat,3);
	BiFT.initialize(latFT,3);
	plan_Bi.initialize(&Bi, &BiFT);
}
if(sim.energy_conservation)
{
	T00_old.initialize(lat,1);
	T00_old.alloc();
	T0i_old.initialize(lat,3);
	T0i_old.alloc();

	T00_new.initialize(lat,1);
	T00_new.alloc();
	T0i_new.initialize(lat,3);
	T0i_new.alloc();
}

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
	double Rbar;
	double Hubble;

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
	tau_Lambda = -1.0;

	if(sim.mg_flag == FOFR)
	{
		if(sim.read_bg_from_file == 1)
			loadBackground(a_spline,H_spline,Rbar_spline,gsl_inpl_type, sim.background_filename);
		if(sim.fofR_type == FOFR_TYPE_HU_SAWICKI)
		{
			rescale_params_Hu_Sawicki(cosmo, fourpiG, &sim);
		}
		if(sim.read_bg_from_file == 1)
		{
			a = gsl_spline_eval(a_spline, tau, gsl_inpl_acc);
			Hubble = gsl_spline_eval(H_spline, tau, gsl_inpl_acc);
			Rbar = gsl_spline_eval(Rbar_spline, tau, gsl_inpl_acc);
		}
		else{
			Rbar = R_initial_fR(a, 2. * fourpiG, cosmo);
			Hubble = H_initial_fR(a,
														Hconf(a, fourpiG, cosmo),
														Rbar,
														F(Rbar, sim.fofR_params, sim.fofR_type),
														FR(Rbar, sim.fofR_params, sim.fofR_type),
														6. * fourpiG * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * FRR(Rbar, sim.fofR_params, sim.fofR_type, 1)/a/a/a // TODO: Check Omega_m term
														);
		}
		if(sim.Cf * dx < sim.steplimit / Hubble) dtau = sim.Cf * dx;
	  else dtau = sim.steplimit / Hubble;
		if(!sim.back_to_GR)
		{
			dtau_osci = sim.fofR_epsilon_fields * sqrt(3.0 * FRR(Rbar, sim.fofR_params, sim.fofR_type, 2))/a;
			if(dtau > dtau_osci) dtau = dtau_osci;
		}
	}
	else
	{
		Hubble = Hconf(a, fourpiG, cosmo);
		Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
		if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
			dtau = sim.Cf * dx;
	  else
		dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
	}

	dtau_old = 0.;

	if(ic.generator == ICGEN_BASIC)
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij); // generates ICs on the fly
	else if(ic.generator == ICGEN_READ_FROM_DISK)
		readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
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

	dtau_old = dtau;
	dtau_new = dtau;

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

	ofstream xi_output("xi_evolution.dat");

	if(sim.follow_xi)
	{
		if(!xi_output)
		{
			cout << "Cannot open file xi_evolution.dat. Closing...\n";
			exit(1);
		}
		xi_output << setprecision(10);
	}
	Site site_xi = check_field(phi, "phi");

	int check_original = sim.check_fields;

	while(true)    // main loop
	{
		// sim.check_fields = check_original * (cycle % CYCLE_INFO_INTERVAL == 0);
		if(sim.check_fields)
		{
			COUT << "\n---------- CYCLE " << cycle << " ----------\n";
			COUT << " Scale factor = " << a << "  Hubble = " << Hubble << "   Rbar = " << Rbar << "   dtau / Hubble time = " << dtau * Hubble;
			if(sim.mg_flag == FOFR) COUT << "   numsteps_bg = " << numsteps_bg;
			COUT << endl;
			cin.get();
		}

		if(sim.energy_conservation)
		{
			write_T00(T00_old, T00_new);
			write_T0i(T0i_old, T0i_new);
			T0i_old.updateHalo();
		}

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor
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

		if(sim.energy_conservation)
		{
			write_T00(T00_new, source, a);
		}

		if(sim.mg_flag == GENREL)
		{
			if(sim.vector_flag == VECTOR_ELLIPTIC)
			{
				projection_init(&Bi);
				projection_T0i_project(&pcls_cdm, &Bi, &phi);
				if(sim.baryon_flag)
				projection_T0i_project(&pcls_b, &Bi, &phi);
				for(i = 0; i < cosmo.num_ncdm; i++)
				{
					if(a >= 1. / (sim.z_switch_Bncdm[i] + 1.))
					projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
				}
				projection_T0i_comm(&Bi);
				if(sim.energy_conservation)
				{
					write_T0i(T0i_new, Bi, a);
					T0i_new.updateHalo();
				}
				phi.updateHalo();
				phidot.updateHalo();
				Bi.updateHalo();
			}
		}
		else if(sim.mg_flag == FOFR) // Write a^4 T^0_i on Si
		{
			projection_init(&Si);
			projection_T0i_project(&pcls_cdm, &Si, &phi);
			if(sim.baryon_flag)
			{
				projection_T0i_project(&pcls_b, &Si, &phi);
			}
			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(a >= 1. / (sim.z_switch_Bncdm[i] + 1.))
				projection_T0i_project(pcls_ncdm+i, &Si, &phi);
			}
			projection_T0i_comm(&Si);
			phi.updateHalo();
			phidot.updateHalo();
			Si.updateHalo();
		}

		projection_init(&Sij); // Write a^3 Tij on Sij
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
			T00hom = 0.;
			for(x.first(); x.test(); x.next())
			{
				T00hom += source(x);
			}

			parallel.sum<Real>(T00hom);
			T00hom /= (Real) numpts3d;

			if(cycle % CYCLE_INFO_INTERVAL == 0)
			{
				//TODO check bg_ncdm() in background.hpp.
				COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
			}

			if(sim.mg_flag == FOFR)
			{
				// Computes 8piG*deltaT = 8piG( (T00 + T11 + T22 + T33) - (-rho_backg + 3P_backg) ), to be stored in {deltaT}
				// At the moment: {source} = -a^3*T00, {Sij} = a^3*Tij
				computeTtrace(deltaT, source, Sij, a, Tbar(a, cosmo), 2.*fourpiG);
				if(!cycle) // If 0-th step, prepare xi and zeta
				{
					// TODO: Verify for hibernation, and see if initial conditions are good
					xi_initial_conditions(xi, deltaT, Rbar, FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type, sim.back_to_GR); // Needs value of deltaT (computeTtrace)
					xi_prev_initial_conditions(xi_prev, xi, sim.back_to_GR);
					xi_prev_initial_conditions(xi_next, xi, sim.back_to_GR); // TODO: change this after debugging
					zeta_initial_conditions(zeta, sim.back_to_GR);
					zeta.updateHalo();
					deltaR_initial_conditions(deltaR, deltaT, zeta); // TODO: might not be necessary

					check_field(zeta, "zeta", "AFTER initial_conditions\n");
					check_field(source, "source");
					check_field(deltaR, "deltaR");
					check_field(deltaT, "deltaT");
					COUT << " Rbar = " << Rbar << endl;

					xi_prev.updateHalo();
					xi.updateHalo();
					xi_next.updateHalo();
				}
				else if(sim.quasi_static)
				{
					// TODO: At the moment, we are fixing xi to the bottom of the potential, i.e. the GR solution deltaR + deltaT = 0
					xi_initial_conditions(xi_next, deltaT, Rbar, FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type, sim.back_to_GR);
					update_xi(xi, xi_next, xi_prev); // Writes: xi -> {xi_old}, xi_next -> {xi}
					xi.updateHalo();

					// TODO: Update value of deltaR with new value of deltaT and OLD value of zeta -- Might be superfluous and/or wrong
					deltaR_initial_conditions(deltaR, deltaT, zeta);
				}

				if(sim.follow_xi && x.setCoord(site_xi.coord(0), site_xi.coord(1), site_xi.coord(2)));
				{
					xi_output << tau << " " << phi(site_xi) << " " << xi(site_xi) << endl;
				}
			}

			if(sim.check_fields)
			{
				check_field(source, "source", "BEFORE prepareFTsource\n");
				check_field(phi, "phi");
				if(sim.mg_flag == FOFR)
				{
					check_field(xi_prev, "xi_prev");
					check_field(xi, "xi");
					check_field(xi_next, "xi_next");
					check_field(zeta, "zeta");
					check_field(deltaR, "deltaR");
					check_field(deltaT, "deltaT");
				}
			}

			if(cycle) // Only after the 0th cycle
			{
				if(sim.mg_flag == FOFR)
				{
					if(sim.read_bg_from_file == 1)
					{
						a = gsl_spline_eval(a_spline, tau, gsl_inpl_acc);
						Hubble = gsl_spline_eval(H_spline, tau, gsl_inpl_acc);
						Rbar = gsl_spline_eval(Rbar_spline, tau, gsl_inpl_acc);
					}

					// Write source term for the 00 equation in {source}
					prepareFTsource_S00<Real>(source, phi, chi, xi, xi_prev, deltaT, zeta, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), dx*dx, dtau, 	Hubble, a, fourpiG / a, Rbar, F(Rbar,sim.fofR_params,sim.fofR_type), FR(Rbar,sim.fofR_params,sim.fofR_type), sim.fofR_params, sim.fofR_type, 	sim.back_to_GR);

					if(sim.check_fields)
					{
						check_field(source, "source", "AFTER prepareFTsource_S00\n");
						check_field(phi, "phi");
						check_field(phidot, "phidot");
						if(sim.mg_flag == FOFR)
						{
							check_field(xi, "xi");
							check_field(deltaR, "deltaR");
							check_field(deltaT, "deltaT");
							check_field(zeta, "zeta");
						}
					}
				}
				else if(sim.gr_flag)
				{
					prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update
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
					if(sim.quasi_static)
					{
						solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Quasi-static, as in Newton
					}
					else
					{
						solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hubble / dtau);
					}
				}
				else  // GR
				{
					solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau);
				}


#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				// go back to position space
				plan_phidot.execute(FFT_BACKWARD);
				update_phi_phidot(phi, phidot, dtau);
				// Update halos for phi and phidot
				phi.updateHalo();
				phidot.updateHalo();

				if(sim.check_fields)
				{
					check_field(source, "source", "AFTER solveModifiedPoissonFT\n");
					check_field(phi, "phi");
					check_field(phidot, "phidot");
					if(sim.mg_flag == FOFR)
					{
						check_field(xi_prev, "xi_prev");
						check_field(xi, "xi");
						check_field(xi_next, "xi_next");
						check_field(zeta, "zeta");
						check_field(deltaR, "deltaR");
						check_field(deltaT, "deltaT");
					}
				}

#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif

				if(sim.mg_flag == FOFR)
				{
					if(sim.back_to_GR)
					{
						for(x.first(); x.test(); x.next()) // TODO: Is this necessary?
						{
							xi(x) = 0.;
						}
						xi.updateHalo();// TODO: Probably unnecessary
					}
					else if(sim.quasi_static)
					{
						xi_initial_conditions(xi_next, deltaT, Rbar, FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type, sim.back_to_GR);
						update_xi(xi, xi_next, xi_prev); // Writes: xi -> {xi_old}, xi_next -> {xi}
						xi.updateHalo();
					}
					else
					{
						// Prepare the source for the trace equation -> {xi_next}, to evolve xi
						prepareFTsource_xi<Real>(deltaT, phi, xi, xi_prev, xi_next, dx * dx, dtau, dtau_old, a * a, Hubble, Rbar, F(Rbar, sim.fofR_params, sim.fofR_type), FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type);

						if(sim.check_fields)
						{
							check_field(source, "source", "AFTER prepareFTsource_xi\n");
							check_field(phi, "phi");
							check_field(phidot, "phidot");
							check_field(xi, "xi");
							check_field(xi_next, "xi_next");
							check_field(xi_prev, "xi_prev");
						}

						plan_xi_next.execute(FFT_FORWARD); // Fourier space, {xi_next} -> {scalarFT}
						// Solve wave-like equation transformed into Poisson-like
						solveTraceXi(scalarFT, scalarFT, 2. * (Hubble + 1./dtau) / dtau_old );
						plan_xi_next.execute(FFT_BACKWARD); // Back to real space, {scalarFT} -> {xi_next}
						if(sim.check_fields)
						{
							check_field(source, "source", "AFTER plan_xi_next\n");
							check_field(phi, "phi");
							check_field(phidot, "phidot");
							check_field(xi, "xi");
							check_field(xi_next, "xi_next");
							check_field(xi_prev, "xi_prev");
						}
						update_xi(xi, xi_next, xi_prev); // Writes: xi -> {xi_old}, xi_next -> {xi}
					}
					xi.updateHalo();
					xi_prev.updateHalo();
					xi_next.updateHalo();// TODO: might not be necessary
					if(sim.check_fields)
					{
						check_field(source, "source", "AFTER update_xi\n");
						check_field(phi, "phi");
						check_field(phidot, "phidot");
						check_field(xi, "xi");
						check_field(xi_next, "xi_next");
						check_field(xi_prev, "xi_prev");
					}
				}
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

			if(sim.check_fields)
			{
				check_field(source, "source", "AFTER solveModifiedPoissonFT\n");
				check_field(phi, "phi");
				check_field(phidot, "phidot");
			}
		}

		// record some background data
		if(kFT.setCoord(0, 0, 0))
		{
			sprintf(filename, "%s%s_background.dat", sim.output_path, sim.basename_generic);
			outfile = fopen(filename, "a");
			if(outfile == NULL)
			{
				cout << " error opening file for background output!" << endl;
			}
			else
			{
				if(cycle == 0)
					fprintf(outfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
					if(sim.mg_flag == FOFR)
					fprintf(outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau, a, Hubble , scalarFT(kFT).real(), T00hom);
					else
					fprintf(outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau, a, Hconf(a, fourpiG, cosmo) / Hconf(1., fourpiG, cosmo), scalarFT(kFT).real(), T00hom);

					fclose(outfile);
			}
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

		if(sim.mg_flag == FOFR && cycle) addXi(chi, xi, sim.back_to_GR); // TODO: check if it should be done only after the 0th time step
		chi.updateHalo();  // communicate halo values

		if(sim.vector_flag == VECTOR_ELLIPTIC)
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
			evolveFTvector(SijFT, BiFT_check, a * a * dtau);
#endif
		}
		else
			evolveFTvector(SijFT, BiFT, a * a * dtau);  // evolve B using vector projection (k-space)

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

		if(sim.mg_flag == FOFR)
		{
			max_FRR = computeDRzeta(deltaT, deltaR, zeta, xi, Rbar, FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type, sim.fofR_target_precision);
			if(sim.check_fields)
			{
				check_xi_constraint(xi, deltaT, zeta, Rbar, FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type);
			}
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
			writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &xi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_xi_next, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &xi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_xi_next, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
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

#ifdef CHECK_B
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &xi_next, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_xi_next, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &xi_next, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_xi_next, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif

			pkcount++;
		}

#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif
if(pkcount >= sim.num_pk && snapcount >= sim.num_snapshot) break; // simulation complete

		// compute number of step subdivisions for particle updates
		numsteps = 1; // TODO: why numsteps = 1? Ask Julian
		for(i = 0; i < cosmo.num_ncdm; i++)
		{
			if(dtau_new * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau_new * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;

			if(numsteps < numsteps_ncdm[i]) numsteps = numsteps_ncdm[i];
		}
		if(numsteps > 1 && numsteps % 2 > 0) numsteps++;   // if >1, make it an even number

		for(i = 0; i < cosmo.num_ncdm; i++)
		{
			if(numsteps / numsteps_ncdm[i] <= 1) numsteps_ncdm[i] = numsteps;
			else if(numsteps_ncdm[i] > 1) numsteps_ncdm[i] = numsteps / 2;
		}

		if(cycle % CYCLE_INFO_INTERVAL == 0)
		{
			if(sim.mg_flag == FOFR && sim.check_fields)
			{
				check_field(xi, "xi", "Checking f(R)\n");
				check_field(zeta, "zeta");
				check_field(deltaR, "deltaR");
				check_field(deltaT, "deltaT");
			}

			numsteps_bg = 1;
			if(sim.mg_flag == FOFR)
			{
				dtau_bg = sim.fofR_epsilon_bg * sqrt(3. * FRR(Rbar, sim.fofR_params, sim.fofR_type, 3)) / a;
				if(dtau_new >= 2*dtau_bg)
				{
					numsteps_bg = (int) (0.5*dtau_new/dtau_bg);
					numsteps_bg *= 2; // So that numsteps_bg is even (dvisible by 2)
					dtau_bg = (double) (dtau_new/numsteps_bg);
				}
				else
				{
					dtau_bg = dtau_new;
				}
			}

			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau_new / dx;
			if(sim.baryon_flag)
			{
				COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau_new / dx;
			}
			if(sim.mg_flag == FOFR)
			{
				COUT << "), time step / Hubble time = " <<  Hubble * dtau_new;
				COUT << ", numsteps_bg = " << numsteps_bg << endl;
			}
			else
			{
				COUT << "), time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau_new;
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

		if(sim.mg_flag == FOFR) tau_temp = tau;

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
					maxvel[0] = pcls_cdm.updateVel(update_q, (dtau_new + dtau) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					if(sim.baryon_flag)
						maxvel[1] = pcls_b.updateVel(update_q, (dtau_new + dtau) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				}
				else
				{
					maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau_new + dtau) / 2., update_cdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
					if(sim.baryon_flag)
						maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau_new + dtau) / 2., update_b_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
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
						maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau_new + dtau) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau_new + dtau) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

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
						pcls_ncdm[i].moveParticles(update_pos, dtau_new / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau_new / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
				}
			}

			if(numsteps == 1)
			{
				if(sim.mg_flag == FOFR)
				{
					if(sim.read_bg_from_file == 1)
					{
						tau_temp += 0.5 * dtau_new;
						a = gsl_spline_eval(a_spline, tau_temp, gsl_inpl_acc);
						Hubble = gsl_spline_eval(H_spline, tau_temp, gsl_inpl_acc);
						Rbar = gsl_spline_eval(Rbar_spline, tau_temp, gsl_inpl_acc);
					}
					else if(numsteps_bg == 1)
					{
						rungekutta_fR(a, Hubble, Rbar, fourpiG, cosmo, 0.5 * dtau_bg, sim.fofR_params, sim.fofR_type);
					}
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR(a, Hubble, Rbar, fourpiG, cosmo, dtau_bg, sim.fofR_params, sim.fofR_type);
						}
					}
				}
				else
				{
					rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau_new);  // evolve background by half a time step
					Hubble = Hconf(a, fourpiG, cosmo);
					Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
				}
			}

			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			if(numsteps == 1 || j == numsteps / 2)
			{
				if(sim.gr_flag > 0)
				{
					pcls_cdm.moveParticles(update_pos, dtau_new, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
					if(sim.baryon_flag)
						pcls_b.moveParticles(update_pos, dtau_new, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
				}
				else
				{
					pcls_cdm.moveParticles(update_pos_Newton, dtau_new, NULL, 0, f_params);
					if(sim.baryon_flag)
						pcls_b.moveParticles(update_pos_Newton, dtau_new, NULL, 0, f_params);
				}

#ifdef BENCHMARK
				moveParts_count++;
				moveParts_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif
			}

			if(numsteps != 1)
			{
				if(sim.mg_flag == FOFR)
				{
					if(sim.read_bg_from_file == 1)
					{
						tau_temp += 0.5 * dtau_new / numsteps;
						a = gsl_spline_eval(a_spline, tau_temp, gsl_inpl_acc);
						Hubble = gsl_spline_eval(H_spline, tau_temp, gsl_inpl_acc);
						Rbar = gsl_spline_eval(Rbar_spline, tau_temp, gsl_inpl_acc);
					}
					else if(numsteps_bg == 1) rungekutta_fR(a, Hubble, Rbar, fourpiG, cosmo, 0.5 * dtau_bg / numsteps, sim.fofR_params, sim.fofR_type);
					else
					{
						for(g=0; g<numsteps_bg/2; g++)
						{
							rungekutta_fR(a, Hubble, Rbar, fourpiG, cosmo, dtau_bg / numsteps, sim.fofR_params, sim.fofR_type);
						}
					}
				}
				else
				{
					rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau_new / numsteps);  // evolve background by half a time step
					Hubble = Hconf(a, fourpiG, cosmo);
					Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
				}
			}

#ifdef OUTPUT_BACKGROUND
					//// write a the f(R) background file
			if(parallel.rank()==0)
			{
				ofstream fileBG;
				fileBG.open("background.txt",std::ifstream::app);
				fileBG <<setprecision(16)<< tau_temp <<" "<< a<< " " << Hconf(a, fourpiG, cosmo)<< " 1.0 1.0 1.0"<<endl;
				fileBG.close();
			}
#endif

			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			for(i = 0; i < cosmo.num_ncdm; i++)
			{
				if(numsteps_ncdm[i] == numsteps)
				{
					if(sim.gr_flag > 0)
						pcls_ncdm[i].moveParticles(update_pos, dtau_new / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau_new / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
				}
			}

			if(sim.mg_flag == FOFR)
			{
				if(sim.read_bg_from_file == 1)
				{
					tau_temp += 0.5 * dtau_new / numsteps;
					a = gsl_spline_eval(a_spline, tau_temp, gsl_inpl_acc);
					Hubble = gsl_spline_eval(H_spline, tau_temp, gsl_inpl_acc);
					Rbar = gsl_spline_eval(Rbar_spline, tau_temp, gsl_inpl_acc);
				}
				else if(numsteps_bg == 1) rungekutta_fR(a, Hubble, Rbar, fourpiG, cosmo, 0.5 * dtau_bg / numsteps, sim.fofR_params, sim.fofR_type);
				else
				{
					for(g=0; g<numsteps_bg/2; g++)
					{
						rungekutta_fR(a, Hubble, Rbar, fourpiG, cosmo, dtau_bg / numsteps, sim.fofR_params, sim.fofR_type);
					}
				}
			}
			else
			{
				rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau_new / numsteps);  // evolve background by half a time step
				Hubble = Hconf(a, fourpiG, cosmo);
				Rbar = 2. * fourpiG * ( (cosmo.Omega_cdm + cosmo.Omega_b) / a / a / a + 4.*cosmo.Omega_Lambda );
			}

			#ifdef OUTPUT_BACKGROUND
					//// write the f(R) background file
					if(parallel.rank()==0)
					{
						ofstream fileBG;
						fileBG.open("background.txt",std::ifstream::app);
						fileBG <<setprecision(16)<< tau_temp <<" "<< a<< " " << Hconf(a, fourpiG, cosmo)<< " 1.0 1.0 1.0"<<endl;
						fileBG.close();
					}
			#endif

		}   // particle update done

		parallel.max<double>(maxvel, numspecies);

		if(sim.gr_flag > 0)
		{
			for(i = 0; i < numspecies; i++)
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		}

		tau += dtau_new;

		if(tau_Lambda < 0. && (cosmo.Omega_m / a / a / a) < cosmo.Omega_Lambda)
		{
			tau_Lambda = tau;
			COUT << "matter-dark energy equality at z=" << ((1./a) - 1.) << endl;
		}

		///////

		if(sim.wallclocklimit > 0.)   // check for wallclock time limit
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
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau_new, cycle);
				}
				else
#endif
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau_new, cycle);
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
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi_check, a, tau, dtau_new, cycle, restartcount);
			}
			else
#endif
			hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, chi, Bi, a, tau, dtau_new, cycle, restartcount);
			restartcount++;
		}


		if(sim.mg_flag == FOFR)
		{
			if(sim.Cf * dx < sim.steplimit / Hubble)
				dtau_new = sim.Cf * dx;
			else
				dtau_new = sim.steplimit / Hubble;
			if(!sim.back_to_GR)
			{
				max_FRR = computeDRzeta(deltaT, deltaR, zeta, xi, Rbar, FR(Rbar, sim.fofR_params, sim.fofR_type), sim.fofR_params, sim.fofR_type, sim.fofR_target_precision);
				dtau_osci = sim.fofR_epsilon_fields * sqrt(3.*max_FRR)/a;
				if(dtau_new > dtau_osci) dtau_new = dtau_osci;
			}
		}
		else
		{
		if(sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
			dtau_new = sim.Cf * dx;
		else
			dtau_new = sim.steplimit / Hconf(a, fourpiG, cosmo);
		}

		dtau_old = dtau;
		dtau = dtau_new;

		cycle++;

#ifdef BENCHMARK
		cycle_time += MPI_Wtime()-cycle_start_time;
#endif
	}

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

	xi_output.close();

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
	COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
#endif

#ifdef EXTERNAL_IO
		ioserver.stop();
	}
#endif

	return 0;
}
