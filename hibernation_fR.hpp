//////////////////////////
// hibernation_fR.hpp
//////////////////////////
//
// Auxiliary functions for hibernation in f(R) gravity
//
//////////////////////////

#ifndef HIBERNATION_FR_HEADER
#define HIBERNATION_FR_HEADER

/////////////////////////////////////
// Restart settings for f(R) gravity
/////////////////////////////////////
/////////////////////////////////////
void writeRestartSettings_fR(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
	string bgfilename,
	const double a,
	const double Hubble,
	const double Rbar,
	const double dot_Rbar,
	const double fbar,
	const double fRbar,
	const double fRRbar,
	const double tau,
	const double dtau,
	const double dtau_old,
	const int cycle,
	const int restartcount = -1
)
{
	if(!parallel.isRoot()) return;

	ofstream outfile_ini;
	ofstream outfile_bin;
	ofstream outfile_bg;
	ifstream bgfile;

	int i;

	string filename;
	filename.reserve(PARAM_MAX_LENGTH);
	filename.assign((string) sim.restart_path + sim.basename_restart + "_");
	if(restartcount > 0)
	{
		for(i=(int)log10(restartcount); i<2; ++i)
		{
			filename += "0";
		}
		filename += to_string(restartcount);
	}
	else if(!restartcount)
	{
		filename += "000";
	}

	////////// Background //////////
	bgfile.open(bgfilename);

	if(bgfile.is_open())
	{
		outfile_bg.open(filename + "_background.dat");

		if(outfile_bg.is_open())
		{
			outfile_bg << bgfile.rdbuf();
			outfile_bg.close();
		}
		else
		{
			cout << "Error writing hibernation background file." << endl;
		}
	}
	else
	{
		cout << "Error reading background file before hibernation." << endl;
	}

	////////// Settings (binary) //////////
	outfile_bin.open(filename + ".ini.bin", ios::out | ios::trunc | ios::binary);

	if(!outfile_bin.is_open())
	{
		cout << "error opening binary restart file" << endl;
	}
	else
	{
		outfile_bin.write((char*) & a, sizeof(double));
		outfile_bin.write((char*) & Hubble, sizeof(double));
		outfile_bin.write((char*) & Rbar, sizeof(double));
		outfile_bin.write((char*) & dot_Rbar, sizeof(double));
		outfile_bin.write((char*) & fbar, sizeof(double));
		outfile_bin.write((char*) & fRbar, sizeof(double));
		outfile_bin.write((char*) & fRRbar, sizeof(double));
		outfile_bin.write((char*) & tau, sizeof(double));
		outfile_bin.write((char*) & dtau, sizeof(double));
		outfile_bin.write((char*) & dtau_old, sizeof(double));
		outfile_bin.write((char*) & sim, sizeof(metadata));
		outfile_bin.write((char*) & cosmo, sizeof(cosmology));
		outfile_bin.close();
	}

	////////// Settings (text) //////////
	outfile_ini.open(filename + ".ini", ios::out | ios::trunc);

	if(!outfile_ini.is_open())
	{
		cout << " error opening file for restart settings!" << endl;
	}
	else
	{
		outfile_ini << "# automatically generated settings for restart after hibernation ";
		if(restartcount < 0)
		{
			outfile_ini << "due to wallclock limit ";
		}
		else
		{
			outfile_ini << "requested ";
		}

		outfile_ini << setprecision(15);

		outfile_ini << "at redshift z = " << (1./a)-1. << endl << endl;
		outfile_ini << "fRevolution version = " << FREVOLUTION_VERSION << endl;
		outfile_ini << endl << "#====== IC generation ======#" << endl;
		outfile_ini << "IC generator        = restart" << endl;
		outfile_ini << "restart count       = " << restartcount << endl;

		if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
		{
			outfile_ini << "particle file       = " << filename << "_cdm.h5";
			if(sim.baryon_flag)
			{
				outfile_ini << ", " << filename << "_b.h5";
			}
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				outfile_ini << ", " << filename << "_ncdm" << i << ".h5";
			}
		}
		else
		{
			outfile_ini << "particle file       = " << filename << "_cdm";
			if(sim.baryon_flag)
			{
				outfile_ini << ", " << filename << "_b";
			}
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				outfile_ini << ", " << filename << "_ncdm" << i;
			}
		}

		outfile_ini << endl;

		outfile_ini << "metric file         = " << filename << "_phi.h5";
		outfile_ini << ", " << filename << "_chi.h5";
		outfile_ini << ", " << filename << "_B.h5";
		outfile_ini << ", " << filename << "_xi.h5";
		outfile_ini << ", " << filename << "_xi_old.h5";
		outfile_ini << ", " << filename << "_deltaR.h5";

		outfile_ini << endl;

		// TODO: check what is necessary
		outfile_ini << "seed                = " << ic.seed << endl;
		outfile_ini << "restart redshift    = " << (1./a) - 1. << endl;
		outfile_ini << "cycle               = " << cycle << endl;
		outfile_ini << "tau                 = " << tau << endl;
		outfile_ini << "dtau                = " << dtau << endl;
		outfile_ini << "dtau_old            = " << dtau_old << endl;
		outfile_ini << "a                 	= " << a << endl;
		outfile_ini << "Hubble            	= " << Hubble << endl;
		outfile_ini << "Rbar                = " << Rbar << endl;
		outfile_ini << "dot_Rbar            = " << dot_Rbar << endl;

		outfile_ini << endl;

		if(ic.flags & ICFLAG_KSPHERE)
		{
			outfile_ini << "k-domain            = sphere" << endl;
		}
		else
		{
			outfile_ini << "k-domain            = cube" << endl;
		}

		outfile_ini << endl << "#====== primordial power spectrum ======#" << endl;
		outfile_ini << "k_pivot = " << ic.k_pivot << endl;
		outfile_ini << "A_s     = " << ic.A_s << endl;
		outfile_ini << "n_s     = " << ic.n_s << endl;

		outfile_ini << endl << "#====== cosmological parameters ======#" << endl;
		outfile_ini << "h         = " << cosmo.h << endl;
		outfile_ini << "Omega_cdm = " << cosmo.Omega_cdm << endl;
		outfile_ini << "Omega_b   = " << cosmo.Omega_b << endl;
		outfile_ini << "Omega_g   = " << cosmo.Omega_g << endl;
		outfile_ini << "Omega_ur  = " << cosmo.Omega_ur << endl;
		outfile_ini << "N_ncdm    = " << cosmo.num_ncdm << endl;
		if(cosmo.num_ncdm > 0)
		{
			outfile_ini << "m_cdm     = ";
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				outfile_ini << cosmo.m_ncdm[i] << ", ";
			}
			outfile_ini << cosmo.m_ncdm[i] << endl;

			outfile_ini << "T_cdm     = ";
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				outfile_ini << cosmo.T_ncdm[i] << ", ";
			}
			outfile_ini << cosmo.T_ncdm[i] << endl;

			outfile_ini << "deg_cdm   = ";
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				outfile_ini << cosmo.deg_ncdm[i] << ", ";
			}
			outfile_ini << cosmo.deg_ncdm[i] << endl;
		}

		outfile_ini << endl << "#====== simulation settings ======#" << endl;
		if(sim.baryon_flag > 0)
		{
			outfile_ini << "baryon treatment    = sample" << endl;
		}

		if(sim.radiation_flag > 0)
		{
			outfile_ini << "radiation treatment = CLASS" << endl;
			outfile_ini << "switch delta_rad    = " << sim.z_switch_deltarad;
			if(cosmo.num_ncdm > 0)
			{
				outfile_ini << "switch delta_ncdm   = ";
				for(i=0; i<cosmo.num_ncdm-1; i++)
				{
					outfile_ini << sim.z_switch_deltancdm[i] << ", ";
				}
				outfile_ini << sim.z_switch_deltancdm[i] << endl;
			}
			outfile_ini << "switch linear chi   = " << sim.z_switch_linearchi << endl;
		}

		if(sim.vector_flag == VECTOR_ELLIPTIC)
		{
			outfile_ini << "vector method       = elliptic" << endl;
		}
		else
		{
			outfile_ini << "vector method       = parabolic" << endl;
		}

		outfile_ini << "initial redshift    = " << sim.z_in << endl;
		outfile_ini << "boxsize             = " << sim.boxsize << endl;
		outfile_ini << "Ngrid               = " << sim.numpts << endl;
		outfile_ini << "tiling factor       = " << ic.numtile[0] << endl;
		outfile_ini << "Courant factor      = " << sim.Cf << endl;
		outfile_ini << "time step limit     = " << sim.steplimit << endl;
		if(cosmo.num_ncdm > 0)
		{
			outfile_ini << "move limit          = " << sim.movelimit << endl;
		}
		outfile_ini << "check fields        = " << sim.check_fields << endl;
		outfile_ini << "CYCLE_INFO_INTERVAL = " << sim.CYCLE_INFO_INTERVAL << endl;
		outfile_ini << "gravity theory      = fr" << endl;

		outfile_ini << endl << "#====== f(R) settings ======#" << endl;
		if(sim.fR_model == FR_MODEL_RN || sim.fR_model == FR_MODEL_R2)
		{
			// TODO: Check that comments are well written -- need something for FR_MODEL_DELTA too?
			if(sim.fR_model == FR_MODEL_R2)
			{
				outfile_ini << "# WARNING: f(R) parameters for R + R^2 model hav already been rescaled. See fR_tools.hpp for additional information." << endl;
			}
			outfile_ini << "f(R) model                  = RN" << endl;
			outfile_ini << "f(R) parameters             = " << sim.fR_params[0] << ", " << sim.fR_params[1];
		}
		else if(sim.fR_model == FR_MODEL_DELTA)
		{
			outfile_ini << "f(R) model                  = DE" << endl;
			outfile_ini << "f(R) parameters             = " << sim.fR_params[0] << ", " << sim.fR_params[1];
		}
		else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
		{
			outfile_ini << "# WARNING: f(R) parameters for Hu-Sawicki model have already been rescaled. See fR_tools.hpp for additional information." << endl;
			outfile_ini << "f(R) model                  = HS" << endl;
			outfile_ini << "f(R) parameters             = " << sim.fR_params[0] << ", " << sim.fR_params[1] << ", " << sim.fR_params[2] << ", " << sim.fR_params[3];
		}
		else
		{
			COUT << " error f(R) model not recognized!";
		}

		 outfile_ini << endl;

		if(sim.relativistic_flag)
		{
			outfile_ini << "Newtonian f(R)              = 0" << endl;
		}
		else
		{
			outfile_ini << "Newtonian f(R)              = 1" << endl;
		}
		outfile_ini << "f(R) epsilon background     = " << sim.fR_epsilon_bg << endl;
		outfile_ini << "f(R) count max              = " << sim.fR_count_max << endl;
		outfile_ini << "f(R) epsilon fields         = " << sim.fR_epsilon_fields << endl;
		outfile_ini << "f(R) target precision       = " << sim.fR_target_precision << endl;
		outfile_ini << "lcdm background             = " << sim.lcdm_background << endl;

		outfile_ini << endl << "#====== Multigrid and relaxation ======#" << endl;
		outfile_ini << "relaxation method          = " << sim.relaxation_method << endl;
		outfile_ini << "relaxation error           = " << sim.relaxation_error << endl;
		outfile_ini << "red black                  = " << sim.red_black << endl;
		outfile_ini << "overrelaxation factor      = " << sim.overrelaxation_factor << endl;
		outfile_ini << "pre-smoothing              = " << sim.pre_smoothing << endl;
		outfile_ini << "post-smoothing             = " << sim.post_smoothing << endl;
		outfile_ini << "multigrid n-grids          = " << sim.multigrid_n_grids << endl;
		outfile_ini << "multigrid n-cycles         = " << sim.multigrid_n_cycles << endl;
		outfile_ini << "multigrid damping          = " << sim.multigrid_damping << endl;
		outfile_ini << "check shape                = " << sim.multigrid_check_shape << endl;
		if(sim.multigrid_shape == MULTIGRID_SHAPE_V)
		{
			outfile_ini << "multigrid shape            = V" << endl;
		}
		else if(sim.multigrid_shape == MULTIGRID_SHAPE_W)
		{
			outfile_ini << "multigrid shape            = W" << endl;
		}
		if(sim.multigrid_restrict_mode == RESTRICT_XI)
		{
			outfile_ini << "restrict mode              = xi" << endl;
		}
		else
		{
			outfile_ini << "restrict mode              = deltaR" << endl;
		}

		outfile_ini << endl << "#====== output ======#" << endl;
		outfile_ini << "output path         = " << sim.output_path << endl;
		outfile_ini << "generic file base   = " << sim.basename_generic << endl;
		outfile_ini << "snapshot file base  = " << sim.basename_snapshot << endl;
		outfile_ini << "Pk file base        = " << sim.basename_pk << endl;
		outfile_ini << "lightcone file base = " << sim.basename_lightcone << endl;
		if(sim.num_snapshot > 0)
		{
			outfile_ini << "snapshot redshifts = ";
			for(i=0; i<sim.num_snapshot - 1; i++)
			{
				outfile_ini << sim.z_snapshot[i] << ", ";
			}
			outfile_ini << sim.z_snapshot[i] << endl;
		}

		if(sim.out_snapshot)
		{
			outfile_ini << "snapshot outputs   = ";
			if(sim.out_snapshot & MASK_PHI)
			{
				outfile_ini << "phi";
				if(sim.out_snapshot > MASK_CHI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_CHI)
			{
				outfile_ini << "chi";
				if(sim.out_snapshot > MASK_POT)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_POT)
			{
				outfile_ini << "psiN";
				if(sim.out_snapshot > MASK_B)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_B)
			{
				outfile_ini << "B";
				if(sim.out_snapshot > MASK_T00)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_T00)
			{
				outfile_ini << "T00";
				if(sim.out_snapshot > MASK_TIJ)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_TIJ)
			{
				outfile_ini << "Tij";
				if(sim.out_snapshot > MASK_RBARE)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_RBARE)
			{
				outfile_ini << "rhoN";
				if(sim.out_snapshot > MASK_HIJ)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_HIJ)
			{
				outfile_ini << "hij";
				if(sim.out_snapshot > MASK_P)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_P)
			{
				outfile_ini << "p";
				if(sim.out_snapshot > MASK_GADGET)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_GADGET)
			{
				outfile_ini << "Gadget2";
				if(sim.out_snapshot > MASK_PCLS)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_PCLS)
			{
				outfile_ini << "particles";
				if(sim.out_snapshot > MASK_DELTA)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_DELTA)
			{
				outfile_ini << "delta";
				if(sim.out_snapshot > MASK_DBARE)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_DBARE)
			{
				outfile_ini << "deltaN";
				if(sim.out_snapshot > MASK_XI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_XI)
			{
				outfile_ini << "xi";
				if(sim.out_snapshot > MASK_ZETA)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_ZETA)
			{
				outfile_ini << "zeta";
				if(sim.out_snapshot > MASK_DELTAT)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_DELTAT)
			{
				outfile_ini << "eightpiG_deltaT";
				if(sim.out_snapshot > MASK_DELTAR)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_DELTAR)
			{
				outfile_ini << "deltaR";
				if(sim.out_snapshot > MASK_LAPLACE_XI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_snapshot & MASK_LAPLACE_XI)
			{
				outfile_ini << "laplace_xi";
			}
			outfile_ini << endl;
		}

		if(sim.out_snapshot & MASK_GADGET)
		{
			outfile_ini << "tracer factor      = ", sim.tracer_factor[0];
			for(i=1; i<=sim.baryon_flag + cosmo.num_ncdm; i++)
			{
				outfile_ini << ", " << sim.tracer_factor[i];
			}
			outfile_ini << endl;
		}

		// OUTPUTS FOR check_fields()
		if(sim.out_check)
		{
			outfile_ini << "check outputs       = ";
			if(sim.out_check & MASK_PHI)
			{
				outfile_ini << "phi";
				if(sim.out_check > MASK_CHI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_CHI)
			{
				outfile_ini << "chi";
				if(sim.out_check > MASK_POT)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_POT)
			{
				outfile_ini << "psiN";
				if(sim.out_check > MASK_B)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_B)
			{
				outfile_ini << "B";
				if(sim.out_check > MASK_T00)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_T00)
			{
				outfile_ini << "T00";
				if(sim.out_check > MASK_TIJ)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_TIJ)
			{
				outfile_ini << "Tij";
				if(sim.out_check > MASK_RBARE)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_RBARE)
			{
				outfile_ini << "rhoN";
				if(sim.out_check > MASK_HIJ)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_HIJ)
			{
				outfile_ini << "hij";
				if(sim.out_check > MASK_P)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_P)
			{
				outfile_ini << "p";
				if(sim.out_check > MASK_GADGET)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_GADGET)
			{
				outfile_ini << "Gadget2";
				if(sim.out_check > MASK_PCLS)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_PCLS)
			{
				outfile_ini << "particles";
				if(sim.out_check > MASK_DELTA)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_DELTA)
			{
				outfile_ini << "delta";
				if(sim.out_check > MASK_DBARE)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_DBARE)
			{
				outfile_ini << "deltaN";
				if(sim.out_check > MASK_XI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_XI)
			{
				outfile_ini << "xi";
				if(sim.out_check > MASK_ZETA)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_ZETA)
			{
				outfile_ini << "zeta";
				if(sim.out_check > MASK_DELTAT)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_DELTAT)
			{
				outfile_ini << "eightpiG_deltaT";
				if(sim.out_check > MASK_DELTAR)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_DELTAR)
			{
				outfile_ini << "deltaR";
				if(sim.out_check > MASK_LAPLACE_XI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_check & MASK_LAPLACE_XI)
			{
				outfile_ini << "laplace_xi";
			}
			outfile_ini << endl;
		}

		if(sim.downgrade_factor > 1)
		{
			outfile_ini << "downgrade factor   = " << sim.downgrade_factor << endl;
		}

		if(sim.num_pk > 0)
		{
			outfile_ini << "Pk redshifts        = ";
			for(i=0; i<sim.num_pk-1; i++)
			{
				outfile_ini <<  sim.z_pk[i] << ", ";
			}
			outfile_ini << sim.z_pk[i] << endl;
		}

		if(sim.out_pk)
		{
			outfile_ini << "Pk outputs          = ";
			if(sim.out_pk & MASK_PHI)
			{
				outfile_ini << "phi";
				if(sim.out_pk > MASK_CHI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_CHI)
			{
				outfile_ini << "chi";
				if(sim.out_pk > MASK_POT)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_POT)
			{
				outfile_ini << "psiN";
				if(sim.out_pk > MASK_B)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_B)
			{
				outfile_ini << "B";
				if(sim.out_pk > MASK_T00)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_T00)
			{
				outfile_ini << "T00";
				if(sim.out_pk > MASK_TIJ)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_TIJ)
			{
				outfile_ini << "Tij";
				if(sim.out_pk > MASK_RBARE)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_RBARE)
			{
				outfile_ini << "rhoN";
				if(sim.out_pk > MASK_HIJ)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_HIJ)
			{
				outfile_ini << "hij";
				if(sim.out_pk > MASK_P)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_P)
			{
				outfile_ini << "p";
				if(sim.out_pk > MASK_XSPEC)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_XSPEC)
			{
				outfile_ini << "X-spectra";
				if(sim.out_pk > MASK_DELTA)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_DELTA)
			{
				outfile_ini << "delta";
				if(sim.out_pk > MASK_DBARE)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_DBARE)
			{
				outfile_ini << "deltaN";
				if(sim.out_pk > MASK_XI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_XI)
			{
				outfile_ini << "xi";
				if(sim.out_pk > MASK_ZETA)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_ZETA)
			{
				outfile_ini << "zeta";
				if(sim.out_pk > MASK_DELTAT)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_DELTAT)
			{
				outfile_ini << "eightpiG_deltaT";
				if(sim.out_pk > MASK_DELTAR)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_DELTAR)
			{
				outfile_ini << "deltaR";
				if(sim.out_pk > MASK_LAPLACE_XI)
				{
					outfile_ini << ", ";
				}
			}

			if(sim.out_pk & MASK_LAPLACE_XI)
			{
				outfile_ini << "laplace_xi";
			}
			outfile_ini << endl;
		}

		outfile_ini << "Pk bins            = " << sim.numbins << endl;

		if(sim.num_lightcone == 1)
		{
			outfile_ini << "lightcone vertex    = " << sim.lightcone[0].vertex[0] << ", " << sim.lightcone[0].vertex[1] << ", " << sim.lightcone[0].vertex[2] << endl;
			outfile_ini << "lightcone outputs   = ";
			if(sim.out_lightcone[0] & MASK_PHI)
			{
				outfile_ini << "phi";
				if(sim.out_lightcone[0] > MASK_CHI)
				{
					outfile_ini << ", ";
				}
			}
			if(sim.out_lightcone[0] & MASK_CHI)
			{
				outfile_ini << "chi";
				if(sim.out_lightcone[0] > MASK_POT)
				{
					outfile_ini << ", ";
				}
			}
			if(sim.out_lightcone[0] & MASK_B)
			{
				outfile_ini << "B";
				if(sim.out_lightcone[0] > MASK_T00)
				{
					outfile_ini << ", ";
				}
			}
			if(sim.out_lightcone[0] & MASK_HIJ)
			{
				outfile_ini << "hij";
				if(sim.out_lightcone[0] > MASK_P)
				{
					outfile_ini << ", ";
				}
			}
			if(sim.out_lightcone[0] & MASK_GADGET)
			{
				outfile_ini << "Gadget2";
			}
			outfile_ini << endl;
			if(sim.lightcone[0].opening > -1.)
			{
				outfile_ini << "lightcone opening half-angle = " << acos(sim.lightcone[0].opening) * 180. / M_PI << endl;
			}
			outfile_ini << "lightcone distance  = " << sim.lightcone[0].distance[1] << ", " << sim.lightcone[0].distance[0] << endl;

			if(sim.lightcone[0].z != 0)
			{
				outfile_ini << "lightcone redshift  = " << sim.lightcone[0].z << endl;
			}
			outfile_ini << "lightcone direction = " << sim.lightcone[0].direction[0] << ", " << sim.lightcone[0].direction[1] << ", " << sim.lightcone[0].direction[2] << endl;
			outfile_ini << "lightcone covering  = " << sim.covering[0] << endl;
			if(sim.Nside[0][0] != sim.Nside[0][1])
			{
				outfile_ini << "lightcone Nside     = " << sim.Nside[0][0] << ", " << sim.Nside[0][1] << endl;
			}
			else
			{
				outfile_ini << "lightcone Nside     = " << sim.Nside[0][0] << endl;
			}
			outfile_ini << "lightcone pixel factor = " << sim.pixelfactor[0] << endl;
			outfile_ini << "lightcone shell factor = " << sim.shellfactor[0] << endl;
		}
		else if(sim.num_lightcone > 1)
		{
			for(i=0; i<sim.num_lightcone; i++)
			{
				outfile_ini << "lightcone " << i << " vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
				outfile_ini << "lightcone " << i << " outputs = ";
				if(sim.out_lightcone[i] & MASK_PHI)
				{
					outfile_ini << "phi";
					if(sim.out_lightcone[i] > MASK_CHI)
					{
						outfile_ini << ", ";
					}
				}
				if(sim.out_lightcone[i] & MASK_CHI)
				{
					outfile_ini << "chi";
					if(sim.out_lightcone[i] > MASK_POT)
					{
						outfile_ini << ", ";
					}
				}
				if(sim.out_lightcone[i] & MASK_B)
				{
					outfile_ini << "B";
					if(sim.out_lightcone[i] > MASK_T00)
					{
						outfile_ini << ", ";
					}
				}
				if(sim.out_lightcone[i] & MASK_HIJ)
				{
					outfile_ini << "hij";
					if(sim.out_lightcone[i] > MASK_P)
					{
						outfile_ini << ", ";
					}
				}
				if(sim.out_lightcone[i] & MASK_GADGET)
				{
					outfile_ini << "Gadget2";
				}
				outfile_ini << "" << endl;
				if(sim.lightcone[i].opening > -1.)
				{
					outfile_ini << "lightcone " << i << " opening half-angle = " << acos(sim.lightcone[i].opening) * 180. / M_PI << endl;
				}
				outfile_ini << "lightcone " << i << " distance  = " << sim.lightcone[i].distance[1] << ", " << sim.lightcone[i].distance[0] << endl;
				if(sim.lightcone[i].z != 0)
				{
					outfile_ini << "lightcone " << i << " redshift  = " << sim.lightcone[i].z << endl;
				}
				outfile_ini << "lightcone " << i << " direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
				outfile_ini << "lightcone " << i << " covering  = " << sim.covering[0] << endl;
				if(sim.Nside[0][0] != sim.Nside[0][1])
				{
					outfile_ini << "lightcone " << i << " Nside     = " << sim.Nside[0][0] << ", " << sim.Nside[0][1] << endl;
				}
				else
				{
					outfile_ini << "lightcone " << i << " Nside     = " << sim.Nside[0][0] << endl;
				}
				outfile_ini << "lightcone " << i << " pixel factor = " << sim.pixelfactor[0] << endl;
				outfile_ini << "lightcone " << i << " shell factor = " << sim.shellfactor[0] << endl;
			}
		}

		outfile_ini << endl << "#====== hibernation ======#" << endl;
		if(sim.num_restart > 0)
		{
			outfile_ini << "hibernation redshifts       = ";
			for(i=0; i<sim.num_restart - 1; i++)
			{
				outfile_ini << sim.z_restart[i] << ", ";
			}
			outfile_ini << sim.z_restart[i] << endl;
		}

		if(sim.wallclocklimit > 0.)
		{
			outfile_ini << "hibernation wallclock limit = " << sim.wallclocklimit << endl;
		}

		if(sim.restart_path[0] != '\0')
		{
			outfile_ini << "hibernation path            = " << sim.restart_path << endl;
		}

		outfile_ini << "hibernation file base       = " << sim.basename_restart << endl;

		if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
		{
			outfile_ini << "particle save mode          = hdf5" << endl;
		}
		else if(sim.hibernation_save_mode == HIB_SAVE_GADGET2)
		{
			outfile_ini << "particle save mode          = gadget2" << endl;
		}

		outfile_ini.close();
	}
}

/////////////////////////////
// Hibernate for f(R) gravity
/////////////////////////////
void hibernate_fR(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
	string bgfilename,
	gadget2_header & hdr,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm,
	Field<Real> & phi,
	Field<Real> & chi,
	Field<Real> & Bi,
	Field<Real> & xi,
	Field<Real> & xi_old,
	Field<Real> & deltaR,
	const double a,
	const double Hubble,
	const double Rbar,
	const double dot_Rbar,
	const double fbar,
	const double fRbar,
	const double fRRbar,
	const double tau,
	const double dtau,
	const double dtau_old,
	const int cycle,
	const int restartcount = -1
)
{
	int i;
	Site x(Bi.lattice());
	string filename;

	filename.reserve(PARAM_MAX_LENGTH);
	filename.assign((string) sim.restart_path + sim.basename_restart + "_");
	if(restartcount > 0)
	{
		for(i=(int)log10(restartcount); i<2; ++i)
		{
			filename += "0";
		}
		filename += to_string(restartcount);
	}
	else if(!restartcount)
	{
		filename += "000";
	}

	writeRestartSettings_fR(sim, ic, cosmo, bgfilename, a, Hubble, Rbar, dot_Rbar, fbar, fRbar, fRRbar, tau, dtau, dtau_old, cycle, restartcount);

	for(x.first(); x.test(); x.next())
	{
		Bi(x,0) /= a * a * sim.numpts;
		Bi(x,1) /= a * a * sim.numpts;
		Bi(x,2) /= a * a * sim.numpts;
	}

#ifdef EXTERNAL_IO
	while(ioserver.openOstream() == OSTREAM_FAIL);

	phi.saveHDF5_server_open(filename + "_phi";
	chi.saveHDF5_server_open(filename + "_chi";
	xi.saveHDF5_server_open(filename + "_xi";
	xi_old.saveHDF5_server_open(filename + "_xi_old";
	deltaR.saveHDF5_server_open(filename + "_deltaR";
	Bi.saveHDF5_server_open(filename + "_B";

	if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
	{
		pcls_cdm->saveHDF5_server_open(filename + "_cdm");
		pcls_cdm->saveHDF5_server_write();

		if(sim.baryon_flag)
		{
			pcls_b->saveHDF5_server_open(filename + "_b");
			pcls_b->saveHDF5_server_write();
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			pcls_ncdm[i].saveHDF5_server_open((string) filename + "_ncdm_" + to_string(i));
			pcls_ncdm[i].saveHDF5_server_write();
		}
	}
	else
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
		pcls_cdm->saveGadget2(filename + "_cdm", hdr, sim.tracer_factor[0]);

		if(sim.baryon_flag)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[1] / sim.tracer_factor[1]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			pcls_b->saveGadget2(filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for (i=0; i<cosmo.num_ncdm; i++)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			pcls_ncdm[i].saveGadget2((string) filename + "_ncdm_" + to_string(i), hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
		}
	}

	phi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	chi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	xi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	xi_old.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	deltaR.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	Bi.saveHDF5_server_write(NUMBER_OF_IO_FILES);

	ioserver.closeOstream();

#else

	if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
	{
		pcls_cdm->saveHDF5(filename + "_cdm", 1);
		if(sim.baryon_flag)
		{
			pcls_b->saveHDF5(filename + "_b", 1);
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			pcls_ncdm[i].saveHDF5((string) filename + "_ncdm_" + to_string(i), 1);
		}
	}
	else
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
		pcls_cdm->saveGadget2(filename + "_cdm", hdr, sim.tracer_factor[0]);

		if(sim.baryon_flag)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[1] / sim.tracer_factor[1]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			pcls_b->saveGadget2(filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			pcls_ncdm[i].saveGadget2((string) filename + "_ncdm_" + to_string(i), hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
		}
	}

	COUT << "Saving fields..." << endl;
	phi.saveHDF5(filename + "_phi.h5");
	chi.saveHDF5(filename + "_chi.h5");
	xi.saveHDF5(filename + "_xi.h5");
	xi_old.saveHDF5(filename + "_xi_old.h5");
	deltaR.saveHDF5(filename + "_deltaR.h5");
	Bi.saveHDF5(filename + "_B.h5");

#endif
}

#endif
