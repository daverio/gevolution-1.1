//////////////////////////
// hibernation.hpp
//////////////////////////
//
// Auxiliary functions for hibernation
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: February 2019
//
//////////////////////////

#ifndef HIBERNATION_HEADER
#define HIBERNATION_HEADER

//////////////////////////
// writeRestartSettings
//////////////////////////
// Description:
//   writes a settings file containing all the relevant metadata for restarting
//   a run from a hibernation point
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   a              scale factor
//   tau            conformal coordinate time
//   dtau           time step
//   cycle          current main control loop cycle count
//   restartcount   restart counter aka number of hibernation point (default -1)
//                  if < 0 no number is associated to the hibernation point
//
// Returns:
//
//////////////////////////

void writeRestartSettings(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
	const double a,
	const double tau,
	const double dtau,
	const int cycle,
	const int restartcount = -1
)
{
	char buffer[2*PARAM_MAX_LENGTH+24];
	FILE * outfile;
	int i;

	if(!parallel.isRoot()) return;

	if(restartcount >= 0)
	{
		sprintf(buffer, "%s%s%03d.ini", sim.restart_path, sim.basename_restart, restartcount);
	}
	else
	{
		sprintf(buffer, "%s%s.ini", sim.restart_path, sim.basename_restart);
	}
	outfile = fopen(buffer, "w");
	if(outfile == NULL)
	{
		cout << " error opening file for restart settings!" << endl;
	}
	else
	{
		fprintf(outfile, "# automatically generated settings for restart after hibernation ");
		if(restartcount < 0)
		{
			fprintf(outfile, "due to wallclock limit ");
		}
		else
		{
			fprintf(outfile, "requested ");
		}

		fprintf(outfile, "at redshift z=%f\n\n", (1./a)-1.);
		fprintf(outfile, "# info related to IC generation\n\n");
		fprintf(outfile, "IC generator       = restart\n");
		if(restartcount >= 0)
		{
			sprintf(buffer, "%03d", restartcount);
		}
		else
		{
			buffer[0] = '\0';
		}

		fprintf(outfile, "particle file      = %s%s%s_cdm.h5", sim.restart_path, sim.basename_restart, buffer);
		if(sim.baryon_flag)
		{
			fprintf(outfile, ", %s%s%s_b.h5", sim.restart_path, sim.basename_restart, buffer);
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			if(sim.numpcl[1+sim.baryon_flag+i] < 1)
			{
				fprintf(outfile, ", /dev/null");
			}
			else
			{
				fprintf(outfile, ", %s%s%s_ncdm%d.h5", sim.restart_path, sim.basename_restart, buffer, i);
			}
		}

		fprintf(outfile, "\n");
		if(sim.relativistic_flag > 0 || sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			fprintf(outfile, "metric file        = %s%s%s_phi.h5", sim.restart_path, sim.basename_restart, buffer);
			fprintf(outfile, ", %s%s%s_chi.h5", sim.restart_path, sim.basename_restart, buffer);
			if(sim.vector_flag == VECTOR_PARABOLIC)
			{
				fprintf(outfile, ", %s%s%s_B.h5\n", sim.restart_path, sim.basename_restart, buffer);
			}
			else
			{
#ifdef CHECK_B
				fprintf(outfile, ", %s%s%s_B_check.h5\n", sim.restart_path, sim.basename_restart, buffer);
#else
				fprintf(outfile, "\n");
#endif
			}
		}
		else if(sim.vector_flag == VECTOR_PARABOLIC)
		{
			fprintf(outfile, "metric file        = %s%s%s_B.h5\n", sim.restart_path, sim.basename_restart, buffer);
		}
#ifdef CHECK_B
		else
		{
			fprintf(outfile, "metric file        = %s%s%s_B_check.h5\n", sim.restart_path, sim.basename_restart, buffer);
		}
#endif

		fprintf(outfile, "restart redshift   = %.15lf\n", (1./a) - 1.);
		fprintf(outfile, "cycle              = %d\n", cycle);
		fprintf(outfile, "tau                = %.15le\n", tau);
		fprintf(outfile, "dtau               = %.15le\n", dtau);
		fprintf(outfile, "fRevolution version = %g\n\n", FREVOLUTION_VERSION);
		fprintf(outfile, "seed               = %d\n", ic.seed);
		if(ic.flags & ICFLAG_KSPHERE)
		{
			fprintf(outfile, "k-domain           = sphere\n");
		}
		else
		{
			fprintf(outfile, "k-domain           = cube\n");
		}
		fprintf(outfile, "\n\n# primordial power spectrum\n\n");
		fprintf(outfile, "k_pivot = %lg\n", ic.k_pivot);
		fprintf(outfile, "A_s     = %lg\n", ic.A_s);
		fprintf(outfile, "n_s     = %lg\n", ic.n_s);
		fprintf(outfile, "\n\n# cosmological parameters\n\n");
		fprintf(outfile, "h         = %lg\n", cosmo.h);
		fprintf(outfile, "Omega_cdm = %.15le\n", cosmo.Omega_cdm);
		fprintf(outfile, "Omega_b   = %.15le\n", cosmo.Omega_b);
		fprintf(outfile, "Omega_g   = %.15le\n", cosmo.Omega_g);
		fprintf(outfile, "Omega_ur  = %.15le\n", cosmo.Omega_ur);
		if(cosmo.Omega_fld > 0.)
		{
			fprintf(outfile, "Omega_fld = %.15le\n", cosmo.Omega_fld);
			fprintf(outfile, "w0_fld    = %lg\n", cosmo.w0_fld);
			fprintf(outfile, "wa_fld    = %lg\n", cosmo.wa_fld);
			if(sim.fluid_flag > 0)
			{
				fprintf(outfile, "cs2_fld   = %lg\n", cosmo.cs2_fld);
			}
		}
		fprintf(outfile, "N_ncdm    = %d\n", cosmo.num_ncdm);
		if(cosmo.num_ncdm > 0)
		{
			fprintf(outfile, "m_cdm     = ");
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				fprintf(outfile, "%9lf, ", cosmo.m_ncdm[i]);
			}
			fprintf(outfile, "%9lf\n", cosmo.m_ncdm[i]);
			fprintf(outfile, "T_cdm     = ");
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				fprintf(outfile, "%9lf, ", cosmo.T_ncdm[i]);
			}
			fprintf(outfile, "%9lf\n", cosmo.T_ncdm[i]);
			fprintf(outfile, "deg_cdm   = ");
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				fprintf(outfile, "%lf, ", cosmo.deg_ncdm[i]);
			}
			fprintf(outfile, "%lg\n", cosmo.deg_ncdm[i]);
		}
		fprintf(outfile, "\n\n# simulation settings\n\n");
		if(sim.baryon_flag > 0)
		{
			fprintf(outfile, "baryon treatment    = sample\n");
		}
		if(sim.radiation_flag > 0)
		{
			fprintf(outfile, "radiation treatment = CLASS\n");
			fprintf(outfile, "switch delta_rad    = %lf\n", sim.z_switch_deltarad);
			if(cosmo.num_ncdm > 0)
			{
				fprintf(outfile, "switch delta_ncdm   = ");
				for(i=0; i<cosmo.num_ncdm - 1; i++)
				{
					fprintf(outfile, "%lf, ", sim.z_switch_deltancdm[i]);
				}
				fprintf(outfile, "%lf\n", sim.z_switch_deltancdm[i]);
			}
			fprintf(outfile, "switch linear chi   = %lf\n", sim.z_switch_linearchi);
		}
		if(sim.fluid_flag > 0)
		{
			fprintf(outfile, "fluid treatment     = CLASS\n");
		}
		if(sim.modified_gravity_flag == MODIFIED_GRAVITY_FLAG_FR)
		{
			fprintf(outfile, "gravity theory      = fR\n");
			if(sim.relativistic_flag > 0)
			{
				fprintf(outfile, "Newtonian f(R)      = 0\n");
			}
			else
			{
				fprintf(outfile, "Newtonian f(R)      = 1\n");
			}
		}
		else if(sim.relativistic_flag > 0)
		{
			fprintf(outfile, "gravity theory      = GR\n");
		}
		else
		{
			fprintf(outfile, "gravity theory      = Newton\n");
		}

		if(sim.vector_flag == VECTOR_ELLIPTIC)
		{
			fprintf(outfile, "vector method       = elliptic\n");
		}
		else
		{
			fprintf(outfile, "vector method       = parabolic\n");
		}
		fprintf(outfile, "\ninitial redshift    = %lg\n", sim.z_in);
		fprintf(outfile, "boxsize             = %lg\n", sim.boxsize);
		fprintf(outfile, "Ngrid               = %d\n", sim.numpts);
		fprintf(outfile, "Courant factor      = %lg\n", sim.Cf);
		fprintf(outfile, "time step limit     = %lg\n", sim.steplimit);
		if(cosmo.num_ncdm > 0)
		{
			fprintf(outfile, "move limit          = %lg\n", sim.movelimit);
		}
		fprintf(outfile, "\n\n# output\n\n");
		fprintf(outfile, "output path         = %s\n", sim.output_path);
		fprintf(outfile, "generic file base   = %s\n", sim.basename_generic);
		fprintf(outfile, "snapshot file base  = %s\n", sim.basename_snapshot);
		fprintf(outfile, "Pk file base        = %s\n", sim.basename_pk);
		fprintf(outfile, "lightcone file base = %s\n", sim.basename_lightcone);
		if(sim.num_snapshot > 0)
		{
			fprintf(outfile, "snapshot redshifts  = ");
			for(i=0; i<sim.num_snapshot - 1; i++)
			{
				fprintf(outfile, "%lg, ", sim.z_snapshot[i]);
			}
			fprintf(outfile, "%lg\n", sim.z_snapshot[i]);
		}
		if(sim.out_snapshot)
		{
			fprintf(outfile, "snapshot outputs    = ");
			if(sim.out_snapshot & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_snapshot > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_snapshot > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if(sim.out_snapshot > MASK_B)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_snapshot > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_T00)
			{
				fprintf(outfile, "T00");
				if(sim.out_snapshot > MASK_TIJ)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if(sim.out_snapshot > MASK_RBARE)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if(sim.out_snapshot > MASK_HIJ)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_snapshot > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_P)
			{
				fprintf(outfile, "p");
				if(sim.out_snapshot > MASK_GADGET)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_GADGET)
			{
				if(sim.out_snapshot & MASK_MULTI)
				{
					fprintf(outfile, "multi-Gadget2");
					if(sim.out_snapshot - MASK_MULTI > MASK_PCLS)
					{
						fprintf(outfile, ", ");
					}
				}
				else
				{
					fprintf(outfile, "Gadget2");
					if(sim.out_snapshot > MASK_PCLS)
					{
						fprintf(outfile, ", ");
					}
				}
			}
			if(sim.out_snapshot & MASK_PCLS)
			{
				fprintf(outfile, "particles");
				if((sim.out_snapshot & MASK_MULTI == 0 && sim.out_snapshot > MASK_DELTA) || sim.out_snapshot - MASK_MULTI > MASK_DELTA)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if((sim.out_snapshot & MASK_MULTI == 0 && sim.out_snapshot > MASK_DBARE) || sim.out_snapshot - MASK_MULTI > MASK_DBARE)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
			}
			fprintf(outfile, "\n");
		}
		if(sim.out_snapshot & MASK_GADGET)
		{
			fprintf(outfile, "tracer factor       = %d", sim.tracer_factor[0]);
			for(i=1; i<=sim.baryon_flag + cosmo.num_ncdm; i++)
			{
				fprintf(outfile, ", %d", sim.tracer_factor[i]);
			}
			fprintf(outfile, "\n");
		}
		if(sim.downgrade_factor > 1)
		{
			fprintf(outfile, "downgrade factor    = %d", sim.downgrade_factor);
		}
		if(sim.num_pk > 0)
		{
			fprintf(outfile, "Pk redshifts        = ");
			for(i=0; i<sim.num_pk - 1; i++)
			{
				fprintf(outfile, "%lg, ", sim.z_pk[i]);
			}
			fprintf(outfile, "%lg\n", sim.z_pk[i]);
		}
		if(sim.out_pk)
		{
			fprintf(outfile, "Pk outputs          = ");
			if(sim.out_pk & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_pk > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_pk > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if(sim.out_pk > MASK_B)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_pk > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_T00)
			{
				fprintf(outfile, "T00");
				if(sim.out_pk > MASK_TIJ)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if(sim.out_pk > MASK_RBARE)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if(sim.out_pk > MASK_HIJ)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_pk > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_P)
			{
				fprintf(outfile, "p");
				if(sim.out_pk > MASK_XSPEC)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_XSPEC)
			{
				fprintf(outfile, "X-spectra");
				if(sim.out_pk > MASK_DELTA)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if(sim.out_pk > MASK_DBARE)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_pk & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
			}
			fprintf(outfile, "\n");
		}
		fprintf(outfile, "Pk bins             = %d\n", sim.numbins);

		if(sim.num_lightcone == 1)
		{
			fprintf(outfile, "lightcone vertex    = %lg, %lg, %lg\n", sim.lightcone[0].vertex[0], sim.lightcone[0].vertex[1], sim.lightcone[0].vertex[2]);
			fprintf(outfile, "lightcone outputs   = ");
			if(sim.out_lightcone[0] & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_lightcone[0] > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_lightcone[0] > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_lightcone[0] > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_lightcone[0] > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_GADGET)
			{
				fprintf(outfile, "Gadget2");
			}
			fprintf(outfile, "\n");
			if(sim.lightcone[0].opening > -1.)
			{
				fprintf(outfile, "lightcone opening half-angle = %lg\n", acos(sim.lightcone[0].opening) * 180. / M_PI);
			}
			fprintf(outfile, "lightcone distance  = %lg, %lg\n", sim.lightcone[0].distance[1], sim.lightcone[0].distance[0]);
			if(sim.lightcone[0].z != 0)
			{
				fprintf(outfile, "lightcone redshift  = %lg\n", sim.lightcone[0].z);
			}
			fprintf(outfile, "lightcone direction = %.15le, %.15le, %.15le\n", sim.lightcone[0].direction[0], sim.lightcone[0].direction[1], sim.lightcone[0].direction[2]);
			fprintf(outfile, "lightcone covering  = %lg\n", sim.covering[0]);
			if(sim.Nside[0][0] != sim.Nside[0][1])
			{
				fprintf(outfile, "lightcone Nside     = %d, %d\n", sim.Nside[0][0], sim.Nside[0][1]);
			}
			else
			{
				fprintf(outfile, "lightcone Nside     = %d\n", sim.Nside[0][0]);
			}
			fprintf(outfile, "lightcone pixel factor = %lg\n", sim.pixelfactor[0]);
			fprintf(outfile, "lightcone shell factor = %lg\n", sim.shellfactor[0]);
		}
		else if(sim.num_lightcone > 1)
		{
			for(i=0; i<sim.num_lightcone; i++)
			{
				fprintf(outfile, "lightcone %d vertex    = %lg, %lg, %lg\n", i, sim.lightcone[0].vertex[0], sim.lightcone[0].vertex[1], sim.lightcone[0].vertex[2]);
				fprintf(outfile, "lightcone %d outputs   = ", i);
				if(sim.out_lightcone[0] & MASK_PHI)
				{
					fprintf(outfile, "phi");
					if(sim.out_lightcone[0] > MASK_CHI)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_CHI)
				{
					fprintf(outfile, "chi");
					if(sim.out_lightcone[0] > MASK_POT)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_B)
				{
					fprintf(outfile, "B");
					if(sim.out_lightcone[0] > MASK_T00)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_HIJ)
				{
					fprintf(outfile, "hij");
					if(sim.out_lightcone[0] > MASK_P)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_GADGET)
				{
					fprintf(outfile, "Gadget2");
				}
				fprintf(outfile, "\n");
				if(sim.lightcone[0].opening > -1.)
				{
					fprintf(outfile, "lightcone %d opening half-angle = %lg\n", i, acos(sim.lightcone[0].opening) * 180. / M_PI);
				}
				fprintf(outfile, "lightcone %d distance  = %lg, %lg\n", i, sim.lightcone[0].distance[1], sim.lightcone[0].distance[0]);
				if(sim.lightcone[0].z != 0)
				{
					fprintf(outfile, "lightcone %d redshift  = %lg\n", i, sim.lightcone[0].z);
				}
				fprintf(outfile, "lightcone %d direction = %.15le, %.15le, %.15le\n", i, sim.lightcone[0].direction[0], sim.lightcone[0].direction[1], sim.lightcone[0].direction[2]);
				fprintf(outfile, "lightcone %d covering  = %lg\n", i, sim.covering[0]);
				if(sim.Nside[0][0] != sim.Nside[0][1])
				{
					fprintf(outfile, "lightcone %d Nside     = %d, %d\n", i, sim.Nside[0][0], sim.Nside[0][1]);
				}
				else
				{
					fprintf(outfile, "lightcone %d Nside     = %d\n", i, sim.Nside[0][0]);
				}
				fprintf(outfile, "lightcone %d pixel factor = %lg\n", i, sim.pixelfactor[0]);
				fprintf(outfile, "lightcone %d shell factor = %lg\n", i, sim.shellfactor[0]);
			}
		}

		fprintf(outfile, "\n\n# hibernations\n\n");
		if(sim.num_restart > 0)
		{
			fprintf(outfile, "hibernation redshifts       = ");
			for(i=0; i<sim.num_restart - 1; i++)
			{
				fprintf(outfile, "%lg, ", sim.z_restart[i]);
			}
			fprintf(outfile, "%lg\n", sim.z_restart[i]);
		}
		if(sim.wallclocklimit > 0.)
		{
			fprintf(outfile, "hibernation wallclock limit = %lg\n", sim.wallclocklimit);
		}
		if(sim.restart_path[0] != '\0')
		{
			fprintf(outfile, "hibernation path            = %s\n", sim.restart_path);
		}
		fprintf(outfile, "hibernation file base       = %s\n", sim.basename_restart);

		fclose(outfile);
	}
}

/////////////////////////////////////
// Restart settings for f(R) gravity
/////////////////////////////////////
/////////////////////////////////////
void writeRestartSettings_fR(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
	const double a,
	const double Hubble,
	const double Rbar,
	const double dot_Rbar,
	const double tau,
	const double dtau,
	const double dtau_old,
	const int cycle,
	const int restartcount = -1
)
{
	char buffer[2*PARAM_MAX_LENGTH+24];
	char buffer_bin[2*PARAM_MAX_LENGTH+24];
	FILE * outfile;
	ofstream outfile_bin;

	int i;

	if(!parallel.isRoot()) return;

	if(restartcount >= 0)
	{
		sprintf(buffer, "%s%s%03d.ini", sim.restart_path, sim.basename_restart, restartcount);
		sprintf(buffer_bin, "%s%s%03d.ini.bin", sim.restart_path, sim.basename_restart, restartcount);
	}
	else
	{
		sprintf(buffer, "%s%s.ini", sim.restart_path, sim.basename_restart);
		sprintf(buffer_bin, "%s%s.ini.bin", sim.restart_path, sim.basename_restart);
	}

	outfile_bin.open(buffer_bin, ios::out | ios::trunc | ios::binary);

	if(outfile_bin.is_open())
	{
		outfile_bin.write((char*) & a, sizeof(double));
		outfile_bin.write((char*) & tau, sizeof(double));
		outfile_bin.write((char*) & dtau, sizeof(double));
		outfile_bin.write((char*) & dtau_old, sizeof(double));
		outfile_bin.write((char*) & Hubble, sizeof(double));
		outfile_bin.write((char*) & Rbar, sizeof(double));
		outfile_bin.write((char*) & dot_Rbar, sizeof(double));
		outfile_bin.write((char*) & sim, sizeof(metadata));
		outfile_bin.write((char*) & cosmo, sizeof(cosmology));
		outfile_bin.close();
	}
	else
	{
		cout << "error opening binary restart file" << endl;
	}

	outfile = fopen(buffer, "w");

	if(outfile == NULL)
	{
		cout << " error opening file for restart settings!" << endl;
	}
	else
	{
		fprintf(outfile, "# automatically generated settings for restart after hibernation ");
		if(restartcount < 0)
		{
			fprintf(outfile, "due to wallclock limit ");
		}
		else
		{
			fprintf(outfile, "requested ");
		}

		fprintf(outfile, "at redshift z=%f\n\n", (1./a)-1.);
		fprintf(outfile, "fRevolution version = %g\n", FREVOLUTION_VERSION); //
		fprintf(outfile, "\n#====== IC generation ======#\n");
		fprintf(outfile, "IC generator        = restart\n");
		if(restartcount >= 0)
		{
			sprintf(buffer, "%03d", restartcount);
		}
		else
		{
			buffer[0] = '\0';
		}

		if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
		{
			fprintf(outfile, "particle file       = %s%s%s_cdm.h5", sim.restart_path, sim.basename_restart, buffer);
			if(sim.baryon_flag)
			{
				fprintf(outfile, ", %s%s%s_b.h5", sim.restart_path, sim.basename_restart, buffer);
			}
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				fprintf(outfile, ", %s%s%s_ncdm%d.h5", sim.restart_path, sim.basename_restart, buffer, i);
			}
		}
		else
		{
			fprintf(outfile, "particle file       = %s%s%s_cdm", sim.restart_path, sim.basename_restart, buffer);
			if(sim.baryon_flag)
			{
				fprintf(outfile, ", %s%s%s_b", sim.restart_path, sim.basename_restart, buffer);
			}
			for(i=0; i<cosmo.num_ncdm; i++)
			{
				fprintf(outfile, ", %s%s%s_ncdm%d", sim.restart_path, sim.basename_restart, buffer, i);
			}
		}

		fprintf(outfile, "\n");
		fprintf(outfile, "metric file         = %s%s%s_phi.h5", sim.restart_path, sim.basename_restart, buffer);
		fprintf(outfile, ", %s%s%s_chi.h5", sim.restart_path, sim.basename_restart, buffer);

		if(sim.vector_flag == VECTOR_PARABOLIC)
		{
			fprintf(outfile, ", %s%s%s_B.h5", sim.restart_path, sim.basename_restart, buffer);
		}
		else
		{
#ifdef CHECK_B
			fprintf(outfile, ", %s%s%s_B_check.h5", sim.restart_path, sim.basename_restart, buffer);
#else
			fprintf(outfile, ", %s%s%s_B.h5", sim.restart_path, sim.basename_restart, buffer);
#endif
		}

		fprintf(outfile, ", %s%s%s_xi.h5", sim.restart_path, sim.basename_restart, buffer);
		fprintf(outfile, ", %s%s%s_xi_old.h5", sim.restart_path, sim.basename_restart, buffer);
		fprintf(outfile, ", %s%s%s_deltaR.h5", sim.restart_path, sim.basename_restart, buffer);
		fprintf(outfile, "\n");
		// TODO: check what is necessary
		fprintf(outfile, "seed                = %d\n", ic.seed);
		fprintf(outfile, "restart redshift    = %f\n", (1./a) - 1.);
		fprintf(outfile, "cycle               = %d\n", cycle);
		fprintf(outfile, "tau                 = %e\n", tau);
		fprintf(outfile, "dtau                = %e\n", dtau);
		fprintf(outfile, "dtau_old            = %e\n", dtau_old);
		fprintf(outfile, "a                 	= %e\n", a);
		fprintf(outfile, "Hubble            	= %e\n", Hubble);
		fprintf(outfile, "Rbar                = %e\n", Rbar);
		fprintf(outfile, "dot_Rbar            = %e\n", dot_Rbar);
		fprintf(outfile, "\n");

		if(ic.flags & ICFLAG_KSPHERE)
		{
			fprintf(outfile, "k-domain            = sphere\n");
		}
		else
		{
			fprintf(outfile, "k-domain            = cube\n");
		}

		fprintf(outfile, "\n#====== primordial power spectrum ======#\n");
		fprintf(outfile, "k_pivot = %lg\n", ic.k_pivot);
		fprintf(outfile, "A_s     = %lg\n", ic.A_s);
		fprintf(outfile, "n_s     = %lg\n", ic.n_s);

		fprintf(outfile, "\n#====== cosmological parameters ======#\n");
		fprintf(outfile, "h         = %lg\n", cosmo.h);
		fprintf(outfile, "Omega_cdm = %.15le\n", cosmo.Omega_cdm);
		fprintf(outfile, "Omega_b   = %.15le\n", cosmo.Omega_b);
		fprintf(outfile, "Omega_g   = %.15le\n", cosmo.Omega_g);
		fprintf(outfile, "Omega_ur  = %.15le\n", cosmo.Omega_ur);
		fprintf(outfile, "N_ncdm    = %d\n", cosmo.num_ncdm);
		if(cosmo.num_ncdm > 0)
		{
			fprintf(outfile, "m_cdm     = ");
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				fprintf(outfile, "%9lf, ", cosmo.m_ncdm[i]);
			}
			fprintf(outfile, "%9lf\n", cosmo.m_ncdm[i]);
			fprintf(outfile, "T_cdm     = ");
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				fprintf(outfile, "%9lf, ", cosmo.T_ncdm[i]);
			}
			fprintf(outfile, "%9lf\n", cosmo.T_ncdm[i]);
			fprintf(outfile, "deg_cdm   = ");
			for(i=0; i<cosmo.num_ncdm - 1; i++)
			{
				fprintf(outfile, "%lf, ", cosmo.deg_ncdm[i]);
			}
			fprintf(outfile, "%lg\n", cosmo.deg_ncdm[i]);
		}

		fprintf(outfile, "\n#====== simulation settings ======#\n");
		if(sim.baryon_flag > 0)
		{
			fprintf(outfile, "baryon treatment    = sample\n");
		}

		if(sim.radiation_flag > 0)
		{
			fprintf(outfile, "radiation treatment = CLASS\n");
			fprintf(outfile, "switch delta_rad    = %lf\n", sim.z_switch_deltarad);
			if(cosmo.num_ncdm > 0)
			{
				fprintf(outfile, "switch delta_ncdm   = ");
				for(i=0; i<cosmo.num_ncdm-1; i++)
				{
					fprintf(outfile, "%lf, ", sim.z_switch_deltancdm[i]);
				}
				fprintf(outfile, "%lf\n", sim.z_switch_deltancdm[i]);
			}
			fprintf(outfile, "switch linear chi   = %lf\n", sim.z_switch_linearchi);
		}

		if(sim.vector_flag == VECTOR_ELLIPTIC)
		{
			fprintf(outfile, "vector method       = elliptic\n");
		}
		else
		{
			fprintf(outfile, "vector method       = parabolic\n");
		}

		fprintf(outfile, "initial redshift    = %lg\n", sim.z_in);
		fprintf(outfile, "boxsize             = %lg\n", sim.boxsize);
		fprintf(outfile, "Ngrid               = %d\n", sim.numpts);
		fprintf(outfile, "tiling factor       = %d\n", ic.numtile[0]);
		fprintf(outfile, "Courant factor      = %lg\n", sim.Cf);
		fprintf(outfile, "time step limit     = %lg\n", sim.steplimit);
		if(cosmo.num_ncdm > 0)
		{
			fprintf(outfile, "move limit          = %lg\n", sim.movelimit);
		}
		fprintf(outfile, "check fields        = %d\n", sim.check_fields);
		fprintf(outfile, "CYCLE_INFO_INTERVAL = %d\n", sim.CYCLE_INFO_INTERVAL);
		fprintf(outfile, "	gravity theory      = fr\n");

		fprintf(outfile, "\n#====== f(R) settings ======#\n");
		if(sim.fR_model == FR_MODEL_RN || sim.fR_model == FR_MODEL_R2)
		{
			// TODO: Check that comments are well written -- need something for FR_MODEL_DELTA too?
			if(sim.fR_model == FR_MODEL_R2)
			{
				fprintf(outfile, "# WARNING: f(R) parameters for R + R^2 model hav already been rescaled. See fR_tools.hpp for additional information.\n");
			}
			fprintf(outfile, "f(R) model                  = RN\n");
			fprintf(outfile, "f(R) parameters             = %.15le, %.15le\n", sim.fR_params[0], sim.fR_params[1]);
		}
		else if(sim.fR_model == FR_MODEL_DELTA)
		{
			fprintf(outfile, "f(R) model                  = DE\n");
			fprintf(outfile, "f(R) parameters             = %.15le, %.15le\n", sim.fR_params[0], sim.fR_params[1]);
		}
		else if(sim.fR_model == FR_MODEL_HU_SAWICKI)
		{
			fprintf(outfile, "# WARNING: f(R) parameters for Hu-Sawicki model have already been rescaled. See fR_tools.hpp for additional information.\n");
			fprintf(outfile, "f(R) model                  = HS\n");
			fprintf(outfile, "f(R) parameters             = %.15le, %.15le, %.15le, %.15le\n",sim.fR_params[0], sim.fR_params[1], sim.fR_params[2], sim.fR_params[3]);
		}
		else
		{
			COUT << " error f(R) model not recognized!" << endl;
		}

		if(sim.relativistic_flag)
		{
			fprintf(outfile, "Newtonian f(R)              = 0\n");
		}
		else
		{
			fprintf(outfile, "Newtonian f(R)              = 1\n");
		}
		fprintf(outfile, "f(R) epsilon background     = %e\n", sim.fR_epsilon_bg);
		fprintf(outfile, "f(R) count max              = %e\n", sim.fR_count_max);
		fprintf(outfile, "f(R) epsilon fields         = %e\n", sim.fR_epsilon_fields);
		fprintf(outfile, "f(R) target precision       = %e\n", sim.fR_target_precision);
		fprintf(outfile, "lcdm background             = %d\n", sim.lcdm_background);

		fprintf(outfile, "\n#====== Multigrid and relaxation ======#\n");
		fprintf(outfile, "relaxation method          = %d\n", sim.relaxation_method);
		fprintf(outfile, "relaxation error           = %e\n", sim.relaxation_error);
		fprintf(outfile, "red black                  = %d\n", sim.red_black);
		fprintf(outfile, "overrelaxation factor      = %e\n", sim.overrelaxation_factor);
		fprintf(outfile, "pre-smoothing              = %d\n", sim.pre_smoothing);
		fprintf(outfile, "post-smoothing             = %d\n", sim.post_smoothing);
		fprintf(outfile, "multigrid n-grids          = %d\n", sim.multigrid_n_grids);
		fprintf(outfile, "multigrid n-cycles         = %d\n", sim.multigrid_n_cycles);
		fprintf(outfile, "multigrid damping          = %d\n", sim.multigrid_damping);
		fprintf(outfile, "check shape                = %d\n", sim.multigrid_check_shape);
		if(sim.multigrid_shape == MULTIGRID_SHAPE_V)
		{
			fprintf(outfile, "multigrid shape            = V\n");
		}
		else if(sim.multigrid_shape == MULTIGRID_SHAPE_W)
		{
			fprintf(outfile, "multigrid shape            = W\n");
		}
		if(sim.multigrid_restrict_mode == RESTRICT_XI)
		{
			fprintf(outfile, "restrict mode              = xi");
		}
		else
		{
			fprintf(outfile, "restrict mode              = deltaR");
		}

		fprintf(outfile, "\n#====== output ======#\n");
		fprintf(outfile, "output path         = %s\n", sim.output_path);
		fprintf(outfile, "generic file base   = %s\n", sim.basename_generic);
		fprintf(outfile, "snapshot file base  = %s\n", sim.basename_snapshot);
		fprintf(outfile, "Pk file base        = %s\n", sim.basename_pk);
		fprintf(outfile, "lightcone file base = %s\n", sim.basename_lightcone);
		if(sim.num_snapshot > 0)
		{
			fprintf(outfile, "snapshot redshifts = ");
			for(i=0; i<sim.num_snapshot - 1; i++)
			{
				fprintf(outfile, "%lg, ", sim.z_snapshot[i]);
			}
			fprintf(outfile, "%lg\n", sim.z_snapshot[i]);
		}

		if(sim.out_snapshot)
		{
			fprintf(outfile, "snapshot outputs   = ");
			if(sim.out_snapshot & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_snapshot > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_snapshot > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if(sim.out_snapshot > MASK_B)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_snapshot > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_T00)
			{
				fprintf(outfile, "T00");
				if(sim.out_snapshot > MASK_TIJ)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if(sim.out_snapshot > MASK_RBARE)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if(sim.out_snapshot > MASK_HIJ)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_snapshot > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_P)
			{
				fprintf(outfile, "p");
				if(sim.out_snapshot > MASK_GADGET)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_GADGET)
			{
				fprintf(outfile, "Gadget2");
				if(sim.out_snapshot > MASK_PCLS)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_PCLS)
			{
				fprintf(outfile, "particles");
				if(sim.out_snapshot > MASK_DELTA)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if(sim.out_snapshot > MASK_DBARE)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
				if(sim.out_snapshot > MASK_XI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_XI)
			{
				fprintf(outfile, "xi");
				if(sim.out_snapshot > MASK_ZETA)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_ZETA)
			{
				fprintf(outfile, "zeta");
				if(sim.out_snapshot > MASK_DELTAT)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_DELTAT)
			{
				fprintf(outfile, "eightpiG_deltaT");
				if(sim.out_snapshot > MASK_DELTAR)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_DELTAR)
			{
				fprintf(outfile, "deltaR");
				if(sim.out_snapshot > MASK_LAPLACE_XI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_snapshot & MASK_LAPLACE_XI)
			{
				fprintf(outfile, "laplace_xi");
			}
			fprintf(outfile, "\n");
		}

		if(sim.out_snapshot & MASK_GADGET)
		{
			fprintf(outfile, "tracer factor      = %d", sim.tracer_factor[0]);
			for(i=1; i <= sim.baryon_flag + cosmo.num_ncdm; i++)
			{
				fprintf(outfile, ", %d", sim.tracer_factor[i]);
			}
			fprintf(outfile, "\n");
		}

		// OUTPUTS FOR check_fields()
		if(sim.out_check)
		{
			fprintf(outfile, "check outputs       = ");
			if(sim.out_check & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_check > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_check > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if(sim.out_check > MASK_B)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_check > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_T00)
			{
				fprintf(outfile, "T00");
				if(sim.out_check > MASK_TIJ)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if(sim.out_check > MASK_RBARE)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if(sim.out_check > MASK_HIJ)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_check > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_P)
			{
				fprintf(outfile, "p");
				if(sim.out_check > MASK_GADGET)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_GADGET)
			{
				fprintf(outfile, "Gadget2");
				if(sim.out_check > MASK_PCLS)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_PCLS)
			{
				fprintf(outfile, "particles");
				if(sim.out_check > MASK_DELTA)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if(sim.out_check > MASK_DBARE)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
				if(sim.out_check > MASK_XI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_XI)
			{
				fprintf(outfile, "xi");
				if(sim.out_check > MASK_ZETA)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_ZETA)
			{
				fprintf(outfile, "zeta");
				if(sim.out_check > MASK_DELTAT)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_DELTAT)
			{
				fprintf(outfile, "eightpiG_deltaT");
				if(sim.out_check > MASK_DELTAR)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_DELTAR)
			{
				fprintf(outfile, "deltaR");
				if(sim.out_check > MASK_LAPLACE_XI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_check & MASK_LAPLACE_XI)
			{
				fprintf(outfile, "laplace_xi");
			}
			fprintf(outfile, "\n");
		}

		if(sim.downgrade_factor > 1)
		{
			fprintf(outfile, "downgrade factor   = %d", sim.downgrade_factor);
		}

		if(sim.num_pk > 0)
		{
			fprintf(outfile, "Pk redshifts        = ");
			for(i=0; i<sim.num_pk-1; i++)
			{
				fprintf(outfile, "%lg, ", sim.z_pk[i]);
			}
			fprintf(outfile, "%lg\n", sim.z_pk[i]);
		}

		if(sim.out_pk)
		{
			fprintf(outfile, "Pk outputs          = ");
			if(sim.out_pk & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_pk > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_pk > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if(sim.out_pk > MASK_B)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_pk > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_T00)
			{
				fprintf(outfile, "T00");
				if(sim.out_pk > MASK_TIJ)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if(sim.out_pk > MASK_RBARE)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if(sim.out_pk > MASK_HIJ)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_pk > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_P)
			{
				fprintf(outfile, "p");
				if(sim.out_pk > MASK_XSPEC)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_XSPEC)
			{
				fprintf(outfile, "X-spectra");
				if(sim.out_pk > MASK_DELTA)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if(sim.out_pk > MASK_DBARE)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
				if(sim.out_pk > MASK_XI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_XI)
			{
				fprintf(outfile, "xi");
				if(sim.out_pk > MASK_ZETA)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_ZETA)
			{
				fprintf(outfile, "zeta");
				if(sim.out_pk > MASK_DELTAT)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_DELTAT)
			{
				fprintf(outfile, "eightpiG_deltaT");
				if(sim.out_pk > MASK_DELTAR)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_DELTAR)
			{
				fprintf(outfile, "deltaR");
				if(sim.out_pk > MASK_LAPLACE_XI)
				{
					fprintf(outfile, ", ");
				}
			}

			if(sim.out_pk & MASK_LAPLACE_XI)
			{
				fprintf(outfile, "laplace_xi");
			}
			fprintf(outfile, "\n");
		}

		fprintf(outfile, "Pk bins            = %d\n", sim.numbins);

		if(sim.num_lightcone == 1)
		{
			fprintf(outfile, "lightcone vertex    = %lg, %lg, %lg\n", sim.lightcone[0].vertex[0], sim.lightcone[0].vertex[1], sim.lightcone[0].vertex[2]);
			fprintf(outfile, "lightcone outputs   = ");
			if(sim.out_lightcone[0] & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if(sim.out_lightcone[0] > MASK_CHI)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if(sim.out_lightcone[0] > MASK_POT)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_B)
			{
				fprintf(outfile, "B");
				if(sim.out_lightcone[0] > MASK_T00)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if(sim.out_lightcone[0] > MASK_P)
				{
					fprintf(outfile, ", ");
				}
			}
			if(sim.out_lightcone[0] & MASK_GADGET)
			{
				fprintf(outfile, "Gadget2");
			}
			fprintf(outfile, "\n");
			if(sim.lightcone[0].opening > -1.)
			{
				fprintf(outfile, "lightcone opening half-angle = %lg\n", acos(sim.lightcone[0].opening) * 180. / M_PI);
			}
			fprintf(outfile, "lightcone distance  = %lg, %lg\n", sim.lightcone[0].distance[1], sim.lightcone[0].distance[0]);
			if(sim.lightcone[0].z != 0)
			{
				fprintf(outfile, "lightcone redshift  = %lg\n", sim.lightcone[0].z);
			}
			fprintf(outfile, "lightcone direction = %.15le, %.15le, %.15le\n", sim.lightcone[0].direction[0], sim.lightcone[0].direction[1], sim.lightcone[0].direction[2]);
			fprintf(outfile, "lightcone covering  = %lg\n", sim.covering[0]);
			if(sim.Nside[0][0] != sim.Nside[0][1])
			{
				fprintf(outfile, "lightcone Nside     = %d, %d\n", sim.Nside[0][0], sim.Nside[0][1]);
			}
			else
			{
				fprintf(outfile, "lightcone Nside     = %d\n", sim.Nside[0][0]);
			}
			fprintf(outfile, "lightcone pixel factor = %lg\n", sim.pixelfactor[0]);
			fprintf(outfile, "lightcone shell factor = %lg\n", sim.shellfactor[0]);
		}
		else if(sim.num_lightcone > 1)
		{
			for(i=0; i<sim.num_lightcone; i++)
			{
				fprintf(outfile, "lightcone %d vertex    = %lg, %lg, %lg\n", i, sim.lightcone[0].vertex[0], sim.lightcone[0].vertex[1], sim.lightcone[0].vertex[2]);
				fprintf(outfile, "lightcone %d outputs   = ", i);
				if(sim.out_lightcone[0] & MASK_PHI)
				{
					fprintf(outfile, "phi");
					if(sim.out_lightcone[0] > MASK_CHI)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_CHI)
				{
					fprintf(outfile, "chi");
					if(sim.out_lightcone[0] > MASK_POT)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_B)
				{
					fprintf(outfile, "B");
					if(sim.out_lightcone[0] > MASK_T00)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_HIJ)
				{
					fprintf(outfile, "hij");
					if(sim.out_lightcone[0] > MASK_P)
					{
						fprintf(outfile, ", ");
					}
				}
				if(sim.out_lightcone[0] & MASK_GADGET)
				{
					fprintf(outfile, "Gadget2");
				}
				fprintf(outfile, "\n");
				if(sim.lightcone[0].opening > -1.)
				{
					fprintf(outfile, "lightcone %d opening half-angle = %lg\n", i, acos(sim.lightcone[0].opening) * 180. / M_PI);
				}
				fprintf(outfile, "lightcone %d distance  = %lg, %lg\n", i, sim.lightcone[0].distance[1], sim.lightcone[0].distance[0]);
				if(sim.lightcone[0].z != 0)
				{
					fprintf(outfile, "lightcone %d redshift  = %lg\n", i, sim.lightcone[0].z);
				}
				fprintf(outfile, "lightcone %d direction = %.15le, %.15le, %.15le\n", i, sim.lightcone[0].direction[0], sim.lightcone[0].direction[1], sim.lightcone[0].direction[2]);
				fprintf(outfile, "lightcone %d covering  = %lg\n", i, sim.covering[0]);
				if(sim.Nside[0][0] != sim.Nside[0][1])
				{
					fprintf(outfile, "lightcone %d Nside     = %d, %d\n", i, sim.Nside[0][0], sim.Nside[0][1]);
				}
				else
				{
					fprintf(outfile, "lightcone %d Nside     = %d\n", i, sim.Nside[0][0]);
				}
				fprintf(outfile, "lightcone %d pixel factor = %lg\n", i, sim.pixelfactor[0]);
				fprintf(outfile, "lightcone %d shell factor = %lg\n", i, sim.shellfactor[0]);
			}
		}

		fprintf(outfile, "\n#====== hibernation ======#\n");
		if(sim.num_restart > 0)
		{
			fprintf(outfile, "hibernation redshifts       = ");
			for(i=0; i<sim.num_restart - 1; i++)
			{
				fprintf(outfile, "%lg, ", sim.z_restart[i]);
			}
			fprintf(outfile, "%lg\n", sim.z_restart[i]);
		}

		if(sim.wallclocklimit > 0.)
		{
			fprintf(outfile, "hibernation wallclock limit = %lg\n", sim.wallclocklimit);
		}

		if(sim.restart_path[0] != '\0')
		{
			fprintf(outfile, "hibernation path            = %s\n", sim.restart_path);
		}

		fprintf(outfile, "hibernation file base       = %s\n", sim.basename_restart);
		if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
		{
			fprintf(outfile, "particle save mode          = hdf5\n");
		}
		else if(sim.hibernation_save_mode == HIB_SAVE_GADGET2)
		{
			fprintf(outfile, "particle save mode          = gadget2\n");
		}

		fclose(outfile);
	}
}


//////////////////////////
// hibernate
//////////////////////////
// Description:
//   creates a hibernation point by writing snapshots of the simulation data and metadata
//
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for non-cold DM
//   phi            reference to field containing first Bardeen potential
//   chi            reference to field containing difference of Bardeen potentials
//   Bi             reference to vector field containing frame-dragging potential
//   a              scale factor
//   tau            conformal coordinate time
//   dtau           time step
//   cycle          current main control loop cycle count
//   restartcount   restart counter aka number of hibernation point (default -1)
//                  if < 0 no number is associated to the hibernation point
//
// Returns:
//
//////////////////////////

void hibernate(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
	Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm,
	Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_b,
	Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm,
	Field<Real> & phi,
	Field<Real> & chi,
	Field<Real> & Bi,
	const double a,
	const double tau,
	const double dtau,
	const int cycle,
	const int restartcount = -1
)
{
	string h5filename;
	char buffer[5];
	int i;
	Site x(Bi.lattice());

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.restart_path);
	h5filename += sim.basename_restart;
	if(restartcount >= 0)
	{
		sprintf(buffer, "%03d", restartcount);
		h5filename += buffer;
	}

	writeRestartSettings(sim, ic, cosmo, a, tau, dtau, cycle, restartcount);

#ifndef CHECK_B
	if(sim.vector_flag == VECTOR_PARABOLIC)
#endif
	for(x.first(); x.test(); x.next())
	{
		Bi(x,0) /= a * a * sim.numpts;
		Bi(x,1) /= a * a * sim.numpts;
		Bi(x,2) /= a * a * sim.numpts;
	}

#ifdef EXTERNAL_IO
	while(ioserver.openOstream()== OSTREAM_FAIL);

	pcls_cdm->saveHDF5_server_open(h5filename + "_cdm");
	if(sim.baryon_flag)
	{
		pcls_b->saveHDF5_server_open(h5filename + "_b");
	}
	for(i=0; i<cosmo.num_ncdm; i++)
	{
		if(sim.numpcl[1+sim.baryon_flag+i] < 1) continue;
		sprintf(buffer, "%d", i);
		pcls_ncdm[i].saveHDF5_server_open(h5filename + "_ncdm" + buffer);
	}

	if(sim.relativistic_flag > 0)
	{
		phi.saveHDF5_server_open(h5filename + "_phi");
		chi.saveHDF5_server_open(h5filename + "_chi");
	}

	if(sim.vector_flag == VECTOR_PARABOLIC)
	{
		Bi.saveHDF5_server_open(h5filename + "_B");
	}
#ifdef CHECK_B
	else
	{
		Bi.saveHDF5_server_open(h5filename + "_B_check");
	}
#endif

	pcls_cdm->saveHDF5_server_write();
	if(sim.baryon_flag)
	{
		pcls_b->saveHDF5_server_write();
	}
	for(i=0; i<cosmo.num_ncdm; i++)
	{
		if(sim.numpcl[1+sim.baryon_flag+i] < 1) continue;
		pcls_ncdm[i].saveHDF5_server_write();
	}

	if(sim.relativistic_flag > 0)
	{
		phi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
		chi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	}

#ifndef CHECK_B
	if(sim.vector_flag == VECTOR_PARABOLIC)
#endif
		Bi.saveHDF5_server_write(NUMBER_OF_IO_FILES);

	ioserver.closeOstream();
#else
	pcls_cdm->saveHDF5(h5filename + "_cdm", 1);
	if(sim.baryon_flag)
	{
		pcls_b->saveHDF5(h5filename + "_b", 1);
	}
	for(i=0; i<cosmo.num_ncdm; i++)
	{
		if(sim.numpcl[1+sim.baryon_flag+i] < 1) continue;
		sprintf(buffer, "%d", i);
		pcls_ncdm[i].saveHDF5(h5filename + "_ncdm" + buffer, 1);
	}

	if(sim.relativistic_flag > 0)
	{
		phi.saveHDF5(h5filename + "_phi.h5");
		chi.saveHDF5(h5filename + "_chi.h5");
	}

	if(sim.vector_flag == VECTOR_PARABOLIC)
	{
		Bi.saveHDF5(h5filename + "_B.h5");
	}
#ifdef CHECK_B
	else
	{
		Bi.saveHDF5(h5filename + "_B_check.h5");
	}
#endif
#endif
}

/////////////////////////////
// Hibernate for f(R) gravity
/////////////////////////////
void hibernate_fR(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
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
	const double tau,
	const double dtau,
	const double dtau_old,
	const int cycle,
	const int restartcount = -1
)
{
	char buffer[12];
	int i;
	Site x(Bi.lattice());
	string h5filename;
	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.restart_path);
	h5filename += sim.basename_restart;
	if(restartcount >= 0)
	{
		sprintf(buffer, "%03d", restartcount);
		h5filename += buffer;
	}

	writeRestartSettings_fR(sim, ic, cosmo, a, Hubble, Rbar, dot_Rbar, tau, dtau, dtau_old, cycle, restartcount);

#ifndef CHECK_B
	if(sim.vector_flag == VECTOR_PARABOLIC)
#endif
	for(x.first(); x.test(); x.next())
	{
		Bi(x,0) /= a * a * sim.numpts;
		Bi(x,1) /= a * a * sim.numpts;
		Bi(x,2) /= a * a * sim.numpts;
	}

#ifdef EXTERNAL_IO
	while (ioserver.openOstream() == OSTREAM_FAIL);
	if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
	{
		pcls_cdm->saveHDF5_server_open(h5filename + "_cdm");
		if(sim.baryon_flag)
		{
			pcls_b->saveHDF5_server_open(h5filename + "_b");
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "%d", i);
			pcls_ncdm[i].saveHDF5_server_open(h5filename + "_ncdm" + buffer);
		}
	}

	phi.saveHDF5_server_open(h5filename + "_phi");
	chi.saveHDF5_server_open(h5filename + "_chi");
	xi.saveHDF5_server_open(h5filename + "_xi");
	xi_old.saveHDF5_server_open(h5filename + "_xi_old");
	deltaR.saveHDF5_server_open(h5filename + "_deltaR");
	Bi.saveHDF5_server_open(h5filename + "_B");

	if(sim.hibernation_save_mode == HIB_SAVE_HDF5)
	{
		pcls_cdm->saveHDF5_server_write();
		if(sim.baryon_flag)
		{
			pcls_b->saveHDF5_server_write();
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
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
		pcls_cdm->saveGadget2(h5filename + "_cdm", hdr, sim.tracer_factor[0]);

		if(sim.baryon_flag)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[1] / sim.tracer_factor[1]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			pcls_b->saveGadget2(h5filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for (i=0; i<cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (unsigned int) (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			pcls_ncdm[i].saveGadget2(h5filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
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
		pcls_cdm->saveHDF5(h5filename + "_cdm", 1);
		if(sim.baryon_flag)
		{
			pcls_b->saveHDF5(h5filename + "_b", 1);
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "%d", i);
			pcls_ncdm[i].saveHDF5(h5filename + "_ncdm" + buffer, 1);
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
		pcls_cdm->saveGadget2(h5filename + "_cdm", hdr, sim.tracer_factor[0]);

		if(sim.baryon_flag)
		{
			hdr.npart[1] = (unsigned int) (sim.numpcl[1] / sim.tracer_factor[1]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			pcls_b->saveGadget2(h5filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for(i=0; i<cosmo.num_ncdm; i++)
		{
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (unsigned int) (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]);
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			pcls_ncdm[i].saveGadget2(h5filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
		}
	}

	COUT << "Saving fields..." << endl;
	phi.saveHDF5(h5filename + "_phi.h5");
	chi.saveHDF5(h5filename + "_chi.h5");
	xi.saveHDF5(h5filename + "_xi.h5");
	xi_old.saveHDF5(h5filename + "_xi_old.h5");
	deltaR.saveHDF5(h5filename + "_deltaR.h5");
	Bi.saveHDF5(h5filename + "_B.h5");

#endif
}

#endif
