//////////////////////////
// ic_read.hpp
//////////////////////////
//
// read initial conditions from disk
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef IC_READ_FR_HEADER
#define IC_READ_FR_HEADER

/////////////////////////////////////////////////
/////////////////////////////////////////////////
void readIC_fR(
	metadata & sim,
	icsettings & ic,
	cosmology & cosmo,
	const double fourpiG,
	double & a,
	double & Hubble,
	double & Rbar,
	double & dot_Rbar,
	double & fbar,
	double & fRbar,
	double & fRRbar,
	double & tau,
	double & dtau,
	double & dtau_old,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b,
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm,
	double * maxvel,
	Field<Real> * phi,
	Field<Real> * chi,
	Field<Real> * Bi,
	Field<Real> * xi,
	Field<Real> * xi_old,
	Field<Real> * deltaR,
	Field<Real> * source,
	Field<Real> * Sij,
	Field<Cplx> * scalarFT,
	Field<Cplx> * BiFT,
	Field<Cplx> * SijFT,
	PlanFFT<Cplx> * plan_phi,
	PlanFFT<Cplx> * plan_chi,
	PlanFFT<Cplx> * plan_Bi,
	PlanFFT<Cplx> * plan_source,
	PlanFFT<Cplx> * plan_Sij,
	int & cycle,
	int & snapcount,
	int & pkcount,
	int & restartcount
)
{
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_b_info;
	part_simple_dataType pcls_b_dataType;
	part_simple_info pcls_ncdm_info[MAX_PCL_SPECIES];
	part_simple_dataType pcls_ncdm_dataType;
	Real boxSize[3] = {1.,1.,1.};
	string filename;
	string temp_filename;
	string buf;
	int i, p;
	char * ext;
	char line[PARAM_MAX_LINESIZE];
	FILE * bgfile;
	struct fileDsc fd;
	gadget2_header hdr;
	long * numpcl;
	Real * dummy1;
	Real * dummy2;
	Site x(Bi->lattice());
	rKSite kFT(scalarFT->lattice());

	filename.reserve(PARAM_MAX_LENGTH);
	temp_filename.reserve(PARAM_MAX_LENGTH);
	hdr.npart[1] = 0;

	filename.assign((string) sim.restart_path + sim.basename_restart);
	if(restartcount > 0)
	{
		filename += "_";
		for(i=(int)log10(restartcount); i<2; ++i)
		{
			filename += "0";
		}
		filename += to_string(restartcount);
	}
	else if(!restartcount)
	{
		filename += "_000";
	}

	temp_filename.assign(filename);

	projection_init(phi);
	metadata mdtemp;
	cosmology cosmotemp;
	double dtemp;
	ifstream file_bin;

	file_bin.open((string) filename + ".ini.bin", ios::in | ios::binary);

	if(file_bin.is_open())
	{
		file_bin.read((char*) & dtemp, sizeof(double));
		a = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		Hubble = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		Rbar = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		dot_Rbar = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		fbar = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		fRbar = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		fRRbar = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		tau = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		dtau = dtemp;
		file_bin.read((char*) & dtemp, sizeof(double));
		dtau_old = dtemp;
		file_bin.read((char*) & mdtemp, sizeof(metadata));
		sim = mdtemp;
		file_bin.read((char*) & cosmotemp, sizeof(cosmology));
		cosmo = cosmotemp;
		// TODO: Add icsettings too?
		file_bin.close();
	}
	else
	{
		cout << "readIC: cannot open binray settings file: " << endl;
	}

	strcpy(pcls_cdm_info.type_name, "part_simple");
	pcls_cdm_info.mass = 0.;
	pcls_cdm_info.relativistic = false;
	pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);

	if((ext = strstr(ic.pclfile[0], ".h5")) != NULL)
	{
		temp_filename.assign(ic.pclfile[0], ext-ic.pclfile[0]);
		get_fileDsc_global(temp_filename + ".h5", fd);
		numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
		dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
		dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
		get_fileDsc_local(temp_filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
		pcls_cdm->loadHDF5(temp_filename, 1);
		free(numpcl);
		free(dummy1);
		free(dummy2);
	}
	else
	{
		i=0;
		do
		{
			temp_filename.assign(ic.pclfile[0]);
			pcls_cdm->loadGadget2(temp_filename, hdr);
			if(hdr.npart[1] == 0) break;
			if(hdr.time / a > 1.001 || hdr.time / a < 0.999)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
			}
			// sim.numpcl[0] += hdr.npart[1]; // TODO: Removed because we already read sim from binary file. Check!
			i++;
			if(hdr.num_files > 1)
			{
				ext = ic.pclfile[0];
				while (strchr(ext, (int) '.') != NULL)
				{
					ext = strchr(ext, (int) '.');
				}
				sprintf(ext+1, "%d", i);
			}
		}
		while(i<hdr.num_files);

		if(sim.baryon_flag == 1)
		{
			pcls_cdm->parts_info()->mass = cosmo.Omega_cdm / (Real) sim.numpcl[0];
		}
		else
		{
			pcls_cdm->parts_info()->mass = (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) sim.numpcl[0];
		}
	}

	COUT << " " << sim.numpcl[0] << " cdm particles read successfully." << endl;

	maxvel[0] = pcls_cdm->updateVel(update_q, 0., &phi, 1, &a);

	if(sim.baryon_flag == 1)
	{
		strcpy(pcls_b_info.type_name, "part_simple");
		pcls_b_info.mass = 0.;
		pcls_b_info.relativistic = false;
		pcls_b->initialize(pcls_b_info, pcls_b_dataType, &(phi->lattice()), boxSize);

		if((ext = strstr(ic.pclfile[1], ".h5")) != NULL)
		{
			temp_filename.assign(ic.pclfile[1], ext-ic.pclfile[1]);
			get_fileDsc_global(temp_filename + ".h5", fd);
			numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
			dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			get_fileDsc_local(temp_filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
			// for(i=0; i<fd.numProcPerFile; i++)// TODO: Removed because we already read sim from binary file. Check!
			// {
				// 	sim.numpcl[1] += numpcl[i];
				// }
				pcls_b->loadHDF5(temp_filename, 1);
				free(numpcl);
				free(dummy1);
				free(dummy2);
			}
			else
			{
				i=0;
				do
				{
					temp_filename.assign(ic.pclfile[1]);
					pcls_b->loadGadget2(temp_filename, hdr);
					if(hdr.npart[1] == 0) break;
					if(hdr.time / a > 1.001 || hdr.time / a < 0.999)
					{
						COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
					}
					// sim.numpcl[1] += hdr.npart[1];// TODO: Removed because we already read sim from binary file. Check!
					i++;
					if(hdr.num_files > 1)
					{
						ext = ic.pclfile[1];
						while (strchr(ext, (int) '.') != NULL)
						{
							ext = strchr(ext, (int) '.');
						}
						sprintf(ext+1, "%d", i);
					}
				}
				while(i<hdr.num_files);

				pcls_b->parts_info()->mass = cosmo.Omega_b / (Real) sim.numpcl[1];
			}

			COUT << " " << sim.numpcl[1] << " baryon particles read successfully." << endl;
			maxvel[1] = pcls_b->updateVel(update_q, 0., &phi, 1, &a);
		}
		else
		{
			sim.baryon_flag = 0;
		}

		for(p=0; p<cosmo.num_ncdm; p++)
		{
			strcpy(pcls_ncdm_info[p].type_name, "part_simple");
			pcls_ncdm_info[p].mass = 0.;
			pcls_ncdm_info[p].relativistic = true;
			pcls_ncdm[p].initialize(pcls_ncdm_info[p], pcls_ncdm_dataType, &(phi->lattice()), boxSize);
			if((ext = strstr(ic.pclfile[sim.baryon_flag+1+p], ".h5")) != NULL)
			{
				temp_filename.assign(ic.pclfile[sim.baryon_flag+1+p], ext-ic.pclfile[sim.baryon_flag+1+p]);
				get_fileDsc_global(temp_filename + ".h5", fd);
				numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
				dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
				dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
				get_fileDsc_local(temp_filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
				// for(i=0; i<fd.numProcPerFile; i++)// TODO: Removed because we already read sim from binary file. Check!
				// {
					// 	sim.numpcl[1] += numpcl[i];
					// }
					pcls_ncdm[p].loadHDF5(temp_filename, 1);
					free(numpcl);
					free(dummy1);
					free(dummy2);
				}
				else
				{
					i=0;
					do
					{
						temp_filename.assign(ic.pclfile[sim.baryon_flag+1+p]);
						pcls_ncdm[p].loadGadget2(temp_filename, hdr);
						if(hdr.npart[1] == 0) break;
						if(hdr.time / a > 1.001 || hdr.time / a < 0.999)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
						}
						// sim.numpcl[sim.baryon_flag+1+p] += hdr.npart[1];// TODO: Removed because we already read sim from binary file. Check!
						i++;
						if(hdr.num_files > 1)
						{
							ext = ic.pclfile[sim.baryon_flag+1+p];
							while(strchr(ext, (int) '.') != NULL)
							{
								ext = strchr(ext, (int) '.');
							}
							sprintf(ext+1, "%d", i);
						}
					}
					while (i<hdr.num_files);

					pcls_ncdm[p].parts_info()->mass = cosmo.Omega_ncdm[p] / (Real) sim.numpcl[sim.baryon_flag+1+p];
				}

				COUT << " " << sim.numpcl[sim.baryon_flag+1+p] << " ncdm particles read successfully." << endl;
				maxvel[sim.baryon_flag+1+p] = pcls_ncdm[p].updateVel(update_q, 0., &phi, 1, &a);
			}

			if(ic.restart_cycle >= 0)
			{
				Bi->loadHDF5(filename + "_B.h5");

				if(sim.vector_flag == VECTOR_PARABOLIC)
				{
					for(x.first(); x.test(); x.next())
					{
						(*Bi)(x,0) *= a * a / (sim.numpts * sim.numpts);
						(*Bi)(x,1) *= a * a / (sim.numpts * sim.numpts);
						(*Bi)(x,2) *= a * a / (sim.numpts * sim.numpts);
					}
					plan_Bi->execute(FFT_FORWARD);
				}

				phi->loadHDF5(filename + "_phi.h5");
				phi->updateHalo();

				chi->loadHDF5(filename + "_chi.h5");
				chi->updateHalo();

				xi->loadHDF5(filename + "_xi.h5");
				xi->updateHalo();

				xi_old->loadHDF5(filename + "_xi_old.h5");
				xi_old->updateHalo();

				deltaR->loadHDF5(filename + "_deltaR.h5");
				deltaR->updateHalo();

				if(parallel.isRoot())
				{
					sprintf(line, "%s%s_background.dat", sim.output_path, sim.basename_generic);
					bgfile = fopen(line, "r");
					if(bgfile == NULL)
					{
						COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to locate file for background output! A new file will be created" << endl;
						bgfile = fopen(line, "w");
						if(bgfile == NULL)
						{
							COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for background output!" << endl;
							parallel.abortForce();
						}
						else
						{
							fprintf(bgfile, "# background statistics\n# cycle   tau/boxsize    a              conformal H    R              phi(k=0)       T00(k=0)\n");
							fclose(bgfile);
						}
					}
					else
					{
						buf.reserve(PARAM_MAX_LINESIZE);
						buf.clear();

						if(fgets(line, PARAM_MAX_LINESIZE, bgfile) == 0)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for background output! A new file will be created" << endl;
						}
						else if(line[0] != '#')
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for background output has unexpected format! Contents will be overwritten!" << endl;
						}
						else
						{
							if(fgets(line, PARAM_MAX_LINESIZE, bgfile) == 0)
							{
								COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for background output! A new file will be created" << endl;
							}
							else if(line[0] != '#')
							{
								COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for background output has unexpected format! Contents will be overwritten!" << endl;
							}
						}

						while(fgets(line, PARAM_MAX_LINESIZE, bgfile) != 0)
						{
							if(sscanf(line, " %d", &i) != 1) break;

							if(i > ic.restart_cycle)
							{
								break;
							}
							else
							{
								buf += line;
							}
						}

						fclose(bgfile);

						sprintf(line, "%s%s_background.dat", sim.output_path, sim.basename_generic);
						bgfile = fopen(line, "w");

						if(bgfile == NULL)
						{
							COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for background output!" << endl;
							parallel.abortForce();
						}
						else
						{
							fprintf(bgfile, "# background statistics\n# cycle   tau/boxsize    a              conformal H    R              phi(k=0)       T00(k=0)\n");
							fwrite((const void *) buf.data(), sizeof(char), buf.length(), bgfile);
							fclose(bgfile);
							buf.clear();
						}
					}
				}
			}
			else
			{
				projection_init(Bi);
				projection_T0i_project(pcls_cdm, Bi, phi);
				if(sim.baryon_flag)
				{
					projection_T0i_project(pcls_b, Bi, phi);
				}
				projection_T0i_comm(Bi);
				plan_Bi->execute(FFT_FORWARD);
				projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);
				plan_Bi->execute(FFT_BACKWARD);
				Bi->updateHalo();

				projection_init(Sij);
				projection_Tij_project(pcls_cdm, Sij, a, phi);
				if(sim.baryon_flag)
				{
					projection_Tij_project(pcls_b, Sij, a, phi);
				}
				projection_Tij_comm(Sij);

				prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);
				plan_Sij->execute(FFT_FORWARD);
				projectFTscalar(*SijFT, *scalarFT);
				plan_chi->execute(FFT_BACKWARD);
				chi->updateHalo();
			}

			if(ic.restart_cycle >= 0)
			{
				cycle = ic.restart_cycle + 1;
			}

			while(snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
			{
				snapcount++;
			}

			while(pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
			{
				pkcount++;
			}

			while(restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
			{
				restartcount++;
			}
		}


		#endif
