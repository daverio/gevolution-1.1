//////////////////////////
// metadata.hpp
//////////////////////////
//
// Constants and metadata structures
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//         f(R):
//				 David Daverio (Cambridge University)
//				 Lorenzo Reverberi (University of Cape Town)
//
// Last modified: December 2016
//
//////////////////////////

#ifndef METADATA_HEADER
#define METADATA_HEADER

#define FREVOLUTION_VERSION 1.0

#ifndef MAX_OUTPUTS
#define MAX_OUTPUTS 32
#endif

#ifndef MAX_FR_PARAMS
#define MAX_FR_PARAMS 8
#endif

#ifndef PARAM_MAX_LENGTH
#define PARAM_MAX_LENGTH 128
#endif

#ifndef PARAM_MAX_LINESIZE
#define PARAM_MAX_LINESIZE 1024
#endif

#ifndef MAX_PCL_SPECIES
#define MAX_PCL_SPECIES 6
#endif

#ifndef METRICFILE_LENGTH
#define METRICFILE_LENGTH 20
#endif

#define MASK_PHI        		1
#define MASK_CHI        		2
#define MASK_POT        		4
#define MASK_B          		8
#define MASK_T00        		16
#define MASK_TIJ        		32
#define MASK_RBARE      		64
#define MASK_HIJ        		128
#define MASK_P          		256
#define MASK_GADGET     		512
#define MASK_PCLS       		1024
#define MASK_XSPEC      		2048
#define MASK_DELTA      		4096
#define MASK_DBARE      		8192
#define MASK_XI  		    		16384
#define MASK_ZETA       		32768
#define MASK_DELTAR	    		65536
#define MASK_DELTAT     		131072
#define MASK_LAPLACE_XI 		262144
#define MASK_PHI_EFFECTIVE  524288

#define ICFLAG_CORRECT_DISPLACEMENT 1
#define ICFLAG_KSPHERE              2

// Identifiers for IC generator modules
#define ICGEN_BASIC                 0
#define ICGEN_READ_FROM_DISK        1
#ifdef ICGEN_PREVOLUTION
#undef ICGEN_PREVOLUTION
#define ICGEN_PREVOLUTION           2
#endif
#ifdef ICGEN_SONG
#undef ICGEN_SONG
#define ICGEN_SONG                  3
#endif
#ifdef ICGEN_FALCONIC
#undef ICGEN_FALCONIC
#define ICGEN_FALCONIC              4
#endif

#define VECTOR_PARABOLIC            0
#define VECTOR_ELLIPTIC             1

// Physical constants
#define C_PLANCK_LAW      4.48147e-7    // omega_g / (T_cmb [K])^4
#define C_BOLTZMANN_CST   8.61733e-5    // Boltzmann constant [eV/K]
#define C_SPEED_OF_LIGHT  2997.92458    // speed of light [100 km/s]
#define C_RHO_CRIT        2.77459457e11 // critical density [M_sun h^2 / Mpc^3]
#define C_FD_NORM         1.80308535    // Integral[q*q/(exp(q)+1), 0, infinity]

// default physical parameters (used in parser.hpp)
#define P_HUBBLE          0.67556       // default value for h
#define P_T_NCDM          0.71611       // default value for T_ncdm
#define P_NCDM_MASS_OMEGA 93.14         // m_ncdm / omega_ncdm [eV]
#define P_N_UR            3.046         // default value for N_ur (= N_eff)
#define P_SPECTRAL_AMP    2.215e-9      // default value for A_s
#define P_SPECTRAL_INDEX  0.9619        // default value for n_s
#define P_PIVOT_SCALE     0.05          // default pivot scale [Mpc^-1]

#ifndef GADGET_LENGTH_CONVERSION
#define GADGET_LENGTH_CONVERSION 0.001  // Gadget length unit in Mpc / h
#endif
#ifndef GADGET_MASS_CONVERSION
#define GADGET_MASS_CONVERSION 1.0e10   // Gadget mass unit in M_sun / h
#endif
#ifndef GADGET_VELOCITY_CONVERSION
#define GADGET_VELOCITY_CONVERSION 3.335640952e-6  // Gadget velocity unit / speed of light
#endif
#ifndef GADGET_ID_BYTES
#define GADGET_ID_BYTES 8
#endif

#ifdef EXTERNAL_IO
#ifndef NUMBER_OF_IO_FILES
#define NUMBER_OF_IO_FILES 4
#endif
#endif

// Hibernation particle save mode
#define HIB_SAVE_HDF5 0
#define HIB_SAVE_GADGET2 1

// Modified gravity
#define FLAG_FR 1   // f(R)

#define FR_TYPE_RN 1
#define FR_TYPE_R2 2
#define FR_TYPE_HU_SAWICKI 3
#define FR_TYPE_DELTA 4

#define FR_EPSILON_BACKGROUND_DEFAULT 0.1
#define FR_EPSILON_FIELDS_DEFAULT 0.1
#define FR_COUNT_MAX_DEFAULT 100
#define FR_TARGET_PRECISION_DEFAULT 1.E-5

#define FR_WRONG 1.E+30
#define FR_WRONG_RETURN 2. * FR_WRONG

// Multigrid
#define MG_SHAPE_V 1
#define MG_SHAPE_W 2
#define PRE_SMOOTHING_DEFAULT 5
#define POST_SMOOTHING_DEFAULT 5

// Relaxation methods
#define METHOD_RELAX 1
#define METHOD_MULTIGRID 2
#define METHOD_FMG 3

#define RELAXATION_VAR_U 1
#define RELAXATION_VAR_XI 2

// Restrict/prolong u or deltaR
#define RESTRICT_SCALARON 1
#define RESTRICT_DELTAR 2


// Relaxation error
#define RELAXATION_ERROR_DEFAULT 1.E-5

// color escape sequences for terminal highlighting (enable with -DCOLORTERMINAL)
#ifdef COLORTERMINAL
#define COLORTEXT_WHITE     "\033[37;1m"
#define COLORTEXT_CYAN      "\033[36;1m"
#define COLORTEXT_GREEN     "\033[32;1m"
#define COLORTEXT_RED       "\033[31;1m"
#define COLORTEXT_YELLOW    "\033[33;1m"
#define COLORTEXT_RESET     "\033[0m"
#else
#define COLORTEXT_WHITE     '\0'
#define COLORTEXT_CYAN      '\0'
#define COLORTEXT_GREEN     '\0'
#define COLORTEXT_RED       '\0'
#define COLORTEXT_YELLOW    '\0'
#define COLORTEXT_RESET     '\0'
#endif

// header structure for GADGET-2 files [V. Springel, N. Yoshida, and S.D. White, New Astron. 6 (2001) 79
// and V. Springel, Mon. Not. R. Astron. Soc. 364 (2005) 1105]

#ifndef GADGET2_HEADER
#define GADGET2_HEADER
struct gadget2_header
{
	uint32_t npart[6];
	double mass[6];
	double time;
	double redshift;
	int32_t flag_sfr;
	int32_t flag_feedback;
	uint32_t npartTotal[6];
	int32_t flag_cooling;
	int32_t num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];   /* fills to 256 Bytes */
};
#endif

struct metadata
{
	int CYCLE_INFO_INTERVAL; // Previously fixed: #define CYCLE_INFO_INTERVAL = 10
	int BACKGROUND_NUMPTS;
	int numpts;
	int downgrade_factor;
	long numpcl[MAX_PCL_SPECIES];
	int tracer_factor[MAX_PCL_SPECIES];
	int baryon_flag;
	int relativistic_flag;
	int modified_gravity_flag;
	int fR_type;
	int vector_flag;
	int radiation_flag;
	int out_pk;
	int out_snapshot;
	int out_check;
	int num_pk;
	int numbins;
	int num_snapshot;
	int num_restart;
	int num_fR_params;
	int quasi_static;
	int background_only;
	double z_fin;
	int lcdm_background;
	double Cf;
	int perturbations_are_GR;
	int check_fields;
	int check_fields_precision;
	int check_pause;
	int hibernation_save_mode;
	double fR_params[MAX_FR_PARAMS];
	double fR_epsilon_bg;
	double fR_epsilon_fields;
	double fR_target_precision;
	double fR_count_max;
	double relaxation_error;
	double overrelaxation_coeff;
	int relaxation_method;
	int relaxation_variable;
	int red_black;
	int pre_smoothing;
	int post_smoothing;
	int multigrid_or_not;
	int multigrid_restrict_mode;
	int multigrid_shape;
	int multigrid_n_grids;
	int multigrid_n_cycles;
	int multigrid_check_shape;
	double movelimit;
	double steplimit;
	double boxsize;
	double wallclocklimit;
	double z_in;
	double z_check; // TODO: Remove after debugging
	double z_snapshot[MAX_OUTPUTS];
	double z_pk[MAX_OUTPUTS];
	double z_restart[MAX_OUTPUTS];
	double z_switch_fR_background;
	double z_switch_deltarad;
	double z_switch_linearchi;
	double z_switch_deltancdm[MAX_PCL_SPECIES-2];
	double z_switch_Bncdm[MAX_PCL_SPECIES-2];
	char basename_snapshot[PARAM_MAX_LENGTH];
	char basename_pk[PARAM_MAX_LENGTH];
	char basename_generic[PARAM_MAX_LENGTH];
	char output_path[PARAM_MAX_LENGTH];
	char restart_path[PARAM_MAX_LENGTH];
	char basename_restart[PARAM_MAX_LENGTH];
};

std::ostream& operator<< (std::ostream& os, const metadata& sim)
{
	os << "===== Metadata structure: =====\n";
	os << "CYCLE_INFO_INTERVAL: " << sim.CYCLE_INFO_INTERVAL << "\n";
	os << "BACKGROUND_NUMPTS: " << sim.BACKGROUND_NUMPTS << "\n";
	os << "numpts: " << sim.numpts << "\n";
	os << "downgrade_factor: " << sim.downgrade_factor << "\n";
	os << "numpcl: " << sim.numpcl[0];
	for(int i=1; i<MAX_PCL_SPECIES; i++)
	{
		os << " , " << sim.numpcl[i];
	}
	os << "\n";

	os << "tracer_factor: " << sim.tracer_factor[0];
	for(int i=1; i<MAX_PCL_SPECIES; i++)
	{
		os << " , " << sim.tracer_factor[i];
	}
	os << "\n";
	os << "baryon_flag: " << sim.baryon_flag << "\n";
	os << "relativistic_flag: " << sim.relativistic_flag << "\n";
	os << "modified_gravity_flag: " << sim.modified_gravity_flag << "\n";
	os << "fR_type: " << sim.fR_type << "\n";
	os << "vector_flag: " << sim.vector_flag << "\n";
	os << "radiation_flag: " << sim.radiation_flag << "\n";
	os << "out_pk: " << sim.out_pk << "\n";
	os << "out_snapshot: " << sim.out_snapshot << "\n";
	os << "num_pk: " << sim.num_pk << "\n";
	os << "numbins: " << sim.numbins << "\n";
	os << "num_snapshot: " << sim.num_snapshot << "\n";
	os << "num_restart: " << sim.num_restart << "\n";
	os << "num_fR_params: " << sim.num_fR_params << "\n";
	os << "quasi_static: " << sim.quasi_static << "\n";
	os << "background_only: " << sim.background_only << "\n";
	os << "background final redshift: " << sim.z_fin << "\n";
	os << "lcdm_background: " << sim.lcdm_background << "\n";

	os << "perturbations_are_GR: " << sim.perturbations_are_GR << "\n";
	os << "check_fields: " << sim.check_fields << "\n";
	os << "check_pause: " << sim.check_pause << "\n";

	os << "Cf: " << sim.Cf << "\n";

	os << "fR_params: " << sim.fR_params[0];
	for(int i=1; i<MAX_FR_PARAMS; i++) os << " , " << sim.fR_params[i];
	os << "\n";

	os << "fR_epsilon_bg: " << sim.fR_epsilon_bg << "\n";
	os << "fR_epsilon_fields: " << sim.fR_epsilon_fields << "\n";
	os << "fR_target_precision: " << sim.fR_target_precision << "\n";
	os << "multigrid_shape: " << sim.multigrid_shape << "\n";
	os << "pre_smoothing: " << sim.pre_smoothing << "\n";
	os << "post_smoothing: " << sim.post_smoothing << "\n";
	os << "relaxation_error: " << sim.relaxation_error << "\n";
	os << "overrelaxation coefficient: " << sim.overrelaxation_coeff << "\n";
	os << "movelimit: " << sim.movelimit << "\n";
	os << "steplimit: " << sim.steplimit << "\n";
	os << "boxsize: " << sim.boxsize << "\n";
	os << "wallclocklimit: " << sim.wallclocklimit << "\n";
	os << "z_in: " << sim.z_in << "\n";
	os << "z_check: " << sim.z_check << "\n";
	os << "z_snapshot: " << sim.z_snapshot[0];
	for(int i=1; i<MAX_OUTPUTS; i++)
	{
		os << " , " << sim.z_snapshot[i];
	}
	os << "\n";

	os << "z_pk: " << sim.z_pk[0];
	for(int i=1; i<MAX_OUTPUTS; i++)
	{
		os << " , " << sim.z_pk[i];
	}
	os << "\n";

	os << "z_restart: " << sim.z_restart[0];
	for(int i=1; i<MAX_OUTPUTS; i++)
	{
		os << " , " << sim.z_restart[i];
	}
	os << "\n";

	os << "z_switch_fR_background: " << sim.z_switch_fR_background << "\n";
	os << "z_switch_deltarad: " << sim.z_switch_deltarad << "\n";
	os << "z_switch_linearchi: " << sim.z_switch_linearchi << "\n";

	os << "z_switch_deltancdm: " << sim.z_switch_deltancdm[0];
	for(int i=1; i<MAX_PCL_SPECIES-2; i++) os << " , " << sim.z_switch_deltancdm[i];
	os << "\n";

	os << "z_switch_Bncdm: " << sim.z_switch_Bncdm[0];
	for(int i=1; i<MAX_PCL_SPECIES-2; i++) os << " , " << sim.z_switch_Bncdm[i];
	os << "\n";

	os << "basename_snapshot: " << sim.basename_snapshot << "\n";
	os << "basename_pk: " << sim.basename_pk << "\n";
	os << "basename_generic: " << sim.basename_generic << "\n";
	os << "output_path: " << sim.output_path << "\n";
	os << "restart_path: " << sim.restart_path << "\n";
	os << "basename_restart: " << sim.basename_restart << "\n";
	os << "---------------------------------\n";

	return os;
}

struct icsettings
{
	int numtile[MAX_PCL_SPECIES];
	int seed;
	int flags;
	int generator;
	int restart_cycle;
	char pclfile[MAX_PCL_SPECIES][PARAM_MAX_LENGTH];
	char pkfile[PARAM_MAX_LENGTH];
	char tkfile[PARAM_MAX_LENGTH];
	char metricfile[METRICFILE_LENGTH][PARAM_MAX_LENGTH];
 	double restart_tau;
	double restart_dtau;
	double restart_version;
	double z_ic;
	double z_relax;
	double Cf;
	double A_s;
	double n_s;
	double k_pivot;

	// f(R) restart
	double restart_dtau_old;
	double restart_dtau_old_2;
	double restart_dtau_osci;
	double restart_dtau_bg;
	double restart_a;
	double restart_Hubble;
	double restart_Rbar;
	double restart_dot_Rbar;
};

struct cosmology
{
	double Omega_cdm;
	double Omega_b;
	double Omega_m;
	double Omega_Lambda;
	double Omega_g;
	double Omega_ur;
	double Omega_rad; //sum og g and ur
	double Omega_ncdm[MAX_PCL_SPECIES-2];
	double h;
	double m_ncdm[MAX_PCL_SPECIES-2];//mass in eV
	double T_ncdm[MAX_PCL_SPECIES-2];
	double deg_ncdm[MAX_PCL_SPECIES-2];
	int num_ncdm;
};

std::ostream& operator<< (std::ostream& os, const cosmology& cosmo)
{
	os << "===== Cosmology structure: =====\n";
	os << "Omega_cdm: " << cosmo.Omega_cdm << "\n";
	os << "Omega_b: " << cosmo.Omega_b << "\n";
	os << "Omega_m: " << cosmo.Omega_m << "\n";
	os << "Omega_Lambda: " << cosmo.Omega_Lambda << "\n";
	os << "Omega_g: " << cosmo.Omega_g << "\n";
	os << "Omega_ur: " << cosmo.Omega_ur << "\n";
	os << "Omega_rad: " << cosmo.Omega_rad << "\n";

	os << "Omega_ncdm: " << cosmo.Omega_ncdm[0];
	for(int i=1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.Omega_ncdm[i];
	os << "\n";

	os << "h: " << cosmo.h << "\n";

	os << "m_ncdm: " << cosmo.m_ncdm[0];
	for(int i=1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.m_ncdm[i];
	os << "\n";

	os << "T_ncdm: " << cosmo.T_ncdm[0];
	for(int i=1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.T_ncdm[i];
	os << "\n";

	os << "deg_ncdm: " << cosmo.deg_ncdm[0];
	for(int i=1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.deg_ncdm[i];
	os << "\n";

	os << "num_ncdm: " << cosmo.num_ncdm << "\n";
	os << "---------------------------------\n";

	return os;
}

#endif
