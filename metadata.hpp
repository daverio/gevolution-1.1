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

#define GEVOLUTION_VERSION 1.1

#ifndef MAX_OUTPUTS
#define MAX_OUTPUTS 32
#endif

#ifndef MAX_FOFR_PARAMS
#define MAX_FOFR_PARAMS 8
#endif

#ifndef PARAM_MAX_LENGTH
#define PARAM_MAX_LENGTH 256
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

#define MASK_PHI       1
#define MASK_CHI       2
#define MASK_POT       4
#define MASK_B         8
#define MASK_T00       16
#define MASK_TIJ       32
#define MASK_RBARE     64
#define MASK_HIJ       128
#define MASK_P         256
#define MASK_GADGET    512
#define MASK_PCLS      1024
#define MASK_XSPEC     2048
#define MASK_DELTA     4096
#define MASK_DBARE     8192
#define MASK_XI  		   16384
#define MASK_ZETA      32768
#define MASK_DELTAR	   65536
#define MASK_DELTAT		 131072

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

//Modified gravity types
#define GENREL 0 // General Relativity
#define FOFR 1   // f(R)

//f(R) models
#define FOFR_TYPE_RN 1
#define FOFR_TYPE_R2 2
#define FOFR_TYPE_HU_SAWICKI 3
#define FOFR_TYPE_DELTA 4


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
	int CYCLE_INFO_INTERVAL = 10; // Previously fixed. #define CYCLE_INFO_INTERVAL = 10
	int BACKGROUND_NUMPTS = 10;
	int numpts;
	int downgrade_factor;
	long numpcl[MAX_PCL_SPECIES];
	int tracer_factor[MAX_PCL_SPECIES];
	int baryon_flag;
	int gr_flag;
	int mg_flag;
	int fofR_type;
	int vector_flag;
	int radiation_flag;
	int out_pk;
	int out_snapshot;
	int num_pk;
	int numbins;
	int num_snapshot;
	int num_restart;
	int num_fofR_params;
	int quasi_static = 0;
	int background_only = 0;
	int background_trace = 0;
	double bg_initial_redshift = 100.;
	double bg_final_redshift = 0.;
	int lcdm_background = 0;

	int S0i_mode;// TODO remove after FULL debugging
	int back_to_GR = 0;
	int check_fields = 0;
	int check_pause = 1;

	double Cf;
	double fofR_params[MAX_FOFR_PARAMS];
	double fofR_epsilon_bg;
	double fofR_epsilon_fields;
	double fofR_target_precision;
	double movelimit;
	double steplimit;
	double boxsize;
	double wallclocklimit;
	double z_in;
	double z_check = -1; // TODO: Remove after debugging
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
	char background_filename[PARAM_MAX_LENGTH];
};

std::ostream& operator<< (std::ostream& os, const metadata& sim)
{
	os << "===== Metadata structure: =====\n";
	os << "CYCLE_INFO_INTERVAL: " << sim.CYCLE_INFO_INTERVAL << "\n";
	os << "BACKGROUND_NUMPTS: " << sim.BACKGROUND_NUMPTS << "\n";
	os << "numpts: " << sim.numpts << "\n";
	os << "downgrade_factor: " << sim.downgrade_factor << "\n";

	os << "numpcl: " << sim.numpcl[0];
	for(int i = 1; i<MAX_PCL_SPECIES; i++) os << " , " << sim.numpcl[i];
	os << "\n";

	os << "tracer_factor: " << sim.tracer_factor[0];
	for(int i = 1; i<MAX_PCL_SPECIES; i++) os << " , " << sim.tracer_factor[i];
	os << "\n";

	os << "baryon_flag: " << sim.baryon_flag << "\n";
	os << "gr_flag: " << sim.gr_flag << "\n";
	os << "mg_flag: " << sim.mg_flag << "\n";
	os << "fofR_type: " << sim.fofR_type << "\n";
	os << "vector_flag: " << sim.vector_flag << "\n";
	os << "radiation_flag: " << sim.radiation_flag << "\n";
	os << "out_pk: " << sim.out_pk << "\n";
	os << "out_snapshot: " << sim.out_snapshot << "\n";
	os << "num_pk: " << sim.num_pk << "\n";
	os << "numbins: " << sim.numbins << "\n";
	os << "num_snapshot: " << sim.num_snapshot << "\n";
	os << "num_restart: " << sim.num_restart << "\n";
	os << "num_fofR_params: " << sim.num_fofR_params << "\n";
	os << "quasi_static: " << sim.quasi_static << "\n";
	os << "background_only: " << sim.background_only << "\n";
	os << "background_trace: " << sim.background_trace << "\n";
	os << "bg_initial_redshift: " << sim.bg_initial_redshift << "\n";
	os << "bg_final_redshift: " << sim.bg_final_redshift << "\n";
	os << "lcdm_background: " << sim.lcdm_background << "\n";

	os << "S0i_mode: " << sim.S0i_mode << "\n";
	os << "back_to_GR: " << sim.back_to_GR << "\n";
	os << "check_fields: " << sim.check_fields << "\n";
	os << "check_pause: " << sim.check_pause << "\n";

	os << "Cf: " << sim.Cf << "\n";

	os << "fofR_params: " << sim.fofR_params[0];
	for(int i = 1; i<MAX_FOFR_PARAMS; i++) os << " , " << sim.fofR_params[i];
	os << "\n";

	os << "fofR_epsilon_bg: " << sim.fofR_epsilon_bg << "\n";
	os << "fofR_epsilon_fields: " << sim.fofR_epsilon_fields << "\n";
	os << "fofR_target_precision: " << sim.fofR_target_precision << "\n";
	os << "movelimit: " << sim.movelimit << "\n";
	os << "steplimit: " << sim.steplimit << "\n";
	os << "boxsize: " << sim.boxsize << "\n";
	os << "wallclocklimit: " << sim.wallclocklimit << "\n";
	os << "z_in: " << sim.z_in << "\n";
	os << "z_check: " << sim.z_check << "\n";

	os << "z_snapshot: " << sim.z_snapshot[0];
	for(int i = 1; i<MAX_OUTPUTS; i++) os << " , " << sim.z_snapshot[i];
	os << "\n";

	os << "z_pk: " << sim.z_pk[0];
	for(int i = 1; i<MAX_OUTPUTS; i++) os << " , " << sim.z_pk[i];
	os << "\n";

	os << "z_restart: " << sim.z_restart[0];
	for(int i = 1; i<MAX_OUTPUTS; i++) os << " , " << sim.z_restart[i];
	os << "\n";

	os << "z_switch_fR_background: " << sim.z_switch_fR_background << "\n";
	os << "z_switch_deltarad: " << sim.z_switch_deltarad << "\n";
	os << "z_switch_linearchi: " << sim.z_switch_linearchi << "\n";

	os << "z_switch_deltancdm: " << sim.z_switch_deltancdm[0];
	for(int i = 1; i<MAX_PCL_SPECIES-2; i++) os << " , " << sim.z_switch_deltancdm[i];
	os << "\n";

	os << "z_switch_Bncdm: " << sim.z_switch_Bncdm[0];
	for(int i = 1; i<MAX_PCL_SPECIES-2; i++) os << " , " << sim.z_switch_Bncdm[i];
	os << "\n";

	os << "basename_snapshot: " << sim.basename_snapshot << "\n";
	os << "basename_pk: " << sim.basename_pk << "\n";
	os << "basename_generic: " << sim.basename_generic << "\n";
	os << "output_path: " << sim.output_path << "\n";
	os << "restart_path: " << sim.restart_path << "\n";
	os << "basename_restart: " << sim.basename_restart << "\n";
	os << "background_filename: " << sim.background_filename << "\n";
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

	///fofR restart
	double restart_dtau_old;
	double restart_dtau_old_2;
	double restart_dtau_osci;
	double restart_dtau_bg;
	double restart_a;
	double restart_Hubble;
	double restart_Rbar;
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
	for(int i = 1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.Omega_ncdm[i];
	os << "\n";

	os << "h: " << cosmo.h << "\n";

	os << "m_ncdm: " << cosmo.m_ncdm[0];
	for(int i = 1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.m_ncdm[i];
	os << "\n";

	os << "T_ncdm: " << cosmo.T_ncdm[0];
	for(int i = 1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.T_ncdm[i];
	os << "\n";

	os << "deg_ncdm: " << cosmo.deg_ncdm[0];
	for(int i = 1; i<MAX_PCL_SPECIES-2; i++) os << " , " << cosmo.deg_ncdm[i];
	os << "\n";

	os << "num_ncdm: " << cosmo.num_ncdm << "\n";
	os << "---------------------------------\n";

	return os;
}

#endif
