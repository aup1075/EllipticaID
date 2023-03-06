/* (c) 2023 Pedro Espino */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <stdio.h>
#include <elliptica_id_reader_lib.h>
//#include <idr_main.h>
//#include <idr_main.c>
#define MAX_NTAB 16001
#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif

/*C*/
/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x > xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,n_nearest_pt);

 k=IMIN(IMAX((*n_nearest_pt)-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}
/*************************************************************************/
/* Load Beta equil file.                                                        */
/*************************************************************************/
void load_beta_equil( const char beta_equil_file[],
               double log_rho0_table[MAX_NTAB],
               double Y_e_table[MAX_NTAB],
               int *n_tab_beta)
{
 int i;                    /* counter */

 /*constants to convert from cgs to cactus units c=G=M_sun=1.0 */
 CCTK_REAL const cactusM= (5.028916268544129e-34);    /*  1/g  */
 CCTK_REAL const cactusL= (6.772400341316594e-06);    /*  1/cm */
 CCTK_REAL const cactusT= (2.0303145448833407e5);    /*  1/s  */
 CCTK_REAL const cactusV= (1.0/(cactusL*cactusL*cactusL));
 double rho0,               /* density */
        ye;                /* electron fraction */

 FILE *f_beta;              /* pointer to beta_equil_file */


    /* OPEN FILE TO READ */

    if((f_beta=fopen(beta_equil_file,"r")) == NULL ) {
       CCTK_VERROR("cannot open beta-equil. file:  %s\n",beta_equil_file);
    }

    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_beta,"%d\n",n_tab_beta);

    /* READ ENERGY DENSITY, P, H, N0 AND CONVERT TO CACTUS UNITS */

    for(i=1;i<=(*n_tab_beta);i++) {
       fscanf(f_beta,"%lf %lf \n",&rho0,&ye) ;
       log_rho0_table[i]=log10(rho0*cactusM*cactusV);     /* multiply by C^2 to get energy density */
       if(ye <= 0.036){
           Y_e_table[i] = 0.036;
       }
       else{
           Y_e_table[i]=ye;
       }
    }
}

//Auxiliary variables for interepolating Beta equilibrium table
int n_tab_beta;
double Y_e_tab[MAX_NTAB], log_rho0_tab_beta[MAX_NTAB];
int n_nearest_beta;

extern "C"
void Elliptica_BHNS_initialize(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up LORENE Bin_NS initial data");
  if(init_real){CCTK_INFO ("(with realistic EOS)"); load_beta_equil(beta_file, log_rho0_tab_beta, Y_e_tab, &n_tab_beta);}

  //TODO: fix these to proper NIST values
  CCTK_REAL const c_light  = 3e8;      // speed of light [m/s]
  CCTK_REAL const nuc_dens = 1e12; // Nuclear density as used in Lorene units [kg/m^3]
  CCTK_REAL const G_grav   = 6e-21;      // gravitational constant [m^3/kg/s^2]
  CCTK_REAL const M_sun    = 2e30;   // solar mass [kg]

  // Cactus units in terms of SI units:
  // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
  // for the magnetic field)
  CCTK_REAL const cactusM = M_sun;
  CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light,2);
  CCTK_REAL const cactusT = cactusL / c_light;

  // Other quantities in terms of Cactus units
  CCTK_REAL const coord_unit = cactusL / 1.0e+3;         // from km (~1.477)
  CCTK_REAL const rho_unit   = cactusM / pow(cactusL,3); // from kg/m^3
  CCTK_INT keyerr = 0, anyerr = 0;

  //Get EOS_Omni handle 
  if (!(*init_eos_key = EOS_Omni_GetHandle(eos_table)))
    CCTK_WARN(0,"Cannot get initial eos handle, aborting...");
  CCTK_VInfo(CCTK_THORNSTRING, "Meudon_Bin_NS will use the %s equation of state.", eos_table);
  CCTK_VInfo(CCTK_THORNSTRING, "Meudon_Bin_NS will use the %d eos handle", *init_eos_key);

  //Set up local coordinate arrays
  CCTK_INFO ("Setting up coordinates");
  int const N_points = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<double> xx(N_points), yy(N_points), zz(N_points);
#pragma omp parallel for
  for (int i=0; i<N_points; ++i) {
    xx[i] = x[i] * coord_unit;
    yy[i] = y[i] * coord_unit;
    zz[i] = z[i] * coord_unit;
  }

  // --------------------------------------------------------------
  //   CHECKING FILE NAME EXISTENCE
  // --------------------------------------------------------------
  FILE *file;
  //THERE SHOULD BE A PARAMETER CALLED filename
  const char *filename="/home/pespino/temp_project/projects/Elliptica_TEST/BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_00/BHNS_BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_3x3x3_00/checkpoint.dat";
  const char *filename="/home/pespino/temp_project/projects/Elliptica_TEST/BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_00/BHNS_BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_3x3x3_00/checkpoint.dat"
  //THERE SHOULD BE A PARAMETER CALLED option
  const char *option = "generic";
  if ((file = fopen(filename, "r")) != NULL) 
     fclose(file);
  else {
     CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                 "File \"%s\" does not exist. ABORTING", filename);
  }
  // when using EOS, check for EOS file.
  //if (strlen(eos_table_filepath) > 0) {
  //  if (setenv("LORENE_TABULATED_EOS_PATH", eos_table_filepath, 1)) {
  //    CCTK_ERROR("Unable to set environment variable LORENE_TABULATED_EOS_PATH");

  //  }
  //}

  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);

    
  try {

  try {
    Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(f,option);
    //Print some scalar values from idr
    CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:       %g", idr->omega);
    CCTK_VInfo (CCTK_THORNSTRING, "dist [km]:           %g", idr->dist);
    CCTK_VInfo (CCTK_THORNSTRING, "dist_mass [km]:      %g", idr->dist_mass);
    CCTK_VInfo (CCTK_THORNSTRING, "mass1_b [M_sun]:     %g", idr->mass1_b);
    CCTK_VInfo (CCTK_THORNSTRING, "mass2_b [M_sun]:     %g", idr->mass2_b);
    CCTK_VInfo (CCTK_THORNSTRING, "mass_ADM [M_sun]:    %g", idr->mass_adm);
    CCTK_VInfo (CCTK_THORNSTRING, "L_tot [G M_sun^2/c]: %g", idr->angu_mom);
    CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_comp [km]:    %g", idr->rad1_x_comp);
    CCTK_VInfo (CCTK_THORNSTRING, "rad1_y [km]:         %g", idr->rad1_y);
    CCTK_VInfo (CCTK_THORNSTRING, "rad1_z [km]:         %g", idr->rad1_z);
    CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_opp [km]:     %g", idr->rad1_x_opp);
    CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_comp [km]:    %g", idr->rad2_x_comp);
    CCTK_VInfo (CCTK_THORNSTRING, "rad2_y [km]:         %g", idr->rad2_y);
    CCTK_VInfo (CCTK_THORNSTRING, "rad2_z [km]:         %g", idr->rad2_z);
    CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_opp [km]:     %g", idr->rad2_x_opp);
    double K = 100.0; // make sure ths is in polytropic units
    idr->ifields = "alpha,adm_gxx,adm_gxy";
    idr->npoints = npoints;
//    double x[2] = {1,2};
    idr->x_coords = xx;
    idr->y_coords = yy;
    idr->z_coords = zz;
    idr->param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
    idr->param("ADM_B1I_form","zero",idr);
    
    
    elliptica_id_reader_interpolate(idr);

//    
//    printf("%g %g\n",
//        idr->field[idr->indx("alpha")][0],
//        idr->field[idr->indx("alpha")][1]);
//
//    printf("%g %g\n",
//        idr->field[idr->indx("adm_gxx")][0],
//        idr->field[idr->indx("adm_gxx")][1]);
//
//    printf("%g %g\n",
//        idr->field[idr->indx("adm_gxy")][0],
//        idr->field[idr->indx("adm_gxy")][1]);
//
//    elliptica_id_reader_free(idr);
}

