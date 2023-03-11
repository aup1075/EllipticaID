/* (c) 2023 Pedro Espino & Alireza Rashti */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <vector>
#include <stdio.h>
#include <cassert>
#include <cstdlib>
#include <elliptica_id_reader_lib.h>

#define MAX_NTAB 16001
#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif

using namespace std;

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
       log_rho0_table[i]=log10(rho0);     /* multiply by C^2 to get energy density */
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

static void set_dt_from_domega (CCTK_ARGUMENTS,
                                CCTK_REAL const* const var,
                                CCTK_REAL      * const dtvar,
                                CCTK_REAL const& omega)
{
  DECLARE_CCTK_ARGUMENTS;

  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  vector<CCTK_REAL> dxvar(npoints), dyvar(npoints);

  Diff_gv (cctkGH, 0, var, &dxvar[0], -1);
  Diff_gv (cctkGH, 1, var, &dyvar[0], -1);

#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    CCTK_REAL const ephix = +y[i];
    CCTK_REAL const ephiy = -x[i];
    CCTK_REAL const dphi_var = ephix * dxvar[i] + ephiy * dyvar[i];
    dtvar[i] = omega * dphi_var;
  }
}

void Elliptica_BHNS_initialize(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Setting up LORENE Bin_NS initial data");
  if(init_real){CCTK_INFO ("(with realistic EOS)"); load_beta_equil(beta_file, log_rho0_tab_beta, Y_e_tab, &n_tab_beta);}

  //TODO: fix these to proper NIST values

  // Other quantities in terms of Cactus units
  CCTK_INT keyerr = 0, anyerr = 0;

  //Get EOS_Omni handle 
  if (!(*init_eos_key = EOS_Omni_GetHandle(eos_table)))
    CCTK_WARN(0,"Cannot get initial eos handle, aborting...");
  CCTK_VInfo(CCTK_THORNSTRING, "Elliptica_BHNS will use the %s equation of state.", eos_table);
  CCTK_VInfo(CCTK_THORNSTRING, "Elliptica_BHNS will use the %d eos handle", *init_eos_key);

  //Set up local coordinate arrays
  CCTK_INFO ("Setting up coordinates");
  int const N_points = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  double *xx = (double *)calloc(N_points,sizeof(*xx)); assert(xx);
  double *yy = (double *)calloc(N_points,sizeof(*yy)); assert(yy);
  double *zz = (double *)calloc(N_points,sizeof(*zz)); assert(zz);
  
#pragma omp parallel for
  for (int i=0; i<N_points; ++i) {
    xx[i] = x[i];
    yy[i] = y[i];
    zz[i] = z[i];
  }

  // --------------------------------------------------------------
  //   CHECKING FILE NAME EXISTENCE
  // --------------------------------------------------------------
  FILE *file;
  // TODO: THERE SHOULD BE A PARAMETER CALLED filename
  const char *filename="/home/pespino/temp_project/projects/Elliptica_TEST/BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_00/BHNS_BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_3x3x3_00/checkpoint.dat";
  //TODO: THERE SHOULD BE A PARAMETER CALLED option
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
    Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(filename,option);
    //Print some scalar values from idr
    CCTK_REAL Omega = 0.011592023119370;
    CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:       %g", Omega);
//    CCTK_VInfo (CCTK_THORNSTRING, "dist [km]:           %g", idr->dist);
//    CCTK_VInfo (CCTK_THORNSTRING, "dist_mass [km]:      %g", idr->dist_mass);
//    CCTK_VInfo (CCTK_THORNSTRING, "mass1_b [M_sun]:     %g", idr->mass1_b);
//    CCTK_VInfo (CCTK_THORNSTRING, "mass2_b [M_sun]:     %g", idr->mass2_b);
//    CCTK_VInfo (CCTK_THORNSTRING, "mass_ADM [M_sun]:    %g", idr->mass_adm);
//    CCTK_VInfo (CCTK_THORNSTRING, "L_tot [G M_sun^2/c]: %g", idr->angu_mom);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_comp [km]:    %g", idr->rad1_x_comp);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad1_y [km]:         %g", idr->rad1_y);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad1_z [km]:         %g", idr->rad1_z);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_opp [km]:     %g", idr->rad1_x_opp);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_comp [km]:    %g", idr->rad2_x_comp);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad2_y [km]:         %g", idr->rad2_y);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad2_z [km]:         %g", idr->rad2_z);
//    CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_opp [km]:     %g", idr->rad2_x_opp);
    double K = 100.0; // make sure ths is in polytropic units
    double Gamma = 2.0; // make sure ths is in polytropic units
    idr->ifields = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,grhd_rho,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
    idr->npoints = N_points;
//    double x[2] = {1,2};
    idr->x_coords = xx;
    idr->y_coords = yy;
    idr->z_coords = zz;
    idr->set_param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
    idr->set_param("ADM_B1I_form","zero",idr);
    
    
    elliptica_id_reader_interpolate(idr);

//    elliptica_id_reader_free(idr);
#pragma omp parallel for
  for (int i=0; i<N_points; ++i) {

    if (CCTK_EQUALS(initial_lapse, "Elliptica_BHNS")) { 
      alp[i] = idr->field[idr->indx("alpha")][i];
    }

    //TODO: this is modified by a negative sign from LORENE ID. Is that needed here?
    if (CCTK_EQUALS(initial_shift, "Elliptica_BHNS")) { 
      betax[i] = idr->field[idr->indx("betax")][i];
      betay[i] = idr->field[idr->indx("betay")][i];
      betaz[i] = idr->field[idr->indx("betaz")][i];
    }

    if (CCTK_EQUALS(initial_data, "Elliptica_BHNS")) {
      gxx[i] = idr->field[idr->indx("adm_gxx")][i];
      gxy[i] = idr->field[idr->indx("adm_gxy")][i];
      gxz[i] = idr->field[idr->indx("adm_gxz")][i];
      gyy[i] = idr->field[idr->indx("adm_gyy")][i];
      gyz[i] = idr->field[idr->indx("adm_gyz")][i];
      gzz[i] = idr->field[idr->indx("adm_gzz")][i];

      kxx[i] = idr->field[idr->indx("adm_Kxx")][i];
      kxy[i] = idr->field[idr->indx("adm_Kxy")][i];
      kxz[i] = idr->field[idr->indx("adm_Kxz")][i];
      kyy[i] = idr->field[idr->indx("adm_Kyy")][i];
      kyz[i] = idr->field[idr->indx("adm_Kyz")][i];
      kzz[i] = idr->field[idr->indx("adm_Kzz")][i];
    }

    if (CCTK_EQUALS(initial_data, "Elliptica_BHNS")) {
      rho[i] = idr->field[idr->indx("grhd_rho")][i];
      if(init_real){
      if(rho[i]>=1e-7){
          double yeres = interp(log_rho0_tab_beta, Y_e_tab, n_tab_beta,log10(rho[i]), &n_nearest_beta);
	  if(yeres <=0.036){yeres = 0.036;} 
	  Y_e[i] = yeres;
          temperature[i] = 0.1;
      }
      else{Y_e[i] = 0.25;}
      temperature[i] = 0.1;
      }
      if (!recalculate_eps){ //we don't know the temperature, so assume epsilon from ID is correct.
        eps[i] = idr->field[idr->indx("grhd_epsl")][i];
      }
      // Pressure from EOS_Omni call 
      if (CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::temperature") > 0 &&
          CCTK_ActiveTimeLevelsVN(cctkGH, "HydroBase::Y_e") > 0)
      {
        EOS_Omni_press(*init_eos_key,recalculate_eps,eos_precision,1,&(rho[i]),&(eps[i]),
                       &(temperature[i]),&(Y_e[i]),&(press[i]),&keyerr,&anyerr);


      }
      else
      {
        EOS_Omni_press(*init_eos_key,recalculate_eps,eos_precision,1,&(rho[i]),&(eps[i]),
                       NULL,NULL,&(press[i]),&keyerr,&anyerr);
      }

      vel[i          ] = idr->field[idr->indx("grhd_vx")][i];
      vel[i+  N_points] = idr->field[idr->indx("grhd_vy")][i];
      vel[i+2*N_points] = idr->field[idr->indx("grhd_vz")][i];

      // Especially the velocity is set to strange values outside of the
      // matter region, so take care of this in the following way
      if (rho[i] < 1.e-15) {
        rho[i          ] = 1.e-15;
        vel[i          ] = 0.0;
        vel[i+  N_points] = 0.0;
        vel[i+2*N_points] = 0.0;
        eps[i          ] = K * pow(rho[i], Gamma-1.) / (Gamma-1.);
        press[i        ] = K * pow(rho[i], Gamma);
      }
    }

  } // for i

  {
    // Angular velocity
    //CCTK_REAL const omega = idr->omega * cactusT;

    // These initial data assume a helical Killing vector field

    if (CCTK_EQUALS(initial_lapse, "Elliptica_BHNS")) { 
      if (CCTK_EQUALS (initial_dtlapse, "Elliptica_BHNS")) {
        CCTK_INFO ("Calculating time derivatives of lapse");
        set_dt_from_domega (CCTK_PASS_CTOC, alp, dtalp, Omega);
      } else if (CCTK_EQUALS (initial_dtlapse, "none") or CCTK_EQUALS(initial_dtlapse,"zero")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error setting dtlapse");
      }
    }

    if (CCTK_EQUALS(initial_shift, "Elliptica_BHNS")) { 
      if (CCTK_EQUALS (initial_dtshift, "Elliptica_BHNS")) {
        CCTK_INFO ("Calculating time derivatives of shift");
        set_dt_from_domega (CCTK_PASS_CTOC, betax, dtbetax, Omega);
        set_dt_from_domega (CCTK_PASS_CTOC, betay, dtbetay, Omega);
        set_dt_from_domega (CCTK_PASS_CTOC, betaz, dtbetaz, Omega);
      } else if (CCTK_EQUALS (initial_dtshift, "none") or CCTK_EQUALS(initial_dtshift,"zero")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error setting dtshift");
      }
    }
  }

  CCTK_INFO ("Done.");
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial data from file '%s': %s", filename, e.what());
  }

  // free  
  if (xx) free(xx);
  if (yy) free(yy);
  if (zz) free(zz);
}

