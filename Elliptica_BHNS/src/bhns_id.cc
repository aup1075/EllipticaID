/* (c) 2023 Pedro Espino & Alireza Rashti */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <vector>
#include <stdio.h>
#include <cassert>
#include <cstdlib>

/* @AVIRAL: isn't it better to include this .h file as part of this thorn ? */
#include <elliptica_id_reader_lib.h>

#define MAX_NTAB (16001)
#define IMAX(a,b) ( a>b ? a : b ) 
#define IMIN(a,b) ( a<b ? a : b ) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#ifndef DBL_EPSILON
// TODO: should this be really this much small?
#define DBL_EPSILON 1e-15
#endif

using namespace std;

//TODO: MOVE THESE FUNCTIONS TO LOCAL SCOPE (same for Meudon)
/*C*/
/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/***************************************************************************/
void huntloc(double xx[], int n, double x, int *jlo)
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
double interploc(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 huntloc(xp,np,xb,n_nearest_pt);

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
void load_beta_equilloc( const char beta_equil_file[],
               double log_rho0_table[MAX_NTAB],
               double Y_e_table[MAX_NTAB],
               int *n_tab_betaloc)
{
 int i;                    /* counter */

 /*constants to convert from cgs to cactus units c=G=M_sun=1.0 */
 CCTK_REAL const cactusM= (5.028916268544129e-34);    /*  1/g  */
 CCTK_REAL const cactusL= (6.772400341316594e-06);    /*  1/cm */
 CCTK_REAL const cactusT= (2.0303145448833407e5);    /*  1/s  */

 /* @AVIRAL: Shouldn't the following be 1./cactusV ? */
 CCTK_REAL const cactusV= (1.0/(cactusL*cactusL*cactusL));

 double rho0,               /* density */
        ye;                /* electron fraction */

 FILE *f_beta;              /* pointer to beta_equil_file */


    /* OPEN FILE TO READ */

    if((f_beta=fopen(beta_equil_file,"r")) == NULL ) {
       CCTK_VERROR("cannot open beta-equil. file:  %s\n",beta_equil_file);
    }

    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_beta,"%d\n",n_tab_betaloc);

    /* READ ENERGY DENSITY, P, H, N0 AND CONVERT TO CACTUS UNITS */

    for(i=1;i<=(*n_tab_betaloc);i++) {
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
int n_tab_betaloc;
double Y_e_tabloc[MAX_NTAB], log_rho0_tab_betaloc[MAX_NTAB];
int n_nearest_betaloc;

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

  if(verbose){
    CCTK_INFO ("Entering Elliptica_BHNS_initialize");
  }

  if(init_real){
    CCTK_INFO ("(with realistic EOS)"); 
    load_beta_equilloc(beta_file, log_rho0_tab_betaloc, Y_e_tabloc, &n_tab_betaloc);
    }

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
  if ((file = fopen(Elliptica_bhns_file, "r")) != NULL) 
     fclose(file);
  else {
     CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                 "File \"%s\" does not exist. ABORTING", Elliptica_bhns_file);
  }
  // when using EOS, check for EOS file.
  //if (strlen(eos_table_filepath) > 0) {
  //  if (setenv("LORENE_TABULATED_EOS_PATH", eos_table_filepath, 1)) {
  //    CCTK_ERROR("Unable to set environment variable LORENE_TABULATED_EOS_PATH");

  //  }
  //}

  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", Elliptica_bhns_file);

    
  try {
    CCTK_VInfo (CCTK_THORNSTRING, "Calling Elliptica_ID_Reader_T");
    Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(Elliptica_bhns_file,Elliptica_bhns_option);
    CCTK_VInfo (CCTK_THORNSTRING, "Here1");
    //Print some scalar values from idr
    //BH properties
//    CCTK_REAL BH_Christodoulou_mass = idr->get_param_dbl("BH_Christodoulou_mass",idr);
//    CCTK_VInfo (CCTK_THORNSTRING, "Here2");
//    CCTK_REAL BH_Komar_mass         = idr->get_param_dbl("BH_Komar_mass",idr);
//    CCTK_REAL BH_irr_mass           = idr->get_param_dbl("BH_irreducible_mass",idr);
//    CCTK_REAL BH_comx               = idr->get_param_dbl("BH_center_x",idr);
//    CCTK_REAL BH_comy               = idr->get_param_dbl("BH_center_y",idr);
//    CCTK_REAL BH_comz               = idr->get_param_dbl("BH_center_z",idr);
//    CCTK_REAL BH_spinx              = idr->get_param_dbl("BH_spin_x",idr);
//    CCTK_REAL BH_spiny              = idr->get_param_dbl("BH_spin_y",idr);
//    CCTK_REAL BH_spinz              = idr->get_param_dbl("BH_spin_z",idr);
//    CCTK_REAL BH_boostx             = idr->get_param_dbl("BH_boost_Vx",idr);
//    CCTK_REAL BH_boosty             = idr->get_param_dbl("BH_boost_Vy",idr);
//    CCTK_REAL BH_boostz             = idr->get_param_dbl("BH_boost_Vz",idr);
//   
//    //NS properties
//    CCTK_REAL NS_ADM_mass           = idr->get_param_dbl("NS_ADM_mass",idr);
//    CCTK_REAL NS_max_radius         = idr->get_param_dbl("NS_max_radius",idr);
//    CCTK_REAL NS_baryonic_mass      = idr->get_param_dbl("NS_baryonic_mass",idr);
//    CCTK_REAL NS_rho0_center        = idr->get_param_dbl("NS_rho0_center",idr);
//    CCTK_REAL NS_comx               = idr->get_param_dbl("NS_x_CM",idr);
//    CCTK_REAL NS_comy               = idr->get_param_dbl("NS_y_CM",idr);
//    CCTK_REAL NS_comz               = idr->get_param_dbl("NS_z_CM",idr);
//    CCTK_REAL NS_spinx              = idr->get_param_dbl("NS_spin_x",idr);
//    CCTK_REAL NS_spiny              = idr->get_param_dbl("NS_spin_y",idr);
//    CCTK_REAL NS_spinz              = idr->get_param_dbl("NS_spin_z",idr);
//    //CCTK_STRING NS_EoS_type         = idr->get_param_dbl("NS_EoS_type",idr);
//    //CCTK_STRING NS_EOS_name         = idr->get_param_dbl("NS_EoS_description",idr);
//    CCTK_REAL NS_K0                 = idr->get_param_dbl("NS_EoS_K0",idr);
//    CCTK_REAL NS_Gamma              = idr->get_param_dbl("NS_EoS_Gamma",idr);
//
//    //System properties
//    CCTK_REAL BHNS_Separation       = idr->get_param_dbl("BHNS_separation",idr);
//    CCTK_REAL BHNS_Mass_ratio       = idr->get_param_dbl("BHNS_mass_ratio (BH/NS)",idr);
//    CCTK_REAL BHNS_ADM_mass         = idr->get_param_dbl("BHNS_ADM_mass",idr);
//    CCTK_REAL BHNS_Komar_mass       = idr->get_param_dbl("BHNS_Komar_mass",idr);
//    CCTK_REAL BHNS_comx             = idr->get_param_dbl("BHNS_x_CM",idr);
//    CCTK_REAL BHNS_comy             = idr->get_param_dbl("BHNS_y_CM",idr);
//    CCTK_REAL BHNS_comz             = idr->get_param_dbl("BHNS_z_CM",idr);
//    CCTK_REAL BHNS_Omega            = idr->get_param_dbl("BHNS_angular_velocity",idr);
//    CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:         %g", BHNS_Omega);
//    CCTK_VInfo (CCTK_THORNSTRING, "dist [km]:             %g", BHNS_Separation);
//    CCTK_VInfo (CCTK_THORNSTRING, "mass ratio:            %g", BHNS_Mass_ratio);
//    CCTK_VInfo (CCTK_THORNSTRING, "Irr. mass BH [M_sun]:  %g", BH_irr_mass);
//    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass NS [M_sun]:   %g", NS_ADM_mass);
//    CCTK_VInfo (CCTK_THORNSTRING, "ADM mass BHNS [M_sun]: %g", BHNS_ADM_mass);
    CCTK_REAL BHNS_Omega            = 0.003678125155199;
    // TODO: reading EOS from Elliptica
    double K = 100.0; // make sure ths is in polytropic units
    double Gamma = 2.0; // make sure ths is in polytropic units
    CCTK_VInfo (CCTK_THORNSTRING, "Here22");
    idr->ifields = "alpha,betax,betay,betaz,adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz,adm_Kxx,adm_Kxy,adm_Kxz,adm_Kyy,adm_Kyz,adm_Kzz,grhd_rho,grhd_epsl,grhd_vx,grhd_vy,grhd_vz";
    idr->npoints = N_points;
    idr->x_coords = xx;
    idr->y_coords = yy;
    idr->z_coords = zz;
    idr->set_param("BHNS_filler_method",BH_filler_method,idr);
    idr->set_param("ADM_B1I_form","zero",idr);
    
    
    elliptica_id_reader_interpolate(idr);

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
      //if using a realistic EOS, set the beta equilibrium conditions
      if(init_real){
      if(rho[i]>=1e-7){
          double yeres = interploc(log_rho0_tab_betaloc, Y_e_tabloc, n_tab_betaloc,log10(rho[i]), &n_nearest_betaloc);
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
    // These initial data assume a helical Killing vector field

    if (CCTK_EQUALS(initial_lapse, "Elliptica_BHNS")) { 
      if (CCTK_EQUALS (initial_dtlapse, "Elliptica_BHNS")) {
        CCTK_INFO ("Calculating time derivatives of lapse");
        set_dt_from_domega (CCTK_PASS_CTOC, alp, dtalp, BHNS_Omega);
      } else if (CCTK_EQUALS (initial_dtlapse, "none") or CCTK_EQUALS(initial_dtlapse,"zero")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error setting dtlapse");
      }
    }

    if (CCTK_EQUALS(initial_shift, "Elliptica_BHNS")) { 
      if (CCTK_EQUALS (initial_dtshift, "Elliptica_BHNS")) {
        CCTK_INFO ("Calculating time derivatives of shift");
        set_dt_from_domega (CCTK_PASS_CTOC, betax, dtbetax, BHNS_Omega);
        set_dt_from_domega (CCTK_PASS_CTOC, betay, dtbetay, BHNS_Omega);
        set_dt_from_domega (CCTK_PASS_CTOC, betaz, dtbetaz, BHNS_Omega);
      } else if (CCTK_EQUALS (initial_dtshift, "none") or CCTK_EQUALS(initial_dtshift,"zero")) {
        // do nothing
      } else {
        CCTK_WARN (CCTK_WARN_ABORT, "internal error setting dtshift");
      }
    }
  }

  // free  
  elliptica_id_reader_free(idr);
  if (xx) free(xx);
  if (yy) free(yy);
  if (zz) free(zz);
  CCTK_INFO ("Done.");
  } catch (ios::failure e) {
    CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not read initial data from file '%s': %s", Elliptica_bhns_file, e.what());
  }

  if(verbose){
    CCTK_INFO("Exiting Elliptica_BHNS_initialize");
  }

}

