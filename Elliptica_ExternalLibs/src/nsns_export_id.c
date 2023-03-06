/*
// Alireza Rashti
// March 2023
*/

/* exporting initial data for evolution codes */

#include "nsns_export_id.h"


/* export ID for a general evolution code */
void nsns_export_id_generic(Elliptica_ID_Reader_T *const idr)
{
  FUNC_TIC
  
  Physics_T *nsns = 0;
  ID_Reader_T *points = idr_init();
  double CM[3] = {0.};
  char fields_name[STR_LEN_MAX] = {'\0'};// elliptica field names
  char **sfield = 0;
  Uint f;
  
  /* adapt fields_notations for Elliptica and ensure given fields are expected.
  // NOTE: the following replaces the name with the same ORDER as the input. */
  assert(sprintf(fields_name,"%s",idr->ifields));
  
  /* metric fields */
  regex_replace(fields_name,"\\balpha\\b",EVO_"alpha",fields_name);
  
  regex_replace(fields_name,"\\bbetax\\b",EVO_"beta_U0",fields_name);
  regex_replace(fields_name,"\\bbetay\\b",EVO_"beta_U1",fields_name);
  regex_replace(fields_name,"\\bbetaz\\b",EVO_"beta_U2",fields_name);
  
  regex_replace(fields_name,"\\badm_gxx\\b",EVO_"adm_g_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_gxy\\b",EVO_"adm_g_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_gxz\\b",EVO_"adm_g_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_gyy\\b",EVO_"adm_g_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_gyz\\b",EVO_"adm_g_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_gzz\\b",EVO_"adm_g_D2D2",fields_name);
  
  regex_replace(fields_name,"\\badm_Kxx\\b",EVO_"adm_Kij_D0D0",fields_name);
  regex_replace(fields_name,"\\badm_Kxy\\b",EVO_"adm_Kij_D0D1",fields_name);
  regex_replace(fields_name,"\\badm_Kxz\\b",EVO_"adm_Kij_D0D2",fields_name);
  regex_replace(fields_name,"\\badm_Kyy\\b",EVO_"adm_Kij_D1D1",fields_name);
  regex_replace(fields_name,"\\badm_Kyz\\b",EVO_"adm_Kij_D1D2",fields_name);
  regex_replace(fields_name,"\\badm_Kzz\\b",EVO_"adm_Kij_D2D2",fields_name);
  
  /* matter fields */
  regex_replace(fields_name,"\\bgrhd_vx\\b",EVO_"grhd_v_U0",fields_name);
  regex_replace(fields_name,"\\bgrhd_vy\\b",EVO_"grhd_v_U1",fields_name);
  regex_replace(fields_name,"\\bgrhd_vz\\b",EVO_"grhd_v_U2",fields_name);
  regex_replace(fields_name,"\\bgrhd_rho\\b",EVO_"grhd_rho",fields_name);
  regex_replace(fields_name,"\\bgrhd_p\\b",EVO_"grhd_p",fields_name);
  regex_replace(fields_name,"\\bgrhd_epsl\\b",EVO_"grhd_epsl",fields_name);
  
  /* check if all fields are expected */
  sfield = read_separated_items_in_string(fields_name,',');
  f = 0;
  while(sfield[f])
  {
    if (!regex_search("^"EVO_,sfield[f]))
    {
      printf("field '%s' is Unexpected!\n",sfield[f]);
      Error1("Unexpected field! Add me!");
    }
    f++;
  }
  free_2d(sfield);
  
  /* set all parameters set in evo code. */
  // NOTE: adding both plain and modify version to make sure things kick in 
  // after loading from checkpoint file.
  for (Uint n = 0; n < idr->nparams; ++n)
  {
    char lv[STR_LEN_MAX];
    sprintf(lv,CHECKPOINT_SET_PARAM_"%s",idr->params_lv[n]);
    
    Psets(idr->params_lv[n],idr->params_rv[n]);
    Psets(lv,idr->params_rv[n]);
  }
  
  /* don't stop for checkpoint */
  Pseti(CHECKPOINT_SET_PARAM_ P_"STOP",0);
  /* read physics from checkpoint */
  Psets("checkpoint_file_path",idr->checkpoint_path);
  nsns = nsns_read_physics_from_checkpoint();
  points->grid = nsns->grid;
  
  /* get (x,y,z) points from evo. NOTE: no allocation done for (x,y,z) */
  CM[0] = Pgetd(P_"x_CM");
  CM[1] = Pgetd(P_"y_CM");
  CM[2] = Pgetd(P_"z_CM");
  idr_find_XYZ_from_xyz(idr,points,CM);
  
  /* go from Omega x r to inertial coords sys asymptotically.
  // EVO needs this!? (this param should be set in the reader) */
  // Psets(CHECKPOINT_SET_PARAM_ "ADM_B1I_form","zero");
  physics(nsns,ADM_UPDATE_Kij);/* need this for nsns_set_evo_fields_generic */
  
  /* set bam fields based on initial data to be usable for evo */
  nsns_set_evo_fields_generic(nsns->grid);
  
  /* write into array */
  idr_interpolate_fields_and_save_in_array(idr,points,fields_name,idr->ifields);
  
  /* finishing up */
  // since no alocation done for (x,y,z):
  points->x = 0;
  points->y = 0;
  points->z = 0;
  idr_free(points);
  free_physics(nsns);
  
  FUNC_TOC  
}

