/*
// Alireza Rashti
// January 2023
*/

/* various functions needed to make a physics object for NS-NS */


#include "nsns_initialize.h"

/* decide how to initialize the new physics */
Physics_T *nsns_initialize_new_physics(Physics_T *const old_phys)
{
  Physics_T *new_phys = 0;
  
  /* if already hit the stop */
  if (Pgeti(P_"STOP") == 1)
  {
    /* update things */
    new_phys = infer_new_physics(old_phys);
    /* save and dump */
    write_checkpoint(new_phys,Pgets(P_"my_directory"));
    free_physics(new_phys);
    
    printf(Pretty0"I'm done!  :)\n");
    return 0;
  }
  
  if (!old_phys)/* if empty, come up with a start off */
  {
    /* if we wanna use a particular checkpoint file, 
    // mostly for debug purposes, set:
    // set:
    // NSNS_start_off       = checkpoint_file
    // checkpoint_file_path = path_to_checkpoint_file */
    if (Pcmps(P_"start_off","checkpoint_file"))
    {
      /* modify output directories.
      // AD-HOC: put everything in top directory. */
      Psets(CHECKPOINT_SET_PARAM_ "top_directory", 
            Pgets("top_directory"));
      Psets(CHECKPOINT_SET_PARAM_ P_"my_directory", 
            Pgets("top_directory"));
      Psets(CHECKPOINT_SET_PARAM_ P_"Diagnostics", 
            Pgets("top_directory"));
      
      new_phys = nsns_read_physics_from_checkpoint();
    }
    
    /* can we resume from a useful checkpoint file */
    else if (can_we_use_checkpoint(Pgets("top_directory")))
      new_phys = nsns_read_physics_from_checkpoint();
      
    else 
      new_phys = guess_new_physics();
  }
  else/* use old physics, tune it and make new physics */
  {
    new_phys = infer_new_physics(old_phys);
  }
  
  return new_phys;
}


/* use old physics to infer the new physics */
static Physics_T *infer_new_physics(Physics_T *const old_nsns)
{
  if (!old_nsns) return 0;
  
  FUNC_TIC
  
  Physics_T *const nsns    = init_physics(0,NSNS);/* the whole system */
  Physics_T *const ns1     = init_physics(nsns,NS1);/* NS1 part */
  Physics_T *const ns2     = init_physics(nsns,NS2);/* NS2 part */
  Physics_T *const old_ns1 = init_physics(old_nsns,NS1);/* NS1 part */
  Physics_T *const old_ns2 = init_physics(old_nsns,NS2);/* NS2 part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  old_ns1->grid_char = grid_char;
  old_ns1->igc       = Ins1;
  old_ns2->grid_char = grid_char;
  old_ns2->igc       = Ins2;
  
  /* update, adjust and tune */
  Psets("NS1_enthalpy_neat","no");
  //physics(old_ns1,STRESS_ENERGY_UPDATE);
  physics(old_ns1,STAR_TUNE_EULER_CONST);
  //physics(old_ns1,STRESS_ENERGY_UPDATE);
  
  Psets("NS2_enthalpy_neat","no");
  //physics(old_ns2,STRESS_ENERGY_UPDATE);
  physics(old_ns2,STAR_TUNE_EULER_CONST);
  //physics(old_ns2,STRESS_ENERGY_UPDATE);
  
  physics(old_nsns,SYS_TUNE_P_ADM);
  
  physics(old_ns1,STRESS_ENERGY_UPDATE);
  physics(old_ns1,STAR_TUNE_FORCE_BALANCE);
  physics(old_ns1,STAR_EXTRAPOLATE_MATTERS);
  physics(old_ns1,STAR_TUNE_CENTER);
  physics(old_ns1,STAR_FIND_SURFACE);
  
  physics(old_ns2,STRESS_ENERGY_UPDATE);
  physics(old_ns2,STAR_TUNE_FORCE_BALANCE);
  physics(old_ns2,STAR_EXTRAPOLATE_MATTERS);
  physics(old_ns2,STAR_TUNE_CENTER);
  physics(old_ns2,STAR_FIND_SURFACE);
  
  /* new grid */
  create_new_grid(grid_char,nsns);
  ns1->grid = nsns->grid;
  ns2->grid = nsns->grid;
  
  /* set and update parameters */
  update_params(nsns);
  physics(nsns,FREE_DATA_SET_PARAMS);
  physics(nsns,ADM_SET_PARAMS);
  physics(nsns,SYS_SET_PARAMS);
  physics(nsns,STRESS_ENERGY_SET_PARAMS);
  physics(nsns,OBSERVE_SET_PARAMS);  
  physics(nsns,STAR_SET_PARAMS);
  
  /* add fields */
  physics(nsns,ADM_ADD_FIELDS);
  physics(nsns,FREE_DATA_ADD_FIELDS);
  physics(nsns,STRESS_ENERGY_ADD_FIELDS);
  physics(nsns,SYS_ADD_FIELDS);
  physics(nsns,OBSERVE_ADD_FIELDS);
  physics(nsns,STAR_ADD_FIELDS);
  
  /* populate fields */
  physics(nsns,FREE_DATA_POPULATE);
  initialize_fields_using_previous_solve(nsns,old_nsns);
  
  /* move Jacobian if possible */
  move_jacobian(nsns,old_nsns);
  
  /* beta = B0+B1 */
  physics(nsns,ADM_UPDATE_B1I);
  update_partial_derivatives(nsns,".*","^dB0_U.+,^ddB0_U.+");
  physics(nsns,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(nsns,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  update_partial_derivatives(ns1,"NS1","^dphi_D.$,^ddphi_D.D.$");
  update_partial_derivatives(ns2,"NS2","^dphi_D.$,^ddphi_D.D.$");
  
  /* update AConf^{ij} */
  physics(nsns,ADM_UPDATE_AConfIJ);
  
  /* update matter fields */
  Psets("NS1_enthalpy_neat","yes");
  physics(ns1,STRESS_ENERGY_UPDATE);
  
  Psets("NS2_enthalpy_neat","yes");
  physics(ns2,STRESS_ENERGY_UPDATE);
  
  /* free */
  free_physics(ns1);
  free_physics(ns2);
  free_physics(old_ns1);
  free_physics(old_ns2);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return nsns;
}  

/* use a known NS1 and NS2 solution to initialize the physics */
static Physics_T *guess_new_physics(void)
{
  FUNC_TIC
  
  Physics_T *const nsns = init_physics(0,NSNS);/* the whole system */
  Physics_T *const ns1  = init_physics(nsns,NS1);/* NS part */
  Physics_T *const ns2  = init_physics(nsns,NS2);/* NS part */
  Grid_Char_T *const grid_char = init_grid_char(0);
  
  /* set parameters */
  physics(nsns,FREE_DATA_SET_PARAMS);
  physics(nsns,ADM_SET_PARAMS);
  physics(nsns,SYS_SET_PARAMS);
  physics(nsns,STRESS_ENERGY_SET_PARAMS);
  physics(nsns,OBSERVE_SET_PARAMS);  
  physics(nsns,STAR_SET_PARAMS);
  
  /* create grid */
  ns1->grid_char = grid_char;
  ns1->igc       = Ins1;
  ns2->grid_char = grid_char;
  ns2->igc       = Ins2;

  physics(ns1,STAR_START);
  physics(ns1,STAR_FIND_SURFACE);
  physics(ns2,STAR_START);
  physics(ns2,STAR_FIND_SURFACE);
  
  create_new_grid(grid_char,nsns);
  ns1->grid = nsns->grid;
  ns2->grid = nsns->grid;
  
  /* add fields */
  physics(nsns,ADM_ADD_FIELDS);
  physics(nsns,FREE_DATA_ADD_FIELDS);
  physics(nsns,STRESS_ENERGY_ADD_FIELDS);
  physics(nsns,SYS_ADD_FIELDS);
  physics(nsns,OBSERVE_ADD_FIELDS);
  physics(nsns,STAR_ADD_FIELDS);
  
  /* populate fields */
  physics(nsns,FREE_DATA_POPULATE);
  physics(nsns,SYS_INITIALIZE_FIELDS);
  /* beta = B0+B1 */
  physics(nsns,ADM_UPDATE_B1I);
  initial_B0I(nsns,".*");
  update_partial_derivatives(nsns,".*","^dB0_U.+,^ddB0_U.+");
  physics(nsns,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(nsns,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  update_partial_derivatives(ns1,"NS1","^dphi_D.$,^ddphi_D.D.$");
  update_partial_derivatives(ns2,"NS2","^dphi_D.$,^ddphi_D.D.$");
  
  /* update AConf^{ij} */
  physics(nsns,ADM_UPDATE_AConfIJ);
  
  /* update stress energy-tensor */
  Psetd("NS1_Euler_equation_constant",
        star_NS_current_Euler_eq_const(ns1));
  Psetd("NS2_Euler_equation_constant",
        star_NS_current_Euler_eq_const(ns2));
        
  Psets("NS1_enthalpy_neat","yes");
  physics(ns1,STRESS_ENERGY_UPDATE);

  Psets("NS2_enthalpy_neat","yes");
  physics(ns2,STRESS_ENERGY_UPDATE);
  
  /* free */
  free_physics(ns1);
  free_physics(ns2);
  free_grid_char(grid_char);
  
  FUNC_TOC
  return nsns;
}

/* based on grid character, make a new grid for the system. */
static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const nsns)
{
  FUNC_TIC
  
  /* a new grid */
  Grid_T *const grid = alloc_grid();
  const double ns1_box_len_ratio = 0.2;/* experimentally */
  const double ns2_box_len_ratio = 0.2;/* experimentally */
  int update_ns1_surface = 1;
  int update_ns2_surface = 1;
  Uint lmax,n;
  
  /* grid for characters */
  grid_char->grid = grid;
  
  if (!Pcmps("grid_kind","SplitCubedSpherical(NS+NS)"))
    Error0(NO_OPTION);
  
  /* set "grid_NS1_central_box_length" and "grid_NS2_central_box_length"
  // automatically, if it is asked. 
  // NOTE: these params are set only for the very first time 
  // and this is important for stability of NS. additionally,
  // one must note that if mass of the object is being iterated, 
  // this auto option is not very ideal, since the final radius is not
  // known yet */
  if (Pcmps("grid_NS1_central_box_length","auto"))
    Psetd("grid_NS1_central_box_length",
          ns1_box_len_ratio*grid_char->params[Ins1]->r_min);
          
  if (Pcmps("grid_NS2_central_box_length","auto"))
    Psetd("grid_NS2_central_box_length",
          ns2_box_len_ratio*grid_char->params[Ins2]->r_min);
  
  /* separation */
  grid_char->S              = Pgetd("NSNS_separation");
  /* NS1 */
  grid_char->params[Ins1]->l = Pgetd("grid_NS1_central_box_length");
  grid_char->params[Ins1]->w = Pgetd("grid_NS1_central_box_length");
  grid_char->params[Ins1]->h = Pgetd("grid_NS1_central_box_length");
  /* NS2 */
  grid_char->params[Ins2]->l = Pgetd("grid_NS2_central_box_length");
  grid_char->params[Ins2]->w = Pgetd("grid_NS2_central_box_length");
  grid_char->params[Ins2]->h = Pgetd("grid_NS2_central_box_length");
    
  /* save the values for a rainy day */
  if (Pgeti("NS1_did_NS_surface_finder_work?"))
  {
    double rel_change = 0.;
    /* change the relative difference using coeffs */
    if (
        /* if prev value exists */
        PgetddEZ("NS1_surface_R|realClm")  && 
        PgetddEZ("NS1_surface_R|imagClm")  &&
        /* if the old and new have the same lmax */
        PgetiEZ("NS1_surface_R|lmax") == (int)grid_char->params[Ins1]->lmax
       )
    {
      lmax = (Uint)Pgeti("NS1_surface_R|lmax");
      n    = Ncoeffs_Ylm(lmax);
      const double *realClm = Pgetdd("NS1_surface_R|realClm");
      const double *imagClm = Pgetdd("NS1_surface_R|imagClm");
      /* diff between old and new */
      double dreal = L2_norm(n,realClm,grid_char->params[Ins1]->relClm);
      double dimag = L2_norm(n,imagClm,grid_char->params[Ins1]->imgClm);
      /* relative change df/f */
      rel_change = (dreal+dimag) /
                   (L2_norm(n,realClm,0)+L2_norm(n,imagClm,0));
    
      /* update if change greater than prescribed */
      if (rel_change > Pgetd("NS1_surface_change_threshold"))
        update_ns1_surface = 1;
      else
        update_ns1_surface = 0;
    }
    /* save new values if ns1 surface must change */
    if (update_ns1_surface)
    {
      n = Ncoeffs_Ylm(grid_char->params[Ins1]->lmax);
      update_parameter_array("NS1_surface_R|realClm",
                             grid_char->params[Ins1]->relClm,n);
      update_parameter_array("NS1_surface_R|imagClm",
                             grid_char->params[Ins1]->imgClm,n);
      Pseti("NS1_surface_R|lmax",(int)grid_char->params[Ins1]->lmax);
      Pseti("NS1_did_NS_surface_change?",1);
    }
    else
    {
      printf(Pretty0"relative change is smaller "
                    "than threshold (%.3e < %.3e).\n",
                    rel_change,Pgetd("NS1_surface_change_threshold"));
      USE_LAST_NS_SURFACE(1);
    }
  }
  /* since surface finder failed use previous value. */
  else
  {
    printf(Pretty0"NS1 surface finder failed.\n");
    USE_LAST_NS_SURFACE(1);
  }
  
  /* save the values for a rainy day */
  if (Pgeti("NS2_did_NS_surface_finder_work?"))
  {
    double rel_change = 0.;
    /* change the relative difference using coeffs */
    if (
        /* if prev value exists */
        PgetddEZ("NS2_surface_R|realClm")  && 
        PgetddEZ("NS2_surface_R|imagClm")  &&
        /* if the old and new have the same lmax */
        PgetiEZ("NS2_surface_R|lmax") == (int)grid_char->params[Ins2]->lmax
       )
    {
      lmax = (Uint)Pgeti("NS2_surface_R|lmax");
      n    = Ncoeffs_Ylm(lmax);
      const double *realClm = Pgetdd("NS2_surface_R|realClm");
      const double *imagClm = Pgetdd("NS2_surface_R|imagClm");
      /* diff between old and new */
      double dreal = L2_norm(n,realClm,grid_char->params[Ins2]->relClm);
      double dimag = L2_norm(n,imagClm,grid_char->params[Ins2]->imgClm);
      /* relative change df/f */
      rel_change = (dreal+dimag) /
                   (L2_norm(n,realClm,0)+L2_norm(n,imagClm,0));
    
      /* update if change greater than prescribed */
      if (rel_change > Pgetd("NS2_surface_change_threshold"))
        update_ns2_surface = 1;
      else
        update_ns2_surface = 0;
    }
    /* save new values if ns2 surface must change */
    if (update_ns2_surface)
    {
      n = Ncoeffs_Ylm(grid_char->params[Ins2]->lmax);
      update_parameter_array("NS2_surface_R|realClm",
                             grid_char->params[Ins2]->relClm,n);
      update_parameter_array("NS2_surface_R|imagClm",
                             grid_char->params[Ins2]->imgClm,n);
      Pseti("NS2_surface_R|lmax",(int)grid_char->params[Ins2]->lmax);
      Pseti("NS2_did_NS_surface_change?",1);
    }
    else
    {
      printf(Pretty0"relative change is smaller "
                    "than threshold (%.3e < %.3e).\n",
                    rel_change,Pgetd("NS2_surface_change_threshold"));
      USE_LAST_NS_SURFACE(2);
    }
  }
  /* since surface finder failed use previous value. */
  else
  {
    printf(Pretty0"NS2 surface finder failed.\n");
    USE_LAST_NS_SURFACE(2);
  }
  
  /* check central box length */
  if (grid_char->params[Ins1]->l > grid_char->params[Ins1]->r_min/2. ||
      grid_char->params[Ins1]->w > grid_char->params[Ins1]->r_min/2. ||
      grid_char->params[Ins1]->h > grid_char->params[Ins1]->r_min/2.)
    Error0("NS1 central box is too big!");
  
  /* check central box length */
  if (grid_char->params[Ins2]->l > grid_char->params[Ins2]->r_min/2. ||
      grid_char->params[Ins2]->w > grid_char->params[Ins2]->r_min/2. ||
      grid_char->params[Ins2]->h > grid_char->params[Ins2]->r_min/2.)
    Error0("NS2 central box is too big!");
  
  set_params_of_split_cubed_spherical_grid(grid_char);
    
  make_patches(grid);
  realize_interfaces(grid);
  
  nsns->grid = grid;
  
  FUNC_TOC
}

/* update partial derivatives of the given field name regex match 
// at the given region */
static void update_partial_derivatives(Physics_T *const phys,
                                       const char *const region,
                                       const char *const regex)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  printf(Pretty0"%s\n",regex);
  fflush(stdout);
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    partial_derivative_regex(patch,regex);
  }
  
  FUNC_TOC
}

/* initial B0^i in beta = B0+B1 */
static void initial_B0I(Physics_T *const phys,
                       const char *const region)
{
  FUNC_TIC
  
  Grid_T *const grid = mygrid(phys,region);
  Uint p;
  
  OpenMP_Patch_Pragma(omp parallel for)
  for (p = 0; p < grid->np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    
    READ_v(beta_U0);
    READ_v(beta_U1);
    READ_v(beta_U2);
    READ_v(B1_U0);
    READ_v(B1_U1);
    READ_v(B1_U2);
    WRITE_v(psi);
    
    REALLOC_v_WRITE_v(B0_U0);
    REALLOC_v_WRITE_v(B0_U1);
    REALLOC_v_WRITE_v(B0_U2);
    
    FOR_ALL_ijk
    {
      double psim4 = pow(psi[ijk],-4.);
      
      B0_U0[ijk] = psim4*beta_U0[ijk]-B1_U0[ijk];
      B0_U1[ijk] = psim4*beta_U1[ijk]-B1_U1[ijk];
      B0_U2[ijk] = psim4*beta_U2[ijk]-B1_U2[ijk];
    }
    
    /* attenuate */
    if (IsItCovering(patch,"outermost"))
    {
      FOR_ALL_ijk
      {
        DEF_RELATIVE_x
        DEF_RELATIVE_y
        DEF_RELATIVE_z
        DEF_RELATIVE_r
        
        psi[ijk]   /= r;
        B0_U0[ijk] /= r;
        B0_U1[ijk] /= r;
        B0_U2[ijk] /= r;
      }
    }
  }
  
  FUNC_TOC
}

/* loading from checkpoint */
Physics_T *nsns_read_physics_from_checkpoint(void)
{
  FUNC_TIC
  Physics_T *const nsns = init_physics(0,NSNS);
  FILE *file = 0;
  
  /* first load grid and parameters */
  file = open_checkpoint_file_then_read_grid_and_params(nsns);
  
  /* it already hit the stop */
  if (Pgeti(P_"STOP") == 1)
  {
    printf(Pretty0" All iterations have already done.\n");
    free_physics(nsns);
    Fclose(file);
    
    FUNC_TOC
    return 0;
  }
  
  Physics_T *const ns1 = init_physics(nsns,NS1);
  Physics_T *const ns2 = init_physics(nsns,NS2);
  
  /* make the patches */
  make_patches(nsns->grid);
  
  /* realizing the geometry */
  realize_interfaces(nsns->grid);
  
  /* set parameters, it's important to add paramters 
  // since these call also reposible to set default functions. */
  physics(nsns,FREE_DATA_SET_PARAMS);
  physics(nsns,ADM_SET_PARAMS);
  physics(nsns,SYS_SET_PARAMS);
  physics(nsns,STRESS_ENERGY_SET_PARAMS);
  physics(nsns,OBSERVE_SET_PARAMS);  
  physics(nsns,STAR_SET_PARAMS);
  
  /* now add fields */
  physics(nsns,ADM_ADD_FIELDS);
  physics(nsns,FREE_DATA_ADD_FIELDS);
  physics(nsns,STRESS_ENERGY_ADD_FIELDS);
  physics(nsns,SYS_ADD_FIELDS);
  physics(nsns,OBSERVE_ADD_FIELDS);
  physics(nsns,STAR_ADD_FIELDS);
  
  /* populate free data fields */
  physics(nsns,FREE_DATA_POPULATE);
  
  /* then read saved fields in param "checkpoint_save" */
  read_fields_from_checkpoint_file(nsns,file);
  Fclose(file);
  
  /* alse we need NS spin vector */
  star_W_spin_vector_idealfluid_update(ns1,"NS1");
  star_W_spin_vector_idealfluid_update(ns2,"NS2");
  
  /* beta = B0+B1 */
  physics(nsns,ADM_UPDATE_B1I);
  update_partial_derivatives(nsns,".*","^dB0_U.+,^ddB0_U.+");
  physics(nsns,ADM_UPDATE_beta);
  
  /* update derivatives */
  update_partial_derivatives(nsns,".*","^dpsi_D.$,^ddpsi_D.D.$,"
                                      "^dalphaPsi_D.$,^ddalphaPsi_D.D.$");
  update_partial_derivatives(ns1,"NS1","^dphi_D.$,^ddphi_D.D.$");
  update_partial_derivatives(ns2,"NS2","^dphi_D.$,^ddphi_D.D.$");
  
  /* update AConf^{ij} */
  physics(nsns,ADM_UPDATE_AConfIJ);
  
  /* update matter fields */
  Psets("NS1_enthalpy_neat","yes");
  physics(ns1,STRESS_ENERGY_UPDATE);

  Psets("NS2_enthalpy_neat","yes");
  physics(ns2,STRESS_ENERGY_UPDATE);
  
  free_physics(ns1);
  free_physics(ns2);

  FUNC_TOC
  return nsns;
}

/* using copy or interpolation from old physics to 
// initialize fields for new physics */
static void initialize_fields_using_previous_solve
            (Physics_T *const new_phys,Physics_T *const old_phys)
{
  FUNC_TIC
  
  Physics_T *const old_ns1 = init_physics(old_phys,NS1);
  Physics_T *const new_ns1 = init_physics(new_phys,NS1);
  Physics_T *const old_ns2 = init_physics(old_phys,NS2);
  Physics_T *const new_ns2 = init_physics(new_phys,NS2);
  
  /* matter fields */
  interpolate_fields_from_old_grid_to_new_grid
    (mygrid(old_ns1,"NS1,NS1_around_IB"),mygrid(new_ns1,"NS1"),"phi,enthalpy",0);
  
  interpolate_fields_from_old_grid_to_new_grid
    (mygrid(old_ns2,"NS2,NS2_around_IB"),mygrid(new_ns2,"NS2"),"phi,enthalpy",0);
  
  /* if resolution changed */
  if(Pgeti(P_"did_resolution_change?"))
  {
    interpolate_fields_from_old_grid_to_new_grid
      (old_phys->grid,new_phys->grid,"psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
  }
  else
  {
    const char *region1 = 0;
    const char *region2 = 0;
    if (new_phys->grid->kind == Grid_SplitCubedSpherical_NSNS)
    {
      /* since filling_box,outermost are fixed, only copy */
      region1 = "filling_box,outermost";
      region2 = "filling_box,outermost";
      interpolate_fields_from_old_grid_to_new_grid
        (mygrid(old_phys,region1),mygrid(new_phys,region2),
         "psi,alphaPsi,B0_U0,B0_U1,B0_U2",1);
      
      region1 = "NS1,NS1_around";
      region2 = "NS1,NS1_around";
      interpolate_fields_from_old_grid_to_new_grid
        (mygrid(old_ns1,region1),mygrid(new_ns1,region2),
         "psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
      
      region1 = "NS2,NS2_around";
      region2 = "NS2,NS2_around";
      interpolate_fields_from_old_grid_to_new_grid
        (mygrid(old_ns2,region1),mygrid(new_ns2,region2),
         "psi,alphaPsi,B0_U0,B0_U1,B0_U2",0);
         
    }
    else
      Error0(NO_OPTION);
  }
  
  /* alse we need NS spin vector */
  star_W_spin_vector_idealfluid_update(new_ns1,"NS1");
  star_W_spin_vector_idealfluid_update(new_ns2,"NS2");
  
  free_physics(old_ns1);
  free_physics(new_ns1);
  free_physics(old_ns2);
  free_physics(new_ns2);
  
  FUNC_TOC
}

/* update some parameters for a new physics */
static void update_params(Physics_T *const phys)
{
  FUNC_TIC  
  UNUSED(phys);
  FUNC_TOC
}

/* move Jacobian if possible */
static void move_jacobian
            (Physics_T *const new_phys,Physics_T *const old_phys)
{
  FUNC_TIC
  
  if(Pgeti(P_"did_resolution_change?") || 
     !new_phys->grid                   || 
     !old_phys->grid)
  {
    FUNC_TOC
    return;
  }
  
  Physics_T *const old_ns1 = init_physics(old_phys,NS1);
  Physics_T *const new_ns1 = init_physics(new_phys,NS1);
  Physics_T *const old_ns2 = init_physics(old_phys,NS2);
  Physics_T *const new_ns2 = init_physics(new_phys,NS2);
  Grid_T *gnew = 0;
  Grid_T *gold = 0;
  const char *name1 = 0;
  const char *name2 = 0;
  
  if(new_phys->grid->kind == Grid_SplitCubedSpherical_NSNS)
  {
    /* move Jacobian of outermost */
    gnew = mygrid(new_phys,"outermost");
    gold = mygrid(old_phys,"outermost");
    
    FOR_ALL_p(gnew->np)
    {
      name1 = strchr(gnew->patch[p]->name,'_');
      for(Uint p2 = 0; p2 < gold->np; ++p2)
      {
        name2 = strchr(gold->patch[p2]->name,'_');
        if(!strcmp(name1,name2))
        {
          break;
        }
      }
    }
    
    /* move Jacobian of NS1 and NS1 around */
    if (!Pgeti("NS1_did_NS_surface_change?"))
    {
      gnew = mygrid(new_ns1,"NS1,NS1_around");
      gold = mygrid(old_ns1,"NS1,NS1_around");
      
      FOR_ALL_p(gnew->np)
      {
        name1 = strchr(gnew->patch[p]->name,'_');
        for(Uint p2 = 0; p2 < gold->np; ++p2)
        {
          name2 = strchr(gold->patch[p2]->name,'_');
          if(!strcmp(name1,name2))
          {
            break;
          }
        }
      }
    }
    /* move Jacobian of NS2 and NS2 around */
    if (!Pgeti("NS2_did_NS_surface_change?"))
    {
      gnew = mygrid(new_ns2,"NS2,NS2_around");
      gold = mygrid(old_ns2,"NS2,NS2_around");
      
      FOR_ALL_p(gnew->np)
      {
        name1 = strchr(gnew->patch[p]->name,'_');
        for(Uint p2 = 0; p2 < gold->np; ++p2)
        {
          name2 = strchr(gold->patch[p2]->name,'_');
          if(!strcmp(name1,name2))
          {
            break;
          }
        }
      }
    }
    
  }/* if(new_phys->grid->kind == Grid_SplitCubedSpherical_NSNS) */
  else
  {
    Error0(NO_OPTION);
  }
  
  free_physics(old_ns1);
  free_physics(new_ns1);
  free_physics(old_ns2);
  free_physics(new_ns2);

  FUNC_TOC
}

