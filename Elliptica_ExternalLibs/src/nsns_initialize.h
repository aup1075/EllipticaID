#include "nsns_header.h"
#include "maths_equation_solvings_lib.h"

/* index of ns1 */
#define Ins1 (0)
/* index of ns2 */
#define Ins2 (1)

/* encapsulate these lines for multiple uses, N arg value can take 1 and 2 */
#define USE_LAST_NS_SURFACE(N) \
  printf(Pretty0"Using the last NS"#N" surface.\n");\
  lmax = (Uint)Pgeti("NS"#N"_surface_R|lmax");\
  n    = Ncoeffs_Ylm(lmax);\
  double *realClm = alloc_ClmYlm(lmax);/* freed in free_grid_char */\
  double *imagClm = alloc_ClmYlm(lmax);/* freed in free_grid_char */\
  double *coeffs  = 0;\
  coeffs = Pgetdd("NS"#N"_surface_R|realClm");\
  for (Uint ij = 0; ij < n; ++ij)\
    realClm[ij] = coeffs[ij];\
  coeffs = Pgetdd("NS"#N"_surface_R|imagClm");\
  for (Uint ij = 0; ij < n; ++ij)\
    imagClm[ij] = coeffs[ij];\
  /* might already have values so free them. */\
  Free(grid_char->params[Ins##N]->relClm);\
  Free(grid_char->params[Ins##N]->imgClm);\
  grid_char->params[Ins##N]->relClm = realClm;\
  grid_char->params[Ins##N]->imgClm = imagClm;\
  grid_char->params[Ins##N]->lmax   = lmax;\
  grid_char->params[Ins##N]->r_min  = Pgetd("NS"#N"_min_radius");\
  grid_char->params[Ins##N]->r_max  = Pgetd("NS"#N"_max_radius");\
  Pseti("NS"#N"_did_NS_surface_change?",0);


Physics_T *nsns_initialize_new_physics(Physics_T *const old_phys);
static Physics_T *guess_new_physics(void);

static void 
  create_new_grid(Grid_Char_T *const grid_char,Physics_T *const nsns);


static void update_partial_derivatives(Physics_T *const phys,
                                       const char *const region,
                                       const char *const regex);

static void initial_B0I(Physics_T *const phys,
                       const char *const region);


Physics_T *nsns_read_physics_from_checkpoint(void);
static Physics_T *infer_new_physics(Physics_T *const old_phys);
static void initialize_fields_using_previous_solve
            (Physics_T *const new_phys,Physics_T *const old_phys);
static void update_params(Physics_T *const phys);
static void move_jacobian
            (Physics_T *const new_phys,Physics_T *const old_phys);
