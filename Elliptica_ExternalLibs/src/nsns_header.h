#ifndef nsns_header_LIB_H
#define nsns_header_LIB_H

#include "core_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "utilities_lib.h"
#include "checkpoint_lib.h"
#include "fields_lib.h"
#include "maths_spectral_methods_lib.h"
#include "maths_calculus_lib.h"
#include "physics_lib.h"
#include "physics_EoS_lib.h"
#include "physics_stress_energy_lib.h"
#include "physics_transformation_lib.h"
#include "physics_observe_lib.h"
#include "physics_system_lib.h"
#include "physics_star_lib.h"

/* prefix internal parameters of this project, 
// PLEASE keep it capitalized.
// NOTE: this prefix must be consistent with physics->ssys
// which initialized in init_physics function, and the reason is
// that some parameters when are read in various parts of physics 
// libraries assume this consistency otherwise it might get 
// error of undefined parameter. */
#define P_ "NSNS_"

Physics_T *nsns_initialize_new_physics(Physics_T *const phys);
void nsns_add_fields(Physics_T *const phys,const char *const region);
void nsns_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen);


Physics_T *nsns_read_physics_from_checkpoint(void);
void nsns_set_bam_fields(Grid_T *const grid);
void nsns_analyze(Physics_T *const phys,const int iteration);
void nsns_solve_equation(Physics_T *const phys);

#endif

