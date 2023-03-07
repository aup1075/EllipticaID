/* (c) 2023 Pedro Espino */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <stdio.h>
#include "elliptica_id_reader_lib.h"
//#include <idr_main.h>
//#include <idr_main.c>

extern "C"
void Elliptica_BHNS_initialize(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    const char *f="/home/pespino/temp_project/projects/Elliptica_TEST/BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_00/BHNS_BH_m7_s0.0_flat--NS_m1.6_O0.0_SLy--d40_full_valgrind_3x3x3_00/checkpoint.dat";
    const char *option = "generic";
    
    Elliptica_ID_Reader_T *idr = elliptica_id_reader_init(f,option);
//    idr->ifields = "alpha,adm_gxx,adm_gxy";
//    idr->npoints = 2;
//    double x[2] = {1,2};
//    idr->x_coords = x;
//    idr->y_coords = x;
//    idr->z_coords = x;
//    idr->param("BHNS_filler_method","ChebTn_Ylm_perfect_s2",idr);
//    idr->param("ADM_B1I_form","zero",idr);
//    
//    elliptica_id_reader_interpolate(idr);
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

