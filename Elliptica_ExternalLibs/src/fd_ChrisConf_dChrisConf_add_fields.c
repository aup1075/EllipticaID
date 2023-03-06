/*
  These C codes generated by Cpi version 2.0
  Copyright (C) 2019-2020 Alireza Rashti.
*/


#include "fd_header.h"


void fd_add_fields_ChrisConf_dChrisConf(Grid_T *const grid);
void fd_add_fields_ChrisConf_dChrisConf(Grid_T *const grid)
{
 Uint p;
 FOR_ALL_PATCHES(p,grid)
 {
 Patch_T *patch = grid->patch[p];


  /* declaring: */
  ADD_AND_ALLOC_FIELD(ChrisConf_U1D0D0)
  ADD_AND_ALLOC_FIELD(ChrisConf_U0D2D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U1D0D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U2D2D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U2D1D1)
  ADD_AND_ALLOC_FIELD(ChrisConf_U2D0D0)
  ADD_AND_ALLOC_FIELD(ChrisConf_U1D1D1)
  ADD_AND_ALLOC_FIELD(ChrisConf_U0D1D1)
  ADD_AND_ALLOC_FIELD(ChrisConf_U0D1D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U1D1D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U0D0D1)
  ADD_AND_ALLOC_FIELD(ChrisConf_U0D0D0)
  ADD_AND_ALLOC_FIELD(ChrisConf_U2D0D1)
  ADD_AND_ALLOC_FIELD(ChrisConf_U0D0D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U2D1D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U1D2D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U2D0D2)
  ADD_AND_ALLOC_FIELD(ChrisConf_U1D0D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D1D1D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D2D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D1D1D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D1D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D1D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D1D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D2D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D2D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D1D1D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D2D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D2D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D1D1D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D1D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D1D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D1D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D1D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D1D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D1D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D2D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D2D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D2D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D0D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D1D1D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D1D1D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D1D1D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D2D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D1D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D0D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D0D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D0D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D0D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D1D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D1D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D0D1D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D1D1D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D0D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D0D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D0D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D0D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D1D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D1D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U1D1D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D1D1D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D2D0)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D2D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U0D0D2D2)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D1D1)
  ADD_AND_ALLOC_FIELD(dChrisConf_U2D0D1D0)


 }
}
