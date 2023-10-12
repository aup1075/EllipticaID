#!/usr/bin/env python

##
## usage: read Elliptica properties file and create ETK param file
## $ ./me -h
##

## ---------------------------------------------------------------------- ##
from math import *
import sys
import re
import argparse
## ---------------------------------------------------------------------- ##

## global vars
g_dtfac      = 0.25
g_regrid     = 64
g_ns_expand  = 1.2 ## 20% bigger than NS radius
g_data       = ''

def sets(rhs,val,comment=''):
  """
  replace the right hand side with the given string
  ex: cmd = @cmd@ -> sets('cmd','git') ==> cmd = git
  """

  global g_data
  if comment != '':
    g_data = g_data.replace(f'@{rhs}@', f'{val} ## {comment}')
  else:
    g_data = g_data.replace(f'@{rhs}@', f'{val}')

def setd(rhs,val,comment=''):
  """
  replace the right hand side with the given double value
  ex: dt = @dt@ -> setd('dt',0.2) ==> dt = 0.2
  """

  global g_data
  if comment != '':
    g_data = g_data.replace(f'@{rhs}@', f'{val:+0.15e} ## {comment}')
  else:
    g_data = g_data.replace(f'@{rhs}@', f'{val:+0.15e}')

def seti(rhs,val,comment=''):
  """
  replace the right hand side with the given int value
  ex: n = @n@ -> seti('n',128) ==> n = 128
  """
  
  global g_data
  if comment != '':
    g_data = g_data.replace(f'@{rhs}@', f'{int(val)} ## {comment}')
  else:
    g_data = g_data.replace(f'@{rhs}@', f'{int(val)}')
  

def parse_cli():
  """
  arg parser
  """

  p = argparse.ArgumentParser(
      description="read Elliptica properties file and create ETK param file")
  p.add_argument("-d", type=str, required=True,
      help="path/to/dummy/par/file")
  p.add_argument("-p", type=str, required=True,
      help="path/to/Elliptica/NSNS_properties.txt")
  p.add_argument("-i", type=str, required=True,
      help="path/to/Elliptica/checkpoint.dat")
  p.add_argument("-r", type=int, required=True,
      help="resolution, ex: 256, 128, etc.")
  p.add_argument("-ft", type=float, default=2500.,
      help="final time")
  p.add_argument("-mt", type=float, default=1000.,
      help="merger time")
  p.add_argument("-wt", type=int, default=2880,
      help="job walltime(minutes)")
      
  args = p.parse_args()

  return args

## ---------------------------------------------------------------------- ##

def main():
  """
  read Elliptica properties file and create ETK param file
  """
  
  args = parse_cli()
  res   = args.r
  
  ## read these params from the properties file
  params = [ "NS1_center_x", "NS1_center_y", "NS1_max_radius",
             "NS2_center_x", "NS2_center_y", "NS2_max_radius",
             "NS1_EoS_K0", "NS1_EoS_Gamma",
             "NS2_EoS_K0", "NS2_EoS_Gamma",
             "NSNS_x_CM", "NSNS_y_CM","NSNS_separation"]
  
  ## read properties file and set the params
  ## NOTE: assumed every param is float
  ## NOTE: assumed single polytrope and both stars using the same eos
  param_dict = dict()
  ## init
  for p in params:
    param_dict[p] = ''

  with open(f'{args.p}', 'r') as fp:
    line = fp.readline()
    while line:
      line = fp.readline()
      ## set
      for p in params:
        if re.search(r'^{}[\W]'.format(p),line):
          ## trim new line and white spaces
          line = line.replace('\n', '')
          line = line.replace(' ', '')
          ## for eos
          line = line.replace('[', '')
          line = line.replace(']', '')
          param_dict[p] = (float(line.split('=')[1]))
  
  if not any(param_dict.values()):
    print(f'some params are empty\n{param_dict}\n');
    exit(-1)
    
  ## adjust coords. according to the CM
  param_dict['NS1_center_x'] -= param_dict['NSNS_x_CM']
  param_dict['NS1_center_y'] -= param_dict['NSNS_y_CM']
  param_dict['NS2_center_x'] -= param_dict['NSNS_x_CM']
  param_dict['NS2_center_y'] -= param_dict['NSNS_y_CM']
  
  ## set g_data and replace
  global g_data
  with open(f'{args.d}', 'r') as file:
    g_data = file.read()

## ---------------------------------------------------------------------- ##
  ## time:
  setd('cctk_final_time',args.ft)
  setd('dtfac',g_dtfac)
  seti('walltime',args.wt)
  
  ## NS positions
  setd('position_x_1',param_dict['NS1_center_x'])
  setd('position_y_1',param_dict['NS1_center_y'])
  setd('position_x_2',param_dict['NS2_center_x'])
  setd('position_y_2',param_dict['NS2_center_y'])
  
  ## regions
  ns_r1 = g_ns_expand*ceil(param_dict['NS1_max_radius'])
  ns_r2 = g_ns_expand*ceil(param_dict['NS2_max_radius'])
  
  ## region1
  setd('2**0*NS_radius1',ns_r1)
  setd('2**1*NS_radius1',2*ns_r1)
  setd('2**2*NS_radius1',2**2*ns_r1)
  setd('2**3*NS_radius1',2**3*ns_r1)
  setd('2**4*NS_radius1',2**4*ns_r1)
  setd('2**5*NS_radius1',2**5*ns_r1)
  
  ## region2
  setd('2**0*NS_radius2',ns_r2)
  setd('2**1*NS_radius2',2*ns_r2)
  setd('2**2*NS_radius2',2**2*ns_r2)
  setd('2**3*NS_radius2',2**3*ns_r2)
  setd('2**4*NS_radius2',2**4*ns_r2)
  setd('2**5*NS_radius2',2**5*ns_r2)

  ## grid:
  min_r = ns_r1 if ns_r1 < ns_r2 else ns_r2
  length = 2**6 * min_r # as biggest level is 2**5
  setd('xmin',-length)
  setd('xmax',length)
  setd('ymin',-length)
  setd('ymax',length)
  setd('zmin',-length)
  setd('zmax',length)
  seti('regrid_every',g_regrid)

  ## resolution:
  seti('ncells_x',args.r)
  seti('ncells_y',args.r)
  seti('ncells_z',args.r)

  ## id path
  sets('id_path',args.i)

  ## eos
  k = param_dict['NS1_EoS_K0']
  g = param_dict['NS1_EoS_Gamma']
  setd('poly_k',k)
  setd('poly_gamma',g)
 
  ## AHF
  setd('tmerger',args.mt)
## ---------------------------------------------------------------------- ##

  ## write back into the parfile
  d = param_dict['NSNS_separation']
  output=f'bns_gamma{g:0.1f}_k{k:0.1f}_n{int(args.r)}.par'
  with open(f'{output}', 'w') as file:
    file.write(f'## final time     = {args.ft:0.4f} M0\n')
    file.write(f'## merger time    = {args.mt:0.4f} M0\n')
    file.write(f'## walltime       = {args.wt/(60.):0.4f} h\n')
    file.write(f'## Length         = {2*length:0.4f} M0\n')
    file.write(f'## dx_reg1_lev6  <= {2*ns_r1/args.r:0.4f} M0\n')
    file.write(f'## dx_reg2_lev6  <= {2*ns_r2/args.r:0.4f} M0\n')
    file.write(f'## dx_NS required<= {min_r/64:0.4f} M0\n')
    file.write(f'## dx_coarsest    = {2*length/args.r:0.4f} M0\n')
    file.write(f'## CFL            = {g_dtfac:0.4f}\n')
    file.write(f'## ncells_xyz     = {args.r}\n')
    file.write(f'## BNS_separation = {d:0.4f} M0\n')
    file.write(f'## ID path        = {args.i}\n')
    file.write(g_data)
  
if __name__ == '__main__': main()
