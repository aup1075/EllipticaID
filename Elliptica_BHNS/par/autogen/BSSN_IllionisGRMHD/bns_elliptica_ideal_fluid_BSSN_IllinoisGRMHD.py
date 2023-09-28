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
g_data = ''

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
  p.add_argument("-r", type=float, required=True,
      help="resolution, ex: 256, 128, etc.")
  p.add_argument("-t", type=float, default=1000.,
      help="merger time")
  p.add_argument("-b", type=str, default='n',
      help="bitant?(y/n)")
      
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
             "NSNS_x_CM", "NSNS_y_CM"]
  
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
    file.close()

  ## set params in the dummy ##
  
  ## id path
  sets('id_path',args.i)

  ## eos
  k = param_dict['NS1_EoS_K0']
  g = param_dict['NS1_EoS_Gamma']
  setd('poly_k',k)
  setd('poly_gamma',g)
  
  ## positions
  setd('position_x_1',param_dict['NS1_center_x'])
  setd('position_y_1',param_dict['NS1_center_y'])
 
  setd('position_x_2',param_dict['NS2_center_x'])
  setd('position_y_2',param_dict['NS2_center_y'])
  
  ## regions
  ns_r1 = 1.15*ceil(param_dict['NS1_max_radius'])
  ns_r2 = 1.15*ceil(param_dict['NS2_max_radius'])
  
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
  
  ## default length
  Max = 1024
  setd('xmin',-Max)
  setd('xmax',Max)
  setd('ymin',-Max)
  setd('ymax',Max)
  setd('zmax',Max)

  ## bitant?
  if args.b == 'y':
    setd('zmin',0)
    sets('reflection_z','yes')
    assert(args.r % 2 == 0)
    seti('ncells_z',args.r/2)
  else:
    setd('zmin',-Max)
    sets('reflection_z','no')
    seti('ncells_z',args.r)

  ## resolution
  seti('ncells_x',args.r)
  seti('ncells_y',args.r)
  
  ## CFL
  setd('dtfac',0.25)
  
  ## merger time
  setd('tmerger',args.t)

  ## checkpoint every
  seti('checkpoint_every',8192)

  ## write back into the parfile
  output=f'bns_gamma{g:0.1f}_k{k:0.1f}_n{int(args.r)}.par'
  with open(f'{output}', 'w') as file:
    file.write(g_data)
    file.close()
  
if __name__ == '__main__': main()

