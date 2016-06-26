import os
import fnmatch
import sys
import subprocess

#-------------------------------------------------------------------------------

# provide reference values for the volume, to test for correct normalization of the weights
translate = { 'prism': 'prism', 'tri': 'simplex2d', 'tet': 'simplex3d',
              'pyr': 'pyramid', 'quad': 'cube2d', 'hex': 'cube3d', 'jacobi1': 'jacobi1',
              'jacobi2': 'jacobi2', 'gauss': 'gauss'}


rules_dir = './rules'
if not os.path.exists(rules_dir):
  sys.stderr.write('ERROR: Given rules directory "' + rules_dir + '" does not exist!')

directory = os.listdir(rules_dir)
subdirs = fnmatch.filter(directory, '*')

n_dirs = len(subdirs)
if n_dirs == 0:
  sys.stderr.write('ERROR: No subdirectories found in rules directory "' + rules_dir + '"!')

for d_ in subdirs:
  d = os.path.join(rules_dir, d_)
  if not os.path.isdir(d):
    continue

  if not d_ in translate:
    sys.stderr.write('ERROR: Unknown geometry type "' + d_ + '"!')

  rule = translate[d_]
  subdir = os.listdir(d)
  filenames = fnmatch.filter(subdir, '*.txt')

  for f_ in filenames:
    f = os.path.join(d, f_)

    order = int(float(f_.split('-')[0]))
    basename = os.path.splitext(f)[0]

    print './convert_rules ' + f + ' ' + str(order) + ' ' + rule + ' > ' + basename + '.csv'
    subprocess.call('./convert_rules ' + f + ' ' + str(order) + ' ' + rule + ' > ' + basename + '.csv', shell=True)
