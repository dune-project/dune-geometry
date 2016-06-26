import os
import csv
import fnmatch
import sys

from mpmath import *
from mako.template import Template

def split_reference(ref):
  words = ref.split()
  lines = []
  line = ""
  for w in words:
    if len(line + w) > 72:
      lines.append(line.strip())
      line = w + " "
    else:
      line += w + " "
  if len(line.strip()) > 0:
    lines.append(line.strip())

  return lines

#-------------------------------------------------------------------------------

# provide reference values for the volume, to test for correct normalization of the weights
reference_volume = { 'prism': 0.5, 'simplex2d': 0.5, 'simplex3d': 1.0/6.0,
                     'pyramid': 1.0/3.0, 'cube2d': 1.0, 'cube3d': 1.0, 'jacobi1': 1.0,
                     'jacobi2': 4.0/3.0, 'gauss': 1.0, 'gausslobatto': 1.0 }

translate = { 'prism': 'prism', 'simplex2d': 'tri', 'simplex3d': 'tet',
              'pyramid': 'pyr', 'cube2d': 'quad', 'cube3d': 'hex', 'jacobi1': 'jacobi1',
              'jacobi2': 'jacobi2', 'gauss': 'gauss', 'gausslobatto': 'gausslobatto'}

geo = { 'prism': 'prism', 'simplex2d': 'simplex', 'simplex3d': 'simplex',
              'pyramid': 'pyramid', 'cube2d': 'cube', 'cube3d': 'cube', 'jacobi1': 'simplex',
              'jacobi2': 'simplex', 'gauss': 'simplex', 'gausslobatto': 'simplex'}


name = 'prism'
geometry = 'prism'

if len(sys.argv) > 1:
  name = sys.argv[1]
if not name in reference_volume:
  sys.stderr.write('ERROR: Given name "' + name + '" is not allowed. Choose one of (prism|simplex2d|simplex3d|pyramid|cube2d|cube3d)!')

geometry = geo[name]

mp.dps = 80
rules_dir = './rules/' + translate[name]
if not os.path.exists(rules_dir):
  sys.stderr.write('ERROR: Given rules directory "' + rules_dir + '" does not exist!')

directory = os.listdir(rules_dir)
filenames = fnmatch.filter(directory, '*.csv')

n_files = len(filenames)
if n_files == 0:
  sys.stderr.write('ERROR: No files found in rules directory "' + rules_dir + '"!')

# extract maximal order from all files.
maxorder = 0
for f_ in filenames:
  f = os.path.join(rules_dir, f_)
  with open(f, 'r') as quadfile:
    quadreader = csv.reader(quadfile, delimiter=' ', skipinitialspace=True)

    header = next(quadreader)
    while header[0][0] == '#':
      header = next(quadreader)

    dim = int(float(header[0]))
    order = int(float(header[1]))
    maxorder = max(maxorder, order)

rules = [None]*(maxorder+1)

for f_ in filenames:
  f = os.path.join(rules_dir, f_)

  rule = {'points' : [], 'weights': [], 'reference': []}
  with open(f, 'r') as quadfile:
    quadreader = csv.reader(quadfile, delimiter=' ', skipinitialspace=True)

    header = next(quadreader)
    reference = ''
    while header[0][0] == '#':
      line = " ".join(header)
      reference += line[1:].strip() + ' '
      header = next(quadreader)

    if len(reference) > 0:
      rule['reference'] = split_reference(reference)

    dim = int(float(header[0]))
    order = int(float(header[1]))
    if dim < 1 or dim > 3:
      sys.stderr.write('ERROR: Wrong dimension given: dim = ' + str(dim) + '!')

    sum_weights = mpf(0)
    for row in quadreader:
      if len(row) < 1:
        break
      if len(row) != dim+1:
        sys.stderr.write('ERROR: Wrong number of columns in csv file in row "' + row + '"!')

      p = []
      for i in range(0,dim):
        p.append(mpf(row[i]))
      w = mpf(row[dim])

      rule['points'].append(p)
      rule['weights'].append(w)
      sum_weights += w

    if abs(sum_weights - reference_volume[name]) < 1.e-10:
      if not rules[order] or not rules[order]['points'] or len(rules[order]['points']) == 0 or \
         len(rules[order]['points']) > len(rule['points']):
        rules[order] = rule
    else:
      sys.stderr.write('ERROR: Weights are not normalized correctly. Expected: ' + str(reference_volume[name]) + ', but given: ' + str(sum_weights) + '! In file: "' + f + '"')

max_i = len(rules)-1
while not rules[max_i] or not rules[max_i]["points"] or len(rules[max_i]["points"]) == 0:
  if max_i == 0:
    sys.stderr.write('ERROR: No valid quadrature rule created!')

  max_i -= 1

if max_i != maxorder:
  sys.stderr.write('ERROR: max_i != maxorder: ' + str(max_i) + ' != ' + str(maxorder))

dim = len(rules[maxorder]["points"][0])
maxp = len(rules[maxorder]["points"])

placeholders = {
  'dim': dim,
  'maxp': maxp,
  'maxorder': maxorder,
  'name': name.capitalize(),
  'geometry': geometry,
  'rules': rules,
  #'cast' : lambda s: s,
  'cast' : lambda s: "cast<ct>(\"" + s + "\")",
  'uppercase': lambda s: s.upper()}

code = Template(filename='rules.templ.hh')
with open('../quadraturerules/' + name + 'quadrature.inc.hh', 'w') as out:
  out.write(code.render(**placeholders))
