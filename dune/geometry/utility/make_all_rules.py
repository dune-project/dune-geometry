#!/usr/bin/python
import subprocess

geometries = ['prism', 'simplex2d', 'simplex3d', 'pyramid', 'cube2d', 'cube3d', 'jacobi1',
              'jacobi2', 'gauss', 'gausslobatto']

for g in geometries:
  print g
  subprocess.call('python make_rules.py ' + g, shell=True)
