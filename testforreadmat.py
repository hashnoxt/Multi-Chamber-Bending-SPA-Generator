# -*- coding: utf-8 -*-
import model_utility as MU
from keyword import kwlist

mat = MU.read_matfile('ogden3.mat')

print mat['D']

print kwlist

a = 5.026
b = round(a, 1)
print b
