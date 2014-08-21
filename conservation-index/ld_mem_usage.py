#! /usr/bin/python3.3
import resource

ITEMS = 1000
l = []
for i in range(1, 100000):
	l.append({'-': 0.0, 'A': 0.0, 'G': 0.0, 'C': 0.0, 'T': 0.0})
	if (len(l) % 10000) == 0:
		print('{:d} dicts needs {:d} KB'.format(len(l), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

