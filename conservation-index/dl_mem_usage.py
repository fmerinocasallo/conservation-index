#! /usr/bin/python3.3
import resource

ITEMS = 1000
elems = ['-', 'A', 'G', 'C', 'T']
dict = {}
for k in range(10, 100):
	size = ITEMS * k
	for key in elems:
		dict[key] = [0.0] * size
	if k % 10 == 0:
		print('{:d} items per list needs {:d} KB'.format(size, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
	dict.clear()

