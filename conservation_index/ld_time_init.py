#! /usr/bin/python3.3
import timeit

for k in range(1, 10):
	size = 10000 * k
	time = timeit.timeit("[{'-': 0.0, 'A': 0.0, 'G': 0.0, 'C': 0.0, 'T': 0.0}] * size",
                             setup="from __main__ import size")
	print('time init {:d}: {:.4f} s'.format(size, time))

