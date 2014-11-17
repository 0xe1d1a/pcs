#!/usr/bin/python
__author__ = 'eldiaman'
import sys, getopt
import numpy as np
import matplotlib.pyplot as pyplot
import pprint


colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

ifile=''
dpath='.'
samples = 5
msg = "Usage: %s -i input (default mresults.log) -d directory (default .)  -s samples (default 5) \n \
		file format: time {space} elements_num {space} threads {newline}" % sys.argv[0]

def plot_lines(times, elements, err):

	for y in times:
		X = elements[y]
		Y = times[y]
		Z = err[y]
		Zup = []
		Zdn = []
		for i in range(len(Y)):
			Zup.append(Y[i] + Z[i])
			Zdn.append(Y[i] - Z[i])
		pyplot.plot( X, Y, '-', marker='o' )
		pyplot.plot (X, Zup, linestyle='None', marker='^')
		pyplot.plot (X, Zdn, linestyle='None', marker='v')
		pyplot.title( 'Mergesort plot' )
		pyplot.xlabel( 'elements' )
		pyplot.ylabel( 'time (seconds)' )
	pyplot.show()

try:
    myopts, args = getopt.getopt(sys.argv[1:],"i:d:s:")
except getopt.GetoptError as e:
    print (str(e))
    print msg
    sys.exit(2)
for o, a in myopts:
    if o == '-i':
        ifile = a
    elif o == '-d':
    	dpath = a
    elif o == '-s':
        samples = int(a)

if ifile == '':
	ifile = 'mresults.log'


f = open(dpath+'/'+ifile, "r")
lines = f.readlines()
i = 0
elements = []
threads = set()
data_times = dict()
data_times_err = dict()
data_elements = dict()
_times = []
for line in lines:
    tmp = line.split(" ")
    print tmp[1]
    _times.append(float(tmp[0])) #time
    i = i + 1
    if i==samples:
        i = 0
        times_avg = (sum(_times)/len(_times))
        times_err = (max(_times)-min(_times))/2
        key = tmp[2].strip('\n')
        if key not in data_times:
        	data_times[key] = []
        	data_elements[key] = []
        	data_times_err[key] = []
    	data_times[key].append(times_avg)
    	data_times_err[key].append(times_err)
    	data_elements[key].append(tmp[1])
    	_times = []

print data_times_err
print data_times
print data_elements

plot_lines(data_times, data_elements, data_times_err)
