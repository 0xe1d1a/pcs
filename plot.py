#!/usr/bin/python
__author__ = 'eldiaman'
import sys, getopt
import numpy as np
import matplotlib.pyplot as plt
import pprint

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']


def plot_speedups(data):
    n_groups = len(data.keys())

    fig, ax = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.2

    opacity = 0.4
    error_config = {'ecolor': '0.3'}

    col_based_values = []
    col_based_errors = []
    threads = None
    i = 0
    for key in data:
        threads = sorted(data[key].keys())
        values = []
        errors = []
        for thread in threads:
            #print data[key][thread] 
            l = data[key][thread]
            values.append(l[0])
        col_based_values.append(values)
    for j in range (0, len(col_based_values[0])):
        tmp_val = []
        for i in range(0,n_groups): 
            tmp_val.append(col_based_values[i][0]/col_based_values[i][j])
            if threads[j] == '1':
                lab = threads[j] + " thread"
            else:
                lab = threads[j] + " threads"
        plt.bar(index + bar_width * j, tmp_val, bar_width,
                         alpha=opacity,
                         color=colors[j],
                         #yerr=tmp_err,
                         error_kw=error_config,
                         label= lab)


    plt.xlabel('Dimensions')
    plt.ylabel('times faster')
    plt.title('Speedup (T1/TN)')
    plt.xticks(index + bar_width, data.keys())
    plt.legend()

    plt.tight_layout()
    plt.grid()
    plt.show()

def plot_flops(data):
    n_groups = len(data.keys())

    fig, ax = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.2

    opacity = 0.4
    error_config = {'ecolor': '0.3'}

    col_based_values = []
    col_based_errors = []
    threads = None
    i = 0
    for key in data:
        threads = sorted(data[key].keys())
        values = []
        errors = []
        for thread in threads:
            #print data[key][thread] 
            l = data[key][thread]
            values.append(l[2])
            errors.append(l[3])
        col_based_values.append(values)
        col_based_errors.append(errors)
    for j in range (0, len(col_based_values[0])):
        tmp_val = []
        tmp_err = []
        for i in range(0,n_groups): 
            tmp_val.append(col_based_values[i][j])
            tmp_err.append(col_based_errors[i][j])
            if threads[j] == '1':
                lab = threads[j] + " thread"
            elif threads[j] == 'seq':
                lab = threads[j] + " version"
            else:
                lab = threads[j] + " threads"
        plt.bar(index + bar_width * j, tmp_val, bar_width,
                         alpha=opacity,
                         color=colors[j],
                         yerr=tmp_err,
                         error_kw=error_config,
                         label= lab)


    plt.xlabel('Dimensions')
    plt.ylabel('gFLOPS')
    plt.title('Effectiveness')
    plt.xticks(index + bar_width, data.keys())
    plt.legend()

    plt.tight_layout()
    plt.grid()
    plt.show()

def plot_times(data):
    n_groups = len(data.keys())

    fig, ax = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.2

    opacity = 0.4
    error_config = {'ecolor': '0.3'}

    col_based_values = []
    col_based_errors = []
    threads = None
    i = 0
    for key in data:
        threads = sorted(data[key].keys())
        values = []
        errors = []
        for thread in threads:
            #print data[key][thread] 
            l = data[key][thread]
            values.append(l[0])
            errors.append(l[1])
        col_based_values.append(values)
        col_based_errors.append(errors)
    for j in range (0, len(col_based_values[0])):
        tmp_val = []
        tmp_err = []
        for i in range(0,n_groups): 
            tmp_val.append(col_based_values[i][j])
            tmp_err.append(col_based_errors[i][j])
            if threads[j] == '1':
                lab = threads[j] + " thread"
            elif threads[j] == 'seq':
                lab = threads[j] + " version"
            else:
                lab = threads[j] + " threads"
        plt.bar(index + bar_width * j, tmp_val, bar_width,
                         alpha=opacity,
                         color=colors[j],
                         yerr=tmp_err,
                         error_kw=error_config,
                         label=lab)


    plt.xlabel('Dimensions')
    plt.ylabel('seconds')
    plt.title('Wall time')
    plt.xticks(index + bar_width, data.keys())
    plt.legend()

    plt.tight_layout()
    plt.grid()
    plt.draw()

ifile=''
dpath='.'
samples = 10
msg = "Usage: %s -i input (default results.log) -d directory (default .)  -s samples (default 10) \n" % sys.argv[0]


try:
    myopts, args = getopt.getopt(sys.argv[1:],"i:d:I:s:")
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
	ifile = 'results.log'


f = open(dpath+'/'+ifile, "r")
lines = f.readlines()

data = dict()
data_flops = []
data_times = []
i = 0
for line in lines:
    tmp = line.split(" ")
    tmp[3] = tmp[3].strip('\n') #cols
    tmp[2] = tmp[2] #rows
    data_flops.append(float(tmp[1])/1000000000) #flops
    data_times.append(float(tmp[0])) #time
    i = i + 1
    if i==samples:
        i = 0
        flops_avg = (sum(data_flops)/len(data_flops))
        flops_err = (max(data_flops)-min(data_flops))/2
        times_avg = (sum(data_times)/len(data_times))
        times_err = (max(data_times)-min(data_times))/2
        key = tmp[2]+'x'+tmp[3]
        if key not in data:
            data[key] = dict()
        threads = tmp[4].strip('\n')
        data[key][threads] = [times_avg, times_err, flops_avg, flops_err]
        data_times = []
        data_flops = []

if i!=0:
    print "Something went wrong..Check results.log or maybe change sample?\n%s" % (msg)
    sys.exit(2)
pprint.pprint(data)

plot_times(data)
plot_flops(data)
plot_speedups(data)
