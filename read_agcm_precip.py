__author__ = 'santiago'


from matplotlib.pyplot import plot, show
from datetime import datetime


def read_agcm_precip(input_file):
    data = open(input_file, 'r')
    lines = data.readlines()
    taxs = []
    prec = []
    prcv = []
    prge = []
    for line in lines:
        taxs.append(datetime.strptime(line.split('\t')[0], '%HZ%d%b%Y'))
        prec.append(float(line.split('\t')[1]))
        prcv.append(float(line.split('\t')[2]))
        prge.append(float(line.split('\t')[3]))
    plot(taxs, prec)
    show()


spdiag3_precip = '/home/santiago/Modelos/spdiag_results/SPDIAG03_precip'
read_agcm_precip(spdiag3_precip)
