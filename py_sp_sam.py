# -*- coding: utf-8 -*-

__name__ = 'Py_SP_SAM'
__author__ = 'Paulo Santiago'
__email__ = 'paulohsm@gmail.com'
__licence__ = 'FY - Fuck You'
__date__ = '2015/09/18'


# Paths for important data:
diagin_ara02 = '/home/santiago/Modelos/SEMIPROG_IN/SEMIPROG_IN-SPDARA02'
diadin_kuo02 = '/home/santiago/Modelos/SEMIPROG_IN/SEMIPROG_IN-SPDKUO02'
diagin_zmc01 = '/home/santiago/Modelos/SEMIPROG_IN/SEMIPROG_IN-SPDKUO02'

diagout_ara02 = '/home/santiago/Modelos/exps_spsam_2015-09-24/ARA02_032x001x4000x20s/SEMIPROG_OUT'
diagout_kuo02 = '/home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT'
diagout_zmc01 = '/home/santiago/Modelos/exps_spsam_2015-09-24/ZMC01_032x001x4000x20s/SEMIPROG_OUT'

agcm1d_ara02 = '/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-ARA02'
agcm1d_kuo02 = '/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-KUO02'
agcm1d_zmc01 = '/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-ZMC01'

trmm200401 = '/home/santiago/Datasets/TRMM-3B42RT/200401'
lonfor = 360.0-39.5
latfor = -4.0


def read_semiprogin(input_file, return_var):
    """
    This Python function simply read a variable from SEMIPROG_IN file given the
    arguments below. SEMIPROG_IN is the file used to store initial conditions
    for SAM superparameterization diagnostic tests.
    :param input_file: input file name, full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """

    import sys

    var_list = ['nt', 'nz', 'tstep', 'lon', 'lat', 'topo', 'lmsk', 'pres',
                'time', 'pslc', 'usst', 'vsst', 'cssf', 'clsf', 'ocis', 'oces',
                'iswf', 'roce', 'olis', 'oles', 'role', 'swtc', 'ocic', 'lwtc',
                'lwbc', 'temp', 'umes', 'liqm', 'icem', 'uvel', 'vvel', 'swrh',
                'lwrh']

    if not return_var in var_list:
        print 'You must provide a valid variable name. Pick one from the list' \
              'below:'
        print var_list
        sys.exit("Exit now. Try again.")

    print 'Loading SEMIPROG_IN (', input_file, '): ', return_var
    var_dict = dict(zip(var_list, range(len(var_list))))
    data = open(input_file, 'r')
    lines = data.readlines()

    nz = int(lines[0].split('\t')[0])
    nt = (len(lines) - 1) / (nz + 2)

    # The following conditional complex structure is necessary in order to take
    # into account the way the data is organized in SEMIPROG_IN file. Every
    # conditional refers to a 'section' in the file.
    if var_dict[return_var] == 0:
        return nt
    elif var_dict[return_var] == 1:
        return nz
    elif var_dict[return_var] in range(2,7):
        var_vals = float(lines[0].split('\t')[var_dict[return_var]-1])
        return var_vals
    elif var_dict[return_var] == 7:
        var_vals = []
        for k in range(nz):
            var_vals.append(float(lines[k+3].split('\t')[0]))
        return var_vals
    elif var_dict[return_var] == 8:
        var_vals = []
        for l in range(nt):
            ll = 1 + l * (nz + 1)
            var_vals.append(lines[ll].split('\t')[0])
        return var_vals
    elif var_dict[return_var] in range(9,17):
        var_vals = []
        for l in range(nt):
            idx = var_dict[return_var] - 8
            ll = 1 + l * (nz + 1)
            var_vals.append(float(lines[ll].split('\t')[idx]))
        return var_vals
    elif var_dict[return_var] in range(17,25):
        var_vals = []
        for l in range(nt):
            idx = var_dict[return_var] - 17
            ll = 2 + l * (nz + 1)
            var_vals.append(float(lines[ll].split('\t')[idx]))
        return var_vals
    else:
        idx = var_dict[return_var] - 24
        var_vals = [[0 for x in range(nz)] for x in range(nt)]
        for l in range(nt):
            for k in range(nz):
                lk = 3 + l * (nz + 1) + k
                var_vals[l][k] = float(lines[lk].split('\t')[idx])
        return var_vals


def read_semiprogout(input_file, return_var):
    """
    Reads data from SAM superparameterization (SP) diagnostic test output
    :param input_file: input file name - full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """
    import sys
    var_list = ['t0', 'nx', 'ny', 'nz', 'nt', 'dt', 'lon', 'lat', 'topo',
                'ocnfrac', 'pmid', 'pdel', 'precc', 'precl', 'precsc', 'precsl',
                'cltot', 'clhgh', 'clmed', 'cllow', 'taux_crm', 'tauy_crm',
                'z0m', 'prectend', 'precstend', 'zmid', 'ultend', 'vltend',
                'qltend', 'qcltend', 'qcitend', 'sltend', 'cld', 'cldr',
                'cldtop', 'gicewp', 'gliqwp', 'mc', 'mcup', 'mcdn', 'mcuup',
                'mcudn', 'crm_qc', 'crm_qi', 'crm_qs', 'crm_qg', 'crm_qr',
                'tkez', 'tkesgsz', 'flux_u', 'flux_v', 'flux_qt', 'fluxsgs_qt',
                'flux_qp', 'pflx', 'qt_ls', 'qt_trans', 'qp_trans', 'qp_fall',
                'qp_evp', 'qp_src', 't_ls']

    if not return_var in var_list:
        print 'You must provide a valid variable name. Pick one from the list' \
              'below:'
        print var_list
        sys.exit("Exit now. Try again.")

    print 'Loading SEMIPROG_OUT (', input_file, '): ', return_var
    var_dict = dict(zip(var_list, range(len(var_list))))
    data = open(input_file, 'r')
    lines = data.readlines()

    nvar = 38
    nrec = len(lines) - 3 / nvar

    nz = int(lines[0].split()[3])
    nt = int(lines[0].split()[4])

    if var_dict[return_var] == 0:
        return lines[0].split()[0]
    elif var_dict[return_var] in range(1,5):
        return int(lines[0].split()[var_dict[return_var]])
    elif var_dict[return_var] in range(5,10):
        return float(lines[0].split()[var_dict[return_var]])
    elif var_dict[return_var] == 10:
        var_vals = []
        for element in lines[1].split():
            var_vals.append(float(element))
        return var_vals
    elif var_dict[return_var] == 11:
        var_vals = []
        for element in lines[2].split():
            var_vals.append(float(element))
        return var_vals
    elif var_dict[return_var] in range(12,26):
        idx = var_dict[return_var] - 12
        var_vals = []
        for l in range(nt):
            var_vals.append(float(lines[3 + l * nvar].split()[idx]))
        return var_vals
    else:
        var_vals = [[0 for x in range(nz)] for x in range(nt)]
        for l in range(nt):
            for k in range(nz):
                n = 4 + l * nt + var_dict[return_var] - 26
                var_vals[l][k] = float(lines[n].split()[k])
                print var_vals[l][k]
        return var_vals


def get_semiprog_time(input_file):
    """
    Given the file path, provides an hourly time axis for plotting SEMIPROG_OUT
    variables
    :param input_file: input file name - full or relative path;
    :return: a list of datetime objects.
    """

    from datetime import datetime

    t0 = read_semiprogout(input_file, 't0')
    nt = read_semiprogout(input_file, 'nt')
    t0py = datetime.strptime(t0, "%HZ%d%b%Y")
    timelist = []
    for tt in range(nt):
        timelist.append


def read_agcmdiagfields1d(input_file, return_var):

    """
    This function reads the variables stored in AGCM's diagnostic surface fields
    given the input file name and the variable.
    :param input_file: input file name, full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """

    import sys

    var_list = ['time', 'uves', 'vves', 'tems', 'umrs', 'agpl', 'tsfc', 't02m',
                 'q02m', 'u10m', 'v10m', 'prec', 'prcv', 'neve', 'rnof', 'evap',
                 'cbnv']

    if not return_var in var_list:
        print 'You must provide a valid variable name. Pick one from the list' \
              'below:'
        print var_list
        sys.exit("Exit now. Try again.")

    print 'Loading AGCM 1D fields (', input_file, '): ', return_var
    var_dict = dict(zip(var_list, range(len(var_list))))
    data = open(input_file, 'r')
    lines = data.readlines()

    var_vals = []
    for ll in range(len(lines)-1):
        l = ll + 1
        if var_dict[return_var] == 0:
            var_vals.append(lines[l].split('\t')[0])
        else:
            var_vals.append(float(lines[l].split('\t')[var_dict[return_var]]))

    return var_vals


def read_agcmdiagfields2d(input_file, return_var):
    """
    This function reads the variables stored in AGCM's diagnostic surface fields
    given the input file name and the variable.
    :param input_file: input file name, full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """

    import sys

    var_list = ['time', 'pres', 'omeg', 'temv', 'zgeo', 'umrl', 'cvlh', 'cvms',
                'lglh', 'lgms', 'acvr', 'ctot', 'cinv', 'csat', 'clwd', 'tke2']

    if not return_var in var_list:
        print 'You must provide a valid variable name. Pick one from the list' \
              'below:'
        print var_list
        sys.exit("Exit now. Try again.")

    print 'Loading AGCM 2D fields (', input_file, '): ', return_var
    var_dict = dict(zip(var_list, range(len(var_list))))
    data = open(input_file, 'r')
    lines = data.readlines()

    nz = 28
    nt = len(lines) / (nz + 1)

    if var_dict[return_var] == 0:
        var_vals = [0 for x in range(nt)]
        for l in range(nt):
            n = l * (nz + 1)
            var_vals[l] = lines[n].split('\t')[0]
        return var_vals
    elif var_dict[return_var] == 1:
        var_vals = [0 for x in range(nz)]
        for k in range(nz):
            var_vals[k] = lines[k+1].split('\t')[0]
        return var_vals
    else:
        var_vals = [[0 for x in range(nz)] for x in range(nt)]
        idx = var_dict[return_var] - 1
        for l in range(nt):
            for k in range(nz):
                nk = l * (nz + 1) + (k + 1)
                var_vals[l][k] = float(lines[nk].split('\t')[idx])
        return var_vals


def accum3(var):
    """
    Computes 3 steps accumulation of any variable (useful to get 3-hour
    accumulated precipitation from hourly measurements)
    :param var: hourly precipitation for a specific location.
    :return: 3 steps accumulated values.
    """
    accum3 = []
    for step in range(2,len(var),3):
        accum3.append(var[step] + var[step-1] + var[step-2])

    return accum3


def get_num(x):
    """
    Returns numbers containted in strings.
    :param x: a string.
    :return: the numbers contained in the string.
    """
    return float(''.join(ele for ele in x if ele.isdigit() or ele=='.'))


def find_nearest_idx(value, array):
    """
    Given a single value, returns the index of the nearest value in an
    array/list.
    :param value: the value one want to check the nearest index in the array.
    :param array: the array where one wants to check the index of the nearest
    value.
    :return: the array index (integer) of the nearest value.
    """
    from numpy import abs
    idx = (abs(array-value)).argmin()
    #print 'Nearest grid point and index:', array[idx], idx
    return idx #return array[idx]


def load_trmm_series(path, long, lati):
    """
    Loads TRMM 3B42RT data as a time series for tha given lon, lat location.
    :param path: the location of TRMM 3B42RT data.
    :param long: longitude in 0:360 range
    :param lati: latitude in -90:90 range.
    :return: a large matrix of precipitation data.
    """

    from pytrmm import TRMM3B42RTFile
    from datetime import datetime
    from glob import glob
    from numpy import linspace

    # loading TRMM data
    pref = '3B42RT.'
    suff = '.7R2.bin.gz'
    time_stamps = []
    trmm_series = []
    print 'Loading a series of TRMM matrices... please wait.'
    filelist = sorted(glob(path + '/' + pref + '*' + suff))
    for filepath in filelist:
        trmmfile = TRMM3B42RTFile(filepath)
        time_string = trmmfile.header()['granule_ID'].split(".")[1]
        print time_string
        time_stamps.append(datetime.strptime(time_string, "%Y%m%d%H"))
        trmm_series.append(trmmfile.precip())

    # variables used to define grids and axes
    nlon = int(trmmfile.header()['number_of_longitude_bins'])
    nlat = int(trmmfile.header()['number_of_latitude_bins'])
    ntim = len(trmm_series)
    # 'first_box_center': '59.875N,0.125E'
    # 'last_box_center': '59.875S,359.875E'
    lon0 = get_num(trmmfile.header()['first_box_center'].split(",")[1])
    lat0 = -get_num(trmmfile.header()['first_box_center'].split(",")[0])
    lon1 = get_num(trmmfile.header()['last_box_center'].split(",")[1])
    lat1 = get_num(trmmfile.header()['last_box_center'].split(",")[0])
    # creates zonal and meridional axes
    lons = linspace(lon0, lon1, nlon)
    lats = linspace(lat0, lat1, nlat)
    lon_idx = find_nearest_idx(lons, long)
    lat_idx = find_nearest_idx(lats, lati)
    trmm_point = []
    for it in range(ntim):
        trmm_point.append(trmm_series[it][lat_idx, lon_idx])

    del(trmm_series)
    return time_stamps, trmm_point


def trmm2npy(path):
    """
    Converts a TRMM 3B42RT time sequence in a Numpy data file.
    :param path: the path where TRMM data is located; also the path where npy
    file will be saved.
    :return: a Numpy file containing TRMM data.
    """

    from pytrmm import TRMM3B42RTFile
    from glob import glob
    from numpy import asarray, save

    pref = '3B42RT'
    suff = '7R2.bin.gz'
    #time_stamps = []
    trmm_series = []

    print 'Loading a series of TRMM matrices... please wait.'
    filelist = sorted(glob(path + '/' + pref + '.*.' + suff))
    for filepath in filelist:
        trmmfile = TRMM3B42RTFile(filepath)
     #   time_string = trmmfile.header()['granule_ID'].split(".")[1]
     #   time_stamps.append(datetime.strptime(time_string, "%Y%m%d%H"))
        trmm_series.append(trmmfile.precip())

    trmmdata = asarray(trmm_series)
    save(path + '/' + pref, trmmdata)


def trmm2file(path, long, lati):
    """
    Saves a timeseries of TRMM data for given lon, lat coordinates
    :param path: Path for TRMM data.
    :param long: Longitude.
    :param lati: Latitude.
    :return: saves the data in file "trmm2file.txt".
    """
    time, prec = load_trmm_series(path, long, lati)
    print type(time)
    fout = open("trmm2file.txt", "w")
    for item in range(len(prec)):
        print type(time[item])
        fout.write(";".join([str(time[item]), str(prec[item])]))
    fout.close()


def trmmfile2py():
    """
    Reads as lists the data stored in "trmm2file.txt" provided by function
    'trmm2file'.
    :return: Time and respective data values.
    """
    trmmfile = open('trmm2file.txt', 'r')
    time = []
    data = []
    for row in trmmfile:
        time.append(row.strip().split(';')[0])
        data.append(row.strip().split(';')[1])

    return time, data



def plot_prec_ara02():
    from matplotlib.pyplot import plot, show
    trmm_time, trmm_prec = load_trmm_series(trmm200401, lonfor, latfor)
    agcm_prec = accum3(read_agcmdiagfields1d(agcm1d_ara02, 'prec'))
    diag_prec = accum3(read_semiprogout(diagout_ara02, 'precc'))
    plot(trmm_time, trmm_prec, trmm_time, agcm_prec, trmm_time, diag_prec)
    show()


def barplot_prec_trmm():
    from matplotlib import pyplot as p

    fig = p.figure()

    width = 0.25
    time, prec = load_trmm_series(trmm200401, lonfor, latfor)
    p.bar(time, prec)
    #p.xticks(time + width / 2)

    fig.autofmt_xdate()
    p.show()


def sbs_prec():
    import matplotlib.pyplot as p
    import seaborn as s

    f = p.plot(figsize=(8,6))
    time, prec = load_trmm_series(trmm200401, lonfor, latfor)
    s.barplot(time, prec, palette="BuGn_d")
    f.set_ylabel("Precipitação")
    p.show()


def plot_trmm_prec():
    from matplotlib.pyplot import plot, legend, show, title
    time, prec = load_trmm_series(trmm200401, lonfor, latfor)
    agcm_prec = accum3(read_agcmdiagfields1d(agcm1d_ara02, 'prec'))
    diag_prec = 1000*3600*accum3(read_semiprogout(diagout_ara02, 'precc'))
    #plot(time, prec, color='black', linewidth=2, linestyle='-.', label="Precipitação TRMM")
    plot(time, prec, color='black', linewidth=2, label="TRMM")
    plot(time, agcm_prec, color='red', linewidth=2, label="MCGA")
    plot(time, diag_prec, color='blue', linewidth=2, label="Superparametrizacao")
    title('Precipitacao Fortaleza')
    legend()
    show()