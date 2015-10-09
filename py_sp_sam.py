# -*- coding: utf-8 -*-

__name__ = 'Py_SP_SAM'
__author__ = 'Paulo Santiago'
__email__ = 'paulohsm@gmail.com'
__licence__ = 'FY - Fuck You'
__date__ = '2015/09/18'


def read_semiprogin(input_file, return_var):
    """
    This Python function simply read a variable from SEMIPROG_IN file given the
    arguments below. SEMIPROG_IN is the file used to store initial conditions
    for SAM superparameterization diagnostic tests.
    :param input_file: input file name, full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """
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
    import sys
    """
    Reads data from SAM superparameterization (SP) diagnostic test output
    :param input_file: input file name - full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """
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

    var_dict = dict(zip(var_list, range(len(var_list))))
    data = open(input_file, 'r')
    lines = data.readlines()

    nrec = len(lines) - 3 / 38
    print nrec

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
            var_vals.append(float(lines[l+3].split()[idx]))
        return var_vals
    else:
        var_vals = [[0 for x in range(nz)] for x in range(nt)]
        for l in range(nt):
            for k in range(nz):
                n = 4 + l * nt + var_dict[return_var] - 26
                var_vals[l][k] = float(lines[n].split()[k])
                print var_vals[l][k]
        return var_vals


def read_agcmdiagfields1d(input_file, return_var):
    import sys
    """
    This function reads the variables stored in AGCM's diagnostic surface fields
    given the input file name and the variable.
    :param input_file: input file name, full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """
    var_list = ['time', 'uves', 'vves', 'tems', 'umrs', 'agpl', 'tsfc', 't02m',
                 'q02m', 'u10m', 'v10m', 'prec', 'prcv', 'neve', 'rnof', 'evap',
                 'cbnv']

    if not return_var in var_list:
        print 'You must provide a valid variable name. Pick one from the list' \
              'below:'
        print var_list
        sys.exit("Exit now. Try again.")

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
    import sys
    """
    This function reads the variables stored in AGCM's diagnostic surface fields
    given the input file name and the variable.
    :param input_file: input file name, full or relative path;
    :param return_var: the name of the variable one want to read in the input.
    :return: a Python list of values relative to the variable name provided.
    """
    var_list = ['time', 'pres', 'omeg', 'temv', 'zgeo', 'umrl', 'cvlh', 'cvms',
                'lglh', 'lgms', 'acvr', 'ctot', 'cinv', 'csat', 'clwd', 'tke2']

    if not return_var in var_list:
        print 'You must provide a valid variable name. Pick one from the list' \
              'below:'
        print var_list
        sys.exit("Exit now. Try again.")

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






time_agcm = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-GRE01', 'time')
prec_agcm = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-GRE01', 'prec')

print accum3(prec_agcm)


#time1d = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-GRE01', 'time')
#uves = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-GRE01', 'uves')

#print len(uves)
#temv = read_agcmdiagfields2d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_2D-GRE01', 'temv')

#temp = read_semiprogin('/home/santiago/Modelos/SEMIPROG_IN/SEMIPROG_IN-SPDARA02', 'temp')
#print temp

#tkez = read_semiprogout('/home/santiago/Modelos/exps_spsam_1200s/ARA02_032x001x4000x20s/SEMIPROG_OUT', 'tkez')
#tkez = read_semiprogout('/home/santiago/Modelos/exps_spsam_1200s/ZMC01_032x001x4000x20s/SEMIPROG_OUT', 'tkez')
#print tkez

# working simulations output:
# /home/santiago/Modelos/exps_spsam_2015-09-24/ARA02_032x001x4000x20s/SEMIPROG_OUT
# /home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT
# /home/santiago/Modelos/exps_spsam_2015-09-24/ZMC01_032x001x4000x20s/SEMIPROG_OUT

