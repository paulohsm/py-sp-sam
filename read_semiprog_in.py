__author__ = 'santiago'
__name__ = 'read_semiprog_in'

def readsemiprogin(input_file, return_var):
    data = open(input_file, 'r')

    lines = data.readlines()

    nz = int(lines[0].split('\t')[0])
    tstep = float(lines[0].split('\t')[1])
    lon = float(lines[0].split('\t')[2])
    lat = float(lines[0].split('\t')[3])
    topo = float(lines[0].split('\t')[4])
    lmsk = float(lines[0].split('\t')[5])

    nrec = (len(lines) - 1) / (nz + 2)

    pres = [0 for x in range(nz)]
    time = [0 for x in range(nrec)]
    pslc = [0 for x in range(nrec)]
    usst = [0 for x in range(nrec)]
    vsst = [0 for x in range(nrec)]
    cssf = [0 for x in range(nrec)]
    clsf = [0 for x in range(nrec)]
    ocis = [0 for x in range(nrec)]
    oces = [0 for x in range(nrec)]
    iswf = [0 for x in range(nrec)]
    roce = [0 for x in range(nrec)]
    olis = [0 for x in range(nrec)]
    oles = [0 for x in range(nrec)]
    role = [0 for x in range(nrec)]
    swtc = [0 for x in range(nrec)]
    ocic = [0 for x in range(nrec)]
    lwtc = [0 for x in range(nrec)]
    lwbc = [0 for x in range(nrec)]
    temp = [[0 for x in range(nz)] for x in range(nrec)]
    umes = [[0 for x in range(nz)] for x in range(nrec)]
    liqm = [[0 for x in range(nz)] for x in range(nrec)]
    icem = [[0 for x in range(nz)] for x in range(nrec)]
    uvel = [[0 for x in range(nz)] for x in range(nrec)]
    vvel = [[0 for x in range(nz)] for x in range(nrec)]
    swrh = [[0 for x in range(nz)] for x in range(nrec)]
    lwrh = [[0 for x in range(nz)] for x in range(nrec)]
    #lwrh2 = [ [0] * nz ] * nrec

    for l in range(nrec):
        n = l * (nz + 2) + 1
        time[l] = lines[n].split('\t')[0]
        pslc[l] = float(lines[n].split('\t')[1])
        usst[l] = float(lines[n].split('\t')[2])
        vsst[l] = float(lines[n].split('\t')[3])
        cssf[l] = float(lines[n].split('\t')[4])
        clsf[l] = float(lines[n].split('\t')[5])
        ocis[l] = float(lines[n].split('\t')[6])
        oces[l] = float(lines[n].split('\t')[7])
        iswf[l] = float(lines[n].split('\t')[8])
        roce[l] = float(lines[n+1].split('\t')[0])
        olis[l] = float(lines[n+1].split('\t')[1])
        oles[l] = float(lines[n+1].split('\t')[2])
        role[l] = float(lines[n+1].split('\t')[3])
        swtc[l] = float(lines[n+1].split('\t')[4])
        ocic[l] = float(lines[n+1].split('\t')[5])
        lwtc[l] = float(lines[n+1].split('\t')[6])
        lwbc[l] = float(lines[n+1].split('\t')[7])
        for k in range(nz):
            nk = n + 2 + k
            pres[k] = float(lines[nk].split('\t')[0])
            temp[l][k] = float(lines[nk].split('\t')[1])
            umes[l][k] = float(lines[nk].split('\t')[2])
            liqm[l][k] = float(lines[nk].split('\t')[3])
            icem[l][k] = float(lines[nk].split('\t')[4])
            uvel[l][k] = float(lines[nk].split('\t')[5])
            vvel[l][k] = float(lines[nk].split('\t')[6])
            swrh[l][k] = float(lines[nk].split('\t')[7])
            lwrh[l][k] = float(lines[nk].split('\t')[8])
            #lwrh2[l][k] = float(lines[nk].split('\t')[8])

    if return_var == 'nz':
        return nz
    elif return_var == 'tstep':
        return tstep
    elif return_var == 'lon':
        return lon
    elif return_var == 'lat':
        return lat
    elif return_var == 'topo':
        return topo
    elif return_var == 'lmsk':
        return lmsk
    elif return_var == 'tstep':
        return tstep
    elif return_var == 'nrec':
        return nrec
    elif return_var == 'pres':
        return pres
    elif return_var == 'time':
        return time
    elif return_var == 'pslc':
        return pslc
    elif return_var == 'usst':
        return usst
    elif return_var == 'vsst':
        return vsst
    elif return_var == 'cssf':
        return cssf
    elif return_var == 'clsf':
        return clsf
    elif return_var == 'ocis':
        return ocis
    elif return_var == 'oces':
        return oces
    elif return_var == 'iswf':
        return oswf
    elif return_var == 'roce':
        return roce
    elif return_var == 'olis':
        return olis
    elif return_var == 'oles':
        return oles
    elif return_var == 'role':
        return role
    elif return_var == 'swtc':
        return swtc
    elif return_var == 'ocic':
        return ocic
    elif return_var == 'lwtc':
        return lwtc
    elif return_var == 'lwbc':
        return lwbc
    elif return_var == 'temp':
        return temp
    elif return_var == 'umes':
        return umes
    elif return_var == 'liqm':
        return liqm
    elif return_var == 'icem':
        return icem
    elif return_var == 'uvel':
        return uvel
    elif return_var == 'vvel':
        return vvel
    elif return_var == 'swrh':
        return swrh
    elif return_var == 'lwrh':
        return lwrh


#from pylab import contourf, show
from numpy import asarray, shape
#f_in = '/home/santiago/Modelos/sp-sam/SEMIPROG_IN'
#temp = map(list, map(None,*readsemiprogin(f_in, 'temp')))
#temp = readsemiprogin(f_in, 'temp')
#ttemp = map(None,*temp)
#nlev = readsemiprogin(f_in, 'nz')
#ntim = len(temp)
#print map(None,*my_list)

#k = range(nlev)
#l = range(ntim)
#print len(l), len(k)
#larr = asarray(l)
#karr = asarray(k)
#temparr = asarray(temp)
#print shape(temparr)
#print len(l), len(k)
#contourf(l,k,ttemp)
#show()


#from pylab import *
#xrange = range(100)
