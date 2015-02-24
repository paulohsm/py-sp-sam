__author__ = 'santiago'
__name__ = 'read_semiprog_in'

def readsemiprogin(input_file):
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


    #print lwrh.__len__()
    #return lwrh
    #print lwrh1
    #print lwrh2

#f_in = '/home/santiago/Modelos/sp-sam/SEMIPROG_IN'
#readsemiprogin(f_in)
#lwrh = readsemiprogin(f_in)
#print lwrh1
#print lwrh2


#nlines = 0
#for line in lines:
#    print line
#    nlines = nlines + 1

#print nlines


#for line in data:

#with open(f_in, 'r') as data:
#    line = data.readlines(1)[0].split('\t')

#nz = line[0]
#tstep = line[1]
#lon = line[2]
#lat = line[3]
#topo = line[4]
#lmsk = line[5]

#print line.

#nrec =

