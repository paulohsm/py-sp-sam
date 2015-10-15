__author__ = 'santiago'


from py_sp_sam import *

#  References

#  Ploting guide
#  https://github.com/jbmouret/matplotlib_for_papers


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

#plot_prec_ara02()
#barplot_prec_trmm()
#sbs_prec()

trmm2file(trmm200401, lonfor, latfor)
trmm_time, trmm_prec = trmmfile2py()
ara2_prec = read_agcmdiagfields1d(agcm1d_ara02, 'prec')
kuo2_prec = read_agcmdiagfields1d(agcm1d_kuo02, 'prec')
ara2_pr3h = accum3(ara2_prec)
kuo2_pr3h = accum3(kuo2_prec)

from matplotlib.pyplot import plot, title, legend, show
plot(trmm_time, trmm_prec, color='black', linewidth=2, label='TRMM')
plot(trmm_time, ara2_pr3h, color='black', linewidth=2, label='MCGA-ARA')
plot(trmm_time, kuo2_pr3h, color='black', linewidth=2, label='MCGA-KUO')
title("Precipitacao Acumulada em 3h")
legend()
show()



#plot_trmm_prec()

quit()



precc1 = read_semiprogout(diagout_ara02, 'precc')
precc2 = read_semiprogout(diagout_kuo02, 'precc')
precc3 = read_semiprogout(diagout_zmc01, 'precc')

nt = read_semiprogout(diagout_ara02, 'nt')
time = range(nt)

from matplotlib.pyplot import plot, show
plot(time, precc1)
show()

quit()


#trmm2npy('/home/santiago/Datasets/TRMM-3B42RT/200401')


t0 = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT', 't0')
nt = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT', 'nt')
dt = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT', 'dt')
time = range(nt)

print t0, nt, dt

#time = get_semiprog_time('/home/santiago/Modelos/exps_spsam_2015-09-24/ARA02_032x001x4000x20s/SEMIPROG_OUT')
precc = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT', 'precc')
precl = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/KUO02_032x001x4000x20s/SEMIPROG_OUT', 'precl')

for tt in range(nt):
    print tt, precc[tt], precl[tt], precc[tt] + precl[tt]
from matplotlib.pyplot import plot, show

plot(time, precc)
plot(time, precl)
show()



quit()

spsam_precc = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/ARA02_032x001x4000x20s/SEMIPROG_OUT', 'precc')
spsam_precl = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/ARA02_032x001x4000x20s/SEMIPROG_OUT', 'precl')
spsam_prec = []
for idx in range(len(spsam_precc)):
    spsam_prec.append(spsam_precc[idx] + spsam_precl[idx])
spsam_prec = accum3(spsam_prec)
agcm_prec = accum3(read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-ARA02', 'prec'))
agcm_time = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-ARA02', 'time')
glob_time = []
for step in range(2,len(agcm_time),3):
    glob_time.append(agcm_time[step])
print glob_time
print len(spsam_prec), len(agcm_prec)

trmm_time, trmm_prec = load_trmm_series('/home/santiago/Datasets/TRMM-3B42RT/200401', 360.0-39.5, -4.0)
print len(spsam_prec), len(agcm_prec), len(trmm_prec), len(trmm_time[1:])
print glob_time[0], trmm_time[0]
print glob_time[-1], trmm_time[-1]
#for tt in range(len(trmm_time)):
#    print trmm_time[tt], agcm_time[tt*3]

from matplotlib.pyplot import plot, plot_date, title, show
#plot(trmm_time[1:], trmm_prec[1:]) #, agcm_prec, spsam_prec)
#plot(trmm_time[1:], agcm_prec)
plot(trmm_time[1:], spsam_prec)
show()

#agcm_time = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-GRE01', 'time')
#agcm_prec = read_agcmdiagfields1d('/home/santiago/Modelos/agcm-diagfields/AGCM_DIAGFIELDS_1D-GRE01', 'prec')
#spsam_time = read_semiprogout('/home/santiago/Modelos/exps_spsam_2015-09-24/ZMC01_032x001x4000x20s/SEMIPROG_OUT')



#print accum3(prec_agcm)
#print load_trmm_series('/home/santiago/Datasets/TRMM-3B42RT/200401', 360.0-39.5, -4.0)


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

