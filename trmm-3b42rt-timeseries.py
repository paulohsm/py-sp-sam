__author__ = 'santiago'

from pytrmm import TRMM3B42RTFile
from numpy import linspace, meshgrid, squeeze, arange, empty
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def plot_trmm(filepath):

    trmmdata =  TRMM3B42RTFile(filepath)
    prec = trmmdata.precip()
    prec_units = 'mm h-1'

    lons = linspace(-179.871, 179.875, prec.shape[1])
    lats = linspace( -59.875,  59.875, prec.shape[0])

    m = Basemap(projection='cyl', llcrnrlat=-60, llcrnrlon=-180, \
                urcrnrlat=60, urcrnrlon=180, resolution='c')

    lon, lat = meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    # Plot data
    cs = m.pcolor(xi, yi, squeeze(prec))

    #  Add grid lines
    parallels = arange(  -60.,  61., 20.)
    meridians = arange( -180., 181., 20.)
    m.drawparallels(parallels, labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=10)

    #  Add coastlines, states, and country boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()

    #  Add colorbar
    cbar = m.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label(prec_units)

    #  Add title
    plt.title('Precipitacao TRMM')

    plt.show()

# plot_trmm('/home/santiago/Datasets/TRMM-3B42RT/200401/3B42RT.2004010103.7R2.bin.gz')

#trmmdata =  TRMM3B42RTFile('/home/santiago/Datasets/TRMM-3B42RT/200401/3B42RT.2004010103.7R2.bin.gz')
#prec = trmmdata.precip()
#datadims = [prec.shape, 8*31]
trmmvals = empty([480, 1440, 248], dtype='float32')
counter = 0
for dday in range(1,32):
    for hour in range(0,24,3):
        filepath = '/home/santiago/Datasets/TRMM-3B42RT/200401/3B42RT.200401'+str(dday).zfill(2)+str(hour).zfill(2)+'.7R2.bin.gz'
        print str(dday).zfill(2)+str(hour).zfill(2), counter, filepath
        trmmvals[:,:,counter] = TRMM3B42RTFile(filepath).precip()
        counter = counter + 1

