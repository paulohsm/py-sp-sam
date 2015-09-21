__author__ = 'santiago'


from pytrmm import TRMM3B42RTFile
from numpy import linspace, meshgrid, squeeze, arange, empty, abs, where, array, save
from matplotlib.pyplot import plot, plot_date, title, show
from matplotlib.dates import date2num
from mpl_toolkits.basemap import Basemap
from datetime import datetime


path = '/home/santiago/Datasets/TRMM-3B42RT/200401'
pref = '3B42RT.200401'
suff = '.7R2.bin.gz'


def get_num(x):
    # returns numbers contained in strings
    return float(''.join(ele for ele in x if ele.isdigit() or ele=='.'))


def load_trmm_series():
    # read a series of TRMM data files, returns a matrix and several metadata
    global nlon, nlat, ntim
    global lon0, lat0, lon1, lat1
    global time_stamps, trmm_series
    # loading the data
    time_stamps = []
    trmm_series = []
    print 'Loading a series of TRMM matrices... please wait.'
    for dday in range(1,32):
        for hour in range(0,24,3):
            filepath = path + '/' + pref + str(dday).zfill(2) + str(hour).zfill(2) + suff
            trmmfile = TRMM3B42RTFile(filepath)
            time_string = trmmfile.header()['granule_ID'].split(".")[1]
            time_stamps.append(datetime.strptime(time_string, "%Y%m%d%H"))
            trmm_series.append(trmmfile.precip())
    print '...created "trmm_series" 3D object containing the data: use it.'
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


def get_coordinates():
    #  creates zonal and meridional axes
    print 'Creating axes for longitude and latitude.'
    lons = linspace(lon0, lon1, nlon)
    lats = linspace(lat0, lat1, nlat)
    return lons, lats


def find_nearest_idx(array, value):
    # given a single value, returns the index os the nearest value in an array/list
    idx = (abs(array-value)).argmin()
    print 'Nearest grid point and index:', array[idx], idx
    return idx #return array[idx]


def plot_trmm(lon, lat):
    trmm_point = []
    lons, lats = get_coordinates()
    lon_idx = find_nearest_idx(lons, lon)
    lat_idx = find_nearest_idx(lats, lat)
    for it in range(ntim):
        trmm_point.append(trmm_series[it][lat_idx, lon_idx])
    print 'Ploting the data...'
    plot(time_stamps, trmm_point)
    show()


def fill_trmm(precip, time_idx):
    field = precip[time_idx]
    prec_units = 'mm h-1'
    lons, lats = get_coordinates()

    print 'Creating a map.'
    shade = Basemap(projection='cyl', llcrnrlat=-60, llcrnrlon=0, \
                    urcrnrlat=60, urcrnrlon=360, resolution='c')

    lon, lat = meshgrid(lons, lats)
    xi, yi = shade(lon, lat)

    colorscale = shade.pcolor(xi, yi, squeeze(field))

    parallels = arange(  -60.,  61., 20.)
    meridians = arange( -180., 181., 20.)
    shade.drawparallels(parallels, labels=[1,0,0,0], fontsize=10)
    shade.drawmeridians(meridians, labels=[0,0,0,1], fontsize=10)

    shade.drawcoastlines()
    shade.drawstates()
    shade.drawcountries()

    cbar = shade.colorbar(colorscale, location='bottom', pad='10%')
    cbar.set_label(prec_units)

    title('TRMM 3B42RT Precip. ' + str(time_stamps[time_idx]))
    show()


#  testing...
load_trmm_series()
flon = 360.0-39.5
flat = -4.0
plot_trmm(flon, flat)
arm9707lon = 360.0-262.5
arm9707lat = 36.5
#plot_trmm(arm9707lon, arm9707lat)
#fill_trmm(trmm_series, 30)








def get_coord_idx(lon, lat):
    # gives the indexes of a lon,lat pair
    lons, lats = get_coordinates()
    nearest_lon = find_nearest(lons, lon)
    nearest_lat = find_nearest(lats, lat)
    lon_idx = where(lons==nearest_lon)
    lat_idx = where(lats==nearest_lat)
    return lon_idx, lat_idx
    #return lat_idx, lon_idx


#def plot_trmm(trmmvals, lon, lat):
def plot_trmm_OLD2(lon, lat):
    # plot a time series given the 3D data matrix and the coordinate to plot
    #trmm_point = empty([ntim], dtype='float')
    idxlon, idxlat = get_coord_idx(lon, lat)
    print idxlon, idxlat
    plot_date(time_stamps, trmm_series[:][idxlon, idxlat])
    #for it in range(ntim):
        #trmm_point(it) = trmm_series[it](idxlat, idxlon)
        #print trmm_series[it][idxlon, idxlat]
    #plot_date(date2num(time_stamps), trmm_point)

    #plot_date(time_stamps, trmmvals(:)[idxlat][idxlon])
    #trmm_series = array(trmmvals[idxlat,idxlon,:]).reshape([ntim])
    #plot_date(time_axis, trmm_series)
    show()


#arm9707lon = 360.0-262.5
#arm9707lat = 36.5
#flon = -39.5
#flat = -4.0
#load_trmm_series()
#plot_trmm(flon, flat)
#print len(trmm_series)
#print type(trmm_series[0])
#print trmm_series[0].shape
#trmm_series0 = trmm_series[0]
#flon, flat = get_coord_idx(360.0-39.5, -4.0)
#print flon, flat
#print trmm_series0[flat, flon]


#plot_trmm(flon, flat)
#plot_date(time_stamps, trmm_series[:](flat,flon))
#prec = load_trmm_series()
#plot_trmm(trmm_series(), arm9707lon, arm9707lat)
#plot_trmm(prec, flon, flat)




def shade_trmm(field):
    prec_units = 'mm h-1'
#   lons, lats = get_coordinates()
    #shade = Basemap(projection='cyl', llcrnrlat=-60, llcrnrlon=-180, \
    #                urcrnrlat=60, urcrnrlon=180, resolution='c')
    shade = Basemap(projection='cyl', llcrnrlat=-60, llcrnrlon=-180, \
                    urcrnrlat=60, urcrnrlon=180, resolution='c')
    lon, lat = meshgrid(lons, lats)
    xi, yi = shade(lon, lat)

    colorscale = shade.pcolor(xi, yi, squeeze(field))

    parallels = arange(  -60.,  61., 20.)
    meridians = arange( -180., 181., 20.)
    shade.drawparallels(parallels, labels=[1,0,0,0], fontsize=10)
    shade.drawmeridians(meridians, labels=[0,0,0,1], fontsize=10)

    shade.drawcoastlines()
    shade.drawstates()
    shade.drawcountries()

    cbar = shade.colorbar(colorscale, location='bottom', pad='10%')
    cbar.set_label(prec_units)

    title('Precipitacao TRMM')
    show()

#print trmm_series[0]
#shade_trmm(trmm_series[0])

#prec = load_trmm_series()
#for it in range(248):
#    shade_trmm(prec[:,:,it])
#save(path+'/TRMM-3B42RT-200401
# .npy', prec)

#arm9707lon = 360.0-262.5
#arm9707lat = 36.5
#prec = load_trmm_series()
#plot_trmm(trmm_mat, arm9707lon, arm9707lat)
#flon = -39.5
#flat =  -4.0
#plot_trmm(prec, flon, flat)
#shade_trmm(prec[:,:,25])






###################################
## old code snippets, tests, etc ##
###################################
def plot_trmm_OLD(filepath):

    trmmdata =  TRMM3B42RTFile(filepath)
    prec = trmmdata.precip()
    prec_units = 'mm h-1'

    lons, lats = get_coordinates()

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
def load_trmm_series_OLD():
    trmmvals = empty([480, 1440, 248], dtype='float32')
    counter = 0
    for dday in range(1,32):
        for hour in range(0,24,3):
            filepath = '/home/santiago/Datasets/TRMM-3B42RT/200401/3B42RT.200401'+str(dday).zfill(2)+str(hour).zfill(2)+'.7R2.bin.gz'
            print str(dday).zfill(2)+str(hour).zfill(2), counter, filepath
            trmmvals[:,:,counter] = TRMM3B42RTFile(filepath).precip()
            #dates =
            counter = counter + 1

            return trmmvals
#           plt.plot(trmmvals[300,500,:])
#           plt.show()


def get_coordinates_OLD():
    lons = linspace(-179.875, 179.875, 1440)
    lats = linspace( -59.875,  59.875, 480)
    return lons, lats


def load_trmm_series_OLD():
    trmmvals = empty([nlat, nlon, ntim], dtype='float32')
    counter = 0
    for dday in range(1,32):
        for hour in range(0,24,3):
            filepath = path + '/' + pref + str(dday).zfill(2) + str(hour).zfill(2) + suff
            trmmvals[:,:,counter] = TRMM3B42RTFile(filepath).precip()
            counter = counter + 1
    return trmmvals