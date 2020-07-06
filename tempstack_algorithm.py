sourcename = 'FRB181017'
ra_hours = '17:05:00'
dec_hours = '68:17:00'
err_ra = 12 #arcminutes
err_dec = 12 #arcminutes
mjd = [58408,58530]
utc = ['23:26:11.8600','15:26:58.029']

blow = 1*30*24*3600
npix = 100
radcut = 3

import os
#assign working directory to where the data is present
working_directory = '/Volumes/Seagate/FRB/lightcurves/'+sourcename+'/'
os.chdir(working_directory)

import gt_apps as gt
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

#parameters
iso_spec = 'iso_P8R3_SOURCE_V2_v1'
ra = np.around((float(ra_hours[:2]) + float(ra_hours[3:5])/60 + float(ra_hours[6:])/3600)*15, 3)
dec = np.around((float(dec_hours[:2]) + float(dec_hours[3:5])/60 + float(dec_hours[6:])/3600),3)
tstart = 239557417
tend = 586266493
emin = 100
emax = 300000
zmax = 100

#burst times
mjd_met = np.subtract(mjd,51910+7.428703703703703e-4)
mjd_met_fixed = list(mjd_met)
mjd_met_fixed.insert(0,0)
utc_met = []
for i in range(len(utc)):
    utc_met.append(float(utc[i][:2])*3600 + float(utc[i][3:5])*60 + float(utc[i][6:]))
utc_met.insert(0,0)
tburst = []
for i in range(len(mjd_met_fixed)):
    tburst.append(mjd_met_fixed[i]*24*3600 + utc_met[i])
met_times = tburst
n = len(tburst)

#max
tmax = [0]
for i in range(1,n):
    tmax.append(tburst[i]+blow)

#dimensions of CCUBE -- npix*binsz < s = rad*sqrt2
binsz = np.around((radcut*np.sqrt(2)/npix), decimals=4)

#make sure there is no overlap

#create lists
tmin = [0]

#add first param
tmin.append(tburst[1]-blow)

#repeat
for i in range(2,n):
    if tburst[i] - tburst[i-1] < 2*blow:
        tmin.append(tmax[i-1])
    else: 
        for j in range(1,i-1):
            if tburst[i] - tburst[i-j-1] < 2*blow:
                tmin.append(tmax[i-j-1])    
        else: 
            tmin.append(tburst[i]-blow)

## 0. No Energy Cut

### Preliminaries
# Filter and GTIs

def regroup():
    os.system("bash -c \"cd "+working_directory+" && ls *_PH* > "+sourcename+"_events.txt\"")
regroup()
def screname():
    os.system("bash -c \"mv "+working_directory+"*_SC* "+working_directory+"spacecraft.fits\"")
screname()

#Filter the data
for i in range (1,n):
    gt.filter['evclass'] = 128
    gt.filter['evtype'] = 3
    gt.filter['ra'] = ra
    gt.filter['dec'] = dec
    gt.filter['rad'] = radcut
    gt.filter['emin'] = emin
    gt.filter['emax'] = emax
    gt.filter['zmax'] = zmax
    gt.filter['tmin'] = tmin[i]
    gt.filter['tmax'] = tmax[i]
    gt.filter['infile'] = '@'+sourcename+'_events.txt'
    gt.filter['outfile'] = sourcename+'_burst'+str(i)+'_filtered.fits'
    gt.filter['chatter'] = 0
    gt.filter.run()

print 'DONE'

#make GTIs
for i in range(1,n):
    gt.maketime['scfile'] = 'spacecraft.fits'
    gt.maketime['filter'] = '(DATA_QUAL>0) && (LAT_CONFIG==1)'
    gt.maketime['roicut'] = 'no'
    gt.maketime['evfile'] = sourcename+'_burst'+str(i)+'_filtered.fits'
    gt.maketime['outfile'] = sourcename+'_burst'+str(i)+'_filtered_gti.fits'
    gt.maketime['chatter'] = 0
    gt.maketime.run()

print 'DONE'

### View Data
# Bin and CMAP

#bin for ccube
for i in range(1,n):
    gt.evtbin['algorithm'] = 'CCUBE'
    gt.evtbin['evfile'] = sourcename+'_burst'+str(i)+'_filtered_gti.fits'
    gt.evtbin['outfile'] = sourcename+'_burst'+str(i)+'_ccube.fits'
    gt.evtbin['scfile'] = 'NONE'
    gt.evtbin['nxpix'] = npix
    gt.evtbin['nypix'] = npix
    gt.evtbin['binsz'] = binsz
    gt.evtbin['coordsys'] = 'CEL'
    gt.evtbin['xref'] = ra
    gt.evtbin['yref'] = dec
    gt.evtbin['axisrot'] = 0
    gt.evtbin['proj'] = 'AIT'
    gt.evtbin['ebinalg'] = 'LOG'
    gt.evtbin['emin'] = emin
    gt.evtbin['emax'] = emax
    gt.evtbin['enumbins'] = 1 #sum over all energies
    gt.evtbin['chatter'] = 0
    gt.evtbin.run()

print 'DONE'

### Stack the data
# Extract the 0th extension (image) from each cmap and add them up

#stacking
summed_image_ = pyfits.getdata(sourcename+'_burst1_ccube.fits', ext=0)
for i in range(2,n):
    summed_image_ += pyfits.getdata(sourcename+'_burst'+str(i)+'_ccube.fits', ext=0) 
summed_image = summed_image_[0,0:npix,0:npix]
summed_image = summed_image.astype(float)

from scipy.ndimage import gaussian_filter
result = gaussian_filter(summed_image, sigma=0.1*200/radcut, mode='constant', cval=np.mean(summed_image))

# with error range
err_ra = float(err_ra)/float(binsz)
err_dec = float(err_dec)/float(binsz)
a=[npix*0.5+(-err_ra/60),npix*0.5+(-err_dec/60)]
b=[npix*0.5+(+err_ra/60),npix*0.5+(-err_dec/60)]
c=[npix*0.5+(+err_ra/60),npix*0.5+(+err_dec/60)]
d=[npix*0.5*(-err_ra/60),npix*0.5+(+err_dec/60)]
width = float(a[0]) - float(b[0])
height = float(a[1]) - float(c[1])

import matplotlib.patches as patches

from scipy.ndimage import gaussian_filter
result = gaussian_filter(summed_image, sigma=0.1*200/radcut, mode='constant', cval=np.mean(summed_image))

ra_counts = np.sum(result,axis=0) #sum down the rows
dec_counts = np.sum(result,axis=1) #sum down the columns

# Set up the axes with gridspec
fig = plt.figure(figsize=(8, 8))
grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
main_ax = fig.add_subplot(grid[:-1, 1:])
y_hist = fig.add_subplot(grid[:-1, 0], sharey=main_ax)
x_hist = fig.add_subplot(grid[-1, 1:], sharex=main_ax)
t = np.linspace(0,100,100)

# scatter points on the main axes
im = main_ax.imshow(result)
main_ax.set_title('Stacked CountCube of '+sourcename+'\'s Bursts \n summed over '+str(emin)+'-'+str(emax)+' MeV (blow up time of '+str(blow/30/24/3600)+' month(s))')
main_ax.set_xticks([])
main_ax.set_yticks([])
main_ax.axvline(npix/2, ls=':', c='white')
main_ax.axhline(npix/2, ls=':', c='white')
# Create a Rectangle patch
rect = patches.Rectangle((a[0],a[1]),abs(width),abs(height),linewidth=1.5,edgecolor='skyblue',facecolor='none',linestyle=':',label='Spatial uncertainty')
# Add the patch to the Axes
main_ax.add_patch(rect)

# histogram on the attached axes
x_hist.plot(t,ra_counts, color='gray')
x_hist.invert_yaxis()
x_hist.set_title('RA ('+str(ra)+') [deg] - '+str(npix)+' pixels across', y=-0.3)

y_hist.plot(dec_counts, t, color='gray')
y_hist.invert_xaxis()
y_hist.invert_yaxis()
y_hist.set_ylabel('DEC ('+str(dec)+') [deg]')


cbaxes = fig.add_axes([0.93, 0.33, 0.03, 0.55]) 
cb = plt.colorbar(im, cax = cbaxes)

plt.show()
