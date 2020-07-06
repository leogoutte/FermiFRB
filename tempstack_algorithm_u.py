## Temporal stacking analysis
# Written by: Leo Goutte
# Created: 21/07/2019
# Last updated: 23/07/2019

# -> To be used within the Fermi environment (Python 2.7.14) in order to analyse the behaviour of 
# potential radio-gamma correlations for FRBs
# -> LAT files should be downloaded from the server before coninuing. Alternatively, all-sky 
# files can be used. This, however, will take much longer. 
# -> This analysis can only be conducted on FRBs with multiple detected bursts (repeaters). Only verified bursts will be included.

import csv
import os
import sys
import gt_apps as gt
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy.time import Time
working_directory = os.getcwd()+"/"

print "\nTemporal stacking analysis"
print '\nThis tool is designed for use within the Fermi environment (Python 2.7.14) in order to analyse the behaviour of potential radio-gamma correlations for FRBs \n'
print "--->>> Before continuing, make sure that all LAT server files and the FRBcat csv are in the working directory. The directory should be called FRB###### <<<--- \n"

continue_prompt = raw_input("Begin? ([y]/n): ")
if continue_prompt == 'n':
    sys.exit()
else:
    pass

#Fixed parameters
print '\n>>> Fixed parameters \n\n**Inputs are case-sensistive and must match the information found in FRBcat (sourcename) and the LAT server**\n'

sourcename_prompt = raw_input("Repeater's ID (# or blank): ")
if sourcename_prompt == '':
    sourcename = os.path.split(os.getcwd())[1]
else:
    sourcename = 'FRB'+sourcename_prompt

rad = float(raw_input("Search radius (deg) (0:180): "))

emin_prompt = raw_input("Lower energy limit (MeV): ")
if emin_prompt == '':
    emin = 30
else:
    emin = float(emin_prompt)

emax_prompt = raw_input("Upper energy limit (MeV): ")
if emax_prompt == '':
    emax = 300000
else:
    emax = float(emax_prompt)

tmin_prompt = raw_input("Start time (MET in s or START): ")
if tmin_prompt == 'START' or 'S' or 's' or '0' or '':
    tstart = 239557417
else:
    tstart = float(tmin_prompt)

tend = float(raw_input("End time (MET in s): "))

zmax_prompt = raw_input("Maximum zenith angle (deg) (0:180): ")
if zmax_prompt == '':
    zmax = 180
else:
    zmax = float(zmax_prompt)

#Variable parameters
print '\n>>> Variable parameters \n'

nxpix = int(raw_input("Image resolution (# of pixels on each axis): "))

blow_prompt = raw_input("Time radius around each burst (s or # of months): ")
if len(blow_prompt) < 3:
    blow = float(blow_prompt)*30*24*3600
else:
    blow = float(blow_prompt)

energy_cut_prompt = raw_input("Lower energy cut (MeV or none): ")
if energy_cut_prompt == 'none' or 'n' or '':
    energy_cut = emin
else: 
    energy_cut = float(energy_cut_prompt)

#dimensions of CCUBE -- nxpix*binsz < s = rad*sqrt2
binsz = np.around(rad*np.sqrt(2)/nxpix, decimals=2)

#open .csv with burst info

path = working_directory
catname = []
for item in os.listdir(path):
    if os.path.isfile(os.path.join(path,item)) and 'frbcat' in item:
        catname.append(item)

with open(catname[0], 'rb') as f:
    reader = csv.reader(f)
    the_list = list(reader)

#ra and dec
for i in range(len(the_list)):
    if the_list[i][0] == sourcename and the_list[i][34] in ['1']:
        ra_hours = the_list[i][7]
        dec_hours = the_list[i][8]

if len(ra_hours) > 7: #if there are seconds
    ra = np.around((float(ra_hours[:2]) + float(ra_hours[3:5])/60 + float(ra_hours[6:])/3600)*15, 3)
elif len(ra_hours) < 7: #if it stops at minutes
    ra = np.around((float(ra_hours[:2]) + float(ra_hours[3:])/60)*15, 3)
    

if dec_hours[0] in ['-']:
    if len(dec_hours) > 8:
        dec = np.around((-(float(dec_hours[1:3]) + float(dec_hours[4:6])/60 + float(dec_hours[7:])/3600)),3)
    elif len(dec_hours) < 8:
        dec = np.around((-(float(dec_hours[1:3]) + float(dec_hours[4:])/60)),3)
elif dec_hours[0] in ['+']:
    if len(dec_hours) > 8:
        dec = np.around((float(dec_hours[1:3]) + float(dec_hours[4:6])/60 + float(dec_hours[7:])/3600),3)
    elif len(dec_hours) < 8:
        dec = np.around((float(dec_hours[1:3]) + float(dec_hours[4:])/60), 3)
else:
    if len(dec_hours) > 7:
        dec = np.around((float(dec_hours[:2]) + float(dec_hours[3:5])/60 + float(dec_hours[6:])/3600),3)
    elif len(dec_hours) < 7:
        dec = np.around((float(dec_hours[:2]) + float(dec_hours[3:])/60),3)

#burst times
rawdates = []
for i in range(1,len(the_list)):
    if the_list[i][0] == sourcename and the_list[i][5] == 't':
        rawdates.append(the_list[i][6])
times = Time(rawdates, format='iso')
met_times = np.subtract(times.unix,978307200)
tburst = [0]
for i in range(len(met_times)): #we want the first time slot to be zero, this will make tburst more readable
    tburst.append(met_times[i]) # in the loops that follow
tburst.sort()
n = len(tburst)

#regroup the photon files under one name
os.system("bash -c \"cd "+working_directory+" && ls *_PH* > "+sourcename+"_events.txt\"")
os.system("bash -c \"mv "+working_directory+"*_SC* "+working_directory+sourcename+"_spacecraft.fits\"")

### Commence the actual stacking plot

#define function to be called later

def stack_plot(nxpixels,radius,zoom):
    #set max times
    tmax = [0]
    for i in range(1,n):
        tmax.append(tburst[i]+zoom)

    #min times
    tmin = [0]
    tmin.append(tburst[1]-zoom)
    for i in range(2,n):
        if tburst[i] - tburst[i-1] < 2*zoom:
            tmin.append(tmax[i-1])
        else: 
            for j in range(1,i-1):
                if tburst[i] - tburst[i-j-1] < 2*zoom:
                    tmin.append(tmax[i-j-1])    
            else: 
                tmin.append(tburst[i]-zoom)

    #Filter the data
    for i in range (1,n):
        gt.filter['evclass'] = 128
        gt.filter['evtype'] = 3
        gt.filter['ra'] = ra
        gt.filter['dec'] = dec
        gt.filter['rad'] = radius
        gt.filter['emin'] = emin
        gt.filter['emax'] = emax
        gt.filter['zmax'] = zmax
        gt.filter['tmin'] = tmin[i]
        gt.filter['tmax'] = tmax[i]
        gt.filter['infile'] = '@'+sourcename+'_events.txt'
        gt.filter['outfile'] = sourcename+'_burst'+str(i)+'_filtered.fits'
        gt.filter['chatter'] = 0
        gt.filter.run()

    print 'Filtered'

    #make GTIs
    for i in range(1,n):
        gt.maketime['scfile'] = sourcename+'_spacecraft.fits'
        gt.maketime['filter'] = '(DATA_QUAL>0) && (LAT_CONFIG==1)'
        gt.maketime['roicut'] = 'no'
        gt.maketime['evfile'] = sourcename+'_burst'+str(i)+'_filtered.fits'
        gt.maketime['outfile'] = sourcename+'_burst'+str(i)+'_filtered_gti.fits'
        gt.maketime['chatter'] = 0
        gt.maketime.run()

    print 'Good time intervals'

    ### View Data
    # Bin and CMAP

    #set binsize
    binsize = np.around(radius*np.sqrt(2)/nxpixels, decimals=2)

    #bin for ccube
    for i in range(1,n):
        gt.evtbin['algorithm'] = 'CCUBE'
        gt.evtbin['evfile'] = sourcename+'_burst'+str(i)+'_filtered_gti.fits'
        gt.evtbin['outfile'] = sourcename+'_burst'+str(i)+'_ccube.fits'
        gt.evtbin['scfile'] = 'NONE'
        gt.evtbin['nxpix'] = nxpixels
        gt.evtbin['nypix'] = nxpixels
        gt.evtbin['binsz'] = binsize
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

    print 'Data cubed'

    ### Stack the data
    # Extract the 0th extension (image) from each cmap and add them up

    #stacking
    summed_image_full = pyfits.getdata(sourcename+'_burst1_ccube.fits', ext=0)
    for i in range(2,n):
        summed_image_full += pyfits.getdata(sourcename+'_burst'+str(i)+'_ccube.fits', ext=0) 
    summed_image = summed_image_full[0,0:nxpix,0:nxpix]

    #plot the resulting map
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
    rect = patches.Rectangle((a[0],a[1]),abs(width),abs(width),linewidth=1.5,edgecolor='skyblue',facecolor='none',linestyle=':',label='Spatial uncertainty')
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
    print '\nPlotted \n\n>>> Close plot window to continue'
    plt.ion()
    plt.show(block=True)
    

    

#run the function
stack_plot(nxpix,rad,blow)

#prompt to re-investigate
keep_going_prompt = raw_input("Repeat with different parameters? (y/[n]: ")
if keep_going_prompt == 'y' or 'yes':
    for i in range(42):
        rad = float(raw_input("Search radius (deg) (0:180): "))
        nxpix = int(raw_input("Image resolution (# of pixels on each axis): "))
        blow_prompt = raw_input("Time radius around each burst (s or # of months): ")
        if len(blow_prompt) < 3:
            blow = float(blow_prompt)*30*24*3600
        else:
            blow = float(blow_prompt)
        stack_plot(nxpix,rad,blow)
        pass
        start_again = raw_input("Repeat with different parameters? (y/[n]): ")
        if start_again == 'y' or 'yes':
            pass
        else:
            break
            print '\n>>> Temporal stacking analysis complete'
            sys.exit()
else:
    print '\n>>> Temporal stacking analysis complete'
print '...Exiting program...'
sys.exit()