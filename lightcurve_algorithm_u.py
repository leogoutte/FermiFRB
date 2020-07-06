## Light Curve Plot - Aperture Photometry Method
# Written by: Leo Goutte
# Created: 21/07/2019
# Last updated: 23/07/2019

# -> To be used within the Fermi environment (Python 2.7.14) in order to quickly plot 
# and investigate lightcurves in the gamma-ray spectrum.
# -> LAT files should be downloaded from the server before coninuing. Alternatively, 
# all-sky files can be used. This, however, will take much longer.

import os
import csv
import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import gt_apps as gt
import numpy as np
from GtApp import GtApp
from astropy.time import Time
from datetime import datetime
import seaborn as sns
sns.set()
working_directory = os.getcwd()+"/"

print '\nLight Curve Plot - Aperture Photometry Method \n'
print 'This tool is designed for use within the Fermi environment (Python 2.7.14) in order to quickly plot and investigate lightcurves \n'
print "--->>> Before continuing, make sure that all LAT server files and the FRBcat csv are in the working directory <<<--- \n"
continue_prompt = raw_input("Begin? ([y]/n): ")
if continue_prompt == 'n':
    sys.exit()
else:
    pass

#Fixed parameters
print '\n>>> Fixed parameters \n\n**Inputs are case-sensistive and must match the information found in FRBcat (sourcename) and the LAT server**\n'
sourcename_prompt = 'FRB'+raw_input("Sourcename's ID (# or blank): ")
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
if tmin_prompt == 'START' or 'S' or 's' or '0' or  '':
    tstart = 239557417
else: 
    tstart = float(tmin_prompt)

tend = float(raw_input("End time (MET in s): "))

zmax_prompt = raw_input("Maximum zenith angle (deg) (0:180): ")
if zmax_prompt == '':
    zmax = 300000
else:
    zmax = float(zmax_prompt)

#Variable parameters
print '\n>>> Variable parameters \n'
binsz = int(raw_input("Number of time bins: "))
bin_type = raw_input("Binning type (LIN/LOG): ")

#Others
trim = 3600
dtime = ((tend-tstart)-2*trim)/binsz #account for trimming
iso_spec = 'iso_P8R3_SOURCE_V2_v1'
index = '2.1' #if more accurate is known, use that instead

#open .csv with burst info
path = working_directory
catname = []
for item in os.listdir(path):
    if os.path.isfile(os.path.join(path,item)) and 'frbcat' in item:
        catname.append(item)

with open(catname[0], 'rb') as f:
    reader = csv.reader(f)
    the_list = list(reader)

ra_hours = []
dec_hours = []
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
    if the_list[i][0] == sourcename and the_list[i][5] in ['t']:
        rawdates.append(the_list[i][6])
utc_times = []
for i in range(len(rawdates)):
    utc_times.append(datetime.strptime(rawdates[i], '%Y-%m-%d  %H:%M:%S.%f').strftime('%Y-%m-%d %H:%M:%S.%f')[:-3])
times = Time(utc_times, format='iso')
met_times = np.subtract(times.unix,978307200)
tburst = [0]
for i in range(len(met_times)): #we want the first time slot to be zero, this will make tburst more readable
    tburst.append(met_times[i]) # in the loops that follow
tburst.sort()
n = len(tburst)

### Regroupements and Definitions

#define gtexposure
gtexposure = GtApp('gtexposure','Likelihood')
gtexposure.command()

#define gtbary
gtbary = GtApp('gtbary','Likelihood')
gtbary.command()

#execute ls *_PH* > events.txt and > spacecraft.fits
#BE SURE TO DELETE ALL OTHER *_PH* AND *_SC* FILES BEFORE CONTINUING
os.system("bash -c \"ls *_PH* > "+sourcename+"_events.txt\"")
os.system("bash -c \"mv "+working_directory+"*_SC* "+working_directory+sourcename+"_spacecraft.fits\"")

### Define the various functions to be called later
#full plot

def full_plot_function():

    #filter
    gt.filter['evclass'] = 128
    gt.filter['evtype'] = 3
    gt.filter['ra'] = ra
    gt.filter['dec'] = dec
    gt.filter['rad'] = rad
    gt.filter['emin'] = emin
    gt.filter['emax'] = emax
    gt.filter['zmax'] = zmax
    gt.filter['tmin'] = tstart + trim
    gt.filter['tmax'] = tend - trim
    gt.filter['infile'] = '@'+sourcename+'_events.txt'
    gt.filter['outfile'] = sourcename+'_filtered.fits'
    gt.filter['chatter'] = 0
    gt.filter.run()

    #maketime
    ##very time consuming
    gt.maketime['scfile'] = sourcename+'_spacecraft.fits'
    gt.maketime['filter'] = '(DATA_QUAL>0) && (LAT_CONFIG==1)'
    gt.maketime['roicut'] = 'no'
    gt.maketime['evfile'] = sourcename+'_filtered.fits'
    gt.maketime['outfile'] = sourcename+'_filtered_gti.fits'
    gt.maketime['chatter'] = 0
    gt.maketime.run()

    #bin
    gt.evtbin['algorithm'] = 'LC'
    gt.evtbin['evfile'] = sourcename+'_filtered_gti.fits'
    gt.evtbin['outfile'] = sourcename+'_lc_'+str(binsz)+'.fits'
    gt.evtbin['scfile'] = sourcename+'_spacecraft.fits'
    gt.evtbin['tbinalg'] = bin_type
    gt.evtbin['tstart'] = tstart + trim
    gt.evtbin['tstop'] = tend - trim
    gt.evtbin['dtime'] = dtime
    gt.evtbin['chatter'] = 0
    gt.evtbin.run()

    #make a copy of file to save in working directory under **_unaltered.fits
    import shutil
    shutil.copy2(working_directory+sourcename+'_lc_'+str(binsz)+'.fits', working_directory+sourcename+'_lc_'+str(binsz)+'_unaltered.fits')

    #run gtexposure
    shape_test = pyfits.getdata(sourcename+'_lc_'+str(binsz)+'.fits',ext=1)
    if len(shape_test[0]) == 5:
        print sourcename+'_lc_'+str(binsz)+'.fits already has exposure column'
        continue_test = raw_input("Do you want to compute exposure anyway? (y/[n]) ")
        if continue_test in ['y']:
            ##very time consuming
            gtexposure['infile'] = sourcename+'_lc_'+str(binsz)+'.fits'
            gtexposure['scfile'] = sourcename+'_spacecraft.fits'
            gtexposure['irfs'] = 'CALDB'
            gtexposure['srcmdl'] = 'none'
            gtexposure['specin'] = index
            gtexposure.run()

    #run gtbary
    #optional step - only required if time bins are small in order to investigate short timescale variability
    gtbary['evfile'] = sourcename+'_lc_'+str(binsz)+'.fits'
    gtbary['outfile'] = sourcename+'_lc_'+str(binsz)+'_bary.fits'
    gtbary['scfile'] = sourcename+'_spacecraft.fits'
    gtbary['ra'] = ra
    gtbary['dec'] = dec
    gtbary['tcorrect'] = 'BARY'
    gtbary.run()

    ### Plot the data

    #take the 1st ext (rate)
    data = pyfits.getdata(sourcename+'_lc_'+str(binsz)+'.fits', ext=1) 

    #extract time and counts
    time = [item[0] for item in data[0:np.size(data)-1]] #get rid of last point
    timedel = [item[1] for item in data[0:np.size(data)-1]]
    counts = [item[2] for item in data[0:np.size(data)-1]]
    error = [item[3] for item in data[0:np.size(data)-1]]
    exposure = [item[4] for item in data[0:np.size(data)-1]]

    #create flux (photon/cm^2)
    np.seterr(divide='ignore', invalid='ignore')
    flux = np.divide(counts,exposure)
    fluxerror = np.divide(error,exposure)

    plt.figure(figsize=(32,8))
    plt.title(sourcename+' Light Curve (in '+str(binsz)+' bins)')
    plt.ylabel("Photon Flux [r'$10^{-6} \cdot \frac{photon}{cm^{2} \cdot s}$]'")
    plt.xlabel('Time [s]')
    plt.errorbar(time,10**6*flux,yerr=10**6*fluxerror,xerr=timedel,ls='')
    plt.scatter(time,10**6*flux,s=40,label='Counts')
    for i in range(1,len(tburst)):
        if tburst[i]>time[0] and tburst[i]<time[-1] :
            plt.axvline(tburst[i], ls=':', c='green', label='Burst'+str(i))
    plt.legend(loc='best')
    plt.show()

#blow up plot

def blow_up_plot(zoom_centre,zoom,binsize,typeofbins,elow):

    tminblow = float(zoom_centre)-float(zoom)
    tmaxblow = float(zoom_centre)+float(zoom)
    dtimeblow = (tmaxblow-tminblow)/int(binsize)
    days = np.around(dtimeblow/60/60/24, 1)

    #filter
    gt.filter['evclass'] = 128
    gt.filter['evtype'] = 3
    gt.filter['ra'] = ra
    gt.filter['dec'] = dec
    gt.filter['rad'] = rad
    gt.filter['emin'] = elow
    gt.filter['emax'] = emax
    gt.filter['zmax'] = zmax
    gt.filter['tmin'] = tminblow
    gt.filter['tmax'] = tmaxblow
    gt.filter['infile'] = sourcename+'_filtered.fits'
    gt.filter['outfile'] = sourcename+'_blow_filtered.fits'
    gt.filter['chatter'] = 0
    gt.filter.run()

    #maketime
    gt.maketime['scfile'] = sourcename+'_spacecraft.fits'
    gt.maketime['filter'] = '(DATA_QUAL>0) && (LAT_CONFIG==1)'
    gt.maketime['roicut'] = 'no'
    gt.maketime['evfile'] = sourcename+'_blow_filtered.fits'
    gt.maketime['outfile'] = sourcename+'_blow_filtered_gti.fits'
    gt.maketime['chatter'] = 0
    gt.maketime.run()

    #bin
    gt.evtbin['algorithm'] = 'LC'
    gt.evtbin['evfile'] = sourcename+'_blow_filtered_gti.fits'
    gt.evtbin['outfile'] = sourcename+'_blow_lc.fits'
    gt.evtbin['scfile'] = sourcename+'_spacecraft.fits'
    gt.evtbin['tbinalg'] = typeofbins
    gt.evtbin['tstart'] = tminblow
    gt.evtbin['tstop'] = tmaxblow
    gt.evtbin['dtime'] = dtimeblow
    gt.evtbin['chatter'] = 0
    gt.evtbin.run()

    #make a copy of file to save in working directory under **_unaltered.fits
    import shutil
    shutil.copy2(working_directory+sourcename+'_blow_lc.fits', working_directory+sourcename+'_blow_lc_unaltered.fits')

    #run gtexposure
    ##very time consuming
    gtexposure['infile'] = sourcename+'_blow_lc_.fits'
    gtexposure['scfile'] = sourcename+'_spacecraft.fits'
    gtexposure['irfs'] = 'CALDB'
    gtexposure['srcmdl'] = 'none'
    gtexposure['specin'] = index
    gtexposure['chatter'] = 0
    gtexposure.run()

    #run gtbary
    #optional step - only required if time bins are small in order to investigate short timescale variability
    gtbary['evfile'] = sourcename+'_blow_lc_.fits'
    gtbary['outfile'] = sourcename+'_blow_lc_bary.fits'
    gtbary['scfile'] = sourcename+'_spacecraft.fits'
    gtbary['ra'] = ra
    gtbary['dec'] = dec
    gtbary['tcorrect'] = 'BARY'
    gtbary['chatter'] = 0
    gtbary.run()

    datablow = pyfits.getdata(sourcename+'_blow_lc_bary.fits', ext=1)

    #extract time and counts
    time_blow = [item[0] for item in datablow[0:np.size(datablow)-1]] #get rid of last point
    timedel_blow = [item[1] for item in datablow[0:np.size(datablow)-1]]
    counts_blow = [item[2] for item in datablow[0:np.size(datablow)-1]]
    error_blow = [item[3] for item in datablow[0:np.size(datablow)-1]]
    exposure_blow = [item[4] for item in datablow[0:np.size(datablow)-1]]

    #create flux (photon/cm^2)
    np.seterr(divide='ignore', invalid='ignore')
    flux_blow = np.divide(counts_blow,exposure_blow)
    fluxerror_blow = np.divide(error_blow,exposure_blow)

    plt.figure(figsize=(32,8))
    plt.title(sourcename+' Light Curve ('+str(days)+' day bins)')
    plt.ylabel("Photon Flux [r'$10^{-6} \cdot \frac{photon}{cm^{2} \cdot s}$]'")
    plt.xlabel('Time [s]')
    plt.errorbar(time_blow,10**6*flux_blow,yerr=10**6*fluxerror_blow,xerr=timedel_blow,ls='')
    plt.scatter(time_blow,10**6*flux_blow,s=40,label='Counts')
    for i in range(1,len(tburst)):
        if tburst[i]>time_blow[0] and tburst[i]<time_blow[-1] :
            plt.axvline(tburst[i], ls=':', c='green', label='Burst'+str(i))
    plt.legend(loc=1)
    plt.show()

#proceed with program

#full plot plot
print '\n>>> Start-End plot \n'
full_plot = raw_input(">>> Do you wish to make a start-end plot? (y/[n]): ")
if full_plot == 'y' or full_plot == 'yes':
    full_plot_function()
    continue_past_full_plot = raw_input("Do you wish to investigate a new search centre? ([y]/n): ")
    if continue_past_full_plot == 'y' or continue_past_full_plot == 'yes':
        pass
    else:
        sys.exit()
else:
    pass

#blow up plot plot
print "\n>>> Investigation of new search centres \n"
for i in range(42):
    blow_centre_prompt = raw_input("New search centre (s or burst #): ")
    if len(blow_centre_prompt) < 5:
        blow_centre = tburst[int(blow_centre_prompt)]
    else:
        blow_centre = float(blow_centre_prompt)
    
    blow_prompt = raw_input("Radius of new search centre (s or # of months): ")
    if len(blow_prompt) < 3:
        blow = float(blow_prompt)*30*24*3600
    else:
        blow = float(blow_prompt)
    
    binsz = raw_input("Number of time bins: ")
    
    bin_type = raw_input("Binning type (LIN/LOG): ")
    
    energy_cut_prompt = raw_input("Lower energy cut (MeV or none): ")
    if energy_cut_prompt == 'none' or 'n' or '':
        energy_cut = emin
    else: 
        energy_cut = float(energy_cut_prompt)
    
    blow_up_plot(blow_centre,blow,binsz,bin_type,energy_cut)
    
    start_again = raw_input("Investigate another search centre? (y/[n]): ")
    if start_again == 'y' or 'yes':
        pass
    else:
        break

print "\n>>> Light Curve Plots complete"

sys.exit()