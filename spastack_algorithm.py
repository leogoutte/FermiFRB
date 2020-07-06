## Temporal stacking analysis
# Written by: Leo Goutte
# Created: 21/07/2019
# Last updated: 23/07/2019

# -> To be used within the Fermi environment (Python 2.7.14) in order to analyse the behaviour of 
# potential radio-gamma correlations for FRBs
# -> LAT all-sky files can and should be used.
# -> This analysis can be conducted on any type of FRBs (repeaters or single-bursts). It is meant as
# an extension of the temporal stacking pipeline. Only verified bursts will be included.

import os 
import csv 
import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits 
import gt_apps as gt
from datetime import datetime
from astropy.time import Time

working_directory = os.getcwd()+"/"

print "\nSpatial stacking analysis"
print '\nThis tool is designed for use within the Fermi environment (Python 2.7.14) in order to analyse the behaviour of potential radio-gamma correlations for FRBs \n'
print "--->>> Before continuing, make sure that all LAT server files and the FRBcat csv are in the working directory. The directory should be called FRB###### <<<--- \n"

continue_prompt = raw_input("Begin? ([y]/n): ")
if continue_prompt == 'n':
    sys.exit()
else:
    pass

#Fixed parameters
print '\n>>> Fixed parameters \n\n**Inputs are case-sensistive and must match the information found in FRBcat (sourcename) and the LAT server**\n'

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

index_cutoff_prompt = raw_input("Burst index cutoff (int or none): ")
if index_cutoff_prompt == '' or 'none' or 'n':
    start_index = 0
else:
    start_index = int(index_cutoff_prompt)

repeater_prompt = raw_input("Include repeaters? ([y]/n): ")
if repeater_prompt == 'n' or 'no':
    repeat = 0
else: 
    repeat = 1

#dimensions of CCUBE -- nxpix*binsz < s = rad*sqrt2
binsz = np.around(rad*np.sqrt(2)/nxpix, decimals=2)

#create a whole lotta different lists to iterate processes

path = working_directory
catname = []
for item in os.listdir(path):
    if os.path.isfile(os.path.join(path,item)) and 'frbcat' in item:
        catname.append(item)

with open(catname[0], 'rb') as f:
    reader = csv.reader(f)
    the_list = list(reader)

sourcename = []
ra = []
dec = []
rawdates = []

if repeat == 0:
    for i in range(1,len(the_list)):
        if the_list[i][34] in ['1'] and the_list[i][5] in ['t']: #rank 1 and verified
            sourcename.append(the_list[i][0])
            ra.append(the_list[i][7])
            dec.append(the_list[i][8])
            rawdates.append(the_list[i][6])
elif repeat == 1:
    for i in range(1,len(the_list)):
        if the_list[i][5] in ['t']: #only verified
            sourcename.append(the_list[i][0])
            ra.append(the_list[i][7])
            dec.append(the_list[i][8])
            rawdates.append(the_list[i][6])

 
#how many?
n = len(sourcename)

#times
times = Time(rawdates, format='iso')
met_times = tburst = np.subtract(times.unix,978307200)

def hours_to_degrees_ra(dates):
    result = []
    for i in range(len(dates)):
        if len(dates[i]) > 7: #if there are seconds
            result.append(np.around((float(dates[i][:2]) + float(dates[i][3:5])/60 + float(dates[i][6:])/3600)*15, 3))
        elif len(dates[i]) < 7: #if it stops at minutes
            result.append(np.around((float(dates[i][:2]) + float(dates[i][3:])/60)*15, 3))
    
    return np.around(result,3)

def hours_to_degrees_dec(dates):
    result = []
    for i in range(len(dates)):
        if dates[i][0] in ['-']:
            if len(dates[i]) > 8:
                result.append(np.around((-(float(dates[i][1:3]) + float(dates[i][4:6])/60 + float(dates[i][7:])/3600)),3))
            elif len(dates[i]) < 8:
                result.append(np.around((-(float(dates[i][1:3]) + float(dates[i][4:])/60)),3))
        elif dates[i][0] in ['+']:
            if len(dates[i]) > 8:
                result.append(np.around((float(dates[i][1:3]) + float(dates[i][4:6])/60 + float(dates[i][7:])/3600),3))
            elif len(dates[i]) < 8:
                result.append(np.around((float(dates[i][1:3]) + float(dates[i][4:])/60), 3))
        else:
            if len(dates[i]) > 7:
                result.append(np.around((float(dates[i][:2]) + float(dates[i][3:5])/60 + float(dates[i][6:])/3600),3))
            elif len(dates[i]) < 7:
                result.append(np.around((float(dates[i][:2]) + float(dates[i][3:])/60),3))
        
    return np.around(result,3)

ra_deg = hours_to_degrees_ra(ra)
dec_deg = hours_to_degrees_dec(dec)

print "Ready for data preparation"


## Prepare

#put all files under one list
os.system("bash -c \"ls *v001_filtered* > filelist.txt \"")

#Regroup all photon files in one fits file
path = working_directory+'lat_alldata.fits'
if os.path.isfile(path):
    print 'Weekly files already assembled'
else:
    gt.filter['evclass'] = 'INDEF'
    gt.filter['evtype'] = 'INDEF'
    gt.filter['ra'] = 0
    gt.filter['dec'] = 0
    gt.filter['rad'] = 180
    gt.filter['emin'] = 30
    gt.filter['emax'] = 1000000
    gt.filter['zmax'] = 180
    gt.filter['tmin'] = 'INDEF'
    gt.filter['tmax'] = 'INDEF'
    gt.filter['infile'] = '@filelist.txt'
    gt.filter['outfile'] = 'lat_alldata.fits'
    gt.filter['chatter'] = 0
    gt.filter.run()
    print '**DONE**'

#filter
print '...Working on filtering...'

for i in range(start_index,n):
    path = sourcename[i]+'_filtered.fits'
    if os.path.isfile(path):
        print sourcename[i]+' already filtered'
    else:
        gt.filter['evclass'] = 128
        gt.filter['evtype'] = 3
        gt.filter['ra'] = ra_deg[i]
        gt.filter['dec'] = dec_deg[i]
        gt.filter['rad'] = radcut
        gt.filter['emin'] = emin
        gt.filter['emax'] = emax
        gt.filter['zmax'] = zmax
        gt.filter['tmin'] = tburst[i]-blow
        gt.filter['tmax'] = tburst[i]+blow
        gt.filter['infile'] = 'lat_alldata.fits'
        gt.filter['outfile'] = sourcename[i]+'_filtered.fits'
        gt.filter['chatter'] = 0
        gt.filter.run()
    
print 'Done filtering'

#make GTIs
print '...Working on GTIs...'

for i in range(start_index,n):
    path = sourcename[i]+'_filtered_gti.fits'
    if os.path.isfile(path):
        print sourcename[i]+' already GTId'
    else:
        gt.maketime['scfile'] = 'lat_spacecraft_merged.fits'
        gt.maketime['filter'] = '(DATA_QUAL>0) && (LAT_CONFIG==1)'
        gt.maketime['roicut'] = 'no'
        gt.maketime['evfile'] = sourcename[i]+'_filtered.fits'
        gt.maketime['outfile'] = sourcename[i]+'_filtered_gti.fits'
        gt.maketime['chatter'] = 0
        gt.maketime.run()

print 'Done GTIs'

#bin the data into energy cubes
print '...Working on CCUBE...'

for i in range(start_index,n):
    path = sourcename[i]+'_ccube.fits'
    if os.path.isfile(path):
        print sourcename[i]+' already cubed'
    else:
        gt.evtbin['algorithm'] = 'CCUBE'
        gt.evtbin['evfile'] = sourcename[i]+'_filtered_gti.fits'
        gt.evtbin['outfile'] = sourcename[i]+'_ccube.fits'
        gt.evtbin['scfile'] = 'NONE'
        gt.evtbin['nxpix'] = nxpix
        gt.evtbin['nypix'] = nxpix
        gt.evtbin['binsz'] = binsz
        gt.evtbin['coordsys'] = 'CEL'
        gt.evtbin['xref'] = ra_deg[i]
        gt.evtbin['yref'] = dec_deg[i]
        gt.evtbin['axisrot'] = 0
        gt.evtbin['proj'] = 'AIT'
        gt.evtbin['ebinalg'] = 'LOG'
        gt.evtbin['emin'] = emin
        gt.evtbin['emax'] = emax
        gt.evtbin['enumbins'] = 1 #sum over all energies
        gt.evtbin['chatter'] = 0
        gt.evtbin.run()

print 'Done CCUBE'

## Stack

#stacking
print '...Working on stacking and plot...'

summed_image_ = pyfits.getdata(sourcename[start_index]+'_ccube.fits', ext=0)
for i in range(start_index+1,n):
    summed_image_ += pyfits.getdata(sourcename[i]+'_ccube.fits', ext=0) 
summed_image = summed_image_[0,0:nxpix,0:nxpix]

#plot the resulting map
plt.figure(figsize=(8,8))
plt.title('Stacked CountCube of multiple FRBs summed over \n '+str(emin)+'-'+str(emax)+' MeV (blow up time of '+str(blow/3600/24)+' days)')
plt.imshow(summed_image, cmap='gnuplot2')
plt.xticks([])
plt.yticks([])
plt.axvline(nxpix/2, ls=':', c='chartreuse')
plt.axhline(nxpix/2, ls=':', c='chartreuse')
plt.xlabel(str(rad)+' degrees across')
plt.ylabel(str(rad)+' degrees across')
plt.colorbar(label='Counts')
plt.show()

print 'Done stackking and plot'
print 'Spatial stacking from '+rawdates[start_index]+' to '+rawdates[n]+' complete'

sys.exit()