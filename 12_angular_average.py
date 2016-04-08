import numpy as np
from nexusformat.nexus import *

import matplotlib as mpl
from matplotlib import pyplot as plt
import argparse
from nxs_handling import *
from angular_integration import *

# analysis tools from scikit-beam
# (https://github.com/scikit-beam/scikit-beam/tree/master/skbeam/core)
import skbeam.core.roi as roi
import skbeam.core.correlation as corr
import skbeam.core.utils as utils

#import xray_vision.mpl_plotting as mpl_plot

# -- parse arguments
parser = argparse.ArgumentParser(description='PLOT ANGULAR')
parser.add_argument('-f', '--data_name', type=str,default =
'fa027_03',
                    help='Which file you want foo?',nargs='+')
parser.add_argument('-T', '--temperature', type=int,
                    help='temperature',default=77,nargs='+')
parser.add_argument('-b', '--batch', type=int,
                    help='Which batch?',default=4,nargs='+')
parser.add_argument('-q', '--gq2', type=int,
                    help='At which q index to plot the g2?',default=40)
args = parser.parse_args()

# -- files and folders
base_dir = '/Users/fivos/Documents/Projects/03 XPCS of Amorhous ices/'
base_dir2 = '/Volumes/Seagate Backup Plus Drive/DESY_july2015/'
data_dir ='01_first_experiment/10_data_analysis/00_data/raw/'
data_dir2 = 'P10_XPCS_amorphous_ice/current/raw/'
mask_dir = '01_first_experiment/10_data_analysis/02_python/masks/'
#dataname= '%s_%d'%(np.array(args.data_name)[0],np.array(args.temperature)[0])#'fa027_03_77'#'fa027_03_77'#'HDA_01_80'
maskname= 'mask_fa027_05_300_l02.npy'

# --  parameters and constants
detector = 'l02'
snapshots_per_batch = 100
nx,ny = 516,1556
center_x,center_y = 276,723
batches = np.array(args.batch)
temperature = np.array(args.temperature)
filename = np.array(args.data_name)
gq2 = args.gq2

inner_radius = 40 # radius of the first ring
width = 2
spacing = 0
num_rings =75
center = (273, 723) # center of the speckle pattern

dpix = 0.055 # um (pixel size)
energy = 8.4 #keV
h = 4.135667516*1e-18 #kev*sec
c = 3*1e8 #m/s
lambda_ = h*c/energy*1e10 #1.5498 Angst # wavelength of the X-rays
Ldet = 5080. # detector to sample distance 

# -- pixels to Q range
edges = roi.ring_edges(inner_radius, width, spacing, num_rings)
two_theta = utils.radius_to_twotheta(Ldet, edges*dpix)
q_val = utils.twotheta_to_q(two_theta, lambda_)
q_ring = np.mean(q_val, axis=1)

#-- load mask
mask = np.load(base_dir+mask_dir+maskname)


# -- loop over filenames and temperatures
plt.figure()

for i_data in range(len(filename)):
    dataname='%s_%d'%(filename[i_data],temperature[i_data])
    # -- import data
    data = np.zeros((len(batches)*snapshots_per_batch,nx,ny))

    # -- color index for ploting
    color_index = float(i_data)/(len(filename)) 
    color = plt.cm.jet(color_index)

    for i_batch in range(len(batches)):
        start_number = batches[i_batch]*snapshots_per_batch+1
        end_number = (batches[i_batch]+1)*snapshots_per_batch 

        # -- load file
        file_path =format_filepath(base_dir2,data_dir2,dataname,detector,start_number,end_number)
        data[i_batch*snapshots_per_batch:(i_batch+1)*snapshots_per_batch,:,:] = np.array(load_nxs_data_series(file_path))

        # -- average and mask data
        #masked_data = data*mask
        avg_data = np.average(data[i_batch*snapshots_per_batch:(i_batch+1)*snapshots_per_batch,:,:],axis=0)#*mask

        # -- circular average
        max_x = num_rings*width+inner_radius
        Iq=roi.circular_average(avg_data*mask,calibrated_center=center,nx=num_rings,min_x=inner_radius,max_x=max_x)

        # -- Plot
        plt.plot(q_ring,Iq[1],label=dataname,linestyle='-',color=color,marker='o',alpha=0.5)

#plt.yscale('log',nonposy='clip')
plt.ylabel('Intensity (photons/pixel)')
plt.xlabel('Q (A^-1)')
plt.legend()
plt.xlim([0.0018,0.005])
plt.show()

'''
# -- rings array
inner_radius = 50 # radius of the first ring
width = 10
spacing = 0
num_rings = 36
center = (273, 723) # center of the speckle pattern
# width of each ring
# no spacing between rings
# number of rings
 #  find the edges of the required rings
edges = roi.ring_edges(inner_radius, width, spacing, num_rings)
#print edges

# -- Convert ring to q --FIX - FIX FIX!
dpix = 0.055 # The physical size of the pixels !! CROSSCHECK
energy = 8.4 #keV
h = 4.135667516*1e-18#kev*sec
c = 3*1e8 #m/s
lambda_ = h*c/energy*1e10#1.5498 # wavelength of the X-rays    !! CROSSCHECK
#print lambda_
Ldet = 5080. # # detector to sample distance 
two_theta = utils.radius_to_twotheta(Ldet, edges*dpix)
q_val = utils.twotheta_to_q(two_theta, lambda_)
q_ring = np.mean(q_val, axis=1)
#print q_ring
rings = roi.rings(edges, center, avg_data.shape)
 
ring_mask = rings*mask

plt.figure()
plt.subplot(2,1,1)
plt.imshow(ring_mask)
plt.subplot(2,1,2)
plt.imshow(np.log10(avg_data*mask+1e-3),vmin=-2,vmax=2)
# plot the figure
#fig, axes = plt.subplots()
#axes.set_title("Ring Mask")
#im = mpl_plot.show_label_array(axes, ring_mask, cmap="Dark2")
#plt.show()

num_levels = 1#7
num_bufs = 100#2
#g2, lag_steps = corr.multi_tau_auto_corr(num_levels, num_bufs, ring_mask, masked_data)
#print g2.shape


plt.figure(figsize=[15,5])
plt.subplot(1,2,1)
for i in range(len(g2[0,:])):
    plt.scatter(lag_steps,g2[:,i])
    plt.plot(lag_steps,g2[:,i],label='q=%.4fA-1'%(q_ring[i]))
plt.xscale('log',nonposy='clip')
plt.ylim([1,2])
plt.xlim([1,200])
#plt.legend()
plt.xlabel('dt (sec)')
plt.ylabel('g2(Q,dt)')

plt.subplot(1,2,2)
plt.plot(q_ring,g2[1,:],'ro-')


num_levels = 1#7
num_bufs = 100#2
num_frames=100
g2_t12=corr.two_time_corr(ring_mask,masked_data,num_frames,num_bufs,num_levels)[0]

plt.figure()
plt.imshow(g2_t12[0,:,:],vmin=1.1,vmax=g2[1,0],interpolation='none')
plt.colorbar()

plt.show()

'''

'''
# -- angular and radial maps
rho = radial_map(center_x,center_y,nx,ny)
phi = angular_map(center_x,center_y,nx,ny)

# -- radial ROI
ROI_radius_min,ROI_radius_max = 40,220
radial_ROI = take_ring(nx,ny,rho,ROI_radius_min,ROI_radius_max)
mask*= radial_ROI

# -- Appply mask
masked_data = avg_data*mask

# -- angular integration - MUST DEBUG
rx,I_r = angular_integrate(masked_data,mask,center_x,center_y)

# -- phi integration - MUST DEBUG
phix,I_phi = radial_integrate(masked_data,mask,center_x,center_y)



# -- calculate g2 -- MUST DEBUG
g2,dg2 = g2_calculate(data,mask,gq2,snapshots_per_batch,center_x,center_y)

dq = 2
g2_radius_min,g2_radius_max = gq2,gq2+dq
g2_ROI = take_ring(nx,ny,rho,g2_radius_min,g2_radius_max)*mask

num_levels = 7
num_bufs = 8
g2, lag_steps = corr.multi_tau_auto_corr(num_levels,num_bufs,g2_ROI,masked_data)#ring_mask,mask_data2)

# -- plots
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.imshow(np.log10(masked_data[center_x-ROI_radius_max:center_x+ROI_radius_max,center_y-ROI_radius_max:center_y+ROI_radius_max]+0.000001),vmax=2,vmin=-1,interpolation='none')
plt.colorbar()
plt.title('%s | %d'%(dataname,batches))

plt.subplot(2,2,2)
plt.plot(rx,I_r,'ro')
plt.yscale('log',nonposy='clip')
plt.xlabel('Radius (pixels)')
plt.ylabel('Intensity')

plt.subplot(2,2,4)
plt.plot(phix,I_phi,'gs-')
plt.xlabel('Angle Phi (degrees)')
plt.ylabel('Intensity')

plt.tight_layout()

plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.imshow(np.log10(masked_data[center_x-ROI_radius_max:center_x+ROI_radius_max,center_y-ROI_radius_max:center_y+ROI_radius_max]*g2_ROI[center_x-ROI_radius_max:center_x+ROI_radius_max,center_y-ROI_radius_max:center_y+ROI_radius_max]+0.000001),vmax=2,vmin=-1,interpolation='none')
plt.colorbar()
plt.title('%s | %d'%(dataname,batches))

plt.subplot(1,2,2)
g2x = np.arange(snapshots_per_batch)
plt.errorbar(lag_steps,g2,fmt='o')
plt.xscale('log',nonposy='clip')

plt.show()
'''
