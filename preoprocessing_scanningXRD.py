"""

 2022 09 07

 Loading and save diffraction data (Merlin) NanoMAX MAy 2022


 -real space shifting implemented
 

 Susanna Hammarberg
 
"""

import numpy as np

import hdf5plugin #for 2020 data
import h5py
import matplotlib.pyplot as plt
import time
date_str = time.strftime("%Y%m%d_%H%M")


nanowire = 'nw1'
#savepath = 'F:/Susanna/NanoMAX_perovNW_May2022/'+nanowire+'/'
savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\diffraction/'

maskpath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\merlin_mask_2020/'
mask_file_name = 'merlin_mask_190222_14keV.h5'

#datapath = 'E:/Jesper_beamtimes/NanoMAX_May2022/20211260/2022051808/raw/sample/'
datapath = 'C:/Users/Sanna/NanoMAX_May2022_Perovskites_raw_selection/s387_411/'

#%%

scans = np.arange(387,411+1).tolist()#[387,411]#,397]#np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ])
scans = [397] 

# center on detector y and x is 286,233 (that is center of delta and gamma 0.)
x_cen = 233 # defined from Si calibration
#y_cen = 92 #defined not from Si calibration cos its too far off. using the center of uncoated NW instead. 
#choose roi on detetor

raw_slicey = slice(50,150)
raw_slicex = slice(150,300)

#raw_slicey = slice(20,180)
#raw_slicex = slice(150,300)
 
position_roi = np.arange(2000,3000,20)    #upp till 2000 väldigt lite från 3000 också lite

data = []
for rot_nbr, scan in enumerate(scans):

    with h5py.File(datapath + '000' + str(scan) +'.h5','r' ) as fp:
        print('loading scan ' + str(scan))

        #data.append(np.array(fp['entry/measurement/merlin/frames'][position_roi,raw_slicey,raw_slicex]))
        data.append(np.array(fp['entry/measurement/merlin/frames'][position_roi]))
        #rocking_curve.append(np.sum(data*mask_Merlin[frames_in_roi,raw_slicey,raw_slicex]))
        #gonphi = np.array(fp['entry/snapshot/gonphi'])
        #gonphis.append(gonphi)



#load mask
with h5py.File(maskpath  + mask_file_name,'r' ) as fp:
    mask_Merlin = np.array(fp['mask'])

#apply mask
#data = data * mask_Merlin[raw_slicey,raw_slicex]        
data = data * mask_Merlin        

#reshape data to [position,angle,det1,det2]
#data= np.swapaxes(data, 0,1)

#Save diffraction data
#np.save(savepath+ date_str, data)
#print(r'Saved diff data to path')
#print(savepath)
print('Shape of data is: ', np.array(data).shape)
plt.figure()
plt.imshow(np.log10(sum((data[rot_nbr]))),cmap='jet', interpolation='none')
plt.title('Summed intensity for all rotations (log)')
plt.colorbar()
#plt.savefig(savepath + 'diffraction')
plt.show()

#%%
#plot all frames in one scan
for pos_nbr, position in enumerate(position_roi.tolist()):
    #plot roi
    plt.figure()
    plt.imshow(np.log10(((data[0][pos_nbr]))),cmap='jet', interpolation='none')
    plt.title('1000 central frames in scan 398. Position %d (log) '%position)
    plt.colorbar()
    plt.savefig(savepath + 'diffraction_position%d'%position)



#%%
#plot all scans
for rot_nbr, scan in enumerate(scans):
    #plot roi
    plt.figure()
    plt.imshow(np.log10(sum((data[rot_nbr]))),cmap='jet', interpolation='none')
    plt.title('Summed intensity forrotation %d (log) '%rot_nbr)
    plt.colorbar()
    plt.savefig(savepath + 'diffraction_scan%d'%scan)
    plt.show()

