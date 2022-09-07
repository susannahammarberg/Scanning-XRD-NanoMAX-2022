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


save_dir = r'C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_dataTEST_'

mask_directory = 'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/'
mask_file_name = 'merlin_mask_200430_8keV.h5'


data_directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'
scans = [429,430,431] #np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ])


#choose roi on detetor
raw_slicey = slice(0,270)
raw_slicex = slice(250,-1)

# choose the number of measuring positions in x and y
Nx_orig = 67
Ny_orig = 13


#shifting vectors compensating for real space sample drift
vertical_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\vertical_shift_vector.npy').tolist()
horizontal_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\horizontal_shift_vector.npy').tolist()                               

# the shape of STXM maps, after shifts are included, is the experimental nbr of positions minus the
# 'range' of the - to + values in the list of shift. 
Nx_new = Nx_orig - (np.max(horizontal_shift) - np.min(horizontal_shift))  
Ny_new = Ny_orig - (np.max(vertical_shift) - np.min(vertical_shift)) 
        
data = []
for scan_nbr, scan in enumerate(scans):
    
    #find the scanning positions that will be used for this scan
    roi_vertically = range(max(vertical_shift) - vertical_shift[scan_nbr], max(vertical_shift) - vertical_shift[scan_nbr] + Ny_new)
    roi_horizontally = range(min(horizontal_shift) + horizontal_shift[scan_nbr], min(horizontal_shift) + horizontal_shift[scan_nbr] + Nx_new)
    
    
    frames_in_roi = [] 
    
    #calculate which frames should be use to real space shift the sample
    #(probably dont need to loop)
    for row in roi_vertically:
        for col in roi_horizontally:                
                    frames_in_roi.append(row * Nx_orig + col)
                    
    
    with h5py.File(data_directory + '000' + str(scan) +'.h5','r' ) as fp:
        print('loading scan ' + str(scan))

        #I0= np.array(fp['entry/measurement/alba/channel1'])

        data.append(np.array(fp['entry/measurement/merlin/frames'][frames_in_roi,raw_slicey,raw_slicex]))
        
        #rocking_curve.append(np.sum(data*mask_Merlin[frames_in_roi,raw_slicey,raw_slicex]))
        #gonphi = np.array(fp['entry/snapshot/gonphi'])
        #gonphis.append(gonphi)

#load mask
with h5py.File(mask_directory  + mask_file_name,'r' ) as fp:
    mask_Merlin = np.array(fp['mask'])

#apply mask
data = data * mask_Merlin[raw_slicey,raw_slicex]        

#reshape data to [position,angle,det1,det2]
data= np.swapaxes(data, 0,1)

#Save diffraction data
np.save(save_dir+ date_str, data)
print(r'Saved diff data to C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_dataTEST_%s.np'%date_str)
print('Shape of data is: ', data.shape)

#plot what you are saving
plt.figure()
plt.imshow(np.log10(sum(sum(data))),cmap='magma', interpolation='none')
plt.title('Summed intensity for all rotations (log)')
plt.colorbar()
plt.savefig('ggg2')
plt.show()

