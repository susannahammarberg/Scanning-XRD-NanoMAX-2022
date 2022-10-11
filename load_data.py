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



#savepath = r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_nparray\s387_411/'
savepath = r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_nparray\s288-308\\'

maskpath = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/'
mask_file_name = r'merlin_mask_200430_8keV.h5'

#datapath = r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_raw_selection\s387_411/'
datapath = r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_raw_selection\s288_308/'

#scans = np.arange(387,411+1).tolist()

scans = np.arange(288,308+1).tolist() # reference nw


det_ceny = 92 #defined not from Si calibration. (Merlin not fliped upsidedown). using the center of uncoated NW  (peak1) instead. 
det_cenx = 233 # defined from Si calibration

#choose roi on detetor
# all peaks
#TODO save with roi in y same as peak1
#raw_slicey = slice(50,150) #(defined for Merlin raw data not flipud)
#raw_slicex = slice(150,300)

# do roi around gamma=0.NW1  all peaks
#raw_slicey = slice(50,150) #(defined for Merlin raw data not flipud)
#raw_slicex = slice(233-100,233+100)


# NW1 peak 1
#raw_slicey = slice(20,160)
#raw_slicex = slice(150,208)

#ref NW
raw_slicey = slice(190-60,190+60)
raw_slicex = slice(175-40,175+40)

#define center of diffrarction roi
y_cen = raw_slicey.start + int((raw_slicey.stop-raw_slicey.start)/2) 
x_cen = raw_slicex.start + int((raw_slicex.stop-raw_slicex.start)/2) 


#TODO insted I should calulcate the offset from the center of the diffraction roi.
# the number that i calculate now is used for defining bragg_theta.

#offset of diffraction centers from hte caliubrated detector centers
#calculate offset in x and y: (in pixels)
#TODO save metadata as file with  save with np data. save metadata. roi, etc etc

#for all peaks analysis #if (0,0) is within roi
y_offset = -( y_cen - det_ceny ) #change sign because merlin is flipped upside down
x_offset = x_cen - det_cenx  

#for peak 1: if (0,0) is to the right of roi NO
#y_offset = -(det_ceny - y_cen ) #change sign because merlin is flipped upside down
#x_offset = (det_cenx - x_cen ) + (det_cenx-raw_slicex.stop)


# there are a lot of scanning thin air, dont load all of that, define a position roi
#position_roi = np.arange(1441,3275,1)    #upp till 2000 väldigt lite från 3000 också lite
# ref NW
#position_roi = np.arange(6000,10000,100) 

# insert the number of original measuring positions in x and y (original as in , be)
#NW1
#Nx_orig = 131
#Ny_orig = 51

#refNW
Nx_orig = 201
Ny_orig = 61


#shifting vectors compensating for real space sample drift
#vertical_shift =  np.loadtxt(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\xrf\y_shift_00387.txt',delimiter =',',dtype=int).tolist()
#horizontal_shift = np.loadtxt(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\xrf\x_shift_00387.txt',delimiter =',',dtype=int).tolist()

# reference NW                             
#horizontal_shift = [0,1,2,3,4,6,7,6,7,7,7,8,8,8,8,10,10,10,10,11,11] # nw 70 - reference
#vertical_shift = [0,4,4,3,3,3,5,6,6,6,9,9,9,8,10,12,13,17,19,20,21] # nw 70 - reference
#horizontal_shift =  [0,-1,-2,-3,-4,-6,-7,-6,-7,-7,-7,-8,-8,-8,-8,-10,-10,-10,-10,-11,-11] # nw 70 - reference
horizontal_shift =  [0,1,2,3,4,6,7,6,7,7,7,8,8,8,8,10,10,10,10,11,11] # nw 70 - reference
vertical_shift = [0,-4,-4,-3,-3,-3,-5,-6,-6,-6,-9,-9,-9,-8,-10,-12,-13,-17,-19,-20,-21] # nw 70 - reference


#%%


# the shape of STXM maps, after shifts are included, is the experimental nbr of positions minus the
# 'range' of the - to + values in the list of shift. 
Nx_new = Nx_orig - (np.max(horizontal_shift) - np.min(horizontal_shift))  
Ny_new = Ny_orig - (np.max(vertical_shift) - np.min(vertical_shift)) 
        
data = []
for scan_nbr, scan in enumerate(scans):
    
    #find the scanning positions that will be used for this scan
    roi_vertically = range(max(vertical_shift) - vertical_shift[scan_nbr], max(vertical_shift) - vertical_shift[scan_nbr] + Ny_new)
    roi_horizontally = range(min(horizontal_shift) + horizontal_shift[scan_nbr], min(horizontal_shift) + horizontal_shift[scan_nbr] + Nx_new)    

    #import pdb; pdb.set_trace()
    frames_in_roi = [] 
    
    #calculate which frames should be use to real space shift the sample
    #(probably dont need to loop)
    for row in roi_vertically:
        for col in roi_horizontally:                
            frames_in_roi.append(row * Nx_orig + col)
              
    #remove some unwanted rows above and below NW
    #TTODO this was not what i meant to do I meant to remove 25 (doesnt matter)
    #NW1
    #frames_in_roi = frames_in_roi[Nx_new*4:-(Nx_new*15)]
    #ref NW
    #does not work when you are aligning in x
    #frames_in_roi = frames_in_roi[Nx_new*5:]
    frames_in_roi = frames_in_roi[Nx_new*25:]
    
    with h5py.File(datapath + '000' + str(scan) +'.h5','r' ) as fp:
        print('loading scan ' + str(scan))
        #I0= np.array(fp['entry/measurement/alba/channel1'])
        data.append(np.array(fp['entry/measurement/merlin/frames'][frames_in_roi,raw_slicey,raw_slicex]))  


Ny_new = int(len(frames_in_roi)/Nx_new)        
#load mask
with h5py.File(maskpath  + mask_file_name,'r' ) as fp:
    mask_Merlin = np.array(fp['mask'])

#apply mask
#this takes some extra memory, can also mask when i load
data = data * mask_Merlin[raw_slicey,raw_slicex]        

#reshape data to [position,angle,det1,det2]
data = np.swapaxes(data, 0, 1)

#plt.figure(); plt.imshow(sum(sum(data)), cmap='jet')

#Flip the merlin images up-side-down
data = np.flip(data, axis = 2)


#%%
#Save diffraction data
np.save(savepath+ date_str, data)
print(r'Saved diff data to file')
print(savepath+date_str)
print('Shape of data is: ', data.shape)
with open(savepath+date_str+'_offset.txt','w') as f:
    f.write(str(y_offset) + ',' + str(x_offset))
print('saved offset')


#plot what you are saving
plt.figure()
plt.imshow(np.log10(sum(sum(data))),cmap='magma', interpolation='none')
plt.title('Summed intensity for all rotations (log)')
plt.colorbar()
plt.savefig(savepath + 'diffraction'+date_str)
plt.show()


rocking_curve = sum(np.sum(np.sum(data, axis=-1),axis=-1))

plt.figure()
plt.plot(scans,rocking_curve,'*-')
#plt.title('Rocking curve for all diffraction S(387-411)')            
plt.title('Rocking curve all diffraction S(387-411)')            
plt.show()


