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
savepath = r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_nparray\s387-411\\'

maskpath = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/'
mask_file_name = r'merlin_mask_200430_8keV.h5'

#datapath = 'E:/Jesper_beamtimes/NanoMAX_May2022/20211260/2022051808/raw/%s/' % sample
datapath = r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_raw_selection\s387_411/'
scans = np.arange(387,411+1).tolist()# , np.arange(491,503) , np.arange(465,491) ])

#choose roi on detetor
#raw_slicey = slice(0,270)
#raw_slicex = slice(250,-1)


#smaller but all peaks
raw_slicey = slice(50,150)
raw_slicex = slice(150,300)

#PEAK 1
#raw_slicey = slice(20,180)
#raw_slicex = slice(150,150+60)


# there are a lot of scanning thin air, dont load all of that, define a position roi
#position_roi = np.arange(1441,3275,1)    #upp till 2000 väldigt lite från 3000 också lite

# insert the number of original measuring positions in x and y (original as in , be)
Nx_orig = 131
Ny_orig = 51

#shifting vectors compensating for real space sample drift
vertical_shift =  np.loadtxt(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\xrf\y_shift_00387.txt',delimiter =',',dtype=int).tolist()
horizontal_shift = np.loadtxt(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\xrf\x_shift_00387.txt',delimiter =',',dtype=int).tolist()
#vertical_shift =  [0]*len(scans)#np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\vertical_shift_vector.npy').tolist()
#horizontal_shift = [0]*len(scans) #np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\horizontal_shift_vector.npy').tolist()                               

#vertical_shift =  np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\vertical_shift_vector.npy').tolist()
#horizontal_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\horizontal_shift_vector.npy').tolist()                               

# the shape of STXM maps, after shifts are included, is the experimental nbr of positions minus the
# 'range' of the - to + values in the list of shift. 
Nx_new = Nx_orig - (np.max(horizontal_shift) - np.min(horizontal_shift))  
Ny_new = Ny_orig - (np.max(vertical_shift) - np.min(vertical_shift)) 
        
data = []
for scan_nbr, scan in enumerate(scans):
    
    #find the scanning positions that will be used for this scan
    roi_vertically = range(max(vertical_shift) - vertical_shift[scan_nbr], max(vertical_shift) - vertical_shift[scan_nbr] + Ny_new)
    roi_horizontally = range(min(horizontal_shift) + horizontal_shift[scan_nbr], min(horizontal_shift) + horizontal_shift[scan_nbr] + Nx_new)    
    #roi_vertically = range(min(vertical_shift) - vertical_shift[scan_nbr], max(vertical_shift) - vertical_shift[scan_nbr] + Ny_new)
    #roi_horizontally = range(min(horizontal_shift) + horizontal_shift[scan_nbr], max(horizontal_shift) + horizontal_shift[scan_nbr] + Nx_new)       
    
    #import pdb; pdb.set_trace()
    frames_in_roi = [] 
    
    #calculate which frames should be use to real space shift the sample
    #(probably dont need to loop)
    for row in roi_vertically:
        for col in roi_horizontally:                
            frames_in_roi.append(row * Nx_orig + col)
              
    #remove some unwanted rows above and below NW
    #TTODO this was not what i meant to do I meant to remove 25 (doesnt matter)
    frames_in_roi = frames_in_roi[Nx_new*4:-(Nx_new*15)]
    #frames_in_roi = frames_in_roi[Nx_new*8:-(Nx_new*23)] #bra roi
    
    with h5py.File(datapath + '000' + str(scan) +'.h5','r' ) as fp:
        print('loading scan ' + str(scan))
        #I0= np.array(fp['entry/measurement/alba/channel1'])
        data.append(np.array(fp['entry/measurement/merlin/frames'][frames_in_roi,raw_slicey,raw_slicex]))    
        
#load mask
with h5py.File(maskpath  + mask_file_name,'r' ) as fp:
    mask_Merlin = np.array(fp['mask'])

#apply mask
#this takes some extra memory, can also mask when i load
data = data * mask_Merlin[raw_slicey,raw_slicex]        

#reshape data to [position,angle,det1,det2]
data= np.swapaxes(data, 0,1)

#%%TEMP REMOVE later
savepath_bf = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\bf\shifted\peak1'
dx = dy = 60e-9
nbr_rows = int(data.shape[0]/Nx_new)#51
nbr_cols = Nx_new#131#128#131

extent_microns = [ 0, dx*(nbr_cols-1)*1E6,0, dy*(nbr_rows-1)*1E6]


print('***BF analysis***')
# function for calculating bright field of data
# input is the data as a 3D array with dimension data[scans_y * scans_x][yPixels][xPixels]
def bright_field(data,y,x):
    index = 0
    photons = np.zeros((y,x)) 
    for row in range(0,y):
        for col in range(0,x):
            photons[row,col] = np.sum(data[index]) #/ max_intensity
            
            #fig, ax = plt.subplots(nrows=1)
            #ax.imshow(np.log10(photons),cmap = 'jet', interpolation='none', extent= [ 0, (x-1)*50E-3,(y-1)*50E-3,0 ])
            #plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]',title=r'Bright field #S440 Position %d'%index)
            ##plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\BF_maps\before_shifting\updated_20210609\S#440_pos%s'%("{0:03}".format(index))) 
            index = index+1
            #print(index)
    return photons

for scan in range(data.shape[1]):
    bf_map = bright_field((data[:,scan]),nbr_rows,nbr_cols)
    #bf_map = bright_field(np.sum((data[:,:]),axis=1),nbr_rows,nbr_cols)
    
    fig, ax = plt.subplots(nrows=1)
    ax.imshow(np.log10(bf_map),cmap = 'jet', interpolation='none', extent= extent_microns)
    plt.hlines(y=0.6, xmin=0, xmax=6.44)#, colors=', linestyles)
    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]')#,title=r'Bright field' 3d?
    plt.show()
    #plt.savefig(savepath_bf+ r'\rot%d'%scan)
#print(r'bf saved to ')
#print(savepath_bf)
#%%
#Save diffraction data
np.save(savepath+ date_str, data)
print(r'Saved diff data to file')
print(savepath+date_str)
print('Shape of data is: ', data.shape)

#plot what you are saving
plt.figure()
plt.imshow(np.log10(sum(sum(data))),cmap='magma', interpolation='none')
plt.title('Summed intensity for all rotations (log)')
plt.colorbar()
plt.savefig(savepath + 'ggg2')
plt.show()


rocking_curve = sum(np.sum(np.sum(data, axis=-1),axis=-1))

plt.figure()
plt.plot(np.arange(387,411+1),rocking_curve,'*-')
#plt.title('Rocking curve for all diffraction S(387-411)')            
plt.title('Rocking curve all diffraction S(387-411)')            
plt.show()
