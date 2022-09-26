"""
 
 Scanning XRD on nanomax 2020 data

 Susanna Hammarberg

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import DivergingNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
date_str = time.strftime("%Y%m%d_%H%M")

print('***Load data np file***')
      
#import matplotlib
#matplotlib.use( 'Qt5agg' )


savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\scanningXRD\unshifted/'
savepath_bf = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_Lucas_PerovNW_may2022\analysis\NW1\bf\shifted'
array = '20220920_1211_allpeaks.npy' # remove ths'20220916_1419.npy' #'20220916_1140.npy'
# load diffraction data
diff_data = np.load(r'C:\Users\Sanna\NanoMAX_May2022_Perovskites_nparray\s387-411/'+array)

#data should be saved as [position,angle,det1,det2]

#extra cropping
#diff_data = diff_data[:,:,0:220]

#plot raw data
plt.figure()
plt.imshow(np.log10(sum(sum(diff_data))),cmap='magma', interpolation='none')
plt.title('All data')
#plt.savefig('all_diffraction_%s'%date_str)
plt.colorbar()
plt.show()



fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True) # figsize=(10,3)
ax[0].imshow(np.log10(sum((diff_data[:,0]))),cmap='jet', interpolation='none')
ax[1].imshow(np.log10(sum((diff_data[:,-1]))),cmap='jet', interpolation='none')
ax[0].set_ylabel('q1 qz ')
ax[0].set_xlabel('q2 qy ')
ax[1].set_xlabel('q2 qy ')
plt.suptitle('First index S#489 (left) and last index S#429 (right) of data array (along qx) ')

#             (qx, qz, qy), or (q3, q1, q2).

rocking_curve = sum(np.sum(np.sum(diff_data, axis=-1),axis=-1))

plt.figure()
plt.plot(rocking_curve)
plt.title('Rocking curve for all diffraction S(387-411)')            
plt.show()
#plt.title('Data is sorted with the highest gonphi first \n (= the highest scan number first) ((why??)')
#NOT NOW, right


## reverse the order of the angles
#diff_data = diff_data[:,::-1]
#
#
#plt.figure(); plt.imshow((((diff_data[100,0]))),cmap='jet', interpolation='none'); plt.title('position 100'); plt.colorbar()
#
#
#rocking_curve = sum(np.sum(np.sum(diff_data, axis=-1),axis=-1))
#
#plt.figure()
#plt.plot(rocking_curve)            
#plt.title('After reversed axis)')
#plt.show()


print('***Define geometry of experiment***')

#define geoemtry class to calculate q1 q2 q3
class Geometry:
  def __init__(self, psize, shape,energy,distance,theta_bragg):
      #Rocking curve step (in degrees) and pixel sizes (in meters)
      #First element is the rocking curve step.
      self.psize = psize
      self.shape = shape                #Nbr of rockingcurve steps and nbr of pixels used on the detector
      self.energy = energy              #energy in units of keV
      self.distance = distance          #sample-bragg detector distance
      self.theta_bragg = theta_bragg    #Bragg_angle (the center of the bragg peak, not gonphi but the calibrated theta (the NW may not be laying flat on the substrate))
      
      # from these inputs i calculate:
      self.wavelength = 1.23984E-9 / energy #wavelength

      self.dq1 = self.psize[1] * 2 * np.pi / self.distance / self.wavelength
      self.dq2 = self.psize[2] * 2 * np.pi / self.distance / self.wavelength
      self.dq3 = np.deg2rad(self.psize[0]) * 4 * np.pi / self.wavelength * self.sintheta()
            
  def sintheta(self):
      return np.sin(np.deg2rad(self.theta_bragg))
    
  def costheta(self):
      return np.cos(np.deg2rad(self.theta_bragg))
    
  def tantheta(self):
      return np.tan(np.deg2rad(self.theta_bragg))

#%% 
##### Input parameters
   
g = Geometry((0.1, 55*1E-6, 55*1E-6), shape=(diff_data.shape[1], diff_data.shape[2], diff_data.shape[3]), energy=15.0, distance=0.30495, theta_bragg=7.9)
dx = dy = 60e-9
nbr_rows = 24#51
nbr_cols = 131#131#128#131

extent_microns = [ 0, dx*(nbr_cols-1)*1E6,0, dy*(nbr_rows-1)*1E6]


#%%
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

#for scan in range(0,diff_data.shape[1],6):
#    bf_map = bright_field((diff_data[:,scan]),nbr_rows,nbr_cols)
#    
#    fig, ax = plt.subplots(nrows=1)
#    ax.imshow(np.log10(bf_map),cmap = 'jet', interpolation='none', extent= extent_microns)
#    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]')#,title=r'Bright field' 3d?
#    plt.show()
    #plt.savefig(savepath_bf+ r'\rot%d'%scan)
#print(r'bf saved to ')
#print(savepath_bf)
    
bf_map = bright_field(np.sum((diff_data[:,:]),axis=1),nbr_rows,nbr_cols)
fig, ax = plt.subplots(nrows=1)
ax.imshow(np.log10(bf_map),cmap = 'jet', interpolation='none', extent= extent_microns)
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]')#,title=r'Bright field' 3d?
plt.show()


#%%

print('***XRD analysis***')

# Function defines q1 q2 q3 + q_abs from the geometry function 
# (See "Bending and tilts in NW..." ppt)
def def_q_vectors():
    global q3, q1, q2, q_abs    
    #  units of reciprocal meters [m-1]
    q_abs = 4 * np.pi / g.wavelength * g.sintheta()
    
    q1 = np.linspace(-g.dq1*g.shape[1]/2.+q_abs/g.costheta(), g.dq1*g.shape[1]/2.+q_abs/g.costheta(), g.shape[1]) # ~z
    # q3 defined as centered around 0, that means adding the component from q1
    q3 = np.linspace(-g.dq3*g.shape[0]/2. + g.sintheta()*q1.min() , g.dq3*g.shape[0]/2.+ g.sintheta()*q1.max(), g.shape[0]) #~x       
    q2 = np.linspace(-g.dq2*g.shape[2]/2., g.dq2*g.shape[2]/2., g.shape[2]) #        ~y
    
def_q_vectors()

#test of orientations
fig, ax = plt.subplots(ncols=1) # figsize=(10,3)
ax.imshow(np.log10(sum((diff_data[:,0]))),cmap='jet', interpolation='none')
tick_interval = 10
plt.yticks(range(0,len(q1),tick_interval), np.round(q1[::tick_interval]*1E-10,2))
plt.xticks(range(0,len(q2),tick_interval), np.round(q2[::tick_interval]*1E-8,2))
ax.set_ylabel('q1 qz [Å-1]')
ax.set_xlabel('q2 qy [*10^8 m-1]')


# --------------------------------------------------------------
# Make a meshgrid of q3 q1 q2 and transform it to qx qz qy.
# Also define the vectors qx qz qy
#----------------------------------------------------------------

# in the transformation is should be input and output: (qx, qz, qy), or (q3, q1, q2).
# make q-vectors into a tuple to transform to the orthogonal system; Large Q means meshgrid, small means vector
Q3,Q1,Q2 = np.meshgrid(q3, q1, q2, indexing='ij') 
  

# COM analysis is more exact if the measured data is transformed to a orthogonal grid
# transform the Q-space grid to from experimental grid to orthoganal
# go from natura to cartesian coordinate system
Qz = g.costheta() * Q1
#obs should have a minus sign
Qy = -Q2
Qx = Q3 - g.sintheta() * Q1
            

# NOTE Q2-Qy should not have changed but the other ones should. note 2, Q3 and Qx should be much larger (never negative).

qx = Qx[:,0,0]
qz = Qz[0,:,0]
qy = Qy[0,0,:]

plt.figure(); plt.plot(qx); plt.title('qx')
plt.figure(); plt.plot(qz); plt.title('qz')
plt.figure(); plt.plot(qy); plt.title('qy')

fig, ax = plt.subplots(ncols=3)
img = ax[0].imshow(Qx[:,:,0]); plt.colorbar(img)
img2 = ax[1].imshow(Qz[:,:,0]); plt.colorbar(img2)


#test of orientations
fig, ax = plt.subplots(ncols=1) # figsize=(10,3)
ax.imshow(np.log10(sum((diff_data[:,0]))),cmap='jet', interpolation='none')
tick_interval = 10
plt.yticks(range(0,len(qz),tick_interval), np.round(qz[::tick_interval]*1E-10,2))
plt.xticks(range(0,len(qy),tick_interval), np.round(qy[::tick_interval]*1E-8,2))
ax.set_ylabel('q1 qz [Å-1]')
ax.set_xlabel('q2 qy [*10^8 m-1]')

#%%

# TODO: remove single photon count, if the COM is calculated for very small values like pixels with 1 photon counts, 
#then the result will be missleading. Set a threshold that keeps the resulting pixel on a mean value, like if sum(sum(sum(diffPattern)))< threshold. sum(sum(sum()))==bkg_value
# input is 3d matrix with [nbr_rotations][nbr_pixels_x][nbr_pixels_y]
def COM_voxels_reciproc(data, Qx, Qz, Qy ):

    # center of mass calculation in reciprocal space with the meshgrids
    COM_qx = np.sum(data* Qx)/np.sum(data)
    COM_qz = np.sum(data* Qz)/np.sum(data)
    COM_qy = np.sum(data* Qy)/np.sum(data)

    #print( 'coordinates in reciprocal space:')
    #print( COM_qx, COM_qz, COM_qy)
    return COM_qx, COM_qz, COM_qy

# loop through all scanning postitions and move the 3D Bragg peak from the 
# natural to the orthogonal coordinate system (to be able to calculate COM)
# Calculate COM for every peak - this gives the XRD matrices
def XRD_analysis():
    position_idx = 0
    XRD_qx = np.zeros((nbr_rows,nbr_cols))
    XRD_qz = np.zeros((nbr_rows,nbr_cols))
    XRD_qy = np.zeros((nbr_rows,nbr_cols))

    for row in range(0,nbr_rows):
        print(row)
        for col in range(0,nbr_cols):
            
            #shifting coordinate system code taken from ptypy
            keep_dims = True
            # create a padded copy of the data array
            shape = diff_data[position_idx].shape
            pad = int(np.ceil(g.sintheta() * shape[1]))
            d = np.pad(diff_data[position_idx], pad_width=(
                (0, pad), (0, 0), (0, 0)), mode='constant')
            # walk along the q1/qz axis and roll the q3/qx axis. the
            # array is padded at the right (high indices) so the
            # actual shifts have to be positive.
            for i in range(shape[1]):

                # roll the q3 axis in the positive direction for more
                # negative q1
                shift = int(round((shape[1] - i) * g.sintheta()))
                d[:, i, :] = np.roll(d[:, i, :], shift, axis=0)
            # optionally crop the new array
            if keep_dims:
                d = d[pad // 2:shape[0] + pad // 2, :, :]
                
            data_orth_coord = d

            # do the 3d COM analysis to find the orthogonal reciprocal space coordinates of each Bragg peak
            #  for this to be correct with the coordinate system the 
            # data should be sorted with higher index = higher theta
            COM_qx, COM_qz, COM_qy = COM_voxels_reciproc(data_orth_coord, Qx, Qz, Qy)
           # insert coordinate in reciprocal space maps 
            XRD_qx[row,col] = COM_qx
            XRD_qz[row,col] = COM_qz
            XRD_qy[row,col] = COM_qy
            
            position_idx += 1
            
            ###plot every other 3d peak and print out the postion of the COM analysis
            #if (position_idx%100==0 and position_idx<501):
            if (position_idx==134 or position_idx==139 or position_idx==161): # 
                #import pdb; pdb.set_trace()
                # TODO very har to say anything about this looking in 2d, need 3d plots!
                #TODO plot the right position in 3D, that is look at the correct slice           
                x_p = np.argwhere(qx>COM_qx)[0][0]
                y_p = np.argwhere(qy<COM_qy)[0][0] #take the first value in qy where
                z_p = np.argwhere(qz>COM_qz)[0][0]  
                #import pdb; pdb.set_trace()
                
                print('figure')
                print(COM_qx*1E-7,np.round(COM_qz*1E-10,3),np.round(COM_qy*1E-8,3))
                print('x_p is', x_p)
                print('z_p',z_p)                
                print('y_p is', y_p)
                
                # it might be plotting like 1 pixel of, but in x that is pretty important. that is why i am plotting 3 rotations. be
                #because the COM is not finding the COM-pixel it is finding just the COM coordinate in the mesh
                for iii in range(-2,3,2):  
                    #import pdb; pdb.set_trace()
#                    fig, ax = plt.subplots(ncols=3)
#                    ax[0].plot(sum(sum(data_orth_coord)))
#                    ax[1].plot(np.sum(sum(data_orth_coord),axis=1))
#                    ax[2].plot(np.sum(np.sum(data_orth_coord,axis=1),axis=1))
  
                    #print(np.sum(data_orth_coord[x_p+iii]))
                    fig, ax = plt.subplots(ncols=1) # figsize=(10,3)
                    plt.imshow(np.log10(sum(data_orth_coord)), cmap='jet')
                    plt.colorbar()
                    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
                    tick_interval = 10
                    plt.yticks(range(0,len(qz),tick_interval), np.round(qz[::tick_interval]*1E-10,3))
                    plt.xticks(range(0,len(qy),tick_interval), np.round(qy[::tick_interval]*1E-8,3))
                    ax.set_ylabel('q1 qz [Å-1]')
                    ax.set_xlabel('q2 qy [*10^8 m-1]')
                    
                    # Find the coordinates of that cell closest to this value:              
                    plt.scatter(y_p, z_p, s=500, c='red', marker='x')
                    #plt.scatter(COM_qy, COM_qz, s=500, c='red', marker='x')
                    #plt.axhline(y=z_p,color='red')
                    #plt.axvline(x=y_p,color='yellow')
                    plt.title('pos%d'%position_idx)

    return XRD_qx, XRD_qz, XRD_qy, data_orth_coord

XRD_qx, XRD_qz, XRD_qy, data_orth_coord = XRD_analysis() # units of 1/m

#plt.figure(); plt.imshow(((sum(diff_data[101]))),cmap='jet', interpolation='none'); plt.title('position 100'); plt.colorbar()

#----------------------------------------------------------
# Convert q-vector from  cartesian coordinates to spherical
# (See "Bending and tilts in NW..." pp)
#----------------------------------------------------------

XRD_absq =  np.sqrt(XRD_qx**2 + XRD_qy**2 + XRD_qz**2)
XRD_alpha = np.arcsin( XRD_qy/ XRD_absq)
XRD_beta  = np.arctan( XRD_qx / XRD_qz)
# detta borde väl definieras som, hur stor är vinkeln från detektorn, när detektorn är på center av denteoretiska bragg vinkeln, alltså

print('OOOOOOOOOOBs, beta definition')
#XRD_beta  = np.arctan( (XRD_qx + qx[int(len(qx)/2)]) / XRD_qz)


#%%
print('***plot XRD maps***')

  
    
fig = plt.figure(figsize=(10,3))
plt.title(r'Full range BF', loc='left', pad =-12, color ='black')
#norm1 = DivergingNorm(vmin=0, vcenter=0.15, vmax = 1)
plt.imshow(bf_map/bf_map.max(), cmap='gray', interpolation='none')#,extent=extent_microns_cut, norm = norm1)
plt.ylabel(r'$y$ [$\mathrm{\mu}$m]') #
plt.xticks([])
#po = plt.colorbar(ticks=(10,20,30,40))#,fraction=0.046, pad=0.04) 
po = plt.colorbar()
tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
# if you want no mask use:
#XRD_mask = np.ones((XRD_mask.shape))


#plt.scatter(13, 2, s=500, c='red', marker='x')
#plt.scatter(18, 2, s=500, c='red', marker='x')
#plt.scatter(40, 2, s=500, c='red', marker='x')

#%%
class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x=%i, y=%i, z=%1.4f' % (x, y, z)
    
  
#def plot_XRD_polar():    
# cut the images in x-range:start from the first pixel: 
start_cutXat = 0

# replace the x-scales end-postion in extent_motorposition. 
extent_motorpos_cut = np.copy(extent_microns)
###extent_motorpos_cut[1] = 2.0194197798363955 segmentedNW
cutXat = -1   
cutYat = -1
origin = np.array([-1.1271288e-06, -1.1271288e-06])
#extent_microns_cut = [origin[0]*1E6, origin[0]*1E6 + dx*(cutXat-start_cutXat)*1E6, origin[1]*1E6, origin[1]*1E6+dy*(cutYat)*1E6]
extent_microns_cut = [0, (0.06)*(cutXat-start_cutXat), 0, 0.06*(cutYat)]

print('mean value BF: ')
print( np.mean(sum(bf_map)))
# create a mask from the BF matrix, for the RSM analysis
XRD_mask = np.copy((bf_map))
XRD_mask[XRD_mask < 6.2e13 ] = np.nan #inP 25000  #ingap 16000   # for homo InP, use 280000
XRD_mask[XRD_mask > 0] = 1       #make binary, all values not none to 1

colormap = 'RdBu_r' #Spectral' #RdYlGn'#coolwarm' # 'bwr' #'PiYG'# #'RdYlBu' # 'PiYG'   # 'PRGn' 'BrBG' 'PuOr'

# plot abs q to select pixels that are 'background', not on wire, and set these pixels to NaN (make them white)
#fig = plt.figure(figsize=(12,3))
fig, ax = plt.subplots(ncols=1, nrows=4)
plt.subplots_adjust(hspace = 0.245)
norm1 = DivergingNorm(vmin=0, vcenter=0.15, vmax = 1)
img0 = ax[0].imshow(bf_map[:cutYat,start_cutXat:cutXat]/bf_map[:,start_cutXat:cutXat].max(), cmap='gray', interpolation='none')#,extent=extent_microns_cut, norm = norm1)
ax[0].set_title(r'Total intensity', loc='left', pad =-5, color ='black')
ax[0].set_ylabel(r'$\mu$m')
ax[0].set_xticks([])
#po = plt.colorbar(ticks=(10,20,30,40))#,fraction=0.046, pad=0.04) 
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(img0, cax=cax, ticks=(0, 0.5, 1))

tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()

# if you want no mask use:
XRD_mask = np.ones((XRD_mask.shape))

#plt.subplot(412)   
#calculate lattice constant c from |q|:       
d =  2*np.pi/  XRD_absq
c_lattice_exp = d * np.sqrt(3)   



#print 'mean lattice constant is %d' %np.mean(a_lattice_exp)
#imshow(a_lattice_exp)
#plt.title('Lattice conastant a [$\AA$]')
mean_a = np.nanmean(XRD_mask[:cutYat,start_cutXat:cutXat]*a_lattice_exp[:cutYat,start_cutXat:cutXat])
#TODO try with reference strain equal to the center of the largest segment (for InP) # tody try with reference from the other NWs
#mean_strain = a_lattice_exp[:cutYat,start_cutXat:cutXat].max() 

strain = 100*XRD_mask[:cutYat,start_cutXat:cutXat]*(a_lattice_exp[:cutYat,start_cutXat:cutXat]-mean_a)/mean_a
#strain = 100*XRD_mask[:cutYat,start_cutXat:cutXat]*(( a_table - a_lattice_exp[:cutYat,start_cutXat:cutXat] )/a_lattice_exp[:cutYat,start_cutXat:cutXat])

norm2 = DivergingNorm( vcenter=0)
img1 = ax[1].imshow(strain, cmap=colormap, interpolation='none')#, extent=extent_microns_cut)#
#plt.imshow(XRD_mask[:cutYat,start_cutXat:cutXat]*a_lattice_exp[:cutYat,start_cutXat:cutXat], cmap='jet',interpolation='none')#,extent=extent_motorpos_cut) 

#plt.imshow(XRD_mask[:cutYat,start_cutXat:cutXat]*XRD_absq[:cutYat,start_cutXat:cutXat]*1E-8, cmap=colormap,interpolation='none',extent=extent_motorpos_cut) 
#plt.title('Length of Q-vector |Q|)', loc='left', pad =-12)
ax[1].set_title(r'd-spacing 004', loc='left', pad =-12)   #
#plt.title('Lattice constant a')
#plt.ylabel(r'$y$ [$\mathrm{\mu}$m]')  
ax[1].set_ylabel(r'$\mu$m')
ax[1].set_xticks([])
divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(img1, cax=cax)#, ticks=(-0.15, 0, 0.15 ,0.30))

tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()

#plt.subplot(413)
norm3 = DivergingNorm( vcenter=0)
img2 = ax[2].imshow(XRD_mask[:cutYat,start_cutXat:cutXat]*1E3*XRD_alpha[:cutYat,start_cutXat:cutXat], cmap='coolwarm', interpolation='none')#,extent=extent_microns_cut, norm=norm3) 
# cut in extent_motorposition. x-pixel nbr 67 is at 2.0194197798363955
ax[2].set_title(r'$\alpha$ (mrad)', loc='left', pad =-12)
ax[2].set_ylabel(r'$\mu$m')
ax[2].set_xticks([])
divider = make_axes_locatable(ax[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(img2, cax=cax)#, ticks=(-1.5, 0, 1.5))

tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()

#po = plt.colorbar(ticks=(0,1,2,3,4))
#po.set_label('Bending around $q_x$ $\degree$')
   
#plt.subplot(414)
norm4 = DivergingNorm( vcenter=0)
img3 = ax[3].imshow(XRD_mask[:cutYat,start_cutXat:cutXat]*1E3*XRD_beta[:cutYat,start_cutXat:cutXat], cmap='bwr',interpolation='none')#,extent=extent_microns_cut, norm=norm4) 
plt.setp(ax[3].xaxis.get_majorticklabels(), rotation=70 )
ax[3].set_title('$\\beta$ (mrad)', loc='left', pad =-12)
#ax[3].set_xticks([])
#plt.ylabel(r'$y$ [$\mathrm{\mu}$m]')
#plt.xlabel(r'$x$ [$\mathrm{\mu m}$]') 
ax[3].set_xlabel(r'$\mu$m')
ax[3].set_ylabel(r'$\mu$m')
divider = make_axes_locatable(ax[3])
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(img3, cax=cax)#, ticks=(-6, -3, 0 ))

#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator; po.update_ticks()

#plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InGaP/%s_xrd'%(date_str))
#plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InGaP/%s_xrd.pdf'%(date_str)) 
#plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InP/%s_xrd'%(date_str))
#plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InP/%s_xrd.pdf'%(date_str)) 
#%%
#plot and save lineout of strain XRD
lineout_strain = strain[2]
lineout_xx = np.linspace(0,extent_microns_cut[1],len(lineout_strain))

plt.figure()
plt.plot(lineout_xx,strain[2])#bf_map[2,start_cutXat:cutXat])#strain[1])
plt.plot(lineout_xx,bf_map[2,start_cutXat:cutXat]/bf_map[2,start_cutXat:cutXat].max())
plt.title('Strain and BF lineout')

#np.save(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\XRD\20220405_lineout_XRD',lineout_strain)
#np.save(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\XRD\20220405_lineout_xx',lineout_xx)  


#%%

def q_abs_map():
    q_map_temp = XRD_mask[:,start_cutXat:cutXat]*XRD_absq[:,start_cutXat:cutXat]
    return q_map_temp
q_map_temp = q_abs_map()    

#mean value of the Q_abs map
mean_absq = np.nanmean(q_map_temp)
print( 'Mean |q| value in the masked map' )
print( mean_absq)

plt.figure()
plt.title('Mean |q| is: %f '%mean_absq )
plt.imshow(q_map_temp); plt.colorbar()
#plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InP/%s_xrd_with_mean_qabs_value.pdf'%(date_str)) 
plt.figure()
plt.title('lattice constant a')
plt.imshow(a_lattice_exp[:,start_cutXat:cutXat]); plt.colorbar()
