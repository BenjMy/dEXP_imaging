import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import utils_dEXP as uEXP

# import exemples.fwd_mag_sphere as magfwd
import examples_in_prep.load_MALM_PM as MALM

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

MainPath= r'E:\Padova\Software\SourceInversion\Potential_field_imaging\dEXP_imaging\examples_in_prep\\'
os.chdir(MainPath)

#%%


## --------- read MALM measured data file (with electrode positions) --------- ##

# RealData = np.loadtxt("./1_Data_2_plot/to_plot.dat",skiprows=0,delimiter='\t') 
out = MALM.load_MALM_Porto_real(MainPath + '/malm_models/',
                          './malm_models/XYObs_real_f_m3.txt',
                          shape=(100,100),
                          radius=200,
                          rcor=10,
                          rot=60,
                          showfig=True)

coord_xyz, coord_xyz_int = out[0:2]
Uload = out[2]
coords_liner = out[3]
shape, max_elevation = out[4:6]

p = out[6]         # line points                                       
# set imaging pseudo-inversion parameters                                                                        
parameters = para.set_par(shape=shape,max_elevation=max_elevation)

scaled = parameters[0]
SI = parameters[1]
zp, qorder, nlay = parameters[2:5]
# minAlt_ridge, maxAlt_ridge = parameters[5:7]

nb_of_peak = 1
# nlay = 50

# max_elevation = 30
minAlt_ridge = max_elevation*0.05
maxAlt_ridge = max_elevation*0.65
    
# nlay = 50

# xp, yp, zp = coord_xyz_int
xp, yp, zp = coord_xyz
# len(xp)
U = Uload[1] # U_raw, Ucor, U_int, Ucor_int
p1 , p2 = p

interp=False
smooth = True
# len(U)
# len(xp)

# plt.figure()
# plt.scatter(xp, yp)
# plt.colorbar()
# plt.plot(coords_liner[:,0],coords_liner[:,1],'-')
# plt.axis('square')
# plt.scatter((Ha[0]+Ha[1])/2,Ha[2])

plt.figure()
plt.scatter(xp, yp, c=U , cmap='viridis',vmax=0.05)
plt.colorbar()
plt.axis('square')
plt.show()

print(U[0:2])


# NameSave = 'Togrid.txt'
# f = open(NameSave,'w')
# np.savetxt(f, np.array([xp, yp, U]).T, delimiter='\t',fmt='%1.6f')   # X is an array
# f.close()



# xx, yy, distance, profile = pEXP.plot_line(xp, yp,U,p1,p2,
#                                             interp=True, smooth=False, Xaxis=x_axis)


# %% change p1p2 axis 

# plt.figure()
# plt.scatter(xp, yp, c=U, label='U_int')
# plt.grid()
# plt.legend()
# plt.xlabel('x (m)')
# plt.ylabel('y (m)')
    
# p1,p2 = uEXP.perp_p1p2(p1,p2, offset=0)
# x_axis = 'y'
x_axis = 'x'
# x_axis = 'dist'

# %% zeros sym

#%% select part of the grid
# uEXP.mirror(p[:,0],p[:,1],data,a,b)
# uEXP.mirror(xp,yp,data,p1,p2)

#%% ------------------------------- smooth the data 

# U = dEXP.smooth2d(xp, yp, U, sigma=5)

#%% ------------------------------- Plot the data 
# U = np.copy(Us)

# _, p1, p2, _ = MALM.definep1p2(path=Main,radius=300)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=True)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=False, smooth=True)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=False, smooth=False)

# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U_f ,p1,p2, interp=interp)


#%%
shape = (100,100)
xp_int,yp_int,U_int = gridder.interp(xp,yp,U,shape,extrapolate=False)


plt.figure()
plt.scatter(xp_int,yp_int, c=U_int, label='U_int')
plt.grid()
plt.legend()
p12x=[p1[0],p2[0]]
p12y=[p1[1],p2[1]]
plt.plot(p12x,p12y,c='red')
p12v = plt.scatter(p12x, p12y,c='red')
txt = ['p1', 'p2']
for i in range(len(p12x)):
    plt.annotate(txt[i], (p12x[i], p12y[i]),c='red')
    
plt.xlabel('x (m)')
plt.ylabel('y (m)')

prl = 60
# shape = shape  (max(xp)-min(xp))/
# shape = (115,115)
xint_scipy, yint_scipy = gridder.regular((min(xp)-prl, max(xp)+prl, 
                          min(yp)-prl, max(yp)+prl),shape=shape)

#%% Solution 1
# extrapolate False and fill with 0 before derivative - mask them later on 
U_int_scipy = gridder.interp_at(xp,yp,U, xint_scipy, yint_scipy, algorithm='cubic', extrapolate=False)
InterpData = np.array([xint_scipy, yint_scipy, U_int_scipy]).T
where_are_NaNs = np.isnan(InterpData)
InterpData[where_are_NaNs] = 0.0074
xint_scipy, yint_scipy, U_int_scipy = InterpData.T

plt.scatter(xint_scipy,yint_scipy)
plt.scatter(xp,yp)


#%% Solution 2
# Extrapolate = True
# U_int_scipy = gridder.interp_at(xp,yp,U, xint_scipy, yint_scipy, algorithm='cubic', extrapolate=True)


# %%
# len(U)
pEXP.plot_field(xint_scipy,yint_scipy,U_int_scipy, shape)

# pEXP.plot_field(xp_int,yp_int,U_int, shape)
# xx, yy, distance, profile = pEXP.plot_line(xp_int,yp_int,U_int,p1,p2,
#                                             interp=True, smooth=False, Xaxis=x_axis)

# %%

plt.figure()
plt.subplot(1,2,1)
plt.scatter(xint_scipy,yint_scipy, c=U_int_scipy, label='U_int')
plt.grid()
plt.legend()
p12x=[p1[0],p2[0]]
p12y=[p1[1],p2[1]]
plt.plot(p12x,p12y,c='red')
p12v = plt.scatter(p12x, p12y,c='red')
txt = ['p1', 'p2']
for i in range(len(p12x)):
    plt.annotate(txt[i], (p12x[i], p12y[i]),c='red')
    
plt.xlabel('x (m)')
plt.ylabel('y (m)')

plt.subplot(1,2,2)
plt.scatter(xp,yp, c=U, label='U_int')
plt.grid()
plt.legend()
p12x=[p1[0],p2[0]]
p12y=[p1[1],p2[1]]
plt.plot(p12x,p12y,c='red')
p12v = plt.scatter(p12x, p12y,c='red')
txt = ['p1', 'p2']
for i in range(len(p12x)):
    plt.annotate(txt[i], (p12x[i], p12y[i]),c='red')
    
plt.xlabel('x (m)')
plt.ylabel('y (m)')


xx, yy, distance, profile, ax, plt = pEXP.plot_line(xint_scipy,yint_scipy,U_int_scipy,p1,p2,
                                            interp=False, smooth=True, Xaxis=x_axis, showfig=True)


# xx, yy, distance, profile = pEXP.plot_line(xint_scipy,yint_scipy,U_int_scipy,p1,p2,
#                                             interp=interp, smooth=smooth, Xaxis=x_axis)

# xx, yy, distance, profile = pEXP.plot_line(xint_scipy,yint_scipy,U_int_scipy,p1,p2,
#                                             interp=True,smooth=True, Xaxis='x_axis')


# # xx, yy, distance, profile = pEXP.plot_line(xg, yg,U_int ,p1,p2, interp=False, smooth=smooth)
# # xx, yy, distance, profile = pEXP.plot_line(xg, yg,U_int ,p1,p2, interp=False, smooth=True)


# grd = uEXP.load_surfer('grid_ascii.grd')
# # xp_int,yp_int,U_int_surf  = grd['y'],grd['x'], grd['data']
# xp, yp, U   = grd['y'],grd['x'], grd['data']
# # xx, yy, distance, profile = pEXP.plot_line(xp_int,yp_int,U_int_surf ,p1,p2, interp=interp, smooth=smooth)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=interp,
#                                            smooth=smooth, Xaxis='dist')

  
#%%

# Select the data that fall inside "section"
# section = [-40, 40, -25, 25]
# section = [p1[0], p2[0],p2[1], p1[1]]

# # Tip: you pass more than one data array as input. Use this to cut multiple
# # data sources (e.g., gravity + height + topography).
# x_sub, y_sub, [data_sub] = gridder.cut(xp, yp, [U], section)

# plt.figure()
# plt.scatter(x_sub, y_sub, c=data_sub, label='U_int')
# plt.grid()
# plt.legend()
# plt.xlabel('x (m)')
# plt.ylabel('y (m)')
# len(data_sub)

# np.sqrt(len(y_sub))

#%% ------------------------------- Pad the edges of grids

xpad,ypad,U_pad, shape_pad = dEXP.pad_edges(xint_scipy,yint_scipy,U_int_scipy,shape,
                                            pad_type=0) # reflexion=5

# xpad,ypad,U_pad, shape_pad = dEXP.pad_edges(xp, yp,U,shape_raw,
#                                             pad_type=0) # reflexion=5


# pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)
# pEXP.plot_line(xpad, ypad,U_pad,p1,p2, interp=interp)

plt.figure()
plt.scatter(xpad,ypad, c=U_pad, cmap='viridis')
cbar = plt.colorbar(orientation='vertical', aspect=50)
cbar.set_label('$u_{pad}$ (V)')
# plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
plt.axis('square')

#%% ------------------------------- Plot the derivatives
# xp, yp, U, shape = xpad,ypad,U_pad,shape_pad
xp, yp, U, shape = xint_scipy,yint_scipy,U_int_scipy, shape
# xp, yp, U, shape = xp, yp, U, (60,60)
# len(xp)
# len(yp)

xderiv = transform.derivx(xp, yp, U, shape,order=qorder)
yderiv = transform.derivy(xp, yp, U, shape,order=qorder)
zderiv = transform.derivz(xp, yp, U, shape,order=qorder)

# interp = True
pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, 
               interp=interp, smooth=smooth,  Xaxis=x_axis,showfig=True)
pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, 
               interp=interp, smooth=smooth,  Xaxis=x_axis,showfig=True)
pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, 
               interp=interp, smooth=smooth,  Xaxis=x_axis,showfig=True)

# # pEXP.plot_line(xp_int,yp_int,U_int ,p1,p2, interp=interp, smooth=smooth)
# pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=True, smooth=smooth,  Xaxis=x_axis)
# pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, interp=True, smooth=smooth,  Xaxis=x_axis)
# pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=True, smooth=smooth,  Xaxis=x_axis)


#%% ------- upward continuation of the field data

zp = np.zeros(len(xp))

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                  zmin=0, zmax=max_elevation, nlayers=nlay, 
                  qorder=qorder)

#%% 

# Xaxis = 'y'
# image.shape
# image = mesh.props[label_prop].reshape(mesh.shape)
# mins, maxs = [image.min(),image.max()]
# fig = plt.subplots()
# ax = plt.gca()


#%%
# pEXP.plot_xy(mesh, label=label_prop) #, ldg=)

# pEXP.plot_xy(mesh, label=label_prop, Xaxis='x') #, ldg=)
pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=p) #, ldg=)

# p1p2 = np.array([1,1,1,1])


# plt, cmap = pEXP.slice_mesh(xp, yp, mesh, label_prop, p1, p2, interp=True, Xaxis='y')
# plt.colorbar(cmap)

# %% ridges identification

dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      fix_peak_nb=nb_of_peak,
                                      interp=interp,smooth=smooth,
                                      method_peak='find_peaks',
                                      showfig=True,
                                      Xaxis=x_axis)  

# or  find_peaks or peakdet or spline_roots
# dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,interp=interp,
#                                       label=label_prop,fix_peak_nb=nb_of_peak,
#                                       smooth=smooth,
#                                       method_peak='find_peaks',
#                                       showfig=True,
#                                       Xaxis=x_axis)  

D = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks',
                                      fix_peak_nb=3,
                                      returnAmp=True,
                                      showfig=True,
                                      Xaxis=x_axis,
                                      interp=interp,
                                      smooth = smooth,
                                      qorder=qorder)  

dfI, dfII, dfIII =  D[0:3]
hI, hII, hIII  = D[3:6]
heights  = D[3:6]
    
# dfI, dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       # minAlt_ridge=minAlt_ridge,
#                                       maxAlt_ridge=maxAlt_ridge,
#                                       fix_peak_nb=None) 
 
#%% ------------------------------- plot ridges over continuated section
    
# fig = plt.figure()
# ax = plt.gca()
# # pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis='y')
# pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=p,ax=ax) #, ldg=)

# # pEXP.slice_mesh(xp, yp, mesh, label_prop, p1, p2, 
# #                             interp=True, ax=ax, Xaxis='dist')
# pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

fig = plt.figure()
ax = plt.gca()

plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=np.array([p1, p2]), ax=ax) #, ldg=)
plt.colorbar(cmap)
pEXP.plot_ridges_harmonic(dfI, dfII, dfIII,ax=ax,label=True)
    

#%% ------------------------------- filter ridges regionally constrainsted)
   

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True)

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            minAlt_ridge,
                                            maxAlt_ridge,
                                            minlength=5,
                                            rmvNaN=True)
                                            # xmin=5.036e6+50)                                            )

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True,
#                                            xmin=284200)

df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

# pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
# pEXP.slice_mesh(xp, yp, mesh, label_prop, p1, p2, 
#                             interp=True, ax=ax, Xaxis='y')
pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=p,ax=ax) #, ldg=)

pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
#                           ridge_type=[0,1,2],ridge_nb=None)

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None,
                          xmin=5.036e6+50, 
                          xmax=5.036e6+400)

#%% DEXP ratio

qratio = [1,0]
mesh_dexp, label_dexp = dEXP.dEXP_ratio(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorders=qratio)
fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
             markerMax=False,qratio=str(qratio),
             p1p2=np.array([p1,p2]), ax=ax, Xaxis='y') #, ldg=)
plt.colorbar(cmap)

qratio = [1,0]
mesh_dexp, label_dexp = dEXP.dEXP_ratio(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorders=qratio)
fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
             markerMax=True,qratio=str(qratio),Vminmax=[0,0.075],
             p1p2=np.array([p1,p2]), ax=ax, Xaxis=x_axis) #, ldg=)
cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('ratio voltage (V)')
# plt.xlim([200,600])
plt.savefig('ratios_field.png', dpi=450)    
    
    
# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
#                           ridge_type=[0,1,2],ridge_nb=None,
#                           xmin=5.036e6+50, 
#                           xmax=5.036e6+400)

# plt.savefig(dataname +  '_fig_DEXP_Ratio_' + x_axis + '.png', r=400)
 
    
# uEXP.multipage('DEXP_PortoM_zerofill.pdf')

