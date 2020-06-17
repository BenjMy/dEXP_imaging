# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:02:59 2020

@author: Benjamin
"""

import numpy as np 
import matplotlib.pyplot as plt

def mirrorImage( a, b, c, x1, y1): 
	temp = -2 * (a * x1 + b * y1 + c) /(a * a + b * b) 
	x = temp * a + x1 
	y = temp * b + y1 
	return (x, y) 

# %% Create matrice 
# X = np.array([[1, 2, 3, 4],  
#         [5, 6, 7, 8], 
#         [9, 10, 11, 12],  
#         [13, 14, 15, 16],
#         [17, 18, 19, 20],
#         [21, 22, 23, 24]]); 

nx, ny = (10,10)
nel = nx*ny
# X = np.logspace(0.1,1,nel); 
X = np.arange(0,nel,1); 
X = X.reshape(nx,ny)

x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)

# %% Create line from two points p1 and p2 

# p1_ind = [1,1]
p2_ind = [9,4]

p1_ind = [0,9]
# p2_ind = [40,20]
# p2_ind = [80,30]


p12_indx=[p1_ind[0],p2_ind[0]]
p12_indy=[p1_ind[1],p2_ind[1]]

p1 = np.array([x[p1_ind[0]],y[p1_ind[1]]])
p2 = np.array([x[p2_ind[0]],y[p2_ind[1]]])

p12x=[p1[0],p2[0]]
p12y=[p1[1],p2[1]]

# %% Calcul slope of the line

slope = (p2[1]-p1[1])/(p2[0]-p1[0]) # (yb-ya)/(xb-xa)
# # convert slope to matrice indexes
# slope_ind = (y[0]-y[1])/(x[0]-x[1]) # (yb-ya)/(xb-xa)
# slp = int(slope_ind/slope)

# %% check position of points against line

isabove = lambda p, p1,p2: np.cross(p-p1, p2-p1) < 0

xvcol = xv.reshape([nx*ny,1])
yvcol = yv.reshape([nx*ny,1])

p = np.hstack([xvcol,yvcol])
c_above=isabove(p,p1,p2)
c_bool = c_above.reshape(xv.shape)

xvcol_f = np.delete(xvcol , np.where(c_above == True))
yvcol_f = np.delete(yvcol , np.where(c_above == True))
X_f = np.delete(X.reshape([nx*ny,1]) , np.where(c_above == True))


X_bool = np.copy(X)

plt.figure()
plt.scatter(xvcol_f, yvcol_f,c=X_f, cmap='viridis')
# plt.scatter(xvcol[c_bool==True], yv[c_bool==True],c='red', cmap='viridis',s=1)
# plt.scatter(xvcol[c_bool==False], yv[c_bool==False],c='black', cmap='viridis',s=1)

plt.colorbar()
plt.plot(p12x,p12y, marker="o", color="k")
plt.axis('square')
plt.show()

plt.figure()
plt.scatter(xv, yv,c=X, cmap='viridis')
plt.scatter(xv[c_bool==True], yv[c_bool==True],c='red', cmap='viridis')
plt.scatter(xv[c_bool==False], yv[c_bool==False],c='black', cmap='viridis')



# %% calcul equation of the line using polyfit

# Calculate the coefficients. This line answers the initial question. 
coefficients = np.polyfit(p12x, p12y, 1)

# Print the findings
print ('a =', coefficients[0])
print ('b =', coefficients[1])

a = -coefficients[0]
b = 1
c = -coefficients[1]

plt.figure()
plt.plot(p12x,p12y)
plt.scatter(xv, yv,c=X_bool, cmap='viridis', vmin=0, vmax=len(X)**2)

ind_dist = []
for i, bool_pi in enumerate(zip(c_above,p)):
    if bool_pi[0] == False:
        xmir, ymir = mirrorImage(a, b, c, bool_pi[1][0], bool_pi[1][1]); 
        plt.scatter(xmir, ymir,c=X_bool.reshape(len(X)**2)[i], cmap='viridis',s=1e2, vmin=0, vmax=len(X)**2)
        plt.annotate(str(i)+ '_m', [xmir, ymir])
        plt.annotate(str(i), [bool_pi[1][0],  bool_pi[1][1]])

plt.axis('square')


# %% re-interpolate on the grid position 

U = gridder.interp_at(x, y, U_raw, xnew, ynew, algorithm='cubic', extrapolate=True)
xp,yp,U_int = gridder.interp(xnew,ynew,U,shape)


#%%

boolmat = np.ones(X.shape)
indlow = np.tril_indices(boolmat.shape[0])

boolmat[indlow] = 0

indlow2x = indlow[0][0:int(len(indlow[0])/2)]
indlow2y = indlow[1][0:int(len(indlow[0])/2)]

boolmat[indlow2x,indlow2y] = 5

plt.figure()
plt.imshow(boolmat,vmin=0,vmax=5)
plt.colorbar()
plt.plot(p12_indx,p12_indy, marker="o", color="k")



from scipy import ndimage, misc
boolmat_angle = ndimage.rotate(boolmat, np.rad2deg(slope), reshape=False)
# full_X_angle = ndimage.rotate(X, 0, reshape=False)

plt.figure()
# for i in range(boolmat.shape[0]):
#     for j in range(i,boolmat.shape[1]):
#         # plt.scatter(x, y, kwargs)
#         plt.annotate(str(boolmat_angle[i][j]), (i, j))

plt.imshow(boolmat_angle,vmin=0,vmax=5)
plt.colorbar()
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.title('rotated matrice')

# plt.figure()
# XX_mir = np.copy(boolmat_angle)
# for i in range(XX_mir.shape[0]):
#     for j in range(i,XX_mir.shape[1]):
#         XX_mir[i][j] = boolmat[j][i]
#         plt.annotate(str(boolmat[j][i]), (i, j))
# plt.imshow(XX_mir,vmin=0,vmax=5)
# plt.plot(p12_indx,p12_indy, marker="o", color="k")
# plt.title('all section mirrored/ diag')



XX_mir_slope = np.copy(boolmat_angle)
plt.figure()
slp = 3

fig, ax = plt.subplots()
ax0 = plt.subplot(1,2,1) 
ax1 = plt.subplot(1,2,2) 

ax0.imshow(boolmat)
ax1.imshow(boolmat_angle)
ax0.plot(p12_indx,p12_indy, marker="o", color="k")
ax1.plot(p12_indx,p12_indy, marker="o", color="k")

for i in range(p1_ind[0],p2_ind[0]+1):
    for j in range(i,p2_ind[0]+1):
        si = (i-slp,j)
        if (i-slp)>0:
            XX_mir_slope[int(i-slp)][j] = boolmat[j][i-slp]
            
            if c_bool[i-slp,j] == True:
                ax1.scatter(i-slp, j, c='red',label='target')
            else:
                ax1.scatter(i-slp, j, c='green')
            ax1.annotate(str(i-slp)+ str(j), (i-slp, j))

            ax0.scatter(j, i-slp, c='black')
            ax0.annotate(str(j)+ str(i-slp), (j, i-slp),label='initial')
        elif (i-slp)<0:
            ax1.scatter(i-slp, j, c='blue')
            pass
        # plt.legend()

plt.figure()
plt.subplot(1,3,1)
plt.imshow(boolmat)
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.subplot(1,3,2)
plt.imshow(boolmat_angle)
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.subplot(1,3,3)
plt.imshow(XX_mir_slope)
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.suptitle('section mirrored against line 2')

# %% Plot matrice

# slope = 3
plt.figure()
for i in range(X.shape[0]):
    for j in range(i,X.shape[1]):
        plt.annotate(str(X[i][j]), (i, j))
        
plt.imshow(X,vmin=0,vmax=len(X)*len(X))
plt.colorbar()
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.title('Initial matrice')

from scipy import ndimage, misc
full_X_angle = ndimage.rotate(X, np.rad2deg(slope), reshape=False)
# full_X_angle = ndimage.rotate(X, 0, reshape=False)

plt.figure()
for i in range(X.shape[0]):
    for j in range(i,X.shape[1]):
        # plt.scatter(x, y, kwargs)
        plt.annotate(str(full_X_angle[i][j]), (i, j))

plt.imshow(full_X_angle,vmin=0,vmax=len(X)*len(X))
plt.colorbar()
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.title('rotated matrice')


# plt.figure()
# plt.imshow(X.T)
# plt.plot(p12_indx,p12_indy, marker="o", color="k")
# plt.title('matrice.T')

# %% Prepare for mirroring

plt.figure()
XX = np.copy(X)
for i in range(X.shape[0]):
    for j in range(i,X.shape[1]):
        XX[i][j] = X[j][i]
        plt.annotate(str(X[j][i]), (i, j))
plt.imshow(XX)
plt.title('all section mirrored/ diag')

plt.figure()
XX_angle = np.copy(full_X_angle)
for i in range(full_X_angle.shape[0]):
    for j in range(i,full_X_angle.shape[1]):
        XX_angle[i][j] = X[j][i]
        plt.annotate(str(X[j][i]), (i, j))
plt.imshow(XX_angle)
plt.title('all section mirrored/ diag')

# %%

# plt.figure()

# angles = np.arange(1,50,10)
# for a, angle in enumerate(angles):
#     XX = np.copy(X)
#     for i in range(p1_ind[0],p2_ind[0]+2):
#         for j in range(i,p2_ind[0]+2):
#             si = (i-angle,j)
#             if (i-angle)<0:
#                 pass
#             else:
#                 print(i,j,i-angle)
#                 XX[i+angle][j-angle] = full_X_angle[j][i]

#         plt.subplot(1,len(angles),a+1)
#         plt.imshow(XX)
#         plt.plot(p12_indx,p12_indy, marker="o", color="k")
#         plt.title('Mirror angle:' + str(angle))



# %%

XX = np.copy(full_X_angle)
indlow = np.tril_indices(XX.shape[0])

fig, ax = plt.subplots()
ax0 = plt.subplot(1,2,1) 
ax1 = plt.subplot(1,2,2) 

ax0.imshow(X)
ax1.imshow(full_X_angle)
ax0.plot(p12_indx,p12_indy, marker="o", color="k")
ax1.plot(p12_indx,p12_indy, marker="o", color="k")

slp = 0
for i in range(p1_ind[0],p2_ind[0]+1):
    for j in range(i,p2_ind[0]+1):
        si = (i-slp,j)
        # print(si)
        # if (indlow[0]- si[0]).any==0 and (indlow[1]- si[1]).any==0:
        if (i-slp)>0:
            XX[int(i-slp)][j] = X[j][i-slp]
            
            if c_bool[i-slp,j] == True:
                ax1.scatter(i-slp, j, c='red',label='target')
            else:
                ax1.scatter(i-slp, j, c='green')
            ax1.annotate(str(i-slp)+ str(j), (i-slp, j))

            ax0.scatter(j, i-slp, c='black')
            ax0.annotate(str(j)+ str(i-slp), (j, i-slp),label='initial')
        elif (i-slp)<0:
            ax1.scatter(i-slp, j, c='blue')
            pass
        # plt.legend()

plt.figure()
plt.subplot(1,3,1)
plt.imshow(X)
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.subplot(1,3,2)
plt.imshow(full_X_angle)
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.subplot(1,3,3)
plt.imshow(XX)
plt.plot(p12_indx,p12_indy, marker="o", color="k")
plt.suptitle('section mirrored against line 2')

# %%

# # X = X + X.T - np.diag(np.diag(X))



def squaremat(r=20):
    # MainPath= 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM'
    # os.chdir(MainPath)
    # grid = pv.read('Ano_field_EA.vtk') 
    # grid.set_active_scalars("_Marker")
    # cpos = grid.plot()
    
    # bodies = grid.split_bodies()
    # for i, body in enumerate(bodies):
    #     print('Body {} volume: {:.3f}'.format(i, body.volume)) 
    
    
    # position of the anomaly
    ya =  5.036227e6
    xa =  (284272.0 + 284250.00)/2
    
    # landfill geometry
    # utm coordinates 
    coords_liner = [(284046.43,	5036328.39),
              (284132.24,	5036277.32),
              (284146,	5036297),
              (284281.07,	5036214.66),
              (284317.46,	5036281.81),
              (284245,	5036313),
              (284097.55,	5036411.08)]
    coords_liner= np.asarray(coords_liner)
    
    slope = (coords_liner[2,1]-coords_liner[3,1])/(coords_liner[2,0]-coords_liner[3,0]) # (yb-ya)/(xb-xa)
    x1 = xa + r*np.cos(slope)
    y1 = ya + r*np.sin(slope)
    
    x2 = xa - r*np.cos(slope)
    y2 = ya - r*np.sin(slope)
    
    plt.figure(figsize=(20,10))
    ax = plt.subplot()
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    for i in range(len(coords_liner[:,0])):
        ax.annotate(str(i), (coords_liner[i,0], coords_liner[i,1]))
        
    plt.scatter(xa,ya,c='red')
    plt.scatter(x1,y1,c='blue')
    plt.scatter(x2,y2,c='blue')
    plt.scatter(x2,y1,c='green')
    plt.scatter(x1,y2,c='green')

    plt.axis('square')

    xp = np.array([x1,x2])
    yp = np.array([y1,y2])
    
    return xp,yp
    
    
    
    
    
    

    return xp, yp

