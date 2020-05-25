from fatiando_dev import dEXP


# ------------------------------- Load data

shape = (50, 50)
x,y,z,gz=np.loadtxt('3d_SENS_FAT.txt', unpack=True)

# remove B return electrode effects 
B = [-104.16666666666667, -104.16666666666667, 0]
gz_cor = dEXP.cor_field_B(x,y,z,gz,B,rho=100)


