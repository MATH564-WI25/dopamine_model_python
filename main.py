import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# b = ld
def VAADC(b):
    k1  = 130 
    V =   10000  
    a = V*b/(k1 + b)
    return a

#catabolism of 5ht in the terminal
# b = extracellular DA
def VCATAB(b):
    k1 = 3  #Gottowik93 and Fowler94
    k2 = 30
    a = (k2*b/(k1 + b))
    return a

#SERTS
# b = eda
def VDAT(b):
    k = .2
    a = 8000*(b/(k + b))
    return a

# b = BH2
# c = NADPH
# d = BH4
# e = NADP
def VDRR(b,c,d,e):
    k1 = 100   #Km for BH2 (BRENDA) (6-650)
    k2 = 75  #Km for NADPH (BRENDA, values 770,110,29) (schumber 70-80)
    V1 = 200 # Vmax forward
    k3 = 10  #Km for BH4 (BRENDA) (1.1 to 17)
    k4 = 75  #Km for NADP (BRENDA)(schumber 70-80)
    V2 = 80 # Vmax backward
    # forward direction from BH2 to BH4
    print(d,e)
    a = V1*b*c/((k1 + b)*(k2 + c)) - V2*d*e/((k3 + d)*(k4 + e))
    return a
#regulatory function in vivo, above 0.03 mM, oxidation of Met146 and Met151
#leads to inactivation of the enzyme due to disruption of the NADH-binding
#site (BRENDA) 1.5.1.34


#vesicular monoamine transporter
# b = cDA
#c= vda
def VMAT(b,c):
    k = 3
    V = 7082
    a = ((V*b/(k + b))-(40)*c)  #40
    return a
# Km for uptake by vesicles = 123 or 252 nM (Slotkin77)
# Km for uptake by vesicles =  198 nM (Rau06)

# b = trp
# c = pool function
def VPOOL(b, c):
    k1 = 9 #to pool
    k2 = .6 # from pool
    a = (k1* b - k2* c)
    return a

# b = tyr
# c = BH4
# d = cda
# e = eda
def VTH(b,c,d,e):
    k1 = 160
    k11 = 110
    ktyr = 46
    kbh4 = 60
    vmax = (2)*125
    n = 8

    # factor for effect of autoreceptors
    auto_factor = ((4.5)/(8*(e/0.002024)**4 + 1) + 0.5)

    a = ((0.56)/(1 + b/k1))*auto_factor*(vmax*b*c)/(b*c + ktyr*c + ktyr*kbh4*(1 + d/k11))
    return a

# tryptophan transporter from blood to cytosol
# b = btyr
def VTYRin(b):
    k1 = 400
    k2 = 64
    a = (k1*b/ (k2 + b))
    return a

"""Modify this function to change firing bahavior"""
def fire(t):
    return 1 if t<0.1 else 0

# ODEs
def rhs(y,t):

    dy=np.zeros(10)

    btyr = 90
    NADP = 1*26
    NADPH = 330
    k1 =   6  #to pool
    k11 = .6  #from pool
    k2 = .2   #cell catab of tyr
    k3 = 10   #cell catab of DA
    k4 = 400  #removal of eda
    k5 = 3.45 #hva catab
    k6 = .2   #tyrpool catab

    dy[1] = VTH(y[3],y[2],y[5],y[7]) - VDRR(y[1],NADPH,y[2],NADP)
    dy[2] = VDRR(y[1],NADPH,y[2],NADP) - VTH(y[3],y[2],y[5],y[7])
    dy[3] = VTYRin(btyr) - VTH(y[3],y[2],y[5],y[7]) - k1*y[3] + k11*y[9] - k2*y[3]
    dy[4] = VTH(y[3],y[2],y[5],y[7]) - VAADC(y[4])
    dy[5] = VAADC(y[4]) - VMAT(y[5],y[6]) + VDAT(y[7]) - k3*y[5]
    dy[6] = VMAT(y[5],y[6]) - fire(t)*y[6]
    dy[7] = fire(t)*y[6] - VDAT(y[7]) -  VCATAB(y[7]) - k4*y[7]
    dy[8] = k3*y[5] + VCATAB(y[7]) - k5*y[8]
    dy[9] = k1*y[3] - k11*y[9] - k6*y[9]

    # y[1] = bh2
    # y[2] = bh4
    # y[3] = tyr
    # y[4] = ld
    # y[5] = cda
    # y[6] = vda
    # y[7] = eda
    # y[8] = hva
    # y[9] = yrppool
    return dy

hours = 5
tspan = np.linspace(0,hours,60*hours) # time span is 5 hours
Y0 = np.array([0, 59.2406, 300.7594, 113.5721, 0.5312, 3.9684, 98.3679, 0.0025, 11.5097, 851.7909])
Y  = odeint(rhs,Y0,tspan)
print(Y)


# save solutions
bh2 = Y[:,1]
bh4 = Y[:,2]
tyr  = Y[:,3]
ld = Y[:,4]
cda = Y[:,5]
vda = Y[:,6]
eda = Y[:,7]
hva = Y[:,8]
tyrpool = Y[:,9]

# plot solutions
fig, axs = plt.subplots(2,3,figsize=(6,6),constrained_layout=True)
labels = ["tyr","ld","cda","vda","eda","hva"]
sols = [tyr, ld, cda, vda, eda, hva]
for i in range(6):
    ax=axs[i//3,i%3]
    ax.plot(tspan,sols[i],lw=2)
    ax.legend('Location','best')
    ax.set_xlabel('minutes')
    ax.set_ylabel('\muM')
    ax.set_ylim([0.75*np.min(sols[i]),1.25*np.max(sols[i])])
    ax.set_title(labels[i])
    print("h")
fig.show()
