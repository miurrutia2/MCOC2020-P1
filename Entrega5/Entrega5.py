
from scipy.integrate import odeint
import numpy as np
from leer_eof import leer_eof
import matplotlib.pyplot as plt
from time import perf_counter


rt = 6371.*10**3
mt = 5.972e24 #kg
omega = -7.2921150e-5 #rad/s
G = 6.67408e-11 #m3 kg-1 s-2
J2 =  (1.75553e10)*(1000**5) #km5⋅s−2
J3 = (-2.61913e11)*(1000**6) #km6⋅s−2
zp = np.zeros(6)

def satelite(z,t):
    c = np.cos(omega*t)
    s = np.sin(omega*t)
    R = np.array([[c,s,0],[-s,c,0],[0,0,1]])
    
    Rp =omega* (np.array([[-s,c,0],[-c,-s,0],[ 0,0,0]]))
    
    Rpp = (omega**2)* (np.array([[-c,-s,0],[s,-c,0],[0,0,0]]))
    
    z1 = z[0:3]
    z2 = z[3:6]
    
    r2 = np.dot(z1,z1)
    r = np.sqrt(r2)
    
    Fg = (-G*mt/r**2)* (R@(z1/r))
    
    zp[0:3] = z2
    zp[3:6] = R.T@(Fg-(2*(Rp@z2) + (Rpp@z1)))
    return zp


#-------p1-------
fname = "S1A_OPER_AUX_POEORB_OPOD_20200813T120829_V20200723T225942_20200725T005942.EOF"

t, x, y, z, vx, vy, vz = leer_eof(fname)  #fname es el nombre de su EOF

z0 = np.array([x[0], y[0], z[0], vx[0], vy[0], vz[0]])

zf = np.array([x[-1], y[-1], z[-1], vx[-1], vy[-1], vz[-1]])

from datetime import datetime

sol = odeint(satelite, z0, t)

plt.figure()
plt.subplot(3,1,1)
plt.title("Posicion Real")
plt.ylabel("X(t)")
plt.plot(t,x)
plt.subplot(3,1,2)
plt.ylabel("Y(t)")
plt.plot(t,y)
plt.subplot(3,1,3)
plt.ylabel("Z(t)")
plt.plot(t,z)
plt.xlabel("tiempo")
plt.savefig("P1.png")
plt.show()

#-------P2-------

def eulerint(zp, z0, t, Nsubdivisiones):
  Nt = len(t)
  Ndim = len(z0)

  z = np.zeros((Nt, Ndim))
  z[0, :] = z0

  for i in range(1, Nt):
    t_anterior = t[i-1]
    dt = (t[i]- t[i-1])/Nsubdivisiones
    z_temp = z[i-1, :].copy()
    for k in range(Nsubdivisiones):
      z_temp += dt*zp(z_temp, t_anterior + k*dt)
    z[i, :] = z_temp
  return z

t1 = perf_counter()
sol_euler = eulerint(satelite, z0, t, 1)
t2 = perf_counter()

t3 = perf_counter()
sol_odeint = odeint(satelite, z0, t)
t4 = perf_counter()

deriva1 = []
for i in range(len(t)):
    dif = sol_odeint[i]-sol_euler[i]    
    mult = np.dot(dif,dif)
    deriva1.append(np.sqrt(mult))

tx = t/3600.

print(f"La solucion de eulerint demora {t2-t1} segundos")
print(f"La solucion de odeint demora {t3-t2} segundos")

plt.figure()
plt.title("Diferencia eulerint y odeint")
plt.xlabel("Tiempo")
plt.ylabel("Deriva")
plt.plot(tx,np.array(deriva1)/1000)
plt.savefig("P2.png")
plt.show()



#---------P3--------
Nsubdivisiones = 10
t1 = perf_counter()
sol_euler = eulerint(satelite, z0, t, Nsubdivisiones)
t2 = perf_counter()

t3 = perf_counter()
sol_odeint = odeint(satelite, z0, t)
t4 = perf_counter()

print(f"La solucion de eulerint para {Nsubdivisiones}  demora {t2-t1} segundos")
print(f"La solucion de odeint demora {t3-t2} segundos")

deriva2 = []
for i in range(len(t)):
    zx = np.array([x[i], y[i], z[i], vx[i], vy[i], vz[i]])
    dif = zx - sol_euler[i]    
    mult = np.dot(dif,dif)
    deriva2.append(np.sqrt(mult))


plt.figure()
plt.title("Distancia entre eulerint y odeint")
plt.xlabel("Tiempo")
plt.ylabel("Deriva")
plt.plot(tx,np.array(deriva2)/1000)
plt.savefig("P2.png")
plt.show()


#--------P4--------

def satelite_J2_J3(z,t):
    c = np.cos(omega*t)
    s = np.sin(omega*t)

    R = np.array([[c,s,0],[-s,c,0],[0,0,1]])
    
    Rp =omega* (np.array([[-s,c,0],[-c,-s,0],[ 0,0,0]]))
    
    Rpp = (omega**2)* (np.array([[-c,-s,0],[s,-c,0],[0,0,0]]))
    
    z1 = z[0:3]
    z2 = z[3:6]
    
    r2 = np.dot(z1,z1)
    r = np.sqrt(r2)

    FxJ2 = (J2*z[0]/(r**7))*(6*z[2]**2 - (3/2)*(z[0]**2 + z[1]**2))
    FyJ2 = (J2*z[1]/(r**7))*(6*z[2]**2 - (3/2)*(z[0]**2 + z[1]**2))
    FzJ2 = (J2*z[2]/(r**7))*(3*z[2]**2 - (9/2)*(z[0]**2 + z[1]**2))

    FJ2 = [FxJ2,FyJ2,FzJ2]
    
    FxJ3 = (J3*z[0]*z[2]/(r**9))*(10*z[2]**2 - (15/2)*(z[0]**2 + z[1]**2))
    FyJ3 = (J3*z[1]*z[2]/(r**9))*(10*z[2]**2 - (15/2)*(z[0]**2 + z[1]**2))
    FzJ3 = (J3/(r**9))*(4*z[2]**2 * (z[2]**2 - 3*(z[0]**2 + z[1]**2) + (3/2)*(z[0]**2 + z[1]**2)**2))

    FJ3 = [FxJ3,FyJ3,FzJ3]

    zp[0:3] = z2
    Fg = (-G*mt/r**2)* (R@(z1/r))
    zp[3:6] = R.T@(Fg - (2 * Rp@z2) + (Rpp@z1)) + FJ2 + FJ3

    return zp



t1 = perf_counter()
sol_euler = eulerint(satelite_J2_J3, z0, t, 1)
t2 = perf_counter()

t3 = perf_counter()
sol_odeint = odeint(satelite_J2_J3, z0, t)
t4 = perf_counter()

deriva1 = []
for i in range(len(t)):
    dif = sol_odeint[i]-sol_euler[i]    
    mult = np.dot(dif,dif)
    deriva1.append(np.sqrt(mult))

tx = t/3600.

print(f"La solucion de eulerint demora {t2-t1} segundos")
print(f"La solucion de odeint demora {t3-t2} segundos")

plt.figure()
plt.title("Diferencia eulerint y odeint J2 y J3")
plt.xlabel("Tiempo")
plt.ylabel("Deriva")
plt.plot(tx,np.array(deriva1)/1000)
plt.savefig("P3.png")
plt.show()




Nsubdivisiones = 10
t1 = perf_counter()
sol_euler = eulerint(satelite_J2_J3, z0, t, Nsubdivisiones)
t2 = perf_counter()

t3 = perf_counter()
sol_odeint = odeint(satelite_J2_J3, z0, t)
t4 = perf_counter()

print(f"La solucion de eulerint para {Nsubdivisiones}  demora {t2-t1} segundos")
print(f"La solucion de odeint demora {t3-t2} segundos")

deriva2 = []
for i in range(len(t)):
    zx = np.array([x[i], y[i], z[i], vx[i], vy[i], vz[i]])
    dif = zx - sol_euler[i]    
    mult = np.dot(dif,dif)
    deriva2.append(np.sqrt(mult))


plt.figure()
plt.title("Distancia entre eulerint y odeint J2 y J3")
plt.xlabel("Tiempo")
plt.ylabel("Deriva")
plt.plot(tx,np.array(deriva2)/1000)
plt.savefig("P4.png")
plt.show()