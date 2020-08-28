
from scipy.integrate import odeint
import matplotlib.pylab as plt
import numpy as np

#Parametros

rt = 6371 #km
omega = (7.27 * 10**-5)*3600 #rad/h
mt = 5.972 * 10**24 #kg
G = (6.6729 * 10**-11)*(3600**2)*(1000**-3) #km^3/(kg^2*h^2)

def R(t):
    r0 = np.array([[np.cos(omega*t), -(np.sin(omega*t)), 0],[np.sin(omega*t), np.cos(omega*t), 0],[0, 0, 1]])
    return r0

def Rp(t):
    r1 = omega*np.array([[-np.sin(omega*t), -(np.cos(omega*t)), 0],[np.cos(omega*t), -np.sin(omega*t), 0],[0, 0, 0]])
    return r1

def Rpp(t):
    r2 = (omega**2)*np.array([[-np.cos(omega*t), np.sin(omega*t), 0],[-np.sin(omega*t), -np.cos(omega*t), 0],[0, 0, 0]])
    return r2

def satelite(z, t):
    zp = np.zeros(6)
    zp[0:3] = z[3:6]
    r = np.sqrt(z[0]**2 + z[1]**2 + z[2]**2)
    z1 = z[0:3]
    z2 = z[3:6]
    zp[3:6] = (-G*mt/r**3)*z1 - R(t).transpose()@(Rpp(t)@z1) - 2*R(t).transpose()@(Rp(t)@z2)

    return zp

vt = 24650 #probando valores

#x,y,z,vx,vy,vz
z0 = np.array([6371 + 700, 0, 0, 0, vt, 0])
t = np.linspace(0, 3.2, 2000) #probando para obtener las 2 vueltas (horas)
sol = odeint(satelite, z0, t)

x = sol[:, 0]
y = sol[:, 1]
z = sol[:, 2]

plt.figure()


plt.subplot(3,1,1)
plt.plot(t,x, color = "pink", label = "x(t)")
plt.ylabel("x(t)")
plt.legend()
plt.subplot(3,1,2)
plt.plot(t,y,color = "lightblue", label = "y(t)")
plt.ylabel("y(t)")
plt.legend()
plt.subplot(3,1,3)
plt.plot(t,z, color = "purple", label = "z(t)")
plt.xlabel("t")
plt.ylabel("z(t)")
plt.legend()
plt.savefig("graficos.png")
plt.show()



plt.plot(x, y, "--", color = "magenta", label = "Satelite")
angulo = np.linspace(0, 2*np.pi, 1000)
xt = rt*np.sin(angulo)
yt = rt*np.cos(angulo)
plt.plot(xt,yt, color = "purple", label = "Tierra")
xa = (rt+80)*np.sin(angulo)
ya = (rt+80)*np.cos(angulo)
plt.plot(xa,ya, color = "lightseagreen", label = "Atmosfera")
plt.xlabel("t")
plt.ylabel("r(t)")
plt.grid()
plt.tight_layout() 
plt.legend()
plt.savefig("orbitas.png")
plt.show()

plt.plot(x, y, "--", color = "magenta", label = "Satelite")
angulo = np.linspace(0, 2*np.pi, 1000)
xt = rt*np.sin(angulo)
yt = rt*np.cos(angulo)
plt.plot(xt,yt, color = "purple", label = "Tierra")
xa = (rt+80)*np.sin(angulo)
ya = (rt+80)*np.cos(angulo)
plt.plot(xa,ya, color = "lightseagreen", label = "Atmosfera")
plt.xlabel("t")
plt.ylabel("r(t)")
plt.grid()
plt.tight_layout() 
plt.legend()
plt.savefig("orbitas.png")
plt.xlim(-7000,-5000)
plt.ylim(-7000,7000)
plt.show()

