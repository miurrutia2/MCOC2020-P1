
from scipy.integrate import odeint
import numpy as np

rt = 6371.*10**3
mt = 5972e24 #kg
omega = -7.2921150e-5 #rad/s
G = 6.67408e-11 #m3 kg-1 s-2

def R(t):
    R = np.array([[np.cos(omega*t), np.sin(omega*t), 0],[-np.sin(omega*t), np.cos(omega*t), 0], [0, 0, 1]])
    return R

def Rp(t):
    Rp = omega*np.array([[-np.sin(omega*t), np.cos(omega*t), 0], [-np.cos(omega*t), -np.sin(omega*t), 0], [0, 0, 0]])
    return Rp

def Rpp(t):
    Rpp = (omega**2)*np.array([[-np.cos(omega*t), -np.sin(omega*t), 0], [np.sin(omega*t), -np.cos(omega*t), 0], [0, 0, 0]])
    return Rpp

def satelite(z, t):
	zp = np.zeros(6)
	z1 = z[0:3]
	z2 = z[3:6]

	r2 = np.dot(z1,z1)
	r = np.sqrt(r2)
	Fg = (-G*mt/r**2) * (R(t)@(z1/r))

	zp[0:3] = z2
	zp[3:6] = R(t).T@(Fg - (2*(Rp(t)@z2) + (Rpp(t)@z1)))

	return zp


from datetime import datetime

ti = "2020-07-23T22:59:42.000000"
ti = ti.split("T")
ti ="{} {}".format(ti[0], ti[1])
ti = datetime.strptime(ti, '%Y-%m-%d %H:%M:%S.%f')

tf = "2020-07-25T00:59:42.000000"
tf = tf.split("T")
tf ="{} {}".format(tf[0], tf[1])
tf = datetime.strptime(tf, '%Y-%m-%d %H:%M:%S.%f')

deltaT = (tf - ti).total_seconds()

# Datos iniciales
# <TAI>TAI=2020-07-23T23:00:19.000000</TAI>
#<UTC>UTC=2020-07-23T22:59:42.000000</UTC>
#<UT1>UT1=2020-07-23T22:59:41.786009</UT1>
#<Absolute_Orbit>+33588</Absolute_Orbit>
#<X unit="m">-1659308.501505</X>
#<Y unit="m">6753819.113038</Y>
#<Z unit="m">-1317158.267893</Z>
#<VX unit="m/s">1875.247563</VX>
#<VY unit="m/s">-954.537857</VY>
#<VZ unit="m/s">-7297.506531</VZ>

xi = -1659308.501505
yi = 6753819.113038
zi = -1317158.26789
vxi = 1875.247563
vyi = -954.537857
vzi = -7297.506531

# Datos finales
#  <TAI>TAI=2020-07-25T01:00:19.000000</TAI>
#  <UTC>UTC=2020-07-25T00:59:42.000000</UTC>
#  <UT1>UT1=2020-07-25T00:59:41.786602</UT1>
#  <Absolute_Orbit>+33604</Absolute_Orbit>
#  <X unit="m">-12440.495142</X>
#  <Y unit="m">3460608.473275</Y>
#  <Z unit="m">6161862.778150</Z>
#  <VX unit="m/s">2412.558073</VX>
#  <VY unit="m/s">6276.959093</VY>
#  <VZ unit="m/s">-3512.365075</VZ>

xf = -12440.495142
yf = 3460608.473275
zf = 6161862.778150
vxf = 2412.558073
vyf = 6276.959093
vzf = -3512.365075




t = np.linspace(0, deltaT, 9361)

z0 = np.array([xi, yi, zi, vxi, vyi, vzi])
sol = odeint(satelite, z0, t)

pos_final = np.array([xf, yf, zf, vxf, vyf, vzf]) - sol[-1]

for i in pos_final:
    print(i)