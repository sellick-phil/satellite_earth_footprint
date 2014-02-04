# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 19:45:54 2014

@author: Simon Oliver

Calculating the effective heading of a descending orbit including satellite skew
Projected position of input point
"""
import math
import ephem
# Calculate heading and skew for an input sensor and time
# variables here are latitude (ephem),oi_deg (tle), orad (tle), av 

#TLE FORMAT######
#LINE 2

#Field 	Columns     Content 	                                  Example
#1 	    01–01 	     Line number 	                               2
#2 	    03–07 	     Satellite number 	                               25544
#3 	    09–16 	     Inclination [Degrees] 	                         51.6416
#4 	    18–25 	     Right Ascension of the Ascending Node [Degrees] 247.4627
#5 	    27–33 	     Eccentricity (decimal point assumed) 	         0006703
#6 	    35–42 	     Argument of Perigee [Degrees] 	               130.5360
#7 	    44–51 	     Mean Anomaly [Degrees] 	                     325.0288
#8 	    53–63 	     Mean Motion [Revs per day] 	                     15.72125391
#9 	    64–68 	     Revolution number at epoch [Revs] 	          56353
#10 	   69–69 	     Checksum (Modulo 10) 	                           7 
#################
longitude = 145
latitude = 0
lat_rad = math.radians(latitude) # test latitude
# Satellite Parameters (Brouwer MOE)					   
								
oi_deg = 98.2034                # Orbital Inclination (OI) [degrees] = 98.2034		
oi_rad = math.radians(oi_deg)   # Orbital Inclination (OI) [radians] = 1.713972666653			
orad = 7077689.2                 # Orbit Radius (R) [m] = 7.08E+06		
np = 5925.816                   # Nodal Period [sec] = 5925.816		
av = 2*math.pi/np          # Angular Velocity (V0) [rad/sec] =	 0.001060307189285 =2*PI()/E8				
sr = 0                          # Sensor Roll (r) [degrees] =	0				
						
# Earth Stuff (WGS84)
oneonf = 298.257223563          # Inverse flattening 1/f = 298.257223563
f = 1/oneonf                    # flattening
r = 6378137                     # Radius (a) [m] =	 6378137
e = 1-math.pow((1-1/oneonf),2)  # Eccentricty (e^2) = 0.00669438 =1-(1-1/I5)^2
wO = 0.000072722052             # rotation (w0) [rad/sec] = 7.2722052E-05

xfac = math.sqrt(1-e*(2-e)*(math.pow(math.sin(math.radians(latitude)),2)))# Xfac = =SQRT(1-$I$7*(2-$I$7)*(SIN(B16*PI()/180)^2))
phi_rad = math.asin((1-e)*math.sin(math.radians(latitude))/xfac) # Phi0' (Geocentric latitude) = =ASIN((1-$I$7)*SIN(B16*PI()/180)/C16)
phi_deg = math.degrees(phi_rad) # Phi0' (Degrees) = =D16*180/PI() 
n = r/math.sqrt(1-e*(math.pow(math.sin(math.radians(latitude)),2)))# N =$I$6/SQRT(1-$I$7*(SIN(B16*PI()/180)^2))
altphi_rad = latitude-180*math.asin(n*e*math.sin(lat_rad)*math.cos(lat_rad)/orad)/math.pi# Alt Phi0'Radians = =B16-180*ASIN(I16*I$7*SIN(PI()*B16/180)*COS(B16*PI()/180)/E$7)/PI() 

rho_rad = math.acos(math.sin(altphi_rad*math.pi/180)/math.sin(oi_rad))# Rho (Radians) = =ACOS(SIN(U16*PI()/180)/SIN($E$6))
print rho_rad
print altphi_rad
print n
beta = -1*(math.atan(1/(math.tan(oi_rad)*math.sin(rho_rad)))*180/math.pi) # Heading Beta (degrees) = =-ATAN(1/(TAN($E$6)*SIN(F16)))*180/PI() 
print beta
xn = n*xfac # Xn = =I16*C16
print xn # Altitude (km) = =($E$7-J16)/1000
altitude = (orad-xn)/1000 # altitude = =($E$7-J16)/1000
print altitude
altitude_ = (orad*math.cos(altphi_rad/180*math.pi)/math.cos(lat_rad)-n)/1000 # =($E$7*COS(U16/180*3.14159)/COS(B16/180*3.14159)-I16)/1000
print altitude_
rotn = math.atan((wO*math.cos(phi_rad)*math.cos(beta*math.pi/180))/(av+wO*math.cos(phi_rad)*math.sin(beta*math.pi/180)))*180/math.pi
# Rotn = =ATAN(($I$8*COS(D16)*COS(G16*PI()/180))/($E$9+$I$8*COS(D16)*SIN(G16*PI()/180)))*180/PI()
print rotn
eh = beta+rotn # Effective Heading = Heading Beta + Rotn
print eh
alpha12 = eh
s = 0.5*185000 #s = distance in metres
