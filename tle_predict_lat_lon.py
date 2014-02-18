# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:21:04 2014

@author: user
"""

# Script to sit and poll for upcoming satellite passes.
# Joseph Armbruster

#import datetime
from datetime import timedelta, datetime
import ephem
import math
import os
import sys
import time
import urllib2
import logging
import wget

from pyorbital.orbital import Orbital
from pyorbital import tlefile

#from datetime import datetime, timedelta
import math
from geographiclib.geodesic import Geodesic
import osgeo.ogr, osgeo.osr
from osgeo import ogr
from osgeo import gdal
import numpy
#from geopy.distance import vincenty

# TODO(joe) : add getopt configuration options for this to main.
#_latlong = ('28.5340','-81.2594') # user lat/long
_latlong = ('-23 42','133 54') # Alice Spring Data Acqusition Facility
_notify = 30 # let us know this many minutes in advance to a pass
_usevoice = True # use voice?
_statussleep = 1 # how many minutes to sleep between status updates

def GetTLEs():
    # GetTLEs(): returns a list of tuples of keps for each satellite.
    # This function currently relies on a url from amsat.org.

    # grab the latest keps
    #tles = urllib2.urlopen('http://www.amsat.org/amsat/ftp/keps/current/nasabare.txt').readlines()
    #tles = urllib2.urlopen('http://www.celestrak.com/norad/elements/resource.txt‎').readlines()
    #urllib2 open doesn't like something about the celestrak text file here's a crumby workaround...revisit
    resource = 'http://www.celestrak.com/norad/elements/resource.txt'
    weather = 'http://www.celestrak.com/norad/elements/weather.txt'
    try:
        os.remove('resource.txt')
        #wget.download(resource)
    except OSError:
        pass
    try:
        os.remove('weather.txt')
        #wget.download(weather)
    except OSError:
        pass
    wget.download(resource)
    wget.download(weather)
    filenames = ['weather.txt', 'resource.txt']
    with open('tles.txt', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    tles = open('tles.txt','r').readlines()
    #tlesweather = open('weather.txt','r').readlines()
    #tlesresource = open('resource.txt','r').readlines()
    #tles = tlesresource,tlesweather

    print "retrieving TLE file.........."
    #print tles
    # strip off the header tokens and newlines
    tles = [item.strip() for item in tles]

    # clean up the lines
    tles = [(tles[i],tles[i+1],tles[i+2]) for i in xrange(0,len(tles)-2,3)]

    return tles
def getVectorFile(attributes,inputpoints,polyorline,ogroutput,ogrformat):
    #getVectorFile(dictionary,list of dicts with lat2 and lon2, 'polygon', SWATH_FILENAME, 'KML')

    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    try:
        os.remove(ogroutput)
    except OSError:
        pass

    if os.path.exists(ogroutput):
        driver.DeleteDataSource(ogroutput)
    ds = driver.CreateDataSource(ogroutput)

    if polyorline is 'polygon':

        geomtype = ogr.wkbPolygon
    if polyorline is 'line':
        geomtype = ogr.wkbLine

    if ds is None:
        print 'Could not create file'
        sys.exit(1)
    layer = ds.CreateLayer(attributes['Satellite name'], geom_type=geomtype)
    #attributes = {'Satellite name': satname,'Orbit height':orad,'Orbit':orb.get_orbit_number(now), 'Current time':now,'Minutes to horizon':minutesaway,'AOS time':rt,'LOS time':st,'Transit time':tt, 'Node':node}
    # create a field for the county name
    fieldDefn = ogr.FieldDefn('Satellite Name               :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn('Orbit height                 :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)
    layer.CreateField(ogr.FieldDefn('Orbit number                 :', ogr.OFTInteger))
    fieldDefn = ogr.FieldDefn('Current UTC time             :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)

    #layer.CreateField(ogr.FieldDefn('Current UTC time             :', ogr.OFTString))
    #layer.CreateField(ogr.FieldDefn('Minutes to horizon           :', ogr.OFTTime))
    #layer.CreateField(ogr.FieldDefn('Acquisition of Signal UTC    :', ogr.OFTDateTime))
    #layer.CreateField(ogr.FieldDefn('Loss of Signal UTC           :',ogr.OFTDateTime))
    #layer.CreateField(ogr.FieldDefn('Transit time                 :', ogr.OFTTime))
    fieldDefn = ogr.FieldDefn('Node                         :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)

    # add the field to the shapefile
    #layer.CreateField(fieldDefn)

    # get the FeatureDefn for the shapefile
    featureDefn = layer.GetLayerDefn()

    if polyorline == 'polygon':
        ring = ogr.Geometry(ogr.wkbLinearRing)
        #ring = ogr.Geometry(ogr.wkbMultiPolygon)
        # add the vertex to the ring
        for x in inputpoints:
            ring.AddPoint(x['lon2'],x['lat2'])

        # now that we've looped through all of the coordinates, create a polygon
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        # create a new feature and set its geometry and attribute
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(poly)
        feature.SetField('Satellite Name               :', attributes['Satellite name'])
        feature.SetField('Orbit height                 :', attributes['Orbit height'])
        feature.SetField('Orbit number                 :', attributes['Orbit'])
        feature.SetField('Current UTC time             :', str(attributes['Current time']))
        #feature.SetField('Minutes to horizon           :', attributes['Minutes to horizon'])
        #feature.SetField('Acquisition of Signal UTC    :', attributes['AOS time'])
        #feature.SetField('Loss of Signal UTC           :', attributes['LOS time'])
        #feature.SetField('Transit time                 :', attributes['Transit time'])
        feature.SetField('Node                         :', attributes['Node'])

        # add the feature to the shapefile
        layer.CreateFeature(feature)

        # destroy the geometries and feature
        ring.Destroy()
        poly.Destroy()

    feature.Destroy()

    # now that we've looped through the entire text file, close the files
    #file.close()
    ds.Destroy()
    ####### Shapefile Creation Complete #########################################
    #build the output vector based on content of an input list or list of dicts

    return ()
def getEffectiveHeading(satellite, oi_deg, latitude, longitude, orad):
    # -*- coding: utf-8 -*-
    """
    Created on Tue Feb  4 19:45:54 2014

    @author: Simon Oliver

    Calculating the effective heading of a descending orbit including satellite skew
    Projected position of input point
    """
#    import math
    #import ephem
    # Calculate heading and skew for an input sensor and time
    # variables here are latitude (ephem),oi_deg (tle), orad (tle)

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
    #longitude = 145
    #latitude = 0
    lat_rad = math.radians(latitude) # test latitude
    # Satellite Parameters (Brouwer MOE)

    #oi_deg = 98.2034                # Orbital Inclination (OI) [degrees] = 98.2034
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

    beta = -1*(math.atan(1/(math.tan(oi_rad)*math.sin(rho_rad)))*180/math.pi) # Heading Beta (degrees) = =-ATAN(1/(TAN($E$6)*SIN(F16)))*180/PI()

    xn = n*xfac # Xn = =I16*C16

    altitude = (orad-xn)/1000 # altitude = =($E$7-J16)/1000

    altitude_ = (orad*math.cos(altphi_rad/180*math.pi)/math.cos(lat_rad)-n)/1000 # =($E$7*COS(U16/180*3.14159)/COS(B16/180*3.14159)-I16)/1000

    rotn = math.atan((wO*math.cos(phi_rad)*math.cos(beta*math.pi/180))/(av+wO*math.cos(phi_rad)*math.sin(beta*math.pi/180)))*180/math.pi

    eh = beta+rotn # Effective Heading = Heading Beta + Rotn

    alpha12 = eh
    s = 0.5*185000 #s = distance in metres
    effective_heading = alpha12
    return effective_heading

if __name__ == '__main__':

    observer = ephem.Observer()
    observer.lat = _latlong[0]
    observer.long = _latlong[1]
    
    tles = GetTLEs()

    while 1:
        now = datetime.utcnow()
        #print now
        print "---------------------------------------"
        # iterate through all the two line element sets
        for tle in tles:
            if tle[0] in ("LANDSAT 8", "LANDSAT 7", "TERRA", "AQUA", "NOAA 18", "NOAA 19", "SUOMI NPP"):   
                sat = ephem.readtle(tle[0],tle[1],tle[2])
                twole = tlefile.read(tle[0],'tles.txt')
                #print twole
                print "---------------------------------------"
                print tle[0]
                oi = float(str.split(tle[2],' ')[3])

                orb = Orbital(tle[0])
                now = datetime.utcnow()
                orad = orb.get_lonlatalt(now)[2]
                attributes = []

                print "Orbit height         = ",orad
                print "Orbit Number         = ",orb.get_orbit_number(now)
                #######
                print "Orbital Inclination  = ",oi
                #print str.split(tle[2],' ')
                rt, ra, tt, ta, st, sa = observer.next_pass(sat)
                observer.date = rt
                sat.compute(observer)
                localrisetime = ephem.localtime(rt)
                
                timeuntilrise = localrisetime-now
                
                minutesaway = timeuntilrise.seconds/60.0

                #attributes = [{'orad':orad,'onum':orb.get_orbit_number(now), 'Minutes to horizon':minutesaway,'AOStime':rt,'LOStime':st,'Transit time':tt}]
                print "Minutes to horizon   = ", minutesaway
                print "AOStime              = ", rt
                print "LOStime              = ", st
                print "Transit time         = ", tt
                sat.compute(rt)
                aos_lat = sat.sublat.real*(180/math.pi)
                sat.compute(st)
                los_lat = sat.sublat.real*(180/math.pi)
                if (aos_lat > los_lat):
                    print "PASS is descending"
                    node = "descending"

                else:
                    print "PASS is ascending"
                    node = "ascending"
                    oi = 360 - oi
                    #print oi
                satname = str(tle[0]).replace(" ","_")
                attributes = {'Satellite name': satname,'Orbit height':orad,'Orbit':orb.get_orbit_number(now), 'Current time':now,'Minutes to horizon':minutesaway,'AOS time':rt,'LOS time':st,'Transit time':tt, 'Node':node}
                ######## Create shapefile for each of the upcoming passes ###################

                ORBIT_FILENAME = satname+"_orbit_centre.kml"

                print "Writing output to    = ", ORBIT_FILENAME
                ogr.UseExceptions()

                driver = ogr.GetDriverByName('KML')
                spatialReference = osgeo.osr.SpatialReference()
                spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
                
                try:
                    os.remove(ORBIT_FILENAME)
                except OSError:
                    pass 
                 
                ds = driver.CreateDataSource(ORBIT_FILENAME)
                layer = ds.CreateLayer(str(tle[0]).replace(" ","_"), spatialReference, ogr.wkbLineString)
                line = ogr.Geometry(type=ogr.wkbLineString)

                # Step from AOS to LOS in 1 second intervals

                delta = timedelta(seconds=100)

                deltatime = datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S") 
                geoeastpoint = []
                geowestpoint = []
                print "Current location     = ", sat.sublong.real*(180/math.pi), sat.sublat.real*(180/math.pi) 

                while deltatime < datetime.strptime(str(st),"%Y/%m/%d %H:%M:%S"):

                    sat.compute(deltatime)

                    line.AddPoint(sat.sublong.real*(180/math.pi), sat.sublat.real*(180/math.pi))

                    eastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad)+90
                    westaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad)+270

                    oneonf = 298.257223563          # Inverse flattening 1/f = 298.257223563
                    f = 1/oneonf                    # flattening
                    r = 6378137
                    #Set ground swath per satellite sensor
                    if tle[0] in ("LANDSAT 8","LANDSAT 7"):
                        swath = 185000/2
                    if tle[0] in ("TERRA","AQUA"):
                        swath = 2330000/2
                    if tle[0] in ("NOAA 18", "NOAA 19"):
                        swath = 1100000/2
                    if tle[0] == "SUOMI NPP":
                        swath = 2200000/2

                    geoeastpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, eastaz, swath))

                    geowestpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, westaz, swath))

                    deltatime = deltatime+delta

                feature = ogr.Feature(feature_def=layer.GetLayerDefn())
                feature.SetGeometryDirectly(line)
                layer.CreateFeature(feature)
                feature.Destroy()

                ds.Destroy()
                ##################################
                # Create polygon
                # get the shapefile driver
                #driver = ogr.GetDriverByName('KML')

                # create a new data source and layer

                SWATH_FILENAME = satname+"_orbit_swath.kml"
                polypoints = []
                print geowestpoint
                for x in geowestpoint:
                    polypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                for x in reversed(geoeastpoint):
                    polypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                #polypoints = (geowestpoint+(geoeastpoint.reverse))
                #print polypoints

                getVectorFile(attributes,polypoints,'polygon', SWATH_FILENAME, 'KML')

                # notify if pass is within 30 minutes of horizon

                if minutesaway <= _notify:
                
                    if _usevoice and sys.platform=='darwin':
                        say = 'say "%s WILL BE MAKING A PASS IN %d MINUTES."' % (tle[0],minutesaway)
                        os.system(say)
                    print tle[0], 'WILL BE MAKING A PASS IN ', minutesaway, " MINUTES"
                    print ' Rise Azimuth: ', ra
                    print ' Transit Time: ', tt
                    print ' Transit Altitude: ', ta
                    print ' Set Time: ', st
                    print ' Set Azimuth: ', sa
                    UPCOMING_SWATH_FILENAME = "upcoming_orbit_swath.kml"
                    getVectorFile(attributes,polypoints,'polygon', UPCOMING_SWATH_FILENAME, 'KML')
                    #UPCOMING_POSITION_NOW_FILENAME = "upcoming_current_position"


        time.sleep(60*_statussleep)