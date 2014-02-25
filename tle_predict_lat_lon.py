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
_latlong = ('-23 42','133 54') # Alice Spring Data Acquisition Facility
_notify = 1440 # let us know this many minutes in advance to a pass
#_notify = 30
_usevoice = True # use voice?
_statussleep = 1 # how many minutes to sleep between status updates

# Earth parameters for heading calculations
oneonf = 298.257223563          # Inverse flattening 1/f = 298.257223563
f = 1/oneonf                    # flattening
r = 6378137
def GetTLEs():
    # GetTLEs(): returns a list of tuples of keps for each satellite.

    resource = 'http://www.celestrak.com/norad/elements/resource.txt'
    weather = 'http://www.celestrak.com/norad/elements/weather.txt'
    try:
        os.remove('resource.txt')
    except OSError:
        pass
    try:
        os.remove('weather.txt')

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

    print "retrieving TLE file.........."
    # strip off the header tokens and newlines
    tles = [item.strip() for item in tles]

    # clean up the lines
    tles = [(tles[i],tles[i+1],tles[i+2]) for i in xrange(0,len(tles)-2,3)]

    return tles
def getVectorFile(attributes,inputpoints,polyorline,ogroutput,ogrformat):
    #example usage: getVectorFile(dictionary,list of dicts with lat2 and lon2, 'polygon', SWATH_FILENAME, 'KML')

    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    # if no points passed for ogr build return
    if len(inputpoints) == 0:
        return ()
    try:
        os.remove(ogroutput)
    except OSError:
        pass
    ogr.UseExceptions()

    driver = ogr.GetDriverByName(ogrformat)

    if os.path.exists(ogroutput):
        driver.DeleteDataSource(ogroutput)
    ds = driver.CreateDataSource(ogroutput)

    if polyorline is 'polygon':
        geomtype = ogr.wkbPolygon
    if polyorline is 'line':
        geomtype = ogr.wkbLineString
    if polyorline is 'point':
        geomtype = ogr.wkbPoint

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
    fieldDefn = ogr.FieldDefn('Minutes to horizon           :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn('Acquisition of Signal UTC    :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn('Loss of Signal UTC           :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn('Transit time                 :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)
    fieldDefn = ogr.FieldDefn('Node                         :', ogr.OFTString)
    fieldDefn.SetWidth(30)
    layer.CreateField(fieldDefn)

    featureDefn = layer.GetLayerDefn()

    feature = ogr.Feature(featureDefn)

    feature.SetField('Satellite Name               :', attributes['Satellite name'])
    feature.SetField('Orbit height                 :', attributes['Orbit height'])
    feature.SetField('Orbit number                 :', attributes['Orbit'])
    feature.SetField('Current UTC time             :', str(attributes['Current time']))
    feature.SetField('Minutes to horizon           :', attributes['Minutes to horizon'])
    feature.SetField('Acquisition of Signal UTC    :', str(attributes['AOS time']))
    feature.SetField('Loss of Signal UTC           :', str(attributes['LOS time']))
    feature.SetField('Transit time                 :', str(attributes['Transit time']))
    feature.SetField('Node                         :', attributes['Node'])

    if polyorline == 'point':
        point = ogr.Geometry(ogr.wkbPoint)
        for x in inputpoints:
            point.AddPoint(x['lon2'],x['lat2'],x['alt2'])


        feature.SetGeometry(point)
        layer.CreateFeature(feature)

        point.Destroy()
    if polyorline == 'line':
        line = ogr.Geometry(type=ogr.wkbLineString)
        for x in inputpoints:
            line.AddPoint(x['lon2'],x['lat2'],x['alt2'])
            #print x
        feature.SetGeometry(line)
        layer.CreateFeature(feature)

        line.Destroy()

    if polyorline == 'polygon':
        ring = ogr.Geometry(ogr.wkbLinearRing)

        for x in inputpoints:
            ring.AddPoint(x['lon2'],x['lat2'])

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        feature.SetGeometry(poly)

        layer.CreateFeature(feature)

        ring.Destroy()
        poly.Destroy()

    feature.Destroy()

    ds.Destroy()

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
def getUpcomingPasses(satellitename,tleinfo,fromtime, totime):
    observer = ephem.Observer()
    observer.lat = _latlong[0]
    observer.long = _latlong[1]
    updatetime=0


    #Get most recent TLE for determining upcoming passes from now

    tles = tleinfo

    # make a list of dicts to hold the upcoming pass information for the selected satellites
    schedule = []
    observer.date = fromtime
    #_notify = totime
    while 1:

        #observer.date = datetime.utcnow()
        #print now
        print "---------------------------------------"
        # iterate through all the two line element sets
        satellites = ("LANDSAT 8", "LANDSAT 7", "TERRA", "AQUA", "NOAA 15", "NOAA 18", "NOAA 19", "SUOMI NPP")
        for tle in tles:
            #for i in satellites:
            #while not predictcomplete:
            index=0
            #if tle[0] in satellites:
            if tle[0] == satellitename:#in ("TERRA"):#, "LANDSAT 7", "TERRA", "AQUA", "NOAA 15", "NOAA 18", "NOAA 19", "SUOMI NPP"):
            #predicted = 0
            #index = 0
            #while predicted==0:
                #if tle[0]==satellites[index]:
                #for i in ("LANDSAT 8", "LANDSAT 7", "TERRA", "AQUA", "NOAA 15", "NOAA 18", "NOAA 19", "SUOMI NPP")

                #TODO clean up the use of ephem versus orbital. Orbital can give a orbit number and does many of the ephem functions
                #TODO enable user entered time range to display
                #TODO better configure the shapefile outputs to include 1.next passes for each satellite 2. current position 3.passes in date range
                #TODO add the individual acquisitions as layers in the same ogr output
                #TODO enable altitude in KML display for current position
                #TODO use an approporiate googleearth icon for satellites at a visible display resolution with a name tag and minutesaway
                satname = str(tle[0]).replace(" ","_")
                sat = ephem.readtle(tle[0],tle[1],tle[2])
                twole = tlefile.read(tle[0],'tles.txt')
                now = datetime.utcnow()
                print "---------------------------------------"
                print tle[0]
                oi = float(str.split(tle[2],' ')[3])

                orb = Orbital(tle[0])

                attributes = []

                #print "Orbit height         = ",orad
                #print "Orbit Number         = ",orb.get_orbit_number(now)
                #######
                #print "Orbital Inclination  = ",oi

                # need a time limit in here for rt i.e. if rt > 24hrs from nowtime overwriting - just update
                ###################
                #try:
                #    rt
                #except NameError:
                #    rt, ra, tt, ta, st, sa = observer.next_pass(sat)
                #else:
                #    print "rt exists"

                #if (updatetime==1):

                #    print schedule[len(schedule)-1]['LOS time']
                #    # use the last AOS time for a given satellite as the start time for the next pass


                #    seq = [x['LOS time'] for x in schedule]
                #    #min(seq)
                #    print max(seq)

                    #for x in reversed(schedule):
                    #    if (x['Satellite name'] == satname):
                    #        print x['Satellite name']
                    #        observer.date = datetime.strptime(x['LOS time'],"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5)
                    #        print "LOS time to use for new compute: ",x['LOS time']
                    #        print "LOS time adjusted for next rise: ",datetime.strptime(x['LOS time'],"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5)
                    #        #observer.date = datetime.strptime(str(schedule[len(schedule)-1]['AOS time']),"%Y/%m/%d %H:%M:%S")
                print "OBSERVER DATE:",observer.date
                #sat.compute(observer)

                    #if any ((x['Satellite name'] == satname) for x in reversed(schedule)):
                    #    observer.date = datetime.strptime(x['LOS time'],"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5)
                    #    #print "AOS time to use for new compute: ",x['AOS time']
                    #    #print "AOS time adjusted for next rise: ",datetime.strptime(x['AOS time'],"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5)
                    #    observer.date = datetime.strptime(str(schedule[len(schedule)-1]['LOS time']),"%Y/%m/%d %H:%M:%S")
                    #    sat.compute(observer)
                #print

                rt, ra, tt, ta, st, sa = observer.next_pass(sat)

                ##################

                #print "observer date ",rt
                ##### turn this off to disable computation of next passes - caused skipping of some passes????
                #observer.date = rt
                #sat.compute(observer)
                #########################
                localrisetime = ephem.localtime(rt)

                timeuntilrise = localrisetime-now
                print timeuntilrise

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
                    print "PASS                 = descending"
                    node = "descending"
                else:
                    print "PASS                 = ascending"
                    node = "ascending"
                    oi = 360 - oi

                orad = orb.get_lonlatalt(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S"))[2]
                attributes = {'Satellite name': satname,'Orbit height':orad,'Orbit':orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")),\
                              'Current time':str(now),'Minutes to horizon':minutesaway,'AOS time':str(rt),\
                              'LOS time':str(st),'Transit time':str(tt), 'Node':node}
                #if not any ((x['Satellite name'] == satname and x['Orbit'] == str(rt))for x in schedule):
                if not any ((x['Satellite name'] == satname and x['Orbit'] == orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))for x in schedule):
                    schedule.append(attributes)

                # Step from AOS to LOS in 100 second intervals
                #sat.compute(observer)
                delta = timedelta(seconds=100)

                deltatime = datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")
                geoeastpoint = []
                geowestpoint = []
                geotrack = []

                print "Current location     = ", sat.sublong.real*(180/math.pi), sat.sublat.real*(180/math.pi)
                print "DELTATIME",deltatime
                print "SETTING TIME",datetime.strptime(str(st),"%Y/%m/%d %H:%M:%S")

                while deltatime < datetime.strptime(str(st),"%Y/%m/%d %H:%M:%S"):

                    sat.compute(deltatime)
                    geotrack.append({'lat2':sat.sublat.real*(180/math.pi),'lon2':sat.sublong.real*(180/math.pi),'alt2':orb.get_lonlatalt(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S"))[2]})

                    eastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad)+90
                    westaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad)+270

                    #Set ground swath per satellite sensor
                    #TODO use view angle check to refine step from satellite track see IFOV
                    if tle[0] in ("LANDSAT 8","LANDSAT 7"):
                        swath = 185000/2
                    if tle[0] in ("TERRA","AQUA"):
                        swath = 2330000/2
                    if tle[0] in ("NOAA 15", "NOAA 18", "NOAA 19"):
                        swath = 1100000/2
                    if tle[0] == "SUOMI NPP":
                        swath = 2200000/2

                    geoeastpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, eastaz, swath))
                    geowestpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, westaz, swath))

                    deltatime = deltatime+delta

                # Create current location ogr output
                nowpoint = [{'lat2':orb.get_lonlatalt(datetime.utcnow())[1],'lon2':orb.get_lonlatalt(datetime.utcnow())[0],'alt2':orb.get_lonlatalt(datetime.utcnow())[2]}]

                CURRENT_POSITION_FILENAME = satname+"_current_position.kml"
                getVectorFile(attributes,nowpoint,'point', CURRENT_POSITION_FILENAME, 'KML')

                polypoints = []

                for x in geowestpoint:
                    polypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                for x in reversed(geoeastpoint):
                    polypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                if len(polypoints)>0:
                    polypoints.append({'lat2':geowestpoint[0]['lat2'],'lon2':geowestpoint[0]['lon2']})

                # Create swath footprint ogr output
                SWATH_FILENAME = satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_swath.kml"
                ORBIT_FILENAME = satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_track.kml"
                TRACKING_SWATH_FILENAME = satname+"_tracking_now.kml"

                # Create currently acquiring polygon
                #TODO def this
                # Step from AOS to current time second intervals

                observer.date=datetime.utcnow()
                sat.compute(observer)
                tkdelta = timedelta(seconds=100)
                tkrt, tkra, tktt, tkta, tkst, tksa = observer.next_pass(sat)
                tkdeltatime = datetime.utcnow()
                tkgeoeastpoint = []
                tkgeowestpoint = []
                tkgeotrack = []

                while tkdeltatime < (datetime.utcnow() or datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")):

                    sat.compute(tkdeltatime)
                    tkgeotrack.append({'lat2':sat.sublat.real*(180/math.pi),'lon2':sat.sublong.real*(180/math.pi),'alt2':orb.get_lonlatalt(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S"))[2]})

                    tkeastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad)+90
                    tkwestaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad)+270


                    tkgeoeastpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, eastaz, swath))
                    tkgeowestpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, westaz, swath))

                    tkdeltatime = tkdeltatime+tkdelta

                tkpolypoints = []

                for x in tkgeowestpoint:
                    tkpolypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                for x in reversed(tkgeoeastpoint):
                    tkpolypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                if len(tkpolypoints)>0:
                    tkpolypoints.append({'lat2':tkgeowestpoint[0]['lat2'],'lon2':tkgeowestpoint[0]['lon2']})

                if not ((attributes['Node']=="ascending")and(satname not in ("AQUA"))):
                    # Create swath ogr output
                    getVectorFile(attributes,polypoints,'polygon', SWATH_FILENAME, 'KML')
                    # Create orbit track ogr output
                    getVectorFile(attributes,geotrack,'line', ORBIT_FILENAME, 'KML')
                    # Create currently acquiring ogr output
                    if ((datetime.utcnow() >= datetime.strptime(str(tkrt),"%Y/%m/%d %H:%M:%S")) and (datetime.utcnow() <= datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S"))):
                        getVectorFile(attributes,tkpolypoints,'polygon', TRACKING_SWATH_FILENAME, 'KML')

                if minutesaway <= _notify:

                    print "---------------------------------------"
                    print tle[0], 'WILL BE MAKING A PASS IN ', minutesaway, " MINUTES"
                    print ' Rise Azimuth: ', ra
                    print ' Transit Time: ', tt
                    print ' Transit Altitude: ', ta
                    print ' Set Time: ', st
                    print ' Set Azimuth: ', sa

                    for x in sorted(schedule, key=lambda k: k['AOS time']):
                        print x

                    print (datetime.strptime(str(schedule[0]['AOS time']),"%Y/%m/%d %H:%M:%S"))

                    if ((datetime.strptime(str(schedule[len(schedule)-1]['AOS time']),"%Y/%m/%d %H:%M:%S")<(datetime.utcnow()+timedelta(minutes=_notify)))):
                        observer.date = (datetime.strptime(str(schedule[len(schedule)-1]['LOS time']),"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5))
                        sat.compute(observer)

                        print "MODIFIED OBSERVER DATE",observer.date
                    else:
                        print "--------------NOTHING TO MODIFY MOVING TO NEXT SATELLITE IN LIST------------"
                        return ()

        time.sleep(1*_statussleep)
    return ()

if __name__ == '__main__':
    tles = GetTLEs()
    # loop through satellite list and execute until whenever
    satellites = ("LANDSAT 8", "LANDSAT 7", "TERRA", "AQUA", "NOAA 15", "NOAA 18", "NOAA 19", "SUOMI NPP")
    while 1:
        for i in satellites:
            getUpcomingPasses(i,tles,datetime.utcnow(),_notify)
