# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:21:04 2014

@author: user
"""

# Script to sit and poll for upcoming satellite passes.
# Joseph Armbruster

# TODO(joe) : clean up these imports

#import datetime
from datetime import timedelta, date, datetime
import ephem
import math
import os
import sys
import time
import urllib2
import logging
import wget

#from datetime import datetime, timedelta
import math
import osgeo.ogr, osgeo.osr
from osgeo import ogr
from osgeo import gdal
import numpy

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
    # urllib2 open doesn't like something about the celestrak text file here's a crumby workaround...revisit
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
     
    #wget.download(resource)
    #wget.download(weather)    
    #filenames = ['weather.txt', 'resource.txt']
    #with open('tles.txt', 'w') as outfile:
    #    for fname in filenames:
    #        with open(fname) as infile:
    #            for line in infile:
    #                outfile.write(line)

    tles = open('tles.txt','r').readlines()
    #tlesweather = open('weather.txt','r').readlines()
    #tles = tlesresource,tlesweather    
    
    print "retrieving TLE file.........."
    #print tles
    # strip off the header tokens and newlines
    tles = [item.strip() for item in tles]
    
    # clean up the lines
    tles = [(tles[i],tles[i+1],tles[i+2]) for i in xrange(0,len(tles)-2,3)]
      
    return tles

if __name__ == '__main__':

    observer = ephem.Observer()
    observer.lat = _latlong[0]
    observer.long = _latlong[1]
    
    tles = GetTLEs()

    while 1:
        now = datetime.now()
        print "---------------------------------------"
        # iterate through all the two line element sets
        for tle in tles:
            if tle[0] in ("LANDSAT 8", "LANDSAT 7", "TERRA", "AQUA", "NOAA 18", "NOAA 19", "SUOMI NPP"):   
                sat = ephem.readtle(tle[0],tle[1],tle[2])
                print "---------------------------------------"
                print tle[0]
                rt, ra, tt, ta, st, sa = observer.next_pass(sat)
                observer.date = rt
                sat.compute(observer)
                localrisetime = ephem.localtime(rt)
                
                timeuntilrise = localrisetime-now
                
                minutesaway = timeuntilrise.seconds/60.0
            
                print "Minutes to horizon   = ", minutesaway
                print "AOStime              = ", rt
                print "LOStime              = ", st
                print "Transit time         = ", tt

                ######## Create shapefile for each of the upcoming passes ###################

                #SHP_FILENAME = str(tle[0]).replace(" ","_")+"_orbit_centre.shp"
                SHP_FILENAME = str(tle[0]).replace(" ","_")+"_orbit_centre.kml"
                print "Writing output to    = ", SHP_FILENAME
                ogr.UseExceptions()
                #driver = ogr.GetDriverByName('ESRI Shapefile')
                driver = ogr.GetDriverByName('KML')
                spatialReference = osgeo.osr.SpatialReference()
                spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
                
                try:
                    os.remove(SHP_FILENAME)
                except OSError:
                    pass 
                 
                ds = driver.CreateDataSource(SHP_FILENAME)
                layer = ds.CreateLayer("data", spatialReference, ogr.wkbLineString)
                
                # Create a new line geometry
                line = ogr.Geometry(type=ogr.wkbLineString)
                
                # Step from AOS to LOS in 5 second intervals
                #deltatime = rt
                delta = timedelta(seconds=1)

                deltatime = datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S") 
                
                print "Current location     = ", sat.sublong.real*(180/math.pi), sat.sublat.real*(180/math.pi) 
                while deltatime < datetime.strptime(str(st),"%Y/%m/%d %H:%M:%S"):

                    sat.compute(deltatime)
                    
                    #print sat.sublong.real*(180/math.pi), sat.sublat.real*(180/math.pi)
                    line.AddPoint(sat.sublong.real*(180/math.pi), sat.sublat.real*(180/math.pi))

                    deltatime = deltatime+delta
     
                # Add line as a new feature to the shapefile
                
             
                feature = ogr.Feature(feature_def=layer.GetLayerDefn())
                feature.SetGeometryDirectly(line)
                layer.CreateFeature(feature)
                # Create a buffer around the path the width of the swath. A quick and dirty to be replace with heading + skew calc - feature.Buffer(line,92500,40)
                feature.Destroy()
                ds.Destroy()
                
                ####### Shapefile Creation Complete #########################################    
                
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
            #else:
            #    print tle[0], " being ignored"
        time.sleep(60*_statussleep)
