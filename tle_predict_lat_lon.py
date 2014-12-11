from datetime import timedelta, datetime
import ephem
import math
import os
import sys
import time
#import urllib2
#import logging
import wget
from pyorbital.orbital import Orbital
from pyorbital import tlefile
from geographiclib.geodesic import Geodesic
import osgeo.ogr
import osgeo.osr
from osgeo import ogr
from osgeo import gdal
import numpy
import shutil
import HTML

#PyKML gumpf
#import lxml
#from lxml import etree
#import pykml


# Check for correct usage
if len(sys.argv)<3:
    print "*--------------------------------------------------------------------*"
    print ""
    print " tle_predict_lat_lon.py computes current position, observer track,    "
    print " and approximate imaging footprint of Earth Observation Satellites    "
    print ""
    print "*--------------------------------------------------------------------*"
    print ""
    print " usage: tle_predict_lat_lon.py <period to predict(mins)> <output path> "
    print ""
    print "*--------------------------------------------------------------------*"
    sys.exit()

# Read arguments

ground_station = ('-23 42', '133 54')   # Alice Spring Data Acquisition Facility
period = int(sys.argv[1])  # Generate passes for this time period from start time
output_path = sys.argv[2]
if not os.path.exists(output_path):
    print "OUTPUT PATH DOESN'T EXIST",output_path
    sys.exit()
sleep_status = 1   # how many minutes to sleep between status updates
schedule = []
# Earth parameters for heading calculations
one_on_f = 298.257223563          # Inverse flattening 1/f = 298.257223563
f = 1/one_on_f                    # flattening
r = 6378137

def get_tles():


    # GetTLEs(): returns a list of tuples of kepler parameters for each satellite.
    resource = 'http://www.celestrak.com/norad/elements/resource.txt'
    weather = 'http://www.celestrak.com/norad/elements/weather.txt'
    # Dove constellations    
    planetlabs = 'http://www.celestrak.com/NORAD/elements/stations.txt'
    try:
	os.remove('resource.txt')
    except OSError:
	pass
    try:
        os.remove('weather.txt')
    except OSError:
        pass
    try:
        os.remove('stations.txt')
    except OSError:
        pass

    try:
        
        wget.download(resource)
    except OSError:
        print "COULD NOT DOWNLOAD resource.txt"
        return ()
    try:
        
        wget.download(weather)
    except OSError:
        print "COULD NOT DOWNLOAD weather.txt"
        return ()
    try:
 
        wget.download(planetlabs)
    except OSError:
        print "COULD NOT DOWNLOAD stations.txt"
        return ()




    file_names = ['weather.txt', 'resource.txt', 'stations.txt']
    with open('tles.txt', 'w') as outfile:
        for fname in file_names:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    tles = open('tles.txt', 'r').readlines()

    print "retrieving TLE file.........."
    # strip off the header tokens and newlines
    tles = [item.strip() for item in tles]

    # clean up the lines
    tles = [(tles[i], tles[i+1], tles[i+2]) for i in xrange(0, len(tles)-2, 3)]

    return tles


def getVectorFile(attributes, input_points, poly_or_line, ogr_output, ogr_format):


    #example usage: getVectorFile(dictionary,list of dicts with lat2 and lon2, 'polygon', SWATH_FILENAME, 'KML')
    spatialReference = osgeo.osr.SpatialReference()
    spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    # if no points passed for ogr build return
    if len(input_points) == 0:
        return ()
    try:
        os.remove(ogr_output)
    except OSError:
        pass
    ogr.UseExceptions()

    driver = ogr.GetDriverByName(ogr_format)

    if os.path.exists(ogr_output):
        driver.DeleteDataSource(ogr_output)
    ds = driver.CreateDataSource(ogr_output)

    if poly_or_line is 'polygon':
        geomtype = ogr.wkbPolygon
    if poly_or_line is 'line':
        geomtype = ogr.wkbLineString
    if poly_or_line is 'point':
        geomtype = ogr.wkbPoint

    if ds is None:
        print 'Could not create file'
        sys.exit(1)
    layer = ds.CreateLayer(attributes['Satellite name'], geom_type=geomtype)

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

    if poly_or_line == 'point':
        point = ogr.Geometry(ogr.wkbPoint)
        for x in input_points:
            point.AddPoint(x['lon2'], x['lat2'], x['alt2'])

        feature.SetGeometry(point)
        layer.CreateFeature(feature)

        point.Destroy()
    if poly_or_line == 'line':
        line = ogr.Geometry(type=ogr.wkbLineString)
        for x in input_points:
            line.AddPoint(x['lon2'], x['lat2'], x['alt2'])
            #print x
        feature.SetGeometry(line)
        layer.CreateFeature(feature)

        line.Destroy()

    if poly_or_line == 'polygon':
        ring = ogr.Geometry(ogr.wkbLinearRing)

        for x in input_points:
            ring.AddPoint(x['lon2'], x['lat2'])

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        feature.SetGeometry(poly)

        layer.CreateFeature(feature)

        ring.Destroy()
        poly.Destroy()

    feature.Destroy()

    ds.Destroy()
    # Add altitude to KML if ogr_format=="KML" and change colour of track to yellow
    if ogr_format=="KML":
        if poly_or_line is 'line':
            replace_string_in_file(ogr_output,'<LineString>', '<LineString><altitudeMode>absolute</altitudeMode>')
            replace_string_in_file(ogr_output,'ff0000ff', 'ffffffff')
        if poly_or_line is 'point':
            replace_string_in_file(ogr_output,'<Point>', '<Point><altitudeMode>absolute</altitudeMode>')
        if poly_or_line is 'polygon':
            replace_string_in_file(ogr_output,'<PolyStyle><fill>0</fill>', '<PolyStyle><color>7f0000ff</color><fill>1</fill>')
        # TODO Group KML for a satellite into folders
            #import lxml
            #from lxml import etree
            #import pykml
            #from pykml.factory import KML_ElementMaker as KML
            #from pykml import parser

            #x = KML.Folder(KML.name("meow"))

            #with open("Scratch Paper.kml", "r+") as f:
            #    doc = parser.parse(f).getroot()
            #    print doc.Document.Folder.Folder[3].name
            #    a = doc.Document.Folder[0]
            #    a.append(x)
            #    finished = (etree.tostring(doc, pretty_print=True))

            #with open("Scratch Paper.kml", "w+") as f:
            #    f.write(finished)

            #print "Done!"
    return ()

def replace_string_in_file(infile,text_to_find,text_to_insert):


    in_file = open(infile, 'r')
    temporary = open(os.path.join(output_path,'tmp.txt'), 'w')
    for line in in_file:
        temporary.write(line.replace(text_to_find, text_to_insert))
    in_file.close()
    temporary.close()
    os.remove(infile)
    shutil.move(os.path.join(output_path,'tmp.txt'),infile)
    return ()


def getEffectiveHeading(satellite, oi_deg, latitude, longitude, tle_orbit_radius, daily_revolutions):

    lat_rad = math.radians(latitude)  # Latitude in radians
    oi_rad = math.radians(oi_deg)   # Orbital Inclination (OI) [radians]
    orbit_radius = tle_orbit_radius*1000.0 # Orbit Radius (R) [m]
    #np = 5925.816                   # Nodal Period [sec] = 5925.816
    np = (24*60*60)/daily_revolutions
    av = 2*math.pi/np          # Angular Velocity (V0) [rad/sec] =	 0.001060307189285 =2*PI()/E8
    sr = 0                          # Sensor Roll (r) [degrees] =	0

    #TODO put earth parameters into a dict and add support for other spheroids GRS1980 etc.
    # Earth Stuff (WGS84)
    one_on_f = 298.257223563          # Inverse flattening 1/f = 298.257223563
    f = 1/one_on_f                    # flattening
    r = 6378137                     # Radius (a) [m] =	 6378137
    e = 1-math.pow((1-1/one_on_f), 2)  # Eccentricity (e^2) = 0.00669438 =1-(1-1/I5)^2
    wO = 0.000072722052             # rotation (w0) [rad/sec] = 7.2722052E-05

    xfac = math.sqrt(1-e*(2-e)*(math.pow(math.sin(math.radians(latitude)), 2)))
    phi_rad = math.asin((1-e)*math.sin(math.radians(latitude))/xfac)  # Phi0' (Geocentric latitude)
    phi_deg = math.degrees(phi_rad)  # Phi0' (Degrees)
    n = r/math.sqrt(1-e*(math.pow(math.sin(math.radians(latitude)), 2))) # N
    altphi_rad = latitude-180*math.asin(n*e*math.sin(lat_rad)*math.cos(lat_rad)/orbit_radius)/math.pi  # Alt Phi0'(Radians)
    rho_rad = math.acos(math.sin(altphi_rad*math.pi/180)/math.sin(oi_rad))  # Rho (Radians)
    beta = -1*(math.atan(1/(math.tan(oi_rad)*math.sin(rho_rad)))*180/math.pi)  # Heading Beta (degrees)
    xn = n*xfac  #  Xn
    altitude = (orbit_radius-xn)/1000  # altitude
    altitude_ = (orbit_radius*math.cos(altphi_rad/180*math.pi)/math.cos(lat_rad)-n)/1000
    rotation = math.atan((wO*math.cos(phi_rad)*math.cos(beta*math.pi/180))/(av+wO*math.cos(phi_rad)*math.sin(beta*math.pi/180)))*180/math.pi
    eh = beta+rotation
    alpha12 = eh
    s = 0.5*185000  # s = distance in metres
    effective_heading = alpha12
    return effective_heading

def getUpcomingPasses(satellite_name, tle_information, passes_begin_time, passes_period):


    observer = ephem.Observer()
    observer.lat = ground_station[0]
    observer.long = ground_station[1]
    observer.horizon = '5:0'
    #updatetime = 0
    period = passes_period
    #Get most recent TLE for determining upcoming passes from now
    tles = tle_information

    # make a list of dicts to hold the upcoming pass information for the selected satellites
    schedule = []
    observer.date = passes_begin_time

    while 1:

        print "---------------------------------------"
        for tle in tles:
            
            if tle[0] == satellite_name:
                #print tle
                #TODO clean up the use of pyephem versus orbital. Orbital can give a orbit number and does many of the pyephem functions
                #TODO add the individual acquisitions as layers in the same ogr output
                #TODO use an appropriate google earth icon for satellites at a visible display resolution with a name tag and minutesaway
                #TODO print output to logging
                satname = str(tle[0]).replace(" ","_")
                # Flock has minus in filename but looks like KML creater doesn't like it                
                satname = satname.replace("-","_")
                #print satname
                
                sat = ephem.readtle(tle[0],tle[1],tle[2])


                twole = tlefile.read(tle[0],'tles.txt')
                now = datetime.utcnow()
                #TODO check age of TLE - if older than x days get_tle()
                print "TLE EPOCH:",twole.epoch
                #if twole.epoch < now - timedelta(days=5):
                #    get_tles()
                #    satname = str(tle[0]).replace(" ","_")
                #    sat = ephem.readtle(tle[0],tle[1],tle[2])
                #    twole = tlefile.read(tle[0],'tles.txt')

                print "---------------------------------------"
                print tle[0]

                oi = float(str.split(tle[2],' ')[3])
                #orb = Orbital(tle[0])
                orb = Orbital(tle[0],"tles.txt", tle[1],tle[2])
                attributes = []

                rt, ra, tt, ta, st, sa = observer.next_pass(sat)
                
                # Determine is pass descending or ascending
                # Confirm that observer details have been computed i.e. are not 'Null'
		
		if rt is None:
		    return ()
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

                AOStime = datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")
                minutesaway = ((AOStime-now).seconds/60.0)+((AOStime-now).days*1440.0)

                print "Minutes to horizon   = ", minutesaway
                print "AOStime              = ", rt
                print "LOStime              = ", st
                print "Transit time         = ", tt

                orad = orb.get_lonlatalt(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S"))[2]
                
		# Create swath footprint ogr output
                SWATH_FILENAME = os.path.join(output_path,satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_swath.kml")
                attributes = {'Satellite name': satname, 'Orbit height': orad, 'Orbit': orb.get_orbit_number(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")), \
                              'Current time': str(now),'Minutes to horizon': minutesaway, 'AOS time': str(rt), \
                              'LOS time': str(st), 'Transit time': str(tt), 'Node': node, \
			      'SWATH_FILENAME': (satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_swath.kml"),'ORBIT_FILENAME': (satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_track.kml")}

                # Append the attributes to the list of acquisitions for the acquisition period
                if not any ((x['Satellite name'] == satname and x['Orbit'] == orb.get_orbit_number(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")))for x in schedule):
                    schedule.append(attributes)
                
                
                # Step from AOS to LOS in 100 second intervals
                delta = timedelta(seconds=100)
                deltatime = datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")

                geoeastpoint = []
                geowestpoint = []
                geotrack = []


                print "DELTATIME", deltatime
                print "SETTING TIME", datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S")

                while deltatime < datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S"):
                    print "delta time is less than satellite LOS time"
                    sat.compute(deltatime)
                   
                    geotrack.append({'lat2': sat.sublat.real*(180/math.pi), \
                                     'lon2': sat.sublong.real*(180/math.pi), \
                                     'alt2': orb.get_lonlatalt(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S"))[2]*1000})
                                        
                    eastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi), orad, sat._n)+90
                    westaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi), orad, sat._n)+270

                    #Set ground swath per satellite sensor
                    #TODO use view angle check to refine step from satellite track see IFOV
                    if tle[0] in ("LANDSAT 8","LANDSAT 7"):
                        swath = 185000/2
                    if tle[0] in ("TERRA","AQUA"):
                        swath = 2330000/2
                    if tle[0] in ("NOAA 15", "NOAA 18", "NOAA 19"):
                        swath = 2399000/2
                    if tle[0] == "SUOMI NPP":
                        swath = 2200000/2
                    else:
                        swath = 14000/2

                    geoeastpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, eastaz, swath))
                    geowestpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, westaz, swath))

                    deltatime = deltatime+delta

                # Create current location ogr output
                                
                nowpoint = [{'lat2':orb.get_lonlatalt(datetime.utcnow())[1],'lon2':orb.get_lonlatalt(datetime.utcnow())[0],'alt2':orb.get_lonlatalt(datetime.utcnow())[2]*1000}]
                               
                #TODO ensure the now attributes are actually attributes for the current position of the satellite and include relevant next pass information...tricky?
                #if ((attributes['Orbit']==orb.get_orbit_number(datetime.utcnow()))and(AOStime<now)):
                now_attributes = {'Satellite name': satname, 'Orbit height': orb.get_lonlatalt(datetime.utcnow())[2], 'Orbit': orb.get_orbit_number(datetime.utcnow()), \
                          'Current time': str(now),'Minutes to horizon': "N/A", 'AOS time': "N/A", \
                          'LOS time': "N/A", 'Transit time': "N/A", 'Node': "N/A"}
                    #now_attributes=attributes
                
                CURRENT_POSITION_FILENAME = os.path.join(output_path,satname+"_current_position.kml")

                #TODO draw the current orbit forward for the passes period time from the satellite position as a long stepped ogr line
                print now_attributes,nowpoint
                getVectorFile(now_attributes,nowpoint,'point', CURRENT_POSITION_FILENAME, 'KML')

                polypoints = []

                for x in geowestpoint:
                    polypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                for x in reversed(geoeastpoint):
                    polypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                if len(polypoints)>0:
                    polypoints.append({'lat2':geowestpoint[0]['lat2'],'lon2':geowestpoint[0]['lon2']})

                # Create swath footprint ogr output
                SWATH_FILENAME = os.path.join(output_path,satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_swath.kml")
                ORBIT_FILENAME = os.path.join(output_path,satname+"."+str(orb.get_orbit_number(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S")))+".ALICE.orbit_track.kml")
                TRACKING_SWATH_FILENAME = os.path.join(output_path,satname+"_tracking_now.kml")

                # Create currently acquiring polygon
                #TODO def this
                # Step from AOS to current time second intervals

                observer.date=now
                sat.compute(observer)
                tkdelta = timedelta(seconds=100)
                #TODO Problem determining rise time if rise time has passed and set time not yet reached
		#solution may be to replace rt with now if rt > st
		tkrt, tkra, tktt, tkta, tkst, tksa = observer.next_pass(sat)
		print tkrt, tkra, tktt, tkta, tkst, tksa
		if tkrt is None:
		    return ()		
		#if datetime.strptime(str(tkrt),"%Y/%m/%d %H:%M:%S") > datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S"):
		#    tkrt = datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S") - datetime.strptime(str(tktt),"%Y/%m/%d %H:%M:%S")
                #    tkrt = str(tkt),"%Y/%m/%d %H:%M:%S"
		#print "NOW: ",now," TKRT: ",datetime.strptime(str(tkrt),"%Y/%m/%d %H:%M:%S")," TKST: ",datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")
		tkdeltatime = datetime.utcnow()
                tkgeoeastpoint = []
                tkgeowestpoint = []
                tkgeotrack = []

                #while tkdeltatime < (datetime.utcnow() or datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")):
		while (tkdeltatime < datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")):

                    sat.compute(tkdeltatime)
                    tkgeotrack.append({'lat2':sat.sublat.real*(180/math.pi),'lon2':sat.sublong.real*(180/math.pi),'alt2':orb.get_lonlatalt(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S"))[2]})

                    tkeastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad,sat._n)+90
                    tkwestaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad,sat._n)+270
                    #TODO use view angle check to refine step from satellite track see IFOV
                    if tle[0] in ("LANDSAT 8","LANDSAT 7"):
                        tkswath = 185000/2
                    if tle[0] in ("TERRA","AQUA"):
                        tkswath = 2330000/2
                    if tle[0] in ("NOAA 15", "NOAA 18", "NOAA 19"):
                        tkswath = 1100000/2
                    if tle[0] == "SUOMI NPP":
                        tkswath = 2200000/2
                    else:
                        tkswath = 14000/2
                    tkgeoeastpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, tkeastaz, tkswath))
                    tkgeowestpoint.append(Geodesic.WGS84.Direct(sat.sublat.real*180/math.pi, sat.sublong.real*180/math.pi, tkwestaz, tkswath))

                    tkdeltatime = tkdeltatime+tkdelta

                tkpolypoints = []

                for x in tkgeowestpoint:
                    tkpolypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                for x in reversed(tkgeoeastpoint):
                    tkpolypoints.append({'lat2':x['lat2'],'lon2':x['lon2']})
                if len(tkpolypoints)>0:
                    tkpolypoints.append({'lat2':tkgeowestpoint[0]['lat2'],'lon2':tkgeowestpoint[0]['lon2']})

                    #if not ((attributes['Node']=="ascending")and(satname not in ("AQUA"))):
                    if (((attributes['Node']=="ascending")and(satname in ("AQUA","SUOMI_NPP")))or ((attributes['Node']=="descending")and(satname not in ("AQUA","SUOMI_NPP")))): 
   			# Create swath ogr output
                    	getVectorFile(attributes,polypoints,'polygon', SWATH_FILENAME, 'KML')
                    	# Create orbit track ogr output
                    	getVectorFile(attributes,geotrack,'line', ORBIT_FILENAME, 'KML')
                    	# Create currently acquiring ogr output
			#print "NOW: ",now," TKRT: ",datetime.strptime(str(tkrt),"%Y/%m/%d %H:%M:%S")," TKST: ",datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")
                    	#if ((datetime.strptime(str(tkrt),"%Y/%m/%d %H:%M:%S")>now) and (datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")<now)):
			if tkrt > tkst:
			    print "Executing tracking swath creation - tkpolypoints = ",tkpolypoints
                       	    getVectorFile(now_attributes,tkpolypoints,'polygon', TRACKING_SWATH_FILENAME, 'KML')




		baseline=1
                if minutesaway <= period:

                    print "---------------------------------------"
                    print tle[0], 'WILL BE MAKING A PASS IN ', minutesaway, " MINUTES"
                    print ' Rise Azimuth: ', ra
                    print ' Transit Time: ', tt
                    print ' Transit Altitude: ', ta
                    print ' Set Time: ', st
                    print ' Set Azimuth: ', sa

                    for x in sorted(schedule, key=lambda k: k['AOS time']):
                        print x
			
                        # For dictionary entries with 'LOS time' older than now time - remove
                        if ((datetime.strptime(str(x['LOS time']),"%Y/%m/%d %H:%M:%S"))<(datetime.utcnow())):
                            # Delete output ogr
                            if os.path.exists(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_swath.kml")):
                                shutil.move(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_swath.kml"),os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_swath.kml.OUTOFDATE"))
                            if os.path.exists(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_track.kml")):
                                shutil.move(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_track.kml"),os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_track.kml.OUTOFDATE"))

                            # Delete dictionary entry for pass
                            schedule.remove(x)

                    # Unlikely - if no entries in the schedule don't try to print it
                    # see if there are any new additions to the schedule
		    		    
		    if len(schedule)>0:
                        print (datetime.strptime(str(schedule[0]['AOS time']),"%Y/%m/%d %H:%M:%S"))
			
                    # If the AOS time is less than now + the time delta, shift the time to the latest recorded pass LOS time
                    if (len(schedule)>0 and ((datetime.strptime(str(schedule[len(schedule)-1]['AOS time']),"%Y/%m/%d %H:%M:%S")<(now+timedelta(minutes=period))))):
    			observer.date = (datetime.strptime(str(schedule[len(schedule)-1]['LOS time']),"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5))

                        # Recompute the satellite position for the update time
                        sat.compute(observer)
                        print "MODIFIED OBSERVER DATE",observer.date

                    else:
                        print "--------NOTHING TO MODIFY MOVING TO NEXT SATELLITE IN LIST------"
                        #TODO - write to html
                        # Exit the def if the schedule isn't able to update because there are no passes in the acquisition window
			
                        output_html = HTML.Table(schedule)
                        t = HTML.Table(header_row=['Transit time', 'Node', 'AOS time', 'Current time', 'Satellite name', 'Minutes to horizon', 'LOS time', 'Orbit height', 'Orbit'])
                        html_output = open(os.path.join(output_path,satname+".schedule.html"),'w')
			html_output.write("<!DOCTYPE html>"+'\n')
                        html_output.write("<html>"+'\n')
                        html_output.write("<head>"+'\n')
                        html_output.write('<meta http-equiv="refresh" content="5">'+'\n')
                        html_output.write("<style>"+'\n')
                        html_output.write("table,th,td"+'\n')
                        html_output.write("{"+'\n')
                        html_output.write("border:1px solid black;"+'\n')
                        html_output.write("padding:15px;"+'\n')
                        html_output.write("}"+'\n')
                        html_output.write("</style>"+'\n')
                        html_output.write("<body>"+'\n')
                        html_output.write('<table style="width:300px">'+'\n')

                        html_output.write("<tr>"+'\n')
                        html_output.write("    <th>Satellite</th>"+'\n')
                        html_output.write("    <th>Orbit</th>"+'\n')
                        html_output.write("    <th>Node</th>"+'\n')
                        html_output.write("    <th>AOS time</th>"+'\n')
                        html_output.write("    <th>LOS time</th>"+'\n')
                        html_output.write("    <th>Minutes to horizon</th>"+'\n')
                        html_output.write("</tr>"+'\n')
                        for x in sorted(schedule):
                            if (((x['Node']=="ascending")and(x['Satellite name'] in ("AQUA","SUOMI_NPP"))) or ((x['Node']=="descending")and(x['Satellite name'] not in ("AQUA","SUOMI_NPP")))):
				html_output.write("<tr>"+'\n')
                            	html_output.write("    <td>"+str(x['Satellite name'])+"</td>"+'\n')
                            	html_output.write('    <td><a href="'+str(x['SWATH_FILENAME'])+'">'+str(x['Orbit'])+"</a></td>"+'\n')
                            	html_output.write('    <td><a href="'+str(x['ORBIT_FILENAME'])+'">'+str(x['Node'])+"</a></td>"+'\n')
                            	html_output.write("    <td>"+str(x['AOS time'])+"</td>"+'\n')
                            	html_output.write("    <td>"+str(x['LOS time'])+"</td>"+'\n')
                            	html_output.write('    <td><a href="'+satname+"_current_position.kml"+'">'+str(x['Minutes to horizon'])+"</a></td>"+'\n')
                            	html_output.write("</tr>"+'\n')
                        html_output.write("</body>"+'\n')
                        html_output.write("</html>"+'\n')
                        return ()

                #    return ()

        time.sleep(1*sleep_status)
    return ()

if __name__ == '__main__':
    tles = get_tles()
    # Loop through satellite list and execute until end of period

    satellites = ("LANDSAT 7", "TERRA", "AQUA", "NOAA 15", "NOAA 18", "NOAA 19", "SUOMI NPP")
    #TODO FLOCKS look to be acquiring on ascending pass but the code labels this descending - a
    #TODO Add support for off nadir acquisitions i.e. standard off nadir Sentinel1
    #satellites = ('SENTINEL-1A','FLOCK 1B-24','FLOCK 1B-23','FLOCK 1B-26','FLOCK 1B-25','FLOCK 1B-15','FLOCK 1B-16','FLOCK 1B-1','FLOCK 1B-2','FLOCK 1B-8','FLOCK 1B-7','FLOCK 1B-18','FLOCK 1B-17')

    while 1:
        for i in satellites:
            print "Looking for ",i
            getUpcomingPasses(i,tles,datetime.utcnow(),period)

