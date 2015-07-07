
#----------------------------------------------------------------------------------
#
# Satellite Foot Print code
# Rev:     June 1, 2014
#
#----------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------
# numpy is not used in this cidei
# It has been temporarily commneted
#----------------------------------------------------------------------------------

import numpy

#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
#   Global Variable and satellite footprint parameters 
#----------------------------------------------------------------------------------

DATA_IN_DIR = "../data/"
OUTPUT_DIR  = "../output/"
LOGS_DIR    = "../logs/"

scp_parameters = "scp_parameters.txt"
scp_log        = "scp_log.txt"

orbit_parameters ="orbit_parameters.txt"

FL_ORB = open(OUTPUT_DIR+orbit_parameters,'a')

file_parameters = DATA_IN_DIR + scp_parameters

LAT_PARAMETERS =''
LONG_PARAMATERS=''

GRAVITATIONAL_CONST=0.0
RADIUS_OF_EARTH=0
GROUND_STATION=''

WO = 0.0


DELTA_TIME_STEP  = 100 #seconds
SLEEP_STATUS     = 1   # how many minutes to sleep between status updates

SCHEDULE         = []
SATELLITE_SWATH  = []
PROCESSED_ORBITS = []

#----------------------------------------------------------------------------------
# Satellite AOS and LOS window
#----------------------------------------------------------------------------------

SAT_AOS_WIN	 = []
#SAT_INTF_LIST    = ["AQUA","TERRA","NOAA_19","SUOMI_NPP"]
SAT_INTF_LIST    = ["NOAA_19","SUOMI_NPP"]
#SAT_INTF_LIST    = ["AQUA","TERRA"]

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

# Check for correct usage

# print len(sys.argv)

if len(sys.argv)<2:
    print "*--------------------------------------------------------------------*"
    print ""
    print " tle_predict_lat_lon.py computes current position, observer track,    "
    print " and approximate imaging footprint of Earth Observation Satellites    "
    print ""
    print "*--------------------------------------------------------------------*"
    print ""
    print " usage: tle_predict_lat_lon.py <period to predict(mins)> "
    print ""
    print "*--------------------------------------------------------------------*"

# Read arguments

period = int(sys.argv[1])  # Generate passes for this time period from start time

if len(sys.argv) == 2:
    output_path=OUTPUT_DIR
else:
    output_path = sys.argv[2]

if not os.path.exists(output_path):
    print "OUTPUT PATH DOESN'T EXIST",output_path
    sys.exit()


#----------------------------------------------------------------------------------
# This routine is for testing Azimith and elevation calculations test
#----------------------------------------------------------------------------------

def test_az_el_cal(orb,observer,rt, ra, tt, ta, st, sa, sat_ht):
	print 'Earth Station  Observer long       = ',observer.long
	print 'Earth Station Observer lat         = ',observer.lat
	print '-------------------get_obser_look method starts ------------------------------------------'
	(Azm, Ele) = orb.get_observer_look(datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S"), 133.9, -23.7,0.545)
	print 'Transit  Time          = ',datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S")
	print "Transit  Azimith       = ",Azm
	print "Transit  Elevation     = ",Ele
	print "------------------ get_observer look finishes -----------------------------------"
		
#		(Azm, Ele) = orb.get_observer_look(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S"), 133.0*3/142857/180,-23.0*3.142857/180,0)
#		print 'Start Time          = ',datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S")
#		print "Start Azimith       = ",Azm
#		print "Start Elevation     = ",Ele
		
#		(Azm, Ele) = orb.get_observer_look(datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S"), observer.long, observer.lat,576)
#		print 'End Time  = ',datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S")
#		print "End Azimith       = ",Azm
#		print "End Elevation     = ",Ele

	(lon,lat,alt) = orb.get_lonlatalt(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S"))
	print 'Sat lon  = ', lon
	print 'sat lat  = ', lat
	print 'Sat alt  = ', alt
	(az,el) = get_site_antenna_az_el(133.9, -23.7, lon, lat,sat_ht)
	print 'Satellite starting az = ',az
	print 'Satellite starting el = ',el

	(lon,lat,alt) = orb.get_lonlatalt(datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S"))
	print 'Sat lon  = ', lon
	print 'Sat lat  = ', lat
	print 'Sat alt  = ', alt
	(az,el) = get_site_antenna_az_el(133.9, -23.7, lon, lat, sat_ht)
	print 'Satellite Ending az = ',az
	print 'Satellite Ending el = ',el
	print 'alt  = ', alt
	print '----------------------------------------------------------'

	


	return

#----------------------------------------------------------------------------------
# Computing Azimiht and Elevation angles t ground station, 
# using longitude, latitude parematers of satellite and earth station
#----------------------------------------------------------------------------------

def get_site_antenna_az_el(site_lon, site_lat, sat_lon, sat_lat, sat_ht):
	
	rad_site_lon = math.radians(site_lon)
	rad_site_lat = math.radians(site_lat)
	rad_sat_lon  = math.radians(sat_lon)
	rad_sat_lat  = math.radians(sat_lat)
	rad_g	     = math.radians(site_lon - sat_lon)

	a =  math.cos(rad_g)
	b =  math.cos(rad_site_lat)
	

#	el = math.atan((a * b - 0.1512)/ \
	el = math.atan((a * b - RADIUS_OF_EARTH/(RADIUS_OF_EARTH + sat_ht))/ \
		math.sqrt(1 - (a*a)*(b*b)))

	el = math.degrees(el)

	az = math.radians(180) + math.atan(math.tan(rad_g)/math.sin(rad_site_lat))
	az = math.degrees(az)	

	return (az,el)

#----------------------------------------------------------------------------------
# Earth parameters for heading calculations
#----------------------------------------------------------------------------------

def get_parameters():

    global GRAVITATIONAL_CONST
    global RADIUS_OF_EARTH
    global LONG_PARAMATERS
    global LAT_PARAMETERS
    global WO

    fl = open(file_parameters)
    count = 0
    for item in fl:
	data_str = []
	if  item[0] <> "#" and len(item) > 5:
		if count == 0:
			data_str = (item.strip()).split(',')
			LAT_PARAMETERS=data_str[0]
			LONG_PARAMATERS=data_str[1]
		if count == 1:
			data_str = (item.strip()).split(',')
			GRAVITATIONAL_CONST = float(data_str[0])
			RADIUS_OF_EARTH = float(data_str[1])	
		if count == 2:
			WO = float(item)
		if count >= 3:
			SATELLITE_SWATH.append(item.strip())
		       	
		count += 1
    fl.close()

    return


def get_time_secs(dt_str):

	str1 = dt_str.split(' ')
	dstr = str1[0].split('/')
	tstr = str1[1].split(':')
	
	tsec = datetime(int(dstr[0]),int(dstr[1]),int(dstr[2]),\
		    int(tstr[0]),int(tstr[1]),int(tstr[2])).strftime('%s')
	return str(tsec)


def set_parameters():

    global GROUND_STATION

    GROUND_STATION = (LAT_PARAMETERS, LONG_PARAMATERS)   # Alice Spring Data Acquisition Facility

    return	

def output_orbit_parameters(orbit_x):

	if not ( orbit_x['Orbit'] in PROCESSED_ORBITS ):
		PROCESSED_ORBITS.append(orbit_x['Orbit'])
#		print orbit_x
		print
		print 'Satellite Name                                  = ', orbit_x['Satellite name']
		print 'Orbit Number                                    = ', orbit_x['Orbit']
		print 'Transit Time UTC                                = ', orbit_x['Transit time'],' Epoch(Secs) : ',get_time_secs(orbit_x['Transit time'])
		print 'Node                                            = ', orbit_x['Node']
		print 'Orbit Height                                    = ', orbit_x['Orbit height']
		print 'AOS ( Aquisition of Signal ) UTC                = ', orbit_x['AOS time'],' Epoch(Secs) : ',get_time_secs(orbit_x['AOS time'])
		print 'LOS time  ( Loss of Signal ) UTC                = ', orbit_x['LOS time'],' Epoch(Secs) : ',get_time_secs(orbit_x['LOS time'])
		print 'Data Aquisition at Ground Station ( Secs )      = ', \
			( int(get_time_secs(orbit_x['LOS time'])) - int(get_time_secs(orbit_x['AOS time'])) )
		print 'Current UTC Time                                = ', orbit_x['Current time']
		print 'Minutes to Horizon                              = ', orbit_x['Minutes to horizon']
		
		orb_str = orbit_x['Transit time'] + ',' + \
			  orbit_x['Node'] + ',' + \
			  orbit_x['AOS time'] + ',' + \
			  orbit_x['Current time']+ ',' + \
			  orbit_x['Satellite name']  + ',' + \
			  str(orbit_x['Minutes to horizon'])  + ',' + \
			  orbit_x['LOS time']  + ',' + \
			  str(orbit_x['Orbit height']) + ',' + \
			  str(orbit_x['Orbit']) + '\n' 
		FL_ORB.write(orb_str)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#	Updating AQUA and NOAAA 91 List
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

		if any(orbit_x['Satellite name'] in s for s in SAT_INTF_LIST):
			str1 = get_time_secs(orbit_x['AOS time']) + ' : ' + get_time_secs(orbit_x['LOS time']) + \
                                ' : ' + str(orbit_x['Orbit']) + ' : ' + get_time_secs(orbit_x['AOS time']) + ' : ' + \
                                str( int(get_time_secs(orbit_x['LOS time'])) - int(get_time_secs(orbit_x['AOS time'])) ) + \
                                ' : ' + orbit_x['Satellite name']
			SAT_AOS_WIN.append(str1.strip())
			
		

	
#	print orbit_x
	return


def get_tles():


    # GetTLEs(): returns a list of tuples of kepler parameters for each satellite.
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
    print  '\n Download Finished'

    cmd = 'mv ' + 'resource.txt' + ' ' + DATA_IN_DIR+'resource.txt'
    os.system(cmd)
    cmd = 'mv ' + 'weather.txt' + ' ' + DATA_IN_DIR+'weather.txt'
    os.system(cmd)

    file_names = [DATA_IN_DIR+'weather.txt', DATA_IN_DIR+'resource.txt']
    with open(DATA_IN_DIR+'tles.txt', 'w') as outfile:
        for fname in file_names:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    tles = open(DATA_IN_DIR+'tles.txt', 'r').readlines()

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
    temporary = open('tmp.txt', 'w')
    for line in in_file:
        temporary.write(line.replace(text_to_find, text_to_insert))
    in_file.close()
    temporary.close()
    os.remove(infile)
    os.rename('tmp.txt',infile)

    return ()


def getEffectiveHeading(satellite, oi_deg, latitude, longitude, tle_orbit_radius, daily_revolutions):

    #print "RADius",orbit_radius

    lat_rad = math.radians(latitude)  # Latitude in radians
    oi_rad = math.radians(oi_deg)   # Orbital Inclination (OI) [radians]
    orbit_radius = tle_orbit_radius*1000.0 # Orbit Radius (R) [m]
    #np = 5925.816                   # Nodal Period [sec] = 5925.816
    np = (24*60*60)/daily_revolutions
    av = 2*math.pi/np          # Angular Velocity (V0) [rad/sec] =	 0.001060307189285 =2*PI()/E8
    sr = 0                          # Sensor Roll (r) [degrees] =	0

    #TODO put earth parameters into a dict and add support for other spheroids GRS1980 etc.
    # Earth Stuff (WGS84)

    one_on_f = GRAVITATIONAL_CONST     # Inverse flattening 1/f = 298.257223563

    f = 1/one_on_f                    # flattening

    r = RADIUS_OF_EARTH	

    e = 1-math.pow((1-1/one_on_f), 2)  # Eccentricity (e^2) = 0.00669438 =1-(1-1/I5)^2

#    WO = 0.000072722052             # rotation (w0) [rad/sec] = 7.2722052E-05

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
    rotation = math.atan((WO*math.cos(phi_rad)*math.cos(beta*math.pi/180))/(av+WO*math.cos(phi_rad)*math.sin(beta*math.pi/180)))*180/math.pi
    eh = beta+rotation
    alpha12 = eh
    s = 0.5*185000  # s = distance in metres
    effective_heading = alpha12
    return effective_heading

def getUpcomingPasses(satellite_name,satellite_swath,tle_information, passes_begin_time, passes_period):


    observer       = ephem.Observer()
    observer.lat   = GROUND_STATION[0]
    observer.long  = GROUND_STATION[1]
    #updatetime = 0
    period = passes_period
    #Get most recent TLE for determining upcoming passes from now
    tles = tle_information

    # make a list of dicts to hold the upcoming pass information for the selected satellites
    SCHEDULE = []
    observer.date = passes_begin_time

    while 1:

        for tle in tles:

            if tle[0].strip()== satellite_name:

                #TODO clean up the use of pyephem versus orbital. Orbital can give a orbit number and does many of the pyephem functions
                #TODO add the individual acquisitions as layers in the same ogr output
                #TODO use an appropriate google earth icon for satellites at a visible display resolution with a name tag and minutesaway
                #TODO print output to logging
                satname = str(tle[0]).replace(" ","_")
                sat = ephem.readtle(tle[0],tle[1],tle[2])

                twole = tlefile.read(tle[0],DATA_IN_DIR+'tles.txt')
                now = datetime.utcnow()
                #TODO check age of TLE - if older than x days get_tle()
#                print "TLE EPOCH:",twole.epoch
                
                oi = float(str.split(tle[2],' ')[3])
                orb = Orbital(tle[0])
                attributes = []
                rt, ra, tt, ta, st, sa = observer.next_pass(sat)

                # Determine is pass descending or ascending
                sat.compute(rt)
                aos_lat = sat.sublat.real*(180/math.pi)
                sat.compute(st)
                los_lat = sat.sublat.real*(180/math.pi)

                if (aos_lat > los_lat):
#                    print "PASS                 = descending"
                    node = "descending"
                else:
#                    print "PASS                 = ascending"
                    node = "ascending"
                    oi = 360 - oi

                AOStime = datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")
                minutesaway = (AOStime-now).seconds/60.0

#		print "Satellie             = ", satname
#               print "Minutes to horizon   = ", minutesaway
#               print "AOStime              = ", rt
#               print "LOStime              = ", st
#               print "Transit time         = ", tt
#	-----------------------------------------------------------------------------
#		This is a test routine for calculating Az, El angles
#	-----------------------------------------------------------------------------



                orad = orb.get_lonlatalt(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S"))[2]

#               print '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
		test_az_el_cal(orb,observer, rt, ra, tt, ta, st, sa,orad)
#               print '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'



#		
                attributes = {'Satellite name': satname, 'Orbit height': orad, 'Orbit': orb.get_orbit_number(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")), \
#                attributes = {'Satellite name': satname, 'Orbit height': orad, 'Orbit': orb.get_orbit_number(datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S")), \
                              'Current time': str(now),'Minutes to horizon': minutesaway, 'AOS time': str(rt), \
                              'LOS time': str(st), 'Transit time': str(tt), 'Node': node}

                # Append the attributes to the list of acquisitions for the acquisition period
                if not any ((x['Satellite name'] == satname and x['Orbit'] == orb.get_orbit_number(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")))for x in SCHEDULE):
#                if not any ((x['Satellite name'] == satname and x['Orbit'] == orb.get_orbit_number(datetime.strptime(str(tt), "%Y/%m/%d %H:%M:%S")))for x in SCHEDULE):
                    SCHEDULE.append(attributes)

                # Step from AOS to LOS in 100 second intervals
#                delta = timedelta(seconds=100)
                delta = timedelta(seconds=DELTA_TIME_STEP)
                deltatime = datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S")
                geoeastpoint = []
                geowestpoint = []
                geotrack = []


#                print "DELTATIME", deltatime
#                print "SETTING TIME", datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S")


#		Tesing for next satellite



#	--------------------------------------------------------------------------------------------
#	--------------------------------------------------------------------------------------------
#		The following set of lines have been for testing while making comparision in seconds 
#		instead of string comparisiom
#	--------------------------------------------------------------------------------------------
#	--------------------------------------------------------------------------------------------

#		print '================ Testing Loop starts ==========================================='
#		print 'deltatime       = ',deltatime
#		print 'Secs Time       = ', get_time_secs(str(deltatime).replace("-","/"))
#		print 'st              = ',str(datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S"))
#		print 'st in Secs Time = ',get_time_secs(str(datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S")).replace('-','/'))
#		print '================ Testing Loop Ends   ==========================================='
#		The following if statement has ben included on the basis of dpoch seconds
#

		if get_time_secs(str(deltatime).replace("-","/")) >= \
			get_time_secs(str(datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S")).replace('-','/')):
			return()

		print 'Delta Time = ',deltatime
		print 'date time  = ',datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S")		
		print '---------------------------'

#		if deltatime >= datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S"):
#			return()

                while deltatime < datetime.strptime(str(st), "%Y/%m/%d %H:%M:%S"):

                    sat.compute(deltatime)
                    geotrack.append({'lat2': sat.sublat.real*(180/math.pi), \
                                     'lon2': sat.sublong.real*(180/math.pi), \
                                     'alt2': orb.get_lonlatalt(datetime.strptime(str(rt), "%Y/%m/%d %H:%M:%S"))[2]*1000})

                    eastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi), orad, sat._n)+90
                    westaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi), orad, sat._n)+270

                    #Set ground swath per satellite sensor
                    #TODO use view angle check to refine step from satellite track see IFOV
                    
		    swath = float(satellite_swath)/2.
		    
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
                #CURRENT_POSITION_FILENAME = satname+"_current_position.kml"

                CURRENT_POSITION_FILENAME = OUTPUT_DIR+satname+"_current_position.kml"

                #TODO draw the current orbit forward for the passes period time from the satellite position as a long stepped ogr line

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

                observer.date=datetime.utcnow()
                sat.compute(observer)

#               tkdelta = timedelta(seconds=100)

                tkdelta = timedelta(seconds=DELTA_TIME_STEP)

                tkrt, tkra, tktt, tkta, tkst, tksa = observer.next_pass(sat)
                tkdeltatime = datetime.utcnow()
                tkgeoeastpoint = []
                tkgeowestpoint = []
                tkgeotrack = []

                while tkdeltatime < (datetime.utcnow() or datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S")):

                    sat.compute(tkdeltatime)
                    tkgeotrack.append({'lat2':sat.sublat.real*(180/math.pi),'lon2':sat.sublong.real*(180/math.pi),'alt2':orb.get_lonlatalt(datetime.strptime(str(rt),"%Y/%m/%d %H:%M:%S"))[2]})

                    tkeastaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad,sat._n)+90
                    tkwestaz = getEffectiveHeading(sat,oi,sat.sublat.real*(180/math.pi), sat.sublong.real*(180/math.pi),orad,sat._n)+270
                    #TODO use view angle check to refine step from satellite track see IFOV

		    tkswath = float(satellite_swath)/2.
		    
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

                if not ((attributes['Node']=="ascending")and(satname not in ("AQUA"))):
                    # Create swath ogr output
                    getVectorFile(attributes,polypoints,'polygon', SWATH_FILENAME, 'KML')
                    # Create orbit track ogr output
                    getVectorFile(attributes,geotrack,'line', ORBIT_FILENAME, 'KML')
                    # Create currently acquiring ogr output
                    if ((now >= datetime.strptime(str(tkrt),"%Y/%m/%d %H:%M:%S")) and (now <= datetime.strptime(str(tkst),"%Y/%m/%d %H:%M:%S"))):
                        getVectorFile(now_attributes,tkpolypoints,'polygon', TRACKING_SWATH_FILENAME, 'KML')

                if minutesaway <= period:

#                    print tle[0], 'WILL BE MAKING A PASS IN ', minutesaway, " MINUTES"
#                    print ' Rise Azimuth: ', ra
#                    print ' Transit Time: ', tt
#                    print ' Transit Altitude: ', ta
#                    print ' Set Time: ', st
#                    print ' Set Azimuth: ', sa
#                    print '================================================='
#		    print 'Satellite Name = ',satellite_name
                    for x in sorted(SCHEDULE, key=lambda k: k['AOS time']):
#			print x
			output_orbit_parameters(x)
                        # For dictionary entries with 'LOS time' older than now time - remove
                        if ((datetime.strptime(str(x['LOS time']),"%Y/%m/%d %H:%M:%S"))<(datetime.utcnow())):
                            # Delete output ogr
                            if os.path.exists(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_swath.kml")):
                                os.remove(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_swath.kml"))
                            if os.path.exists(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_track.kml")):
                                os.remove(os.path.join(output_path,satname+"."+str(x['Orbit'])+".ALICE.orbit_track.kml"))
                            # Delete dictionary entry for pass
                            SCHEDULE.remove(x)

                    # Unlikely - if no entries in the SCHEDULE don't try to print it

                    if len(SCHEDULE)>0:
			print (datetime.strptime(str(SCHEDULE[0]['AOS time']),"%Y/%m/%d %H:%M:%S"))

                    # If the AOS time is less than now + the time delta, shift the time to the latest recorded pass LOS time

                    if ((datetime.strptime(str(SCHEDULE[len(SCHEDULE)-1]['AOS time']),"%Y/%m/%d %H:%M:%S")<(datetime.utcnow()+timedelta(minutes=period)))):
                        observer.date = (datetime.strptime(str(SCHEDULE[len(SCHEDULE)-1]['LOS time']),"%Y/%m/%d %H:%M:%S")+timedelta(minutes=5))
                        # Recompute the satellite position for the update time
                        sat.compute(observer)
#                        print "MODIFIED OBSERVER DATE",observer.date
                    else:
#                       print "--------NOTHING TO MODIFY MOVING TO NEXT SATELLITE IN LIST------"
                        #TODO - write to html

                        # Exit the def if the SCHEDULE isn't able to update because there are no passes in the acquisition window
                        return ()

#	print 'Before Time Sleep ......'
#	print 'Loop for While .........'
	print '============================================================================='

        time.sleep(1*SLEEP_STATUS)
    return ()

def SAT_Interference():

	print '='*130
	SAT_AOS_WIN.sort()
	for index in range(len(SAT_AOS_WIN)-1):

		val1 = SAT_AOS_WIN[index].split(':')
		val2 = SAT_AOS_WIN[index+1].split(':')
		if val1[5].strip() == val2[5].strip():
			continue
		if ( ( int(val2[0].strip()) < int(val1[0].strip()) ) and \
			( int(val2[0].strip()) > int(val1[1].strip()) ) ) or \
			( ( int(val2[1].strip()) > int(val1[0].strip()) ) and \
			( int(val2[1].strip()) < int(val1[1].strip()) ) ):
				print '%36s'%' Orbit Interference Detected for ', '%9s'%val1[5].strip(),'%8s'%'Orbit : ',int(val1[2].strip()) \
				, ' : AOS Time = ',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(int(val1[3]))) \
				, ' : Earth station Data Collection Interval ',val1[4].rjust(4),' secs '
				print '%36s'%' Orbit Interference Detected for ','%9s'%val2[5].strip(),'%8s'%'Orbit : ',int(val2[2].strip()) \
				, ' : AOS Time = ',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(int(val2[3]))) \
				, ' : Earth station Data Collection Interval ',val2[4].rjust(4),' secs '
				print '='*130
		else:
				print '%36s'%' No Orbit Interference Detected for ','%9s'%val1[5].strip(),'%8s'%'Orbit : ',int(val1[2].strip()) \
				, ' : AOS Time = ',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(int(val1[3]))) \
				, ' : Earth station Data Collection Interval ',val1[4].rjust(4),' secs '
				print '%36s'%' No Orbit Interference Detected for ','%9s'%val2[5].strip(),'%8s'%'Orbit : ',int(val2[2].strip()) \
				, ' : AOS Time = ',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(int(val2[3]))) \
				, ' : Earth station Data Collection Interval ',val2[4].rjust(4),' secs '
				print '='*130
	
			
if __name__ == '__main__':

    get_parameters()
    set_parameters()
    tles = get_tles()

# Loop through satellite list and execute until end of period

# Commnented from the original code
#    while 1:
#        for item in SATELLITE_SWATH:
#	    sl = item.split(',')
#            getUpcomingPasses(sl[0],sl[1],tles,datetime.utcnow(),period)
#	sys.exit()
for item in SATELLITE_SWATH:
	sl = item.split(',')
	getUpcomingPasses(sl[0].strip(),sl[1].strip(),tles,datetime.utcnow(),period)


SAT_Interference()
