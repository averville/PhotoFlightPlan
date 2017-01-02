#! /usr/bin/python3
#############################################################
# photoflightplan1.01.py by Kildir Technologies, GPL license
# André Verville 2016.12.20 and +
# Photographic flight plan generation for manned aircrafts
# produces a flight plan (.pfp) data file for Collimator,
# a Raspberry Pi based flight navigation system
# and exports a kml output to see the results in Google Earth
#############################################################

# -----------------------------------
# Point class for computing functions
# -----------------------------------
class Point(object):
    global x, y, Point
    def __init__(self, x, y):
        self.x = x
        self.y = y

# ------------------------------------------------------------------------------
# Function to create a list from a string of values read in the parameters table
# ------------------------------------------------------------------------------
def extractlist(string, separator):
    lststring = string.split(separator)
    ndx = 0
    lst = []
    for x in lststring:
        lst = lst + [float(lststring[ndx])]
        ndx += 1
    return lst

# ----------------------------------------------------------------------
# Function to determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This function
# returns True or False. The algorithm is called the "Ray Casting Method
# In our case, the x and y are latitudes and longitudes
# which are angular values of different scale on the ground
# but approximating these values to projected coordinates is acceptable
# in the context of an accuracy within a few meters: we stay in lat-lon
# Credit: Joel Lawhead, geospatialpython.com
# -----------------------------------------------------------------------
def pointinpolygon(x, y, poly):
    inside = False
    p1x, p1y = poly[0], poly[1]
    for index in range(0, len(poly), 2):
        p2x = poly[index]
        p2y = poly[index + 1]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y: xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints: inside = not inside
        p1x, p1y = p2x, p2y
    return inside

# Function to end program prematurely in case of error
# ----------------------------------------------------
def abnormalend():
    global logfile
    msghd = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
    logfile.write(msghd + "Abnormal End\n")
    print(msghd + "Stop - Please check log file")
    exit()

# Vincenty Direct geodetic computation function (Andre Verville)
# generates latitude, longitude and inverse bearing
# from a given point in lat lon, bearing and distance
# Formulas from my Bachelor Degrees thesis, Laval University, 1978
# Credits: Thaddeus Vincenty, Survey Review, Vol.XXIII, No. 176, April 1975
# https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
# ----------------------------------------------------------------------------

def vincentydirect(lat1deg, lon1deg, brg12deg, dist):
    import math
    if dist < 0:  # if the distance is negative, we reverse the bearing (can be useful)
        dist = abs(dist)
        brg12deg = brg12deg + 180
    brg12deg = brg12deg % 360  # and make sure it does not exceed 360 degrees
    a, b = 6378137.0, 6356752.3141 # GRS80 ellipsoid, edit as required

    lat1rad = math.radians(lat1deg)
    lon1rad = math.radians(lon1deg)
    brg12rad = math.radians(brg12deg)

    e2 = (a * a - b * b) / (a * a)
    N1 = a / ((1.0 - e2 * (math.sin(lat1rad) ** 2.0))) ** 0.5
    M1 = (a * (1.0 - e2)) / ((1 - e2 * math.sin(lat1rad) ** 2.0) ** 1.5)
    R = (M1 * N1) / ((M1 * math.sin(brg12rad) ** 2.0) + (N1 * math.cos(brg12rad) ** 2.0))
    d12 = dist - (dist ** 3.0) / (24.0 * (R ** 2.0))
    X1 = N1 * math.cos(lat1rad)
    Z1 = N1 * (1.0 - e2) * math.sin(lat1rad)
    mu = math.asin(d12 / (2.0 * R))

    while 1: # iterative section to compute X2, Y2 and Z2 while minimizing h
        X2 = X1 - d12 * (math.cos(lat1rad) * math.sin(mu) + math.sin(lat1rad)\
         * math.cos(brg12rad) * math.cos(mu))
        Y2 = d12 * math.sin(brg12rad) * math.cos(mu)
        Z2 = Z1 + d12 * (math.cos(lat1rad) * math.cos(brg12rad) * math.cos(mu)\
         - math.sin(lat1rad) * math.sin(mu))
        h = ((X2 ** 2.0) + (Y2 ** 2.0) + ((Z2 ** 2.0) / (1.0 - e2))) ** 0.5 - a
        sigmamu = 180.0 * h / (math.pi * d12 * math.cos(mu))
        mu = mu + sigmamu
        if h < 0.0001:
            break

    lat2rad = math.atan(Z2 / ((1.0 - e2) * ((X2 ** 2.0) + (Y2 ** 2.0)) ** 0.5))
    deltalon = math.atan(Y2 / X2)
    lon2rad = lon1rad + deltalon
    N2 = a / ((1 - e2 * math.sin(lat2rad) ** 2) ** 0.5)
    temp = (Z2 - Z1) * math.cos (lat2rad) - (N2 * math.cos(lat2rad)+
        - X1 * math.cos(deltalon)) * math.sin(lat2rad)
    brg21rad = math.atan(X1 * math.sin(deltalon) / temp)

    lat2deg = math.degrees(lat2rad)
    lon2deg = math.degrees(lon2rad)
    brg21deg = math.degrees(brg21rad)
    if abs(brg21deg - brg12deg) < 10:  # make sure the reverse bearing is at 180 degrees
        brg21deg = brg21deg + 180
        brg21deg = brg21deg % 360
    return lat2deg, lon2deg, brg21deg

# Function to read parameters file data
# gets parameters and polygon
# but also returns end of file status if no more data to read
# -----------------------------------------------------------
def readparams():
    global parfile, aoiname, aoiextent, aoikml, altitudeaglft, lensfocal, resolution
    global sensorxpix, sensorxmm, sensorypix, sensorymm, winddir
    global overlap, sidelap, flightblock, polygon
    polygon = []  # initial state: empty list for the polygon data

    while 1:
        paramdata = parfile.readline() # read parameters line by line
        if not paramdata: # when no more lines to read
            parfile.close()
            return False  # then return False, saying there is no more AOI to process
        if paramdata == "": continue  # go read another line if nothing to read
        if "aoiname =" in paramdata: aoiname = paramdata.split('"')[1]
        if "aoiextent =" in paramdata: aoiextent = float(paramdata.split('"')[1])
        if "aoikml =" in paramdata: aoikml = paramdata.split('"')[1]
        if "altitudeaglft =" in paramdata: altitudeaglft = paramdata.split('"')[1]
        if "lensfocal =" in paramdata: lensfocal = paramdata.split('"')[1]
        if "resolution =" in paramdata: resolution = float(paramdata.split('"')[1])
        if "sensor_xpix =" in paramdata: sensorxpix = float(paramdata.split('"')[1])
        if "sensor_xmm =" in paramdata: sensorxmm = float(paramdata.split('"')[1])
        if "sensor_ypix =" in paramdata: sensorypix = float(paramdata.split('"')[1])
        if "sensor_ymm =" in paramdata: sensorymm = float(paramdata.split('"')[1])
        if "wind_direction =" in paramdata: winddir = float(paramdata.split('"')[1])
        if "overlap =" in paramdata: overlap = float(paramdata.split('"')[1])
        if "sidelap =" in paramdata: sidelap = float(paramdata.split('"')[1])
        if "flightblock =" in paramdata:
            flightblock = paramdata.split('"')[1]
            return True  # getting out of the function with a positive feedback

        # Case where we have Area of Interest (AOI) polygon points to read
        # within the parameters file itself, storing them into a list named polygon
        # -------------------------------------------------------------------------
        if "polygon =" in paramdata and aoikml == "none":
            poly = extractlist(paramdata.split('"')[1], " ")
            polygon = polygon + [poly[1], poly[2]]  # add lat and lon to polygon list
        continue

# ========================================
# Main execution and block processing loop
# ========================================
import os, sys, time
global parfile, aoiname, aoiextent, aoikml, altitudeaglft, lensfocal, resolution
global sensorxpix, sensorxmm, sensorypix, sensorymm, winddir
global overlap, sidelap, flightblock, polygon
kmlheader = False  # flag that keeps track if the kml header has been written

# Opening log file for execution time and error tracking
# ------------------------------------------------------
logfile = open("photoflightplan.log", "a")
loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
logfile.write("\n" + loghead + "Startup\n")

# Reading provided arguments from the command line
# Defaulting to "flightplan" if no file name provided
# ---------------------------------------------------
args, nbargs = sys.argv, len(sys.argv)
for index in range(1, nbargs, 1):
    if ".par" in sys.argv[index]: parname = sys.argv[index]
    if ".pfp" in sys.argv[index]: pfpname = sys.argv[index]
    if ".kml" in sys.argv[index]: kmlname = sys.argv[index]
try: parname
except NameError: parname = "photoflightplan.par"
try: pfpname
except NameError: pfpname = "photoflightplan.pfp"
try: kmlname
except NameError: kmlname = "photoflightplan.kml"

# Opening the flight plan file for flight data and report
# This file is opened by Collimator to get its project info
# ---------------------------------------------------------
planfile = open(pfpname, "w")

# Parameters file opening and execution parameters reading
# --------------------------------------------------------
try:
    parfile = open(parname, "r")  # the parameters file reading buffer
except:
    loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
    logfile.write(loghead + "Error file " + parname + " not found\n")
    abnormalend()

# Printing startup message
# ------------------------
print("------------------------------------------------------------")
print(time.strftime("%Y.%m.%d %H:%M:%S"), " PhotoFlightPlan V1.01 by André Verville")
print("------------------------------------------------------------")

while readparams():  # process if this is a dataset
    # Checking at least if the parameters have been provided
    # ------------------------------------------------------
    paramerror = ""  # initialize paramerror to nothing
    try: aoiname
    except NameError: paramerror = 'No parameter provided for "aoiname"\n'
    try: aoiextent
    except NameError: paramerror = 'No parameter provided for "aoiname"\n'
    try: aoikml
    except NameError: paramerror = 'No parameter provided for "aoikml"\n'
    try: altitudeaglft
    except NameError: paramerror = 'No parameter provided for "altitudeaglft"\n'
    try: lensfocal
    except NameError: paramerror = 'No parameter provided for "lensfocal"\n'
    try: resolution
    except NameError: paramerror = 'No parameter provided for "resolution"\n'
    try: sensorxpix
    except NameError: paramerror = 'No parameter provided for "sensor_xpix"\n'
    try: sensorxmm
    except NameError: paramerror = 'No parameter provided for "sensor_xmm"\n'
    try: sensorypix
    except NameError: paramerror = 'No parameter provided for "sensor_ypix"\n'
    try: sensorymm
    except NameError: paramerror = 'No parameter provided for "sensor_ymm"\n'
    try: winddir
    except NameError: paramerror = 'No parameter provided for "wind_direction"\n'
    try: overlap
    except NameError: paramerror = 'No parameter provided for "overlap"\n'
    try: sidelap
    except NameError: paramerror = 'No parameter provided for "sidelap"\n'
    try: polygon
    except NameError: paramerror = 'No parameter provided for "area of interest"\n'
    if paramerror != "":  # if paramerror contains something, something is wrong
        loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
        logfile.write(loghead + paramerror)
        logfile.write(loghead + "Abnormal End\n")
        print(loghead, "Abnormal End - see .log file for error status")
        exit()

    # Parameters input done, starting processing
    # ------------------------------------------
    loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
    logfile.write(loghead + "Processing block" + flightblock + "\n")

    # If the user provides an AOI polygon within a kml file
    # we extract the AOI polygon from the provided file
    # ----------------------------------------------------
    if aoikml != "none":
        try:
            aoikmlfile = open(aoikml, "r")
            count = False
            index = 0
            for item in aoikmlfile:
                if "<Polygon>" in item: count = True # the coords are 6 lines after
                if count == True: index += 1
                if index == 6: polygonstring = item
            aoikmlfile.close

        except:  # catch all for any type of error but most probably "file not found"
            loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
            logfile.write(loghead + "Error file " + aoikml + " not found\n")
            abnormalend()

        # Parsing the extracted coordinates string to create the polygon list
        # -------------------------------------------------------------------
        polygonstring = polygonstring.replace("\t", "")
        polygonstring = polygonstring.replace(",0 ", ",")
        polygonstring = polygonstring.replace("\n", "")
        polygonstring = polygonstring[:-1]
        lonlatpoly = extractlist(polygonstring, ",")  # list cleaned up from garbage

        polygon = []  # initial state: empty list
        for index in range(0, len(lonlatpoly), 2):  # inversing pairs for lat before lon
            polygon = polygon + [lonlatpoly[index + 1]]
            polygon = polygon + [lonlatpoly[index]]
        loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
        logfile.write(loghead + "Using AOI polygon from " + aoikml + " \n")
    else:
        loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
        logfile.write(loghead + "Using AOI polygon provided within parameters file\n")

    # Writing down the block header in the flight plan file (.pfp)
    # ------------------------------------------------------------
    planfile.write("================================================\n")
    planfile.write(time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01 output\n")
    planfile.write("================================================\n")

    # Initial computations
    # --------------------
    aoinbpoints = int(len(polygon) / 2) - 1
    imagewidth = sensorxpix * resolution / 100  # in meters at ground level
    imageheight = sensorypix * resolution / 100  # in meters at ground level
    triggerdist = (1 - overlap / 100) * imageheight
    corridorwidth = round((1 - sidelap / 100) * imagewidth)

    # Computing flight altitude from imposed lens focal
    # or computing lens focal from flight altitude
    # -------------------------------------------------
    if altitudeaglft == "compute" and lensfocal == "compute":
        loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
        msg = "Error - cannot compute altitudeaglft AND lensfocal at the same time\n"
        logfile.write(loghead + msg)
        abnormalend()
    if altitudeaglft == "compute":
        lensfocal = float(lensfocal)
        altitudeaglm = imagewidth * lensfocal / sensorxmm
        altmode = " (computed)\n"
        loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
        logfile.write(loghead + "Given lensfocal, flight altitude computed\n")
    else:
        altitudeaglft = float(altitudeaglft)
        altmode = " (imposed)\n"
    if lensfocal == "compute":
        altitudeaglm = float(altitudeaglft) * 0.3048
        lensfocal = altitudeaglm * sensorxmm / imagewidth
        lensmode = " (computed)\n"
        loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
        logfile.write(loghead + "Given flight altitude, lens focal computed\n")
    else:
        lensfocal = float(lensfocal)
        lensmode = " (imposed)\n"
    altitudeaglft = altitudeaglm / 0.3048

    # Write down prameters in the flight plan file (.pfp)
    # ---------------------------------------------------
    planfile.write("parameters file:   " + parname + "\n")
    planfile.write("flight plan file:  " + pfpname + " (this file)\n")
    planfile.write("Google Earth file: " + kmlname + "\n")
    planfile.write("flightblock:       " + flightblock + "\n")
    planfile.write("aoiname:           " + str(aoiname) + "\n")
    planfile.write("aoiextent:         " + str(round(aoiextent)) + " m\n")
    planfile.write("lensfocal:         " + str(round(lensfocal, 1)) + " mm" + lensmode)
    planfile.write("resolution:        " + str(round(resolution, 1)) + " cm\n")
    planfile.write("sensor_xpix:       " + str(round(sensorxpix)) + " pixels\n")
    planfile.write("sensor_xmm:        " + str(round(sensorxmm, 1)) + " mm\n")
    planfile.write("sensor_ypix:       " + str(round(sensorypix)) + " pixels\n")
    planfile.write("sensor_ymm:        " + str(round(sensorymm, 1)) + " mm\n")
    planfile.write("wind_direction:    " + str(round(winddir)) + " degrees\n")
    planfile.write("overlap:           " + str(round(overlap)) + "%\n")
    planfile.write("sidelap:           " + str(round(sidelap)) + "%\n")
    planfile.write("imagewidth:        " + str(round(imagewidth)) + " m\n")
    planfile.write("imageheight:       " + str(round(imageheight)) + " m\n")
    planfile.write("altitudeaglm:      " + str(round(altitudeaglm)) + " m" + altmode)
    planfile.write("altitudeaglft:     " + str(round(altitudeaglft)) + " ft" + altmode)
    planfile.write("triggerdist:       " + str(round(triggerdist)) + " m\n")
    planfile.write("corridorwidth:     " + str(round(corridorwidth)) + " m\n")

    # Write down the AOI polygon data and where it comes from
    # -------------------------------------------------------
    planfile.write("\n--------------------------------\n")
    planfile.write("AOI Polygon points (" + str(aoinbpoints) + "):\n")
    if aoikml != "none": planfile.write("From file: " + aoikml + "\n")
    else: planfile.write("From parameters file\n")
    planfile.write("--------------------------------\n")

    for index in range(0, len(polygon), 2):
        if index == len(polygon) - 2: planfile.write("polygon: 01 ")  # Last repeats 1st
        else: planfile.write("polygon: " + str(int((index/2 + 1))).zfill(2) + " ")
        planfile.write(str(round(polygon[index],6)).ljust(9, "0") + " ")
        planfile.write(str(round(polygon[index + 1],6)).ljust(10, "0") + "\n")
    planfile.write("--------------------------------\n")

    # Find the approximate gravity center of the polygon
    # flight lines will be spread over this central point
    # ---------------------------------------------------
    latmin, latmax, lonmin, lonmax = 90, -90, 180, -180
    for index in range(0, len(polygon) - 2, 2):  # remove the last pair (same as first)
        if polygon[index] < latmin: latmin = polygon[index]
        if polygon[index] > latmax: latmax = polygon[index]
        if polygon[index + 1] < lonmin: lonmin = polygon[index + 1]
        if polygon[index + 1] > lonmax: lonmax = polygon[index + 1]
    gclat = (latmin + latmax) / 2  # gravity center latitude
    gclon = (lonmin + lonmax) / 2  # gravity center longitude

    # Compute the scan width and angles in all directions
    # --------------------------------------------------
    scanwidth = corridorwidth * 50  # see comment below (outer loop comment)
    northacross = winddir + 90
    if northacross > 90 and northacross < 270:  # make sure we are going North-something
        northacross = northacross + 180
    northacross = northacross % 360  # and that it does not exceed 360 degrees

    # Move along the lines to find where we enter and exit the polygon
    # ---------------------------------------------------------------
    isaflightline = False # beginning status
    flightlines = []

    # Outer loop where we move by corridor width from south to north (+/-)
    # starting well off the polygon, finishing well off also
    # this routine good for AOIs not bigger than 50 corridor widths
    # if by any chance it is the case, we can change scanwidth to something
    # a greater multiple of corridorwidth, will just take longer to process
    # --------------------------------------------------------------------
    start = -int((scanwidth / 2))  # start with a negative scan distance
    end = -start + corridorwidth  # do the last in range also
    for distacross in range(start, end, corridorwidth):
        if distacross == 0:  # at central point we cannot compute with distance zero
            reflat = gclat
            reflon = gclon
        else:
            reflat, reflon, brg = vincentydirect(gclat, gclon, northacross, distacross)

        # Scan process indication - increasing percentage on the same line
        # tricky assembly to satisfy IOS Pythonista app that is a bit fussy
        # formatting the dynamic one-liner overwrites
        # ----------------------------------------------------------------
        msghd = '\rScanning across AOI "'+ aoiname + '" for block ' + flightblock
        percent = round(distacross + scanwidth / 2) / scanwidth
        msg = msghd + " ({:2.0%}) done".format(percent)  # one-liner display for % done
        print(msg, end = "")
        
        # Initial status for inner loop
        # -----------------------------
        justinlat, justinlon, justoutlat, justoutlon = 0, 0, 0, 0
        distin, distout = 0, 0
        distincrement = int(round(corridorwidth / 100))  # increment for the scan along

        # Inner loop where we move by 1/100 corridor width increments along a potential
        # flight line at the beginning and at the end, we may not cross the perimeter
        # but when we do, the positive returns from the "pointinpolygon" function
        # will confirm we are crossing it and therefore record a flight line
        # -----------------------------------------------------------------------------
        inpolygon = False  # inside or outside status for "memory" effect to flag changes
        for distalong in range(start, end, distincrement):
            if distalong == 0:  # at central point, we cannot compute with distance zero
                chklat, chklon = reflat, reflon
            else:
                chklat, chklon, brg = vincentydirect(reflat, reflon, winddir, distalong)
            inside = pointinpolygon(chklat, chklon, polygon)  # check if within polygon
            if inside:
                inpolygon = True  # true until we get out
                if justinlat == 0:  # if we just entering into the polygon
                    justinlat, justinlon = chklat, chklon  # we remember the lat-lon
                    distin = distalong  # if required for debug purposes, the distance
            else:
                if inpolygon:  # if previously in, we are getting out of the perimeter
                    justoutlat, justoutlon = chklat, chklon  # we remember the lat-lon
                    distout = distalong  # distance at which we got out of the polygon
                    inpolygon = False  # confirm we are out. Note: we might enter again!

        if distout != 0:  # if we got out at some distance, we have crossed the perimeter
            point1 = Point(justinlon, justinlat)  # encode first point in flight line
            point2 = Point(justoutlon, justoutlat)  # encode second point in flight line
            flightlines = flightlines + [justinlat] + [justinlon]
            flightlines = flightlines + [justoutlat] + [justoutlon] 
    nbflightlines = int(len(flightlines) / 4)

    # Format and write the flight lines to the flight plan file
    # ---------------------------------------------------------
    planfile.write("\nFlight lines (" + str(nbflightlines) + "):\n")
    planfile.write("-----------------------------------------------------------\n")

    for index in range(0, len(flightlines), 4):
        planfile.write("flightline: " + flightblock + " ")
        planfile.write(str(int((index / 4 + 1))).zfill(2) + " ")
        planfile.write(str(round(flightlines[index],6)).ljust(9, "0") + " ")
        planfile.write(str(round(flightlines[index+1],6)).ljust(10, "0") + " ")
        planfile.write(str(round(flightlines[index+2],6)).ljust(9, "0") + " ")
        planfile.write(str(round(flightlines[index+3],6)).ljust(10, "0") + "\n")
    planfile.write("-----------------------------------------------------------\n\n")

    # Generate the kml file for Google Earth
    # First, we put the kml header and style data
    # warning: erases any previous kml file without notice
    # note: header data generated only once
    # ----------------------------------------------------
    while not kmlheader:
        kmlfile = open(kmlname, "w")
        kmlfile.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        kmlfile.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http:')
        kmlfile.write('//www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net')
        kmlfile.write('/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n')
        kmlfile.write("<Document>\n")
        kmlfile.write("<name>flightplan.kml</name>\n")
        kmlfile.write('<StyleMap id="kildirtech">\n')
        kmlfile.write("</StyleMap>\n")
        kmlfile.write('<Style id="flightplan1">\n')
        kmlfile.write("<LineStyle>\n")
        kmlfile.write("<color>ffff00ff</color>\n")  # opacity/blue/green/red in hex
        kmlfile.write("<width>3</width>\n")
        kmlfile.write("</LineStyle>\n")
        kmlfile.write("</Style>\n")
        kmlfile.write('<Style id="flightplan2">\n')
        kmlfile.write("<LineStyle>\n")
        kmlfile.write("<color>ffff0000</color>\n")  # opacity/blue/green/red in hex
        kmlfile.write("<width>5</width>\n")
        kmlfile.write("</LineStyle>\n")
        kmlfile.write("<PolyStyle>\n")
        kmlfile.write("<color>33ff0000</color>\n")  # opacity/blue/green/red in hex
        kmlfile.write("</PolyStyle>\n")
        kmlfile.write("</Style>\n")
        kmlfile.write("<Folder>\n")
        kmlfile.write("<name>Collimator flight plan</name>\n")
        kmlfile.write("<open>1</open>\n")
        kmlheader = True

    # Insert the AOI polygon
    # ----------------------
    kmlfile.write("<Placemark>\n")
    kmlfile.write("<name>AOI " + aoiname + "</name>\n")
    kmlfile.write("<styleUrl>#flightplan2</styleUrl>\n")
    kmlfile.write("<Polygon>\n")
    kmlfile.write("<tessellate>1</tessellate>\n")
    kmlfile.write("<outerBoundaryIs>\n")
    kmlfile.write("<LinearRing>\n")
    kmlfile.write("<coordinates>\n")

    for index in range(0, len(polygon), 2):
        kmlfile.write(str(round(polygon[index+1],6)).ljust(10, "0") + ",")
        kmlfile.write(str(round(polygon[index],6)).ljust(9, "0") + ",0 ")

    kmlfile.write("</coordinates>\n")
    kmlfile.write("</LinearRing>\n")
    kmlfile.write("</outerBoundaryIs>\n")
    kmlfile.write("</Polygon>\n")
    kmlfile.write("</Placemark>\n")

    # Insert flight lines, one folder per block
    # -----------------------------------------
    kmlfile.write("<Folder>\n")
    kmlfile.write("<name>" + "Block " + flightblock + "</name>")
    kmlfile.write("<open>1</open>")
    for index in range(0, len(flightlines), 4):
        kmlfile.write("<Placemark>\n")
        kmlfile.write("<name>Line " + str(int((index/4+1))).zfill(2) + "</name>\n")
        kmlfile.write("<styleUrl>#flightplan1</styleUrl>\n")
        kmlfile.write("<LineString>\n")
        kmlfile.write("<tessellate>1</tessellate>\n")
        kmlfile.write("<coordinates>\n")
        kmlfile.write(str(round(flightlines[index + 1],6)).ljust(10, "0"))
        kmlfile.write("," + str(round(flightlines[index],6)).ljust(9, "0"))
        kmlfile.write(",0 " + str(round(flightlines[index + 3],6)).ljust(10, "0"))
        kmlfile.write("," + str(round(flightlines[index + 2],6)).ljust(9, "0") + ",0 \n")
        kmlfile.write("</coordinates>\n")
        kmlfile.write("</LineString>\n")
        kmlfile.write("</Placemark>\n")
    kmlfile.write("</Folder>")

    # Print summary for this block
    # ----------------------------
    print("\nFound and recorded", nbflightlines, "flight lines for block", flightblock)
    loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
    logfile.write(loghead + str(nbflightlines) + " flight lines computed for block")
    logfile.write(flightblock + "\n")

    # While loop ending statements
    # go read another block as long as the readparams() function returns True
    # get out of the loop when readparams() returns False
    # -----------------------------------------------------------------------
    continue

# Case where there is not even one "flightblock" in the parameters file
# PhotoFlightPlan would just get through without doing anything
# Here, we catch this singular type of error to warn the user
# ---------------------------------------------------------------------
try: flightblock
except NameError: paramerror = 'No parameter provided for "flightblock"\n'
if paramerror != "":  # if paramerror contains something, meaning something is wrong
    loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
    logfile.write(loghead + paramerror)
    logfile.write(loghead + "Abnormal End\n")
    print(loghead, "Abnormal End - see .log file for error status")
    exit()

# Insert the kml file tail and close it
# -------------------------------------
kmlfile.write("</Folder>\n")
kmlfile.write("</Document>\n")
kmlfile.write("</kml>\n")
kmlfile.close()

# Cleanup, say goodbye and exit
# -----------------------------
msg = 'Flight plan(s) generated in "' + pfpname + '" and "' + kmlname + '"\n'
print(msg)
loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.01: "
logfile.write(loghead + msg)
logfile.write(loghead + "Normal  End\n")
logfile.close()
planfile.close()
exit()
