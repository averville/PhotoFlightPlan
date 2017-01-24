#! /usr/bin/python3
############################################################
# photoflightplan1.07.py by Kildir Technologies, GPL license
# André Verville 2016.12.20 and +
# Photographic flight plan generation for manned aircrafts
# produces flight plan (.pfp) data files for Collimator,
# a Raspberry Pi based flight navigation system
# and exports kml outputs to see the results in Google Earth
# Note: output files formatted for Windows (cd/lf line ends)
############################################################

# -----------------------------------
# Point class for computing functions
# -----------------------------------
class Point(object):
    global x, y, Point
    def __init__(self, x, y):
        self.x = x
        self.y = y

# -------------------------------------------------------------------------
# Function to create a list from a string of values in the parameters table
# -------------------------------------------------------------------------
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
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

# ----------------------------------------------------
# Function to end program prematurely in case of error
# ----------------------------------------------------
def abnormalend():
    global logfile
    msghd = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.07: "
    logfile.write(msghd + "Abnormal End\r\n\r\n")
    print(msghd + "Stop - Please check log file")
    exit()

# -------------------------------------------------------------------------
# Vincenty Direct geodetic computation function (Andre Verville)
# generates latitude, longitude and inverse bearing
# from a given point in lat lon, bearing and distance
# Formulas from my Bachelor Degrees thesis, Laval University, 1978
# Credits: Thaddeus Vincenty, Survey Review, Vol.XXIII, No. 176, April 1975
# https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
# -------------------------------------------------------------------------
def vincentydirect(lat1deg, lon1deg, brg12deg, dist):
    import math
    if dist < 0:  # if distance negative, we reverse bearing (can be useful)
        dist = abs(dist)
        brg12deg = brg12deg + 180
    brg12deg = brg12deg % 360  # make sure it does not exceed 360 degrees
    a, b = 6378137.0, 6356752.3141 # GRS80 ellipsoid, edit as required

    lat1rad = math.radians(lat1deg)
    lon1rad = math.radians(lon1deg)
    brg12rad = math.radians(brg12deg)

    e2 = ((a * a) - (b * b)) / (a * a)
    sinlat1 = math.sin(lat1rad)
    coslat1 = math.cos(lat1rad)
    sinbrg12 = math.sin(brg12rad)
    cosbrg12 = math.cos(brg12rad)
    N1 = a / ((1.0 - e2 * (sinlat1 * sinlat1))) ** 0.5
    M1 = (a * (1.0 - e2)) / ((1 - (e2 * sinlat1 * sinlat1)) ** 1.5)
    R = (M1 * N1) / ((M1 * sinbrg12 * sinbrg12) + (N1 * cosbrg12 * cosbrg12))
    d12 = dist - (dist * dist * dist) / (24.0 * R * R)
    X1 = N1 * math.cos(lat1rad)
    Z1 = N1 * (1.0 - e2) * sinlat1
    mu = math.asin(d12 / (2.0 * R))
    sinmu = math.sin(mu)
    cosmu = math.cos(mu)

    while 1: # iterative section to compute X2, Y2 and Z2 while minimizing h
        X2 = X1 - d12 * (coslat1 * math.sin(mu) + sinlat1 * cosbrg12 * cosmu)
        Y2 = d12 * sinbrg12 * math.cos(mu)
        Z2 = Z1 + d12 * (coslat1 * cosbrg12 * math.cos(mu) - sinlat1 * sinmu)
        h = (((X2 * X2) + (Y2 * Y2) + ((Z2 * Z2) / (1.0 - e2))) ** 0.5) - a
        sigmamu = h / (d12 * math.cos(mu))
        mu = mu + sigmamu
        if abs(h) < 0.0001:
            break

    lat2rad = math.atan2(Z2, ((1.0 - e2) * ((X2 * X2) + (Y2 * Y2)) ** 0.5))
    sinlat2 = math.sin(lat2rad)
    coslat2 = math.cos(lat2rad)
    deltalon = math.atan2(Y2, X2)
    cosdeltalon = math.cos(deltalon)
    lon2rad = lon1rad + deltalon
    N2 = a / ((1 - e2 * math.sin(lat2rad) ** 2) ** 0.5)
    temp = (Z2 - Z1) * coslat2 - (N2 * coslat2 - X1 * cosdeltalon) * sinlat2
    brg21rad = math.atan(X1 * math.sin(deltalon) / temp)

    lat2deg = math.degrees(lat2rad)
    lon2deg = math.degrees(lon2rad)
    brg21deg = math.degrees(brg21rad)
    if abs(brg21deg - brg12deg) < 10:  # reverse bearing if at 180 degrees
        brg21deg = brg21deg + 180
        brg21deg = brg21deg % 360
    return lat2deg, lon2deg, brg21deg

# -------------------------------------------------------------------------
# Vincenty Inverse geodetic computation function (Andre Verville)
# generates bearing and distance from two given points in lat lon
# Formulas from my Bachelor Degrees thesis, Laval University, 1978
# Credits: Thaddeus Vincenty, Survey Review, Vol.XXIII, No. 176, April 1975
# https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
# -------------------------------------------------------------------------
def vincentyinverse(lat1deg, lon1deg, lat2deg, lon2deg):
    import math
    a, b = 6378137.0, 6356752.3141 # GRS80 ellipsoid, edit as required

    lat1rad = math.radians(lat1deg)
    lon1rad = math.radians(lon1deg)
    lat2rad = math.radians(lat2deg)
    lon2rad = math.radians(lon2deg)

    e2 = ((a * a) - (b * b)) / (a * a)
    sinlat1 = math.sin(lat1rad)
    coslat1 = math.cos(lat1rad)
    sinlat2 = math.sin(lat2rad)
    coslat2 = math.cos(lat2rad)

    deltalon = lon2rad - lon1rad
    N1 = a / ((1.0 - e2 * (sinlat1 * sinlat1))) ** 0.5
    M1 = (a * (1.0 - e2)) / ((1 - (e2 * sinlat1 * sinlat1)) ** 1.5)
    N2 = a / ((1.0 - e2 * (sinlat2 * sinlat2))) ** 0.5
    M2 = (a * (1.0 - e2)) / ((1 - (e2 * sinlat2 * sinlat2)) ** 1.5)
    X1 = N1 * coslat1
    Z1 = N1 * (1.0 - e2) * sinlat1
    X2 = N2 * coslat2
    Z2 = N2 * (1.0 - e2) * sinlat2
    Y = math.sin(deltalon) * X2
    temp = (N1 * coslat1) - (math.cos(deltalon) * X2)
    X = (temp * sinlat1) - ((Z1 - Z2) * coslat1)
    brg12 = math.atan2(Y, X)
    Y = math.sin(deltalon) * X1
    temp = (N2 * coslat2) - (math.cos(deltalon) * X1)
    X = ((Z2 - Z1) * coslat2) - (temp * sinlat2)
    brg21 = math.atan2(Y, X)
    X2P = math.cos(deltalon) * X2
    Y2 = math.sin(deltalon) * X2
    d12 = (((X2P - X1) * (X2P - X1)) + (Y2 * Y2) + ((Z2 - Z1) * (Z2 - Z1)))
    d12 = d12 ** 0.5
    brgb12 = (brg12 + brg21) / 2.0
    sinsqb = (math.sin(brgb12)) ** 2.0
    cossqb = (math.cos(brgb12)) ** 2.0
    Y = (N1 + N2) * (M1 + M2)
    X = (((M1 + M2) * sinsqb) + ((N1 + N2) * cossqb)) * 2.0
    RB = Y / X
    dist = d12 + ((d12 * d12 * d12) / (24.0 * RB * RB))
    brg12deg = math.degrees(brg12)
    brg21deg = math.degrees(brg21)
    if brg12deg < 0.0: brg12deg = brg12deg + 360.0  # adjust bearing (0-360)
    if brg21deg < 0.0: brg21deg = brg21deg + 360.0  # adjust bearing (0-360)
    brg21deg = (brg21deg + 180.0) % 360.0  # inverse bearing at 180 degrees
    return brg12deg, brg21deg, dist  # return 2 bearings and distance

# -----------------------------------------------------------
# Function to read parameters file data
# gets parameters and polygon
# but also returns end of file status if no more data to read
# -----------------------------------------------------------
def readparams():
    global parfile, aoiname, aoiextent, aoikml
    global altitudeaglft, lensfocal, resolution
    global sensorxpix, sensorxmm, sensorypix, sensorymm, winddir
    global overlap, sidelap, flightblock, polygon
    polygon = []  # initial state: empty list for the polygon data

    while 1:
        paramdata = parfile.readline() # read parameters line by line
        if not paramdata: # when no more lines to read
            parfile.close()
            return False  # return False if no more AOI to process
        if paramdata == "": continue  # go read another line if nothing to read
        if "aoiname =" in paramdata:
            aoiname = paramdata.split('"')[1]
        if "aoiextent =" in paramdata:
            aoiextent = float(paramdata.split('"')[1])
        if "aoikml =" in paramdata:
            aoikml = paramdata.split('"')[1]
        if "altitudeaglft =" in paramdata:
            altitudeaglft = paramdata.split('"')[1]
        if "lensfocal =" in paramdata:
            lensfocal = paramdata.split('"')[1]
        if "resolution =" in paramdata:
            resolution = float(paramdata.split('"')[1])
        if "sensor_xpix =" in paramdata:
            sensorxpix = float(paramdata.split('"')[1])
        if "sensor_xmm =" in paramdata:
            sensorxmm = float(paramdata.split('"')[1])
        if "sensor_ypix =" in paramdata:
            sensorypix = float(paramdata.split('"')[1])
        if "sensor_ymm =" in paramdata:
            sensorymm = float(paramdata.split('"')[1])
        if "wind_direction =" in paramdata:
            winddir = float(paramdata.split('"')[1])
        if "overlap =" in paramdata:
            overlap = float(paramdata.split('"')[1])
        if "sidelap =" in paramdata:
            sidelap = float(paramdata.split('"')[1])
        if "flightblock =" in paramdata:
            flightblock = paramdata.split('"')[1]
            return True  # getting out of the function with a positive feedback

        # Case where we have Area of Interest (AOI) polygon points to read
        # within the parameters file itself, storing them to list named polygon
        # ---------------------------------------------------------------------
        if "polygon =" in paramdata:
            if aoikml == "none":
                poly = extractlist(paramdata.split('"')[1], " ")
                polygon = polygon + [poly[1], poly[2]]  # add lat/lon to list
        continue

# ========================================
# Main execution and block processing loop
# ========================================
import os, sys, time
global parfile, aoiname, aoiextent, aoikml
global altitudeaglft, lensfocal, resolution
global sensorxpix, sensorxmm, sensorypix, sensorymm, winddir
global overlap, sidelap, flightblock, polygon
linelength = 0  # total flight line length to compute number of photo frames

# Opening log file for execution time and error tracking
# ------------------------------------------------------
logfile = open("photoflightplan.log", "a")
loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.07: "
logfile.write(loghead + "Startup\r\n")

# Reading optional argument from the command line
# Defaulting to "photoflightplan" if no parameters file name provided
# -------------------------------------------------------------------
args, nbargs = sys.argv, len(sys.argv)
for index in range(1, nbargs, 1):
    if ".par" in sys.argv[index]: parname = sys.argv[index]

try: parname
except NameError: parname = "photoflightplan.par"

# Parameters file opening and execution parameters reading
# --------------------------------------------------------
try:
    parfile = open(parname, "r")  # the parameters file reading buffer
except:
    loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.07: "
    logfile.write(loghead + "Error file " + parname + " not found\r\n")
    abnormalend()

# Printing startup message
# ------------------------
print("------------------------------------------------------------")
print(time.strftime("%Y.%m.%d %H:%M:%S"),\
     " PhotoFlightPlan V1.07 by André Verville")
print("------------------------------------------------------------")

while readparams():  # process if this is a dataset
    # Checking at least if the parameters have been provided
    # ------------------------------------------------------
    paramerror = ""  # initialize paramerror to nothing
    try: aoiname
    except NameError: paramerror = 'No parameter provided for "aoiname"\r\n'
    try: aoiextent
    except NameError: paramerror = 'No parameter provided for "aoiname"\r\n'
    try: aoikml
    except NameError: paramerror = 'No parameter provided for "aoikml"\r\n'
    try: altitudeaglft
    except NameError:
        paramerror = 'No parameter provided for "altitudeaglft"\r\n'
    try: lensfocal
    except NameError: paramerror = 'No parameter provided for "lensfocal"\r\n'
    try: resolution
    except NameError: paramerror = 'No parameter provided for "resolution"\r\n'
    try: sensorxpix
    except NameError:
        paramerror = 'No parameter provided for "sensor_xpix"\r\n'
    try: sensorxmm
    except NameError: paramerror = 'No parameter provided for "sensor_xmm"\r\n'
    try: sensorypix
    except NameError:
        paramerror = 'No parameter provided for "sensor_ypix"\r\n'
    try: sensorymm
    except NameError: paramerror = 'No parameter provided for "sensor_ymm"\r\n'
    try: winddir
    except NameError:
        paramerror = 'No parameter provided for "wind_direction"\r\n'
    try: overlap
    except NameError: paramerror = 'No parameter provided for "overlap"\r\n'
    try: sidelap
    except NameError: paramerror = 'No parameter provided for "sidelap"\r\n'
    try: polygon
    except NameError:
        paramerror = 'No parameter provided for "area of interest"\r\n'
    if paramerror != "":  # if paramerror not empty, something is wrong
        loghead = time.strftime("%Y.%m.%d %H:%M:%S")
        loghead = loghead + " PhotoFlightPlan V1.07: "
        logfile.write(loghead + paramerror)
        logfile.write(loghead + "Abnormal End\r\n\r\n")
        print(loghead, "Abnormal End - see .log file for error status")
        exit()

    # Parameters input done, starting processing
    # ------------------------------------------
    logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
    logfile.write(" PhotoFlightPlan V1.07: ")
    logfile.write("Processing block" + flightblock + "\r\n")

    # If the user provides an AOI polygon within a kml file
    # we extract the AOI polygon from the provided file
    # ----------------------------------------------------
    if aoikml != "none":
        try:
            aoikmlfile = open(aoikml, "r")
            count = False
            index = 0
            for item in aoikmlfile:
                if "<Polygon>" in item: count = True # coords are 6 lines after
                if count == True: index += 1
                if index == 6: polygonstring = item
            aoikmlfile.close

        except:  # catch all errors but most probably "file not found"
            logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
            logfile.write(" PhotoFlightPlan V1.07: ")
            logfile.write("Error file " + aoikml + " not found\r\n")
            abnormalend()

        # Parsing the extracted coordinates string to create the polygon list
        # -------------------------------------------------------------------
        polygonstring = polygonstring.replace("\t", "")
        polygonstring = polygonstring.replace(",0 ", ",")
        polygonstring = polygonstring.replace("\n", "")
        polygonstring = polygonstring[:-1]
        lonlatpoly = extractlist(polygonstring, ",")  # cleanup from garbage

        polygon = []  # initial state: empty list
        for index in range(0, len(lonlatpoly), 2):  # inversing lat/lon pairs
            polygon = polygon + [lonlatpoly[index + 1]]
            polygon = polygon + [lonlatpoly[index]]
        logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
        logfile.write(" PhotoFlightPlan V1.07: ")
        logfile.write("Using AOI polygon from " + aoikml + " \r\n")
    else:
        logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
        logfile.write(" PhotoFlightPlan V1.07: ")
        logfile.write("Using AOI polygon provided within parameters file\r\n")

    # Check if polygon list is empty, meaning we found no polygon points
    # ------------------------------------------------------------------
    if polygon == []:
        logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
        logfile.write(" PhotoFlightPlan V1.07: ")
        logfile.write("No AOI polygon points found, please check input\r\n")
        abnormalend()

    # Opening the flight plan file for flight data and report
    # This file is opened by Collimator to get its project info
    # ---------------------------------------------------------
    pfpname = parname[:-4] + "-" + flightblock + ".pfp"
    kmloutputname = parname[:-4] + "-" + flightblock + ".kml"
    planfile = open(pfpname, "w")

    # Writing down the block header in the flight plan file (.pfp)
    # ------------------------------------------------------------
    planfile.write("================================================\r\n")
    planfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
    planfile.write(" PhotoFlightPlan V1.07 output\r\n")
    planfile.write("================================================\r\n")

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
        loghead = time.strftime("%Y.%m.%d %H:%M:%S")
        loghead = loghead + " PhotoFlightPlan V1.07: "
        msg = "Error - cannot compute altitudeaglft"
        msg = msg + " AND lensfocal at the same time\r\n"
        logfile.write(loghead + msg)
        abnormalend()
    if altitudeaglft == "compute":
        lensfocal = float(lensfocal)
        altitudeaglm = imagewidth * lensfocal / sensorxmm
        altmode = " (computed)\r\n"
        loghead = time.strftime("%Y.%m.%d %H:%M:%S")
        loghead = loghead + " PhotoFlightPlan V1.07: "
        logfile.write(loghead)
        logfile.write("Given lensfocal, flight altitude computed\r\n")
    else:
        altitudeaglft = float(altitudeaglft)
        altmode = " (imposed)\r\n"
    if lensfocal == "compute":
        altitudeaglm = float(altitudeaglft) * 0.3048
        lensfocal = altitudeaglm * sensorxmm / imagewidth
        lensmode = " (computed)\r\n"
        loghead = time.strftime("%Y.%m.%d %H:%M:%S")
        loghead = loghead + " PhotoFlightPlan V1.07: "
        logfile.write(loghead)
        logfile.write("Given flight altitude, lens focal computed\r\n")
    else:
        lensfocal = float(lensfocal)
        lensmode = " (imposed)\r\n"
    altitudeaglft = altitudeaglm / 0.3048

    # Write down prameters in the flight plan file (.pfp)
    # ---------------------------------------------------
    planfile.write("parameters file:   " + parname + "\r\n")
    planfile.write("flight plan file:  " + pfpname + " (this file)\r\n")
    planfile.write("Google Earth file: " + kmloutputname + "\r\n")
    planfile.write("flightblock:       " + flightblock + "\r\n")
    planfile.write("aoiname:           " + str(aoiname) + "\r\n")
    planfile.write("aoiextent:         " + str(round(aoiextent)) + " m\r\n")
    planfile.write("lensfocal:         ")
    planfile.write(str(round(lensfocal, 1)) + " mm" + lensmode)
    planfile.write("resolution:        ")
    planfile.write(str(round(resolution, 1)) + " cm\r\n")
    planfile.write("sensor_xpix:       ")
    planfile.write(str(round(sensorxpix)) + " pixels\r\n")
    planfile.write("sensor_xmm:        ")
    planfile.write(str(round(sensorxmm, 1)) + " mm\r\n")
    planfile.write("sensor_ypix:       ")
    planfile.write(str(round(sensorypix)) + " pixels\r\n")
    planfile.write("sensor_ymm:        ")
    planfile.write(str(round(sensorymm, 1)) + " mm\r\n")
    planfile.write("wind_direction:    ")
    planfile.write(str(round(winddir)) + " degrees\r\n")
    planfile.write("overlap:           " + str(round(overlap)) + "%\r\n")
    planfile.write("sidelap:           " + str(round(sidelap)) + "%\r\n")
    planfile.write("imagewidth:        " + str(round(imagewidth)) + " m\r\n")
    planfile.write("imageheight:       " + str(round(imageheight)) + " m\r\n")
    planfile.write("altitudeaglm:      ")
    planfile.write(str(round(altitudeaglm)) + " m" + altmode)
    planfile.write("altitudeaglft:     ")
    planfile.write(str(round(altitudeaglft)) + " ft" + altmode)
    planfile.write("triggerdist:       " + str(round(triggerdist)) + " m\r\n")
    planfile.write("corridorwidth:     ")
    planfile.write(str(round(corridorwidth)) + " m\r\n")

    # Find the approximate gravity center of the polygon
    # flight lines will be spread over this central point
    # ---------------------------------------------------
    latmin, latmax, lonmin, lonmax = 90, -90, 180, -180
    for index in range(0, len(polygon) - 2, 2):  # last pair same as first
        if polygon[index] < latmin: latmin = polygon[index]
        if polygon[index] > latmax: latmax = polygon[index]
        if polygon[index + 1] < lonmin: lonmin = polygon[index + 1]
        if polygon[index + 1] > lonmax: lonmax = polygon[index + 1]
    gclat = (latmin + latmax) / 2  # gravity center latitude
    gclon = (lonmin + lonmax) / 2  # gravity center longitude

    # Compute the scan width and angles in all directions
    # Using AOI maximum size measured from lat lon extents
    # ---------------------------------------------------
    brg12, brg21, maxdist = vincentyinverse(latmin, lonmin, latmax, lonmax)
    scanwidth = (int((maxdist * 2.0) / corridorwidth) + 1) * corridorwidth
    scanstr = str(scanwidth)
    print("Scanning width: " + scanstr + "m by " + scanstr + "m")
    logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
    logfile.write(" PhotoFlightPlan V1.07: ")
    logfile.write("Scanning width: " + scanstr + "m by " + scanstr + "m\r\n")

    northacross = winddir + 90
    if northacross > 90 and northacross < 270:  # make sure going Northwards
        northacross = northacross + 180
    northacross = northacross % 360  # and that it does not exceed 360 degrees

    # Move along the lines to find where we enter and exit the polygon
    # ---------------------------------------------------------------
    flightlines = []
    linelength = 0

    # Outer loop where we move by corridor width from south to north (+/-)
    # starting well off the polygon, finishing well off also
    # this routine good for AOIs not bigger than 100 corridor widths
    # if by any chance it is the case, we can change scanwidth to something
    # a greater multiple of corridorwidth, will just take longer to process
    # --------------------------------------------------------------------
    start = -int((scanwidth / 2))  # start with a negative scan distance
    end = -start + corridorwidth  # do the last in range also
    for distacross in range(start, end, corridorwidth):
        if distacross == 0:  # at central pt, cannot compute at distance zero
            reflat = gclat
            reflon = gclon
        else:
            reflat, reflon, brg = vincentydirect(gclat, gclon,\
              northacross, distacross)

        # Scan process indication - increasing percentage on the same line
        # tricky assembly to satisfy IOS Pythonista app that is a bit fussy
        # formatting the dynamic one-liner overwrites
        # ----------------------------------------------------------------
        msghd = '\rScanning across AOI "'+ aoiname
        msghd = msghd + '" for block ' + flightblock
        percent = round(distacross + scanwidth / 2) / scanwidth
        msg = msghd + " ({:2.0%}) done".format(percent)  # display for % done
        print(msg, end = "")

        # Initial status for inner loop
        # -----------------------------
        justinlat, justinlon, justoutlat, justoutlon = 0, 0, 0, 0
        distin, distout = 0, 0
        distincrement = 10  # scan along at 10 meters increments

        # Inner loop where we move by 1/100 corridor width increments
        # along a potential flight line at the beginning and at the end,
        # we may not cross the perimeter but when we do, the positive returns
        # from the "pointinpolygon" function will confirm we are crossing it
        # and therefore record a flight line
        # -----------------------------------------------------------------------------
        inpolygon = False  # flag if inside or outside polygon - initial status
        for distalong in range(start, end, distincrement):
            if distalong == 0:  # at central pt, cannot compute with dist. zero
                chklat, chklon = reflat, reflon
            else:
                chklat, chklon, brg = vincentydirect(reflat, reflon,\
                  winddir, distalong)
            inside = pointinpolygon(chklat, chklon, polygon)  # check if within
            if inside:
                inpolygon = True  # true until we get out
                if justinlat == 0:  # if we just entering into the polygon
                    justinlat, justinlon = chklat, chklon  # remember lat-lon
                    distin = distalong  # distance where we got in
            else:
                if inpolygon:  # if previously in, getting out of the perimeter
                    justoutlat, justoutlon = chklat, chklon  # remember lat-lon
                    distout = distalong  # distance at which we got out
                    inpolygon = False  # confirm out but we might enter again!

        if distout != 0:  # if we got out somehow, we crossed the perimeter
            point1 = Point(justinlon, justinlat)  # encode 1st flight line pt
            point2 = Point(justoutlon, justoutlat)  # encode second
            flightlines = flightlines + [justinlat] + [justinlon]
            flightlines = flightlines + [justoutlat] + [justoutlon]
            linelength = linelength + distout - distin + 2.0 * aoiextent

    # Finish initial computations and add to report
    # ---------------------------------------------
    planfile.write("linear coverage:   " + str(round(linelength)) + " m\r\n")
    frames = linelength / triggerdist
    planfile.write("frame count:       " + str(round(frames)))
    planfile.write(" photos (approx)\r\n")
    nbflightlines = int(len(flightlines) / 4)

    # Write down the AOI polygon data and where it comes from
    # -------------------------------------------------------
    planfile.write("\r\n--------------------------------\r\n")
    planfile.write("AOI Polygon points (" + str(aoinbpoints) + "):\r\n")
    if aoikml != "none": planfile.write("From file: " + aoikml + "\r\n")
    else: planfile.write("From parameters file\r\n")
    planfile.write("--------------------------------\r\n")

    for index in range(0, len(polygon), 2):
        if index == len(polygon) - 2:  # note: last repeats first
            planfile.write("polygon: 01 ")
        else:
            planfile.write("polygon: ")
            planfile.write(str(int((index/2 + 1))).zfill(2) + " ")
        planfile.write(str(round(polygon[index],6)).ljust(9, "0") + " ")
        planfile.write(str(round(polygon[index + 1], 6)).ljust(10, "0"))
        planfile.write("\r\n")
    planfile.write("--------------------------------\r\n")

    # Format and write the flight lines to the flight plan file
    # ---------------------------------------------------------
    planfile.write("\r\nFlight lines (" + str(nbflightlines) + "):\r\n")
    planfile.write("-------------------------------")
    planfile.write("----------------------------\r\n")

    for index in range(0, len(flightlines), 4):
        planfile.write("flightline: " + flightblock + " ")
        planfile.write(str(int((index / 4 + 1))).zfill(2) + " ")
        planfile.write(str(round(flightlines[index],6)).ljust(9, "0"))
        planfile.write(" ")
        planfile.write(str(round(flightlines[index+1],6)).ljust(10, "0") + " ")
        planfile.write(str(round(flightlines[index+2],6)).ljust(9, "0") + " ")
        planfile.write(str(round(flightlines[index+3],6)).ljust(10, "0"))
        planfile.write("\r\n")
    planfile.write("-------------------------------")
    planfile.write("----------------------------\r\n\r\n")
    planfile.close()

    # Generate the kml file for Google Earth
    # First, we put the kml header and style data
    # warning: erases any previous kml file without notice
    # ----------------------------------------------------
    kmlfile = open(kmloutputname, "w")
    kmlfile.write('<?xml version="1.0" encoding="UTF-8"?>\r\n')
    kmlfile.write('<kml xmlns="http://www.opengis.net/kml/2.2"')
    kmlfile.write(' xmlns:gx="http:')
    kmlfile.write('//www.google.com/kml/ext/2.2" xmlns:')
    kmlfile.write('kml="http://www.opengis.net')
    kmlfile.write('/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">')
    kmlfile.write('\r\n')
    kmlfile.write("<Document>\r\n")
    kmlfile.write("<name>" + kmloutputname + "</name>\r\n")
    kmlfile.write('<StyleMap id="kildirtech">\r\n')
    kmlfile.write("</StyleMap>\r\n")
    kmlfile.write('<Style id="flightplan1">\r\n')
    kmlfile.write("<LineStyle>\r\n")
    kmlfile.write("<color>ffff00ff</color>\r\n")  # opacity/blue/green/red
    kmlfile.write("<width>3</width>\r\n")
    kmlfile.write("</LineStyle>\r\n")
    kmlfile.write("</Style>\r\n")
    kmlfile.write('<Style id="flightplan2">\r\n')
    kmlfile.write("<LineStyle>\r\n")
    kmlfile.write("<color>ffff0000</color>\r\n")  # opacity/blue/green/red
    kmlfile.write("<width>5</width>\r\n")
    kmlfile.write("</LineStyle>\r\n")
    kmlfile.write("<PolyStyle>\r\n")
    kmlfile.write("<color>33ff0000</color>\r\n")  # opacity/blue/green/red
    kmlfile.write("</PolyStyle>\r\n")
    kmlfile.write("</Style>\r\n")
    kmlfile.write("<Folder>\r\n")
    kmlfile.write("<name>Collimator flight plan</name>\r\n")
    kmlfile.write("<open>1</open>\r\n")

    # Insert the AOI polygon
    # ----------------------
    kmlfile.write("<Placemark>\r\n")
    kmlfile.write("<name>AOI " + aoiname + "</name>\r\n")
    kmlfile.write("<styleUrl>#flightplan2</styleUrl>\r\n")
    kmlfile.write("<Polygon>\r\n")
    kmlfile.write("<tessellate>1</tessellate>\r\n")
    kmlfile.write("<outerBoundaryIs>\r\n")
    kmlfile.write("<LinearRing>\r\n")
    kmlfile.write("<coordinates>\r\n")

    for index in range(0, len(polygon), 2):
        kmlfile.write(str(round(polygon[index+1],6)).ljust(10, "0") + ",")
        kmlfile.write(str(round(polygon[index],6)).ljust(9, "0") + ",0 ")

    kmlfile.write("</coordinates>\r\n")
    kmlfile.write("</LinearRing>\r\n")
    kmlfile.write("</outerBoundaryIs>\r\n")
    kmlfile.write("</Polygon>\r\n")
    kmlfile.write("</Placemark>\r\n")

    # Insert flight lines, one folder per block
    # -----------------------------------------
    kmlfile.write("<Folder>\r\n")
    kmlfile.write("<name>" + "Block " + flightblock + "</name>\r\n")
    kmlfile.write("<open>1</open>\r\n")
    for index in range(0, len(flightlines), 4):
        kmlfile.write("<Placemark>\r\n")
        kmlfile.write("<name>Line " + str(int((index/4+1))).zfill(2))
        kmlfile.write("</name>\r\n")
        kmlfile.write("<styleUrl>#flightplan1</styleUrl>\r\n")
        kmlfile.write("<LineString>\r\n")
        kmlfile.write("<tessellate>1</tessellate>\r\n")
        kmlfile.write("<coordinates>\r\n")
        kmlfile.write(str(round(flightlines[index + 1],6)).ljust(10, "0"))
        kmlfile.write("," + str(round(flightlines[index],6)).ljust(9, "0"))
        kmlfile.write(",0 ")
        kmlfile.write(str(round(flightlines[index + 3],6)).ljust(10, "0"))
        kmlfile.write("," + str(round(flightlines[index + 2],6)).ljust(9, "0"))
        kmlfile.write(",0 \r\n")
        kmlfile.write("</coordinates>\r\n")
        kmlfile.write("</LineString>\r\n")
        kmlfile.write("</Placemark>\r\n")
    kmlfile.write("</Folder>\r\n")

    # Insert the kml file tail and close it
    # -------------------------------------
    kmlfile.write("</Folder>\r\n")
    kmlfile.write("</Document>\r\n")
    kmlfile.write("</kml>\r\n")
    kmlfile.close()

    # Print and log summary for this block
    # ------------------------------------
    msg = "\nFound and recorded " + str(nbflightlines)
    msg = msg + " flight lines for block " + str(flightblock)
    print(msg)
    logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
    logfile.write(" PhotoFlightPlan V1.07: " + str(nbflightlines))
    logfile.write(" flight lines computed for block " + flightblock + "\r\n")
    logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
    logfile.write(" PhotoFlightPlan V1.07: ")
    logfile.write("Produced Flight Plan in " + pfpname + "\r\n")
    logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
    logfile.write(" PhotoFlightPlan V1.07: ")
    logfile.write("Produced Google kml in " + kmloutputname + "\r\n")

    # While loop end
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
if paramerror != "":  # if paramerror not empty, something is wrong
    loghead = time.strftime("%Y.%m.%d %H:%M:%S") + " PhotoFlightPlan V1.07: "
    logfile.write(loghead + paramerror)
    logfile.write(loghead + "Abnormal End\r\n\e\n")
    print(loghead, "Abnormal End - see .log file for error status\r\n\r\n")
    exit()

# Cleanup, say goodbye and exit
# -----------------------------
msg = 'Flight plan(s) generated in "' + parname[:-4] + "-xx.pfp"
msg = msg + '" and "' + parname[:-4] + '-xx.kml"\r\n'
print(msg)
logfile.write(time.strftime("%Y.%m.%d %H:%M:%S"))
logfile.write(" PhotoFlightPlan V1.07: Normal  End\r\n\r\n")
logfile.close()
exit()
