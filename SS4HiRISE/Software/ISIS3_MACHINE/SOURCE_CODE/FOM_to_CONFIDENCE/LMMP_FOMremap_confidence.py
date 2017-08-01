#!/usr/bin/env python
###############################################################################
# $Id$
#
# lmmp_remap_confidence.py
#  Description: simple remapping program to LMMP confidence values
#  from a socet set/GXP FOM file
#
# SS FOM         LMMP  
# 0-1          0 = NoDATA , outside boundary (e.g. out of stereo pair overlap, not all groups use)  
# 2            1 = shadowed   (if group has capability to assign)
# 21           2 = saturated
# 3,5-20,28,31-39
#              3 = suspicious (edge, corner, did not correlate, other bad value, don't use)  
# 4,30         4 = interpolated / extrapolated (e.g. from neighbor pixels)  
# 40-99          10 - 14 = value range of success; linear fit => poor(10) to best(14)  
# 22-24,25-27,29 15 = manually interpolated - mass-point edit tools in SS interpolate
# n/a            17 = seed points (e.g. tied to LOLA Shot)
#
# For Socet Set a potential mapping for success 40-99:
# 40 - 59 => 10 poor(low) correlation
# 60 - 69 => 11
# 70 - 79 => 12 medium correlation
# 80 - 89 => 13
# 90 - 99 => 14 best (max) correlation
#
#
# Author: Trent Hare, Oct 2010
#   April 4, 2011 - Trent changed FOM 25 to manually interpolated - TH
#   July 19, 2012 - Jay Laura changed numpy.choose to a numpy.iterate (2x as fast)
# based on GDAL/numpy sample by Andrey Kiselev, dron@remotesensing.org
###############################################################################
# Copyright (c) 2003, Andrey Kiselev <dron@remotesensing.org>
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
###############################################################################

try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    gdal.TermProgress = gdal.TermProgress_nocb
except ImportError:
    import gdal
    from gdalconst import *

try:
    import numpy
except ImportError:
    import Numeric as numpy

#to help make script compatible with Python2
try:
    range
except NameError:
    range = xrange

import sys

# =============================================================================
def Usage():
    print('Usage: lmmp_FOMremap_confidence.py infile_FOM outfile.tif\n')
    sys.exit( 1 )

# =============================================================================
def ParseType(type):
    if type == 'Byte':
        return GDT_Byte
    elif type == 'Int16':
        return GDT_Int16
    elif type == 'UInt16':
        return GDT_UInt16
    elif type == 'Int32':
        return GDT_Int32
    elif type == 'UInt32':
        return GDT_UInt32
    elif type == 'Float32':
        return GDT_Float32
    elif type == 'Float64':
        return GDT_Float64
    elif type == 'CInt16':
        return GDT_CInt16
    elif type == 'CInt32':
        return GDT_CInt32
    elif type == 'CFloat32':
        return GDT_CFloat32
    elif type == 'CFloat64':
        return GDT_CFloat64
    else:
        return GDT_Byte
# =============================================================================

infile = None
outfile = None
format = 'GTiff'
type = GDT_Byte
quiet = False

# Parse command line arguments.
i = 1
while i < len(sys.argv):
    arg = sys.argv[i]

    if arg == '-innd':
        i = i + 1
        inNoData = float(sys.argv[i])
    elif arg == '-outnd':
        i = i + 1
        outNoData = float(sys.argv[i])
    elif arg == '-of':
        i = i + 1
        format = sys.argv[i]
    elif arg == '-ot':
        i = i + 1
        type = ParseType(sys.argv[i])
    elif infile is None:
        infile = arg
    elif outfile is None:
        outfile = arg
    else:
        Usage()
    i = i + 1

if infile is None:
    Usage()
if  outfile is None:
    Usage()

indataset = gdal.Open( infile, GA_ReadOnly )

#define output format, name, size, type and set projection
out_driver = gdal.GetDriverByName(format)
outdataset = out_driver.Create(outfile, indataset.RasterXSize, indataset.RasterYSize, indataset.RasterCount, type)

outdataset.SetProjection(indataset.GetProjection())
outdataset.SetGeoTransform(indataset.GetGeoTransform())

#Map the 99+ Socet Set FOM values into LMMP Confidence map using a "map to" indexed array.
#So if FOM pixel=2 (thus index=2) then remap=1; if FOM pixel=31 then remap=3
#Probably not the best way to map ranges but good for singletons...

#These index lines just help with mapping into array
# index #  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
FOMmap =  [0,0,1,3,4,3,3,3,3,3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2,15,15,15,15,15,15, 3,15, 4,
# index # 31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
           3, 3, 3, 3, 3, 3, 3, 3, 3,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
# index # 58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85
          10,10,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,
# index # 86,87,88,89,90,91,92,93,94,95,96,97,98,99
          13,13,13,13,14,14,14,14,14,14,14,14,14,14,
# index # Any values beyond the valid 99 FOM mapped to 0 but they should not exist.
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

#Make the values list into a numpy array
FOMarray = numpy.array(FOMmap)

for iBand in range(1, indataset.RasterCount + 1):
    inband = indataset.GetRasterBand(iBand)
    outband = outdataset.GetRasterBand(iBand)

    for i in range(inband.YSize - 1, -1, -1):
        scanline = inband.ReadAsArray(0, i, inband.XSize, 1, inband.XSize, 1)

        #Numpy.choose bug (some versions) which only allows 32 elements - so changed to iteration
        #scanline = numpy.choose(scanline, FOMmap, mode='clip')
        #Let numpy iterate and replace...
        scanline[:] = FOMarray[scanline]  
        
        outband.WriteArray(scanline, 0, i)
        
        #update progress line
        if not quiet:
            gdal.TermProgress_nocb( (float(inband.YSize - i + 1) / inband.YSize ) )
