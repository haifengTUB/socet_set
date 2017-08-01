# Astrogeology Socet Set Source Code

Requires Python 2 or Python 3 and GDAL/Python. Recommendation is to 
install Anaconda Python and once installed, type "conda instal gdal".

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

Usage: 
lmmp_FOMremap_confidence.py infile_FOM.cub outfile_CONF.tif

Note the PERL script gdal_FOM_to_CONFIDENCE.pl just runs this Python script and gdaldem to colorize.