#!/usr/bin/perl -s
###############################################################################
#
#_TITLE  gdal_FOM_to_CONFIDENCE.pl - uses gdal Python and gdaldem  
#              to convert a SOCET FOM to LMMP confidence map(s) 
#              Needs gdal 1.7.x binaries with python and numpy
#              written for Mars2020 (although same process used for LMMP).
#
#_ARGS  
#  input image (most GDAL supported formats)
#  output grayscale GeoTiff and colorized GeoTiff
#
#_USER  Command line entry options [optional parameters]:
#
#   gdal_FOM_to_CONFIDENCE.pl input_FOM.cub
#
# Requirements:
#       GDAL Library (gdaldem)
#
#_DESC make CONFIDENCE and colorized confidence with copied legend.
#
#_CALL  List of calls:
#
#_HIST
#       July 21 2017 - Trent Hare - original version
#   
#_END
#######################################################################
######################################################################
# For help - user can enter this perl script and return
######################################################################
   if ($#ARGV < 0) {
      print " \n\n          *** HELP ***\n\n";
      print "gdal_FOM_to_CONFIDENCE.pl -  Create 8bit confidence and colorized confidence with legend from SS FOM\n\n";
      print "Command line: \n";
      print "gdal_FOM_to_CONFIDENCE.pl input_FOM.cub\n";
      exit 1;
    }

    $input = $ARGV[0];
    chomp $input;
    @fname = split('\.',$input);
    $root = $fname[0];
    $root =~ s/_isis3//;
    $root =~ s/_FOM//;
    $root =~ s/FOM_//;

    ###############################################
    # process  CONFIDENCE
    ###############################################
    $outConf  =  "CONF_".$root.".tif";
    $outClrConf  =  "ClrCONF_".$root.".tif";

    #if (-e $output) {
       #print "\nOutput file $output Exists! Please remove file and run again.\n\n";
       #exit -1;
    #}
    #calculate grayscale confidence
    $cmd = "/usgs/dev/contrib/bin/LMMP_FOMremap_confidence.py $input $outConf";
    print $cmd."\n";
    system($cmd);

    #calculate colrized confidence
    $cmd = "gdaldem color-relief $outConf /usgs/dev/contrib/bin/LMMP_color_confidence.lut $outClrConf";
    print $cmd."\n";
    system($cmd);
  
    #add existing PNG legend to top left corner
    $legend = "ClrCONF_".$root."_legend.png";
    $cmd = "cp /usgs/cdev/contrib/bin/LMMP_Confidence_legend.png $legend";
    print $cmd."\n";
    system($cmd);
  
    print "***  - files $outConf, $outClrConf, and $legend created\n\n";
