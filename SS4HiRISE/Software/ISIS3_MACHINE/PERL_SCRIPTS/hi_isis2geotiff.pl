#!/usr/bin/perl

################################################################################
# NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE
#
#      This script is not supported by ISIS.
#      If you have problems please contact the Astrogeology Photogrammetry group
#      at PlanetaryPhotogrammetry@usgs.gov
#
# NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE NOTICE
################################################################################

use Getopt::Long;
use FileHandle;         # Buffer Flushing

my ($progname) = ($0 =~ m#([^/]+)$#);  # get the name of this program

my $usage = "
Command:  $progname -ortho <ortho.cub> -dem <dem.cub> -drad <deltaradius.cub> -h
where
   ortho = orthoimage (isis2 or isis3 cube)
   dem = digital elevation model (isis2 or isis3 cube)
   drad = delta_radius (isis2 or isis3 cube)
   h = help 

Description:

       $progname converts isis2 or isis3 Mars orthoimages,
       DEMS and/or delta_radius cubes to geotiff images.

       Output geotiff names will be constructed for you by stripping the
       isis-version designator (i.e. _isis3 or _isis2, case-insensitive)
       from the end of the base name (if it exists), and replacing the
       extension with .tif.  Also created are tiff world files with the
       .tfw extension.

       A report of errors encountered in the processing goes to file:
       \"hi_isis2geotiff.err\".

**************************************************************************
**************************************************************************
NOTICE:
       This script is not supported by ISIS.
       If you have problems please contact the Astrogeology Photogrammetry group
       at PlanetaryPhotogrammetry\@usgs.gov
**************************************************************************
**************************************************************************
";

#####################################################################
#  MAIN APPLICATION SECTION
#  Author: Elpitha Howington-Kraus
#  Date:   January 22 2009
#  Version: 1.1
#  History: Jan 22 2009 - E Howington-Kraus, USGS, Flagstaff Original Version
#           Feb  2 2009 - EHK, modified to strip isis version designator
#                         from output file name
#           Jul 15 2015 - EHK, updated contact information
#####################################################################

#--------------------------------------------------------------------
# Forces a buffer flush after every print, printf, and write on the
# currently selected output handle.  Let's you see output as it's
# happening.
#---------------------------------------------------------------------
#   $| = 1;

#---------------------------------------------------------------------
# Check the argument list
#---------------------------------------------------------------------

   if (@ARGV == 0)
     {
       print "$usage\n";
       exit;
     }

   my $help = '';          # Help option

   # Get options
   my $opt     = GetOptions ( "h"        => \$help,
                              "ortho=s"  => \$orthoCub,
                              "dem=s"    => \$demCub,
                              "drad=s"   => \$dradCub);

   if (!$opt)
     {
       print "$usage\n"; 
       exit;
     }

   if ($help)
     {
       print "$usage\n";
       exit;
     } 

#---------------------------------------------------------------------
# If the "hi_isis2geotiff.err" file exist, delete it
#---------------------------------------------------------------------

   if (-e "hi_isis2geotiff.err") {unlink("hi_isis2geotiff.err");}

#---------------------------------------------------------------------
# Open LOG file
#---------------------------------------------------------------------

   $log = "dem2deltaradii.err";
   open (LOG,">$log") or die "\n Cannot open $log\n";

#---------------------------------------------------------------------
# If an orthoimage was supplied, convert it to a geotiff
#---------------------------------------------------------------------

   if (defined $orthoCub)
     {
       $orthoCub =~ s/\s+//;
       chomp $orthoCub;

       $cubExt = index($orthoCub,".cub");
       if ($cubExt > 0)
         { $core_name = substr($orthoCub,0,$cubExt); }
       else
         {
           $core_name = $orthoCub;
           $orthoCub = $orthoCub . ".cub";
         }

       if (!-e $orthoCub)
         {
	   ReportErrAndDie ("[ERROR] $orthoCub does not exist\n");
           exit;
         }

       print "\nconverting orthoimage: $orthoCub....\n";

       $uc_core_name = uc($core_name);
       $isisDesignator = index($uc_core_name,"_ISIS");

       if ($isisDesignator > 0)
         {$orthoTif = substr($core_name,0,$isisDesignator) . ".tif";}
       else
         {$orthoTif = $core_name . ".tif";}

       $cmd = "gdal_translate -of Gtiff -co \"tfw=YES\" $orthoCub $orthoTif";
       system($cmd) == 0 || ReportErrAndDie ("gdal_translate failed on $orthoCub");

       print "completed conversion to geotiff: $orthoTif\n";
     }

#---------------------------------------------------------------------
# If a dem was supplied, convert it to a geotiff
#---------------------------------------------------------------------

   if (defined $demCub)
     {
       $demCub =~ s/\s+//;
       chomp $demCub;

       $cubExt = index($demCub,".cub");
       if ($cubExt > 0)
         { $core_name = substr($demCub,0,$cubExt); }
       else
         {
           $core_name = $demCub;
           $demCub = $demCub . ".cub";
         }

       if (!-e $demCub)
         {
	   ReportErrAndDie ("[ERROR] $demCub does not exist\n");
           exit;
         }

       print "\nconverting DEM: $demCub....\n";

       $uc_core_name = uc($core_name);
       $isisDesignator = index($uc_core_name,"_ISIS");

       if ($isisDesignator > 0)
         {$demTif = substr($core_name,0,$isisDesignator) . ".tif";}
       else
         {$demTif = $core_name . ".tif";}

       $cmd = "gdal_translate -of Gtiff -co \"tfw=YES\" $demCub $demTif";
       system($cmd) == 0 || ReportErrAndDie ("gdal_translate failed on $demCub");
       print "completed conversion to geotiff: $demTif\n";

     }

#---------------------------------------------------------------------
# If a delta-radius cube was supplied, convert it to a geotiff
#---------------------------------------------------------------------

   if (defined $dradCub)
     {
       $dradCub =~ s/\s+//;
       chomp $dradCub;

       $cubExt = index($dradCub,".cub");
       if ($cubExt > 0)
         { $core_name = substr($dradCub,0,$cubExt); } 
       else
         {
           $core_name = $dradCub;
           $dradCub = $dradCub . ".cub";
         }

       if (!-e $dradCub)
         {
           ReportErrAndDie ("[ERROR] $dradCub does not exist\n");
           exit;
         }

       print "\nconverting delta-radius cube: $dradCub....\n";

       $uc_core_name = uc($core_name);
       $isisDesignator = index($uc_core_name,"_ISIS");

       if ($isisDesignator > 0)
         {
           $dradVrt = substr($core_name,0,$isisDesignator) . ".vrt";
           $dradTif = substr($core_name,0,$isisDesignator) . ".tif";
         }
       else
         {
           $dradVrt = $core_name . ".vrt";
           $dradTif = $core_name . ".tif";
         }

       $cmd = "gdal_translate -of VRT $dradCub $dradVrt";
       system($cmd) == 0 || ReportErrAndDie ("gdal_translate failed on $dradCub");
       #---------------------------------------------------------------
       #Get number of lines in dradVrt by loading the file into an
       #array
       #---------------------------------------------------------------

       open (IN,$dradVrt) || ReportErrAndDie ("[Error] Problem opening input file: $dradVrt, $!\n");
       @vrt_lines = <IN>;
       close IN;

       #---------------------------------------------------------------
       # Insert Scale and Offset values into the VRT file using pop and
       # push to add these lines to the array @vrt_lines, and then
       # outputing the revised array
       #---------------------------------------------------------------

       $last_line = pop(@vrt_lines);
       $next2last_line = pop(@vrt_lines);
       push(@vrt_lines, "    <Scale>1.0</Scale>\n");
       push(@vrt_lines, "    <Offset>3396000.0</Offset>\n");
       push(@vrt_lines, $next2last_line);
       $vrt_nl = push(@vrt_lines, $last_line);  #push returns the size of the
                                                #array....THIS IS ZERO-BASED!

       open (OUT,">$dradVrt") || ReportErrAndDie ("[Error] Problem opening output file: $dradVrt, $!\n");
       for ($i=0; $i<=$vrt_nl; $i++)
           { print OUT $vrt_lines[$i]; }

       $cmd = "gdal_translate -of Gtiff -co \"tfw=YES\" $dradVrt $dradTif";
       system($cmd) == 0 || ReportErrAndDie ("gdal_translate failed on $dradVrt");
       print "completed conversion to geotiff: $dradTif\n";
     }

#---------------------------------------------------------------------
# Close the LOG file.
# If an error was detected, print out the log file
#---------------------------------------------------------------------

   close (LOG);

   @lines = `cat $log`;
   if (scalar(@lines) > 0)
      {
      print "\n*** Errors detected in processing ***\n\n";
      print @lines;
      print "\n";
      print "\n*** See hinoproj.prt and hi4socet.prt for details ***\n\n";
      }
   else
      {
      unlink ($log);
      }

#--------------------------------
# Delete temporary files and exit
#--------------------------------

   unlink ($dradVrt);

   exit;

##############################################################################
#  Error Handling Subroutine
##############################################################################
sub ReportErrAndDie
    {
    my $ERROR=shift;

    print "$ERROR\n";
    print "$progname aborted\n";

    print LOG "$ERROR\n";
    close(LOG);

    exit 1;

    }

