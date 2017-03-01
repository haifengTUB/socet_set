#!/usr/bin/perl -s
################################
#  program: pedrTAB2SHP.pl
#
#  Perl script to add field headers to PEDR2TAB file and create a OGR CSV and
#  VRT files and then call ogr2ogr to create a 3D Shapefile
#  for OGR see: http://www.gdal.org/ogr/
#
#  Written to help support Socet Set ingestion of MOLA PEDR2TAB points
#
#  Trent Hare, USGS
#  Mar 2008
#  Mar 2017, add in Mars 2000 projection and binary field mapping for ogr2ogr (*.csvt)
#

$filelist1= @ARGV[0];
$skipline = @ARGV[1];
#$outcsv = @ARGV[2];
#$outvrt = @ARGV[3];

#Create tempory name for conversion
$root = "xxtempZ";
$outcsv = $root.".csv";
$outcsvt= $root.".csvt";
$outvrt = $root.".vrt";
$outprj = $root.".prj";
$outshp = $root.".shp";
$outshx = $root.".shx";
$outdbf = $root.".dbf";

chomp $skipline;
$skipline = $skipline * 1;
if ($#ARGV != 1) {
  print "usage: pedrTAB2SHP.pl file1.tab skip_line\n";
  print "\nDescription: Convert MOLA PEDR Tab file to OGR Vrt and then Shapefile \n";
  print "     The output file name will be the input name appended with a 'Z'\n";
  print "     The created shapefile will have a *.shp, *.shx, *.dbf and *.prj filename.\n";
  print "\nexamples:\n";
  print "pedrTAB2SHP.pl input.tab 2\n";
  print "     For default PEDR2TAB output TAB file use a skip_line = 2\n\n";
  exit;
}

chomp $filelist1;
#chomp $outcsv;
#chomp $outvrt;
open FILESIN1, $filelist1;
open OUTCSV, "> $outcsv ";
open OUTCSVT, "> $outcsvt ";
open OUTVRT, "> $outvrt ";

@infname = split('\.',$filelist1);
$inroot = $infname[0];

$iflag = 1;

while (<FILESIN1>) {
  if ($iflag != $skipline) {
    #check to see if comma or other common delimited, replace with space
    $_ =~ s/,/ /gi;
    $_ =~ s/\t/ /gi;
    $_ =~ s/;/ /gi;
    $_ =~ s/:/ /gi;
    @line = split(/\s+/,$_,2000);
    $cnt = @line * 1;
    $flg=1;
    #write out ID
    #print OUTCSV "$iflag ";
    foreach $string (@line) {
        if (length($string) > 0) {
          if ($string eq "long_East") {
              $lonCol = $flg;
          }
          #switch longitudes to -180 to 180 from 0 to 360
          if (($iflag != 1) && ($flg == $lonCol)) {
            $lon = $string * 1;
            if ($lon > 180) {
               $lon = $lon - 360;
           }
            print OUTCSV "$lon";
          } else {
            print OUTCSV "$string";
          }
          $flg = $flg + 1;
          if ($flg < $cnt) {
            print OUTCSV ","
          }
        } else {
          $cnt = $cnt - 1;
        }
    }
    print OUTCSV "\n";
  }
  $iflag = $iflag + 1;
}
close OUTCSV;

#Create CSVT field mapping. Assume input for tab is
#long_East lat_North topography MOLArange  planet_rad c A  offndr  EphemerisTime  areod_lat areoid_rad shot  pkt   orbit gm
$fieldMapping = "\"Real\",\"Real\",\"Real\",\"Real\",\"Real\",\"Integer\",\"Integer\",\"Real\",\"String\",\"Real\",\"Real\",\"Integer\",\"Integer\",\"Integer\",\"Integer\"";
print OUTCSVT $fieldMapping;
close OUTCSVT;

# Create ogr2ogr virtual header (*.vrt)
print OUTVRT  "<OGRVRTDataSource>\n";
print OUTVRT  "    <OGRVRTLayer name=\"$root\">\n";
#print OUTVRT  "        <SrcDataSource relativeToVRT=\"1\">.</SrcDataSource>\n";
print OUTVRT  "        <SrcDataSource>$outcsv</SrcDataSource>\n";
#print OUTVRT  "        <SrcDataSource relativeToVRT=\"1\">.</SrcDataSource>\n";
#print OUTVRT  "        <SrcLayer>$root</SrcLayer>\n";
#print OUTVRT  "        <GeometryType>wkbPoint</GeometryType>\n";
#print OUTVRT "        <LayerSRS>WGS84</LayerSRS>\n";
print OUTVRT  "        <GeometryField encoding=\"PointFromColumns\" x=\"long_East\" y=\"areod_lat\" z=\"topography\"/>\n";
print OUTVRT  "    </OGRVRTLayer>\n";
print OUTVRT  "</OGRVRTDataSource>\n";

close OUTVRT;

#Run OGR2OGR conversion
#@args = ("ogr2ogr", "-f \"ESRI Shapefile\"", ".", "$outcsv");
#$args = "ogr2ogr -f \"ESRI Shapefile\" . $outcsv";
#if (system($args) != 0) {
     #print "system $args failed: $?\n";
#}

#@args = ("ogr2ogr", "-lco \"SHPT=POINTZ\"", "-f \"ESRI Shapefile\"", ".", "$outvrt");
$args = "ogr2ogr -lco SHPT=POINTZ -f \"ESRI Shapefile\" . $outvrt";
if (system($args) != 0) {
     print "system $args failed: $?\n";
}


# create projection file
unlink ("$outprj");
open  OUTPROJ, "> $outprj ";
print OUTPROJ "GEOGCS[\"Mars 2000\",DATUM[\"D_Mars_2000\",SPHEROID[\"Mars_2000_IAU_IAG\",3396190.0,0.0]],PRIMEM[\"Reference_Meridian\",0.0],UNIT[\"Degree\",0.0174532925199433]]\n";
close OUTPROJ;

if (-e $outshp) {
   print " -renaming output files\n";
   rename($outshp, $inroot."Z.shp");
   rename($outshx, $inroot."Z.shx");
   rename($outdbf, $inroot."Z.dbf");
   rename($outprj, $inroot."Z.prj");
   print " -deleting temporary files\n";
   unlink $outcsv;
   unlink $outcsvt;
   unlink $outvrt;
   print "\n Output shapefile file generated: $final\n\n";
} else {
   print "\n Shapefile not generated...something's wrong\n\n";
}

