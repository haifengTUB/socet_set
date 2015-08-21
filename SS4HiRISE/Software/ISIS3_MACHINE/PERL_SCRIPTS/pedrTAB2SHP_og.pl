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

$filelist1= @ARGV[0];
$skipline = @ARGV[1];
#$outtxt = @ARGV[2];
#$outvrt = @ARGV[3];
$outtxt = "xxtemp.csv";
$outvrt = "xxtemp.vrt";


chomp $skipline;
$skipline = $skipline * 1;
if ($#ARGV != 1) {
  print "usage: pedrTAB2SHP.pl file1.tab skip_line\n";
  print "\nDescription: Convert MOLA PEDR Tab file to OGR Vrt and then Shapefile \n";
  print "     The output file name will be the input name appended with a 'Z'\n";
  print "     The created shapefile will have a *.shp, *.shx, and *.dbf filename.\n";
  print "\nexamples:\n";
  print "pedrTAB2SHP.pl input.tab 2\n";
  print "     For default PEDR2TAB output TAB file use a skip_line = 2\n\n";
  exit;
}

chomp $filelist1;
chomp $outtxt;
chomp $outvrt;
open FILESIN1, $filelist1;
open OUTTXT, "> $outtxt ";
open OUTVRT, "> $outvrt ";

@fname = split('\.',$outtxt);
$root = $fname[0];
$outshp = $fname[0]."Z";
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
    #print OUTTXT "$iflag ";
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
            print OUTTXT "$lon";
          } else {
            print OUTTXT "$string";
          }
          $flg = $flg + 1;
          if ($flg < $cnt) {
            print OUTTXT ","
          }
        } else {
          $cnt = $cnt - 1;
        }
    }
    print OUTTXT "\n";
  }
  $iflag = $iflag + 1;
}

print OUTVRT  "<OGRVRTDataSource>\n";
print OUTVRT  "    <OGRVRTLayer name=\"$outshp\">\n";
print OUTVRT  "        <SrcDataSource relativeToVRT=\"1\">.</SrcDataSource>\n";
print OUTVRT  "        <SrcLayer>$root</SrcLayer>\n";
print OUTVRT  "        <GeometryType>wkbPoint</GeometryType>\n";
#print OUTVRT "        <LayerSRS>WGS84</LayerSRS>\n";
print OUTVRT  "        <GeometryField encoding=\"PointFromColumns\" x=\"long_East\" y=\"areod_lat\" z=\"topography\"/>\n";
print OUTVRT  "    </OGRVRTLayer>\n";
print OUTVRT  "</OGRVRTDataSource>\n";

close OUTTXT;
close OUTVRT;

#Run OGR2OGR conversion
#@args = ("ogr2ogr", "-f \"ESRI Shapefile\"", ".", "$outtxt");
$args = "ogr2ogr -f \"ESRI Shapefile\" . $outtxt";
if (system($args) != 0) {
     print "system $args failed: $?\n";
}

#@args = ("ogr2ogr", "-lco \"SHPT=POINTZ\"", "-f \"ESRI Shapefile\"", ".", "$outvrt");
$args = "ogr2ogr -lco SHPT=POINTZ -f \"ESRI Shapefile\" . $outvrt";
if (system($args) != 0) {
     print "system $args failed: $?\n";
}

$outshp2 = $root."Z.shp";
$outshx2 = $root."Z.shx";
$outdbf2 = $root."Z.dbf";
$outdbf = $root.".dbf";
if (-e $outshp2) {
   $final = $inroot."Z.shp";
   print " -renaming output files\n";
   rename($outshp2, $final);
   rename($outshx2, $inroot."Z.shx");
   rename($outdbf2, $inroot."Z.dbf");
   print " -deleting temporary files\n";
   unlink $outtxt;
   unlink $outvrt;
   unlink $outdbf;
   print "\n Output shapefile file generated: $final\n\n";
} else {
   print "\n Shapefile not generated...something's wrong\n\n";
}

