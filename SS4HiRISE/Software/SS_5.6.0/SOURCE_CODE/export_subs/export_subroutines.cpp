////////////////////////////////////////////////////////////////////////////////
//
//_Title export_subroutines
//
//_Desc  Subroutines in common with dem2isis3 and ortho2isis3
//
//_Hist Sep 20 2010 EHK  Orig Version extracted from dem2isis3 and ortho2isis3
//                       -Added notes in isis_script describing sections
//                       -Updated ISIS commands for isis3.2.0
//      Sep 28 2010 EHK  When creating FOM filename, added search for DTM
//                       along with DEM (to be replaced by FOM)
//      Feb 11 2011 EHK  Added triaxial ellipsoids that are treated
//                       cartographically as spheres (as defined by IAU in 2006)
//                       (see http://planetarynames.wr.usgs.gov/TargetCoordinates).
//                       For portability to other institutions, also changed approach
//                       in using a local pck file for storing spherical radius, to
//                       explicitly suppling the radii values when running
//                       maptemplate
//
//                       Also updated setisis command to isis3.2.1 and added logic that
//                       setisis command only be run on flagstaff machines
//      Oct 17 2011 EHK  Added SS_ prefix to output *.raw files to avoid confusion with
//                       the 'standard' ISIS cubes having the same core name, but different
//                       number of lines and samples
//
//      Jan 10 2012 EHK  Modified to allow export of Polar Stereographic products of
//                       ellipsoidal bodies, and updated setisis command to isis3.3.0
//      Jul 09 2012 EHK  Updated setisis command to isis3.4.0
//      Mar 27 2013 EHK  Updated setisis command to isis3.4.3
//      Nov 21 2013 EHK  Added Sinusodial Projection, and updated setisis command to isis3.4.4
//      Jan 28 2014 EHK, Updated setisis command to isis3.4.5 and added MERCURY_MSGR
//      Jan 29 2014 EHK, Debug problem convering east180 logitudes that cross 180/-180
//                                TEMPORARILY set Dione, Rhea, Enceladus and Tethys to pos East lon
//      Mar 17 2014 EHK, Corrected bug when exporting Geographic products of a Martian sphere.  The
//                                 clon/clat were incorrectly set to the "standard" values used for Mars,
//                                 rather than the clon/clat used in the project.
//      Mar 19 2014 EHK, To be compatible with ISIS3 and ARCMAP defaults, this code will output
//                                 all bodies, EXCEPT TITAN, as positive East longitudes.  TITAN will remain
//                                 positive West longitude because the CASSINI RADAR images are
//                                 map-projected images with postitive West lons.
//                                 NOTE: this means IAU standards are no longer being met, generally, w.r.t.
//                                 positive longitude direction.
//      Apr 09 2014 EHK,  To aid non-Astro users of Socet Set, (1) added support for Jupiter, Saturn, Neptune, Pluto
//                                 and their moons,  and (2) generallized getTargetInfo (details below).  This modification
//                                 includes updating the Geodetic database files so that they currently conform to the IAU 2006
//                                 report and are consistent with ISIS3.  Also note that for triaxial bodies, we are using the
//                                 IAU 2006 Mean Radius for all three axes in Socet Set.  See the
//                                 <ss_install>\internal_dbs\GEODETIC\geodetic.doc file for additional information.  (Although
//                                 this file has a *.doc extention, it is a text file, not an MS Word document.)
//
//                                 In order to generallize getTargetInfo,  I replaced hard coded values of the polar radius
//                                 with a computation of the polar radius based on the values of Socet Set Project parameters
//                                 A_EARTH and E_EARTH.  Also, getTargetInfo has been modified so that if a user adds a
//                                 new  Ellipsoid definition for one of the existing ellipsoids in the Geodetic database files,
//                                 *and keeps the body's name in the new ellipsoid name*, getTargetInfo will not require modification.
//                                 So, for example, for current ellipsoid Rhea2006 in the geodetic database, a user can add
//                                 the new IAU 2009 Rhea definition as Rhea2009 to the database files (or replace Rhea2006
//                                 with Rhea2009), and dem2isis3 & ortho2isis3 will continue to function without the need to modify code.
//      Apr 11 2014 EHK,  Updated logic that setisis also be run on the astro-guest virtual machine.   (The GROUP
//                                 envirnonment for the dpw-user guest account is dpw-user, rather than flagstaff for
//                                 Astro employees.)
//
//                                  Also, to keep this routine as generic as possible, removed use of ISIS default target definitions
//                                  and in all cases, explicitly supply the equatorial and polar radii for the body when running maptemplate.
// 
//                                  To maintain compatability with ARC, and keep IAU pos East/West designations, output ISIS3 cubes
//                                   with a clon of 0, and a londom of 180.  (Trent recommends clon=0 and londom=180.  The alternative is
//                                   clon=180, londom=360.)  For polar stereographic, the Socet Set grids are set up with
//                                   clon = 0.0, so we are ARC compatible, otherwisie resampling would be required to set clon=0.0
//     May 30 2014 EHK, Updated Socet Set Native geographic DEMs to also set clon=0, and a londom of 180 for ARC compatibility,
//                                and added the ARC compatibility details and IAU postive lon direction to the commnets in the output script,
//                                and updated setisis to isis3.4.6
//_End
//
////////////////////////////////////////////////////////////////////////////////


#include <system_includes.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

//SOCET SET
#include <key/handle_key.h>
#include <dtm/dtm.h>
#include <dtmUtil/dtm_util.h>
#include <dtmAccess/DtmGrid.h>
#include <dtmAccess/DtmHeader.h>
#include <dtmAccess/dtm_edit_util.h>
#include <dtmAccess/fom_defs.h>
#include <project/proj.h>
#include <arith/mm_matrix.h>
#include <util/string_dpw.h>
#include <util/init_socet_app.h>
#include <ground_point.h>
#include <image_point.h>

#define FILELEN 512
#define NO_ERRS 0
#define PARINV_ERR -1


// For the record, these are the ISIS NULL values I tried to use under
// windows, but what should have resulted as NULL pixels were LRS
// (round-off???):
//#define NULL3 -0.3402822655089E+39       (value from dem2isis2)
#define NULL3 -3.4028226550889044521e+38 // (modified value Trent uses in gdal) 
//Set ISIS3 NULL value
//const int INULL4 = 0xFF7FFFFB;
//const float NULL4 =(*((const float *) &INULL4));

// prototypes
int parse_label(char *file, char *keyword, char *value);
int getTargetInfo (char *uc_ellipsoid, char *isisTargName,
	                char *ographicPosLonDir, char *ocentricPosLonDir);
int writeToScript(char *isis_script, char *command);

/***********  generate_ss2isis_script **************
*                                                  *
*  This routine generates the ss->isis script for  *
*  dem/fom/ortho export                            *
*                                                  *
*  by E. Howington-Kraus                           *
*                                                  *
****************************************************/
int generate_ss2isis_script (char *isis_script,
                               char *prj, //SS project with full path and extension
                               char *productType,
                               char *byteOrder,
                               char *outcub_name, //DEM or ORTHO output cube name
                               char *layout_flag,
                               int lines,
                               int samples,
                               double x_realspacing,
                               double y_realspacing,
                               double ulcenter_Xlon,  // X (meters) or lon (radians) at center of upperleft pixel
                               double ulcenter_Ylat   // Y (meters) or lat (radians) at center of upperleft pixel 
                               )

{
   // ISIS variables
   double minlat, minlon, maxlat, maxlon, clat, clon, sscale, km_per_dg;
   double east180_minlon, east180_maxlon, west180_minlon, west180_maxlon;
   double east360_minlon, east360_maxlon, west360_minlon, west360_maxlon;
   double east180_clon, west180_clon, east360_clon, west360_clon;
   double ogRef_lat, ogRef_lon, ocRef_lat, ocRef_lon, avglat, temp;
   double xSpacingDG, ySpacingDG; //real spacing of a DEM converted to degrees
   float null;
   char outcub[FILELEN];
   char rawDEM[FILELEN];
   char rawFOM[FILELEN];
   char rawORTHO[FILELEN];
   char socetset_map[FILELEN];
   char standard_map[FILELEN];
   char tempcub[FILELEN];
   char input_cub[FILELEN];
   char input_fomcub[FILELEN];
   char sqrcub[FILELEN];
   char FOM_sqrcub[FILELEN];
   char outcub_tiled[FILELEN];
   char FOM_outcub_name[FILELEN];
   char FOM_outcub[FILELEN];
   char FOM_outcub_tiled[FILELEN];
   char SS_outcub[FILELEN];
   char SS_FOM_outcub[FILELEN];
   char layout_cub[FILELEN];
   
   char note[256];
    
   // ISIS declarations
   char isisTargName[FILELEN];
   char ographicPosLonDir[13];
   char ocentricPosLonDir[13];
   char lattype[15];
   int londom;
   double eqradius =0.0;
   double polradius =0.0;
   
   // Project File Variables
   int coord_sys, xy_units, z_units;
   double ecc;
   char ellipsoid[FILELEN];
   char uc_ellipsoid[FILELEN];  // Socet Set ellipsoid name all in upper case
   char projection[FILELEN];
   char polarAspect[FILELEN];
   double prj_scale;
   int zone;
   
   //Misc
   int ret;
   char value[FILELEN];
   char command[512];
   double rad2deg = 180.0 / M_PI;  //convert deg to radians and back
    
  /////////////////////////////////////////////////////////////////////////////
  // Get Project file parameters that pertain to both
  // geographic and map projected projects...
  /////////////////////////////////////////////////////////////////////////////

   ret = parse_label(prj, "COORD_SYS", value);
   coord_sys = atoi(value);

   ret = parse_label(prj, "XY_UNITS", value);
   xy_units = atoi(value);

   ret = parse_label(prj, "Z_UNITS", value);
   z_units = atoi(value);;
   if (z_units != 1) {
     printf("WARNING: Z units of project must be meters\n");
     printf("         Please make appropriate changes\n");
     printf("         and try again\n");
     exit(1);
   }

   ret = parse_label(prj, "ELLIPSOID", ellipsoid);

   strcpy(uc_ellipsoid,ellipsoid);
   upper_case(uc_ellipsoid);

   ret = parse_label(prj, "A_EARTH", value);
   eqradius = atof(value);

   ret = parse_label(prj, "E_EARTH", value);
   ecc = atof(value);

   // Compute polar radius (needed for ISIS)
   if (ecc == 0.0)
     polradius = eqradius;
   else
     polradius = sqrt(eqradius*eqradius*(1.0-ecc*ecc));
   
   // Get ISIS target name and direction of positive longitude for this planet (ellipsoid)
   ret = getTargetInfo(uc_ellipsoid, isisTargName, ographicPosLonDir,
                       ocentricPosLonDir);
   
   if (ret != 0) {
     printf("No ISIS Target Info for SOCET SET ellipsoid: %s\n",ellipsoid);
      exit (-1);
   }
   
   /////////////////////////////////////////////////////////////////////////////
   //output setisis command compatible with output ISIS script
   /////////////////////////////////////////////////////////////////////////////

   sprintf(command,"set testgroup = `printenv GROUP`");
   writeToScript(isis_script,command);
   sprintf(command,"if ($testgroup == \"flagstaf\" || $testgroup == \"dpw-user\") then ");
   writeToScript(isis_script,command);
   
   sprintf(command,"   setisis isis3.4.6");
   writeToScript(isis_script,command);
   
   sprintf(command,"endif\n");
   writeToScript(isis_script,command);

   ////////////////////////////////////////////////////////////////////////////
   // Form output and temporary file names for all types of products for now
   // (DEM/FOM/ORTHO)
   /////////////////////////////////////////////////////////////////////////////

   // Get outcub (DEM or ORTHO)
   strcpy(outcub,concat(outcub_name,".cub"));

   // Form output FOM file name
   strcpy (FOM_outcub_name,outcub_name);
   upper_case(FOM_outcub_name);
   char* pos=NULL;
   if (pos=strstr(FOM_outcub_name,"DEM")) {
      strcpy (FOM_outcub_name,outcub_name); // Reset original case of filename
      strncpy (pos,"FOM",3);               // Replace "DEM" with FOM, if DEM is in filename
   }
   else if (pos=strstr(FOM_outcub_name,"DTM")) {
      strcpy (FOM_outcub_name,outcub_name); // Reset original case of filename
      strncpy (pos,"FOM",3);             // Replace "DTM" with FOM, if DTM is in filename
   }
   else {       
      strcpy (FOM_outcub_name,"FOM_");   // else, add FOM_ prefix
      strcat (FOM_outcub_name,outcub_name);
   }
   
   // Get raw binary input files
   strcpy(rawDEM,"SS_");
   strcat(rawDEM,outcub_name);
   strcat(rawDEM,".raw");

   strcpy(rawFOM,"SS_");
   strcat (rawFOM,FOM_outcub_name);
   strcat (rawFOM,".raw");
   
   strcpy(rawORTHO,"SS_");
   strcat(rawORTHO,outcub_name);
   strcat(rawORTHO,".raw");

   // Get output FOM filename
   strcpy(FOM_outcub,concat(FOM_outcub_name,".cub"));

   // Add "SS_" prefix to output cube names for the native SS version of the DEM or ORTHO
   strcpy (SS_outcub,"SS_");
   strcat (SS_outcub,outcub);

   strcpy (SS_FOM_outcub,"SS_");
   strcat (SS_FOM_outcub,FOM_outcub);

   // form socetset_map, standard_map, tempcub, input_cub, input_fomcub, sqrcub,
   // FOM_sqrcub, and layout_cub file names

   strcpy(socetset_map,outcub_name);
   strcat(socetset_map,"_socetset.map");

   strcpy(standard_map,outcub_name);
   strcat(standard_map,"_standard.map");

   strcpy(tempcub,outcub_name);
   strcat(tempcub,".temp.cub");

   strcpy(input_cub,"input_");  // input DEM or ORTHO
   strcat(input_cub,outcub);

   strcpy(input_fomcub,"input_");
   strcat(input_fomcub,FOM_outcub);
   
   strcpy(sqrcub,outcub_name);
   strcat(sqrcub,".sqr.cub");

   strcpy(outcub_tiled,outcub_name);
   strcat(outcub_tiled,".tiled.cub");

   strcpy(FOM_sqrcub,FOM_outcub_name);
   strcat(FOM_sqrcub,".sqr.cub");

   strcpy(FOM_outcub_tiled,FOM_outcub_name);
   strcat(FOM_outcub_tiled,".tiled.cub");

   strcpy(layout_cub,outcub_name);
   strcat(layout_cub,"_layout.cub");

   ///////////////////////////////////////////////////////////////////////
   // Generate basic ISIS cubes without mapping labels.  This will be the
   // temporary input cube for the native/standard output cubes
   ///////////////////////////////////////////////////////////////////////
   
   // Output description of this section to script file
   sprintf(command,"######################################################");
   writeToScript(isis_script,command);
   sprintf(command,"## Generate basic ISIS cube(s) without mapping labels");
   writeToScript(isis_script,command);
   sprintf(command,"######################################################\n");
   writeToScript(isis_script,command);
      
   if (strstr(productType,"DEM")) {
      // Convert the raw binary DEM file to a basic ISIS cube
      sprintf(command,"raw2isis from=%s to=%s samples=%d lines=%d bands=1 bittype=real byteorder=%s\n",
              rawDEM,tempcub,samples,lines,byteOrder);
      writeToScript(isis_script,command);

      // Stretch any LRS pixels to NULL in the DEM (this is a workaround
      // until we find out why NULL pixels aren't coming over
      // from WINDOWS)
      sprintf(command,"stretch from=%s to=%s lrs=NULL\n",
              tempcub,input_cub);
      writeToScript(isis_script,command);

      // Convert the raw binary FOM file to a basic ISIS cube
      sprintf(command,"raw2isis from=%s to=%s samples=%d lines=%d bands=1\n",
              rawFOM,input_fomcub,samples,lines);
      writeToScript(isis_script,command);
   }
   else { // productType=ORTHO
      // Convert the raw binary ortho file to a basic ISIS cube
      sprintf(command,"raw2isis from=%s to=%s samples=%d lines=%d bands=1\n",
                 rawORTHO,tempcub,samples,lines);
      writeToScript(isis_script,command);

      // Strech any special pixels (LIS,LRS,HRS,HIS) to 1 or 254 as appropriate
      sprintf(command,"stretch from=%s to=%s+8bit+1:254 pairs=\"0:0 1:1 254:254\" lis=2.0 lrs=2.0 his=253 hrs=253\n",
              tempcub,input_cub);
      writeToScript(isis_script,command);
   }
   
   // Remove tempcub
   sprintf(command,"/bin/rm -f %s\n",tempcub);
   writeToScript(isis_script,command);
   
   /////////////////////////////////////////////////////////////////////////////
   // Based on coordinate system, generate SS_ and standard ISIS cubes
   /////////////////////////////////////////////////////////////////////////////

   switch (coord_sys) {
      case 1: // Geographic Coordinates

         // Geographic DEMs/ORTHOs are stored in EQUI projection with 
         // clat/clon equal to the geographic reference point 
         // established by the user during SOCET Set Project
         // creation.  XY_UNITS of Geographic DEMs/ORTHOs are radians,
         // so set projection to EQUI and get CLAT and CLON from
         // project file and convert to decimal degrees for ISIS

         strcpy(projection,"equirectangular");

         ret = parse_label(prj, "GP_ORIGIN_Y", value);
         clat = atof(value) * rad2deg;
  
         ret = parse_label(prj, "GP_ORIGIN_X", value);
         east180_clon = atof(value) * rad2deg;

         // Convert real spacing from radians to decimal degrees -- for
         // x-spacing, get spacing at the equator for ISIS (SOCET stores
         // the x-spacing at the clat)

         xSpacingDG = x_realspacing * rad2deg * cos(clat/rad2deg);
         ySpacingDG = y_realspacing * rad2deg;

         // Get Lat/Lon boundary of DEM/ORTHO
         // -------------------------------
         // (NOTE: Lons out of SOCET SET are +East, so we are working
         // with +East longitudes here)
         
         // From the (ulcenter_Xlon, ulcenter_Ylat) coordinate, and number of
         // lines/samps of the DEM/FOM/ORTHO, calculate the lat/lon range of
         // the product.  Old-ISIS (version 1) expects the lat/lon range to
         // cover from "edge-to-edge" of the image space, so offset the 
         // lat/lon range by half the  x/y spacing to get an edge-to-edge
         // lat/lon range, and convert to decimal degrees.
         //
         // (ulcenter_Xlon, ulcenter_Ylat) is a Socet Set coordinate, so
         // longitudes are +East +/- 180, latitudes are ographic

         minlat = (ulcenter_Ylat-(y_realspacing*(lines-1))-y_realspacing/2.0 ) * rad2deg; //minlat
         east180_minlon = (ulcenter_Xlon - x_realspacing/2.0 ) * rad2deg; //minlon
         maxlat = (ulcenter_Ylat+y_realspacing/2.0)*rad2deg; //maxlat
         east180_maxlon = (ulcenter_Xlon+(x_realspacing*(samples-1))+x_realspacing/2.0)*rad2deg; //maxlon
 
         // Make sure east180 lons are in +/- 180 system, and minlon is less than maxlon
         if (east180_minlon > 180.0) east180_minlon = east180_minlon - 360.0;
         if (east180_maxlon > 180.0) east180_maxlon = east180_maxlon - 360.0;
         if (east180_clon > 180.0) east180_clon = east180_clon - 360.0;
         if (east180_minlon > east180_maxlon) {
            temp = east180_minlon;
            east180_minlon = east180_maxlon;
            east180_maxlon = temp;
         } 
	 
         // Convert Lon to 0->360 if cross the 180 | -180 border
         if ((east180_minlon > -180) && ( east180_maxlon > 180))
            east180_minlon = 360 + east180_minlon;

         // Now get lons in all possible combinations (we may not need them now, but keep as
         // place holders for future bodies and standard products:
         
         // +East, 0->360 system
         east360_minlon = east180_minlon;
         east360_maxlon = east180_maxlon;
         east360_clon = east180_clon;

         if (east360_minlon < 0) east360_minlon = east360_minlon +360;
         if (east360_maxlon < 0) east360_maxlon = east360_maxlon +360;
         if (east360_clon < 0) east360_clon = east360_clon +360;
	 
	 // Handle Lon range that crosses the 0 | 360 border
         if ((east360_minlon > 180) && ( east360_maxlon < 180))
            east360_minlon = east360_minlon - 360;

         // +West, +/- 180
         // note we switch min/max from +East coordinates to insure minlon < maxlon for ISIS
         west180_minlon = -1 * east180_maxlon;
         west180_maxlon = -1 * east180_minlon;
         west180_clon = -1 * east180_clon;
         // Convert Lon to 0->360 if cross the 180 | -180 border
         if ((west180_minlon < 0) && ( west180_maxlon > 0))
            west180_minlon = 360 + west180_minlon;   

         // +West, 0->360 system
         west360_minlon = west180_minlon;
         west360_maxlon = west180_maxlon;
         west360_clon = west180_clon;
         if (west360_minlon < 0) west360_minlon = west360_minlon +360;
         if (west360_maxlon < 0) west360_maxlon = west360_maxlon +360;
         if (west360_clon < 0) west360_clon = west360_clon +360;
         
//************************************************************************
// Generate the native Socet Set version of the ISIS cube IF the planet
// is treated as an ellipsoid
// NOTE: By default, spherical bodies have square pixels out of
// SOCET Set 
//************************************************************************

         if (ecc != 0) {

            // Output description of this section to script file
            sprintf(command,"#####################################################################");
            writeToScript(isis_script,command);
            sprintf(command,"## Generate SS_ cube(s) (for non-spherical bodies)");
            writeToScript(isis_script,command);
            sprintf(command,"##\n## These cubes are in the native Socet Set format and not resampled,");
	    writeToScript(isis_script,command);
            sprintf(command,"## however the positive longitude direction is set as per IAU standards.");
	    writeToScript(isis_script,command);
            sprintf(command,"##\n## Note that Pixel Scale in not compatible with ISIS!");
            writeToScript(isis_script,command);
	    sprintf(command,"##\n## For compatability with ArcMap, set a londom=180 and clon=0 so that");
	    writeToScript(isis_script,command);
	    sprintf(command,"## ArcMap imports bodies with Positive West Longitudes correctly.");
            writeToScript(isis_script,command);
            sprintf(command,"#####################################################################\n");
            writeToScript(isis_script,command);
         
            // Set minlon, maxlon, clon, ogRef_lat, ogRef_lon to East/West values,
            // as needed for SS native cube
            // NOTE: If longitudes are positive west, maxlat,maxlon refer
            //       to the upper-left-corner of the upper-left pixel,
            //       which is line=0.5, samp=0.5 in ISIS
            //       If longitudes are positive east, we want maxlat, minlon
             
            strcpy(lattype,"planetographic");  // All Socet Set projects are planetographic
            londom=180;                             // All Socet Set projects are in the +/- 180 lon system
            if (strcmp (ographicPosLonDir,"positiveWest") == 0) {
               minlon = west180_minlon;
               maxlon = west180_maxlon;
               clon = west180_clon;
               ogRef_lat = maxlat;
               ogRef_lon = maxlon;
            }
            else {
               minlon = east180_minlon;
               maxlon = east180_maxlon;
               clon = east180_clon;
               ogRef_lat = maxlat;
               ogRef_lon = minlon;
            }
            
            ///////////////////////////////////////////////////////////////
            // Run maptemplate to set up a map file for the Socet Set equi
            // projection.
	    //
	    //  NOTE: for ArcMap to behave properly with positive West longitude
	    //  systems, we either need londom=180 with clon=0, or londom=360 with clon=180.
	    //  Trent Hare recommends we use londom=180, clon=0, so that is
	    //   what we will do.  For Equi Rectangular projection, changing
	    //  the clon does not require resampling of data,  so override all the previous
	    //  clon conversions here. (I've kept all the clon conversions in the code, and commented
	    //  out the maptemplate command with a variable clon, in case we have a reason to revert back.)
            //////////////////////////////////////////////////////////////
	    //sprintf(command,"maptemplate map=%s projection=%s clon=%.8f clat=%.8f targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=%s londom=%d rngopt=user minlat=%.8f maxlat=%.8f minlon=%.8f maxlon=%.8f resopt=ppd resolution=%.8f\n",socetset_map,projection,clon,clat,isisTargName,eqradius,polradius,lattype,ographicPosLonDir,londom,minlat,maxlat,minlon,maxlon,1.0/ySpacingDG);
	    sprintf(command,"maptemplate map=%s projection=%s clon=0.0 clat=%.8f targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=%s londom=%d rngopt=user minlat=%.8f maxlat=%.8f minlon=%.8f maxlon=%.8f resopt=ppd resolution=%.8f\n",socetset_map,projection,clat,isisTargName,eqradius,polradius,lattype,ographicPosLonDir,londom,minlat,maxlat,minlon,maxlon,1.0/ySpacingDG);
	    writeToScript(isis_script,command);

            // For portability, convert tiled input cubes to BSQ format.
            // following the SS_ naming convention for the native formated files
            sprintf(command,"cubeatt from=%s to=%s+BandSequential+Lsb+Attached\n",
                    input_cub,SS_outcub);
            writeToScript(isis_script,command);

            if (strstr(productType,"DEM")) {
               sprintf(command,"cubeatt from=%s to=%s+BandSequential+Lsb+Attached\n",
                      input_fomcub,SS_FOM_outcub);
               writeToScript(isis_script,command);
            }

            // Run maplab to add the Socet Set mapping keywords to the DEM/FOM/ORTHO labels
            sprintf(command,"maplab from=%s map=%s sample=0.5 line=0.5 coordinates=latlon lat=%.8f lon=%.8f\n",
                  SS_outcub,socetset_map,ogRef_lat,ogRef_lon);
            writeToScript(isis_script,command);

            if (strstr(productType,"DEM")) {
               sprintf(command,"maplab from=%s map=%s sample=0.5 line=0.5 coordinates=latlon lat=%.8f lon=%.8f\n",
                                SS_FOM_outcub,socetset_map,ogRef_lat,ogRef_lon);
               writeToScript(isis_script,command);
            }

            // For ellipsoidal bodies, the pixels out of SOCET Set are not
            // scaled the same in the x-direction as in ISIS, so add comment to
            // ISIS labels and issue warning
            sprintf(command,"editlab from=%s options=addg grpname=SS2ISIS_IMPORT_NOTES\n",SS_outcub);
            writeToScript(isis_script,command);
            
            if (strstr(productType,"DEM")) {
               sprintf(command,"editlab from=%s options=addg grpname=SS2ISIS_IMPORT_NOTES\n",SS_FOM_outcub);
               writeToScript(isis_script,command);
            }

            sprintf(note,"PIXEL SCALE NOT ISIS COMPATIBLE, X-dg/px: %.14f, Y-dg/px: %.14f at equator", xSpacingDG, ySpacingDG);
            sprintf(command,"editlab from=%s grpname=SS2ISIS_IMPORT_NOTES keyword=NOTE1 value=\"%s\"\n",SS_outcub,note);
            writeToScript(isis_script,command);
            
            if (strstr(productType,"DEM")) {
               sprintf(command,"editlab from=%s grpname=SS2ISIS_IMPORT_NOTES keyword=NOTE1 value=\"%s\"\n",SS_FOM_outcub,note);
               writeToScript(isis_script,command);
            }
         }
 
//************************************************************************
// Generate standard distribution version of the ISIS cube
// NOTE: By default, spherical bodies have compatable pixels out of
// SOCET Set with ISIS 
//************************************************************************

         sprintf(command,"####################################################################");
         writeToScript(isis_script,command);
         if (ecc != 0.0) {
            sprintf(command,"## Generate standard ISIS cube(s)");
            writeToScript(isis_script,command);
            sprintf(command,"##\n## For ellipsoidal bodies, resample for pixel scale, ocentric");
            writeToScript(isis_script,command);
            sprintf(command,"## latitudes, and when needed, map projection.");
            writeToScript(isis_script,command);
            sprintf(command,"## Also set positive longitude direction as per IAU standards.");
            writeToScript(isis_script,command);
            if (strstr(uc_ellipsoid,"MARS") == NULL) {
               sprintf(command,"## For compatability with ArcMap, set londom=180 and clon=0");
               writeToScript(isis_script,command);
               sprintf(command,"## so that ArcMap imports bodies with Positive West Longitudes correctly.");
               writeToScript(isis_script,command);
            }
           else {
               sprintf(command,"##\n## For MARS, standard products between +/-65 degrees latitude");
               writeToScript(isis_script,command);
               sprintf(command,"## are in equi-rectangular projection; and from 65 degrees to the");
               writeToScript(isis_script,command);
               sprintf(command,"## pole are in polar stereographic projection");
               writeToScript(isis_script,command);
               sprintf(command,"## Also set londom=360, clon=180 and clat=0 in accordance with HiRISE conventions.");
               writeToScript(isis_script,command);
            }
         }
         if (ecc == 0.0) {
            sprintf(command,"## Generate ISIS cube(s)");
            writeToScript(isis_script,command);
            sprintf(command,"##\n## For spherical bodies, no resampling for pixel scale, latitude type,");
            writeToScript(isis_script,command);
            sprintf(command,"## or map projection is required.  Simply set lattype to");
            writeToScript(isis_script,command);
            sprintf(command,"## planetocentric because ocentric==ographic latitudes for spheroids; and");
            writeToScript(isis_script,command);
            sprintf(command,"## set positive longitude direction as per IAU standards.");
            writeToScript(isis_script,command);
            sprintf(command,"##\n## For compatability with ArcMap, set londom=180 and clon=0");
            writeToScript(isis_script,command);
            sprintf(command,"## so that ArcMap imports bodies with Positive West Longitudes correctly.");
            writeToScript(isis_script,command);
            }
         ////////////////}
         sprintf(command,"####################################################################\n");
         writeToScript(isis_script,command);

         // Set minlon, maxlon, clon, ocRef_lat, ocRef_lon to East/West values,
         // as needed for Sstandard ISIS cube
         // NOTE: If longitudes are positive west, maxlat,maxlon refer
         //       to the upper-left-corner of the upper-left pixel,
         //       which is line=0.5, samp=0.5 in ISIS
         //       If longitudes are positive east, we want maxlat, minlon
             
         strcpy(lattype,"planetocentric");  // All default ISIS cubes are planetocentric with londom=360
	                                              //  HOWEVER, for ArcMap to behave properly, we either need
	                                              //  londom=180 with clon=0, or londom=360 with clon=180.
						      //  Trent Hare recommends we use londom=180, clon=0, so that is
	                                              //   what we will do.  For Equi Rectangular projection, changing
	                                              //   the clon or londom does not require resampling of data.
						      //   Note: for Spheroids, ographic==ocentric latitudes, so don't
						      //    have to resample for lattype
         londom=180;
	 clon = 0;
         
         if (strcmp (ocentricPosLonDir,"positiveEast") == 0) {		 
            minlon = east180_minlon;
            maxlon = east180_maxlon;
            ocRef_lat = maxlat;
            ocRef_lon = minlon;
         }
         else {
            minlon = west180_minlon;
            maxlon = west180_maxlon;
            ocRef_lat = maxlat;
            ocRef_lon = maxlon;
         }

        if (strstr(uc_ellipsoid,"MARS") != NULL && (ecc != 0)) { // For standard Mars equirectangular products, londom=360, clon=180, clat=0 (==simp)
            londom=360;
            clon = 180;
            clat = 0.0;
            minlon = east360_minlon;
            maxlon = east360_maxlon;
            ocRef_lat = maxlat;
            ocRef_lon = minlon;
         }
            
         if (ecc != 0.0) {
            // Calculate scale factor to apply in sample direction
            // Use dg/px rather than m/px
            sscale = xSpacingDG / ySpacingDG;

            // run ISIS routine enlarge or reduce (as appropriate), to
            // scale the pixels in the DEM/FOM/ORTHO cubes to fit the ISIS
            // definition of radius of curvature
            // For the FOM cube, use nearest neighbor interpolation to
            // preserve original FOM values as much as possible
            if (sscale > 1.0) {
              sprintf(command,"enlarge from=%s to=%s sscale=%.10f lscale=1.0\n",
                           input_cub, sqrcub, sscale);
              writeToScript(isis_script,command);

              if (strstr(productType,"DEM")) {
                 sprintf(command,"enlarge from=%s to=%s sscale=%.10f lscale=1.0 interp=nearestneighbor\n",
                         input_fomcub, FOM_sqrcub, sscale);
                writeToScript(isis_script,command);
              }
             }
             else {
               //NOTE: scale factor must be greater than 1.0, so pass inverse
               //      of sscale to the reduce program
               sprintf(command,"reduce from=%s to=%s sscale=%.10f lscale=1.0\n",
                            input_cub, sqrcub, 1.0/sscale);
               writeToScript(isis_script,command);

               if (strstr(productType,"DEM")) {
                  sprintf(command,"reduce from=%s to=%s sscale=%.10f lscale=1.0 algorithm=nearest\n",
                         input_fomcub, FOM_sqrcub, 1.0/sscale);
                  writeToScript(isis_script,command);
               }
             }
            
             // Run maplab to add the socetset map projection to the labels
             // and convert from +W <-> +E longitudes if necessary
             // NOTE: If longitudes are positive west, maxlat,maxlon refer
             //       to the upper-left-corner of the upper-left pixel,
             //       which is line=0.5, samp=0.5 in ISIS
             //       If longitudes are positive east, we want maxlat, minlon
           
             sprintf(command,"maplab from=%s map=%s sample=0.5 line=0.5 coordinates=latlon lat=%.8f lon=%.8f\n",
                     sqrcub,socetset_map,ogRef_lat,ogRef_lon);
             writeToScript(isis_script,command);
             
             if (strstr(productType,"DEM")) {
                sprintf(command,"maplab from=%s map=%s sample=0.5 line=0.5 coordinates=latlon lat=%.8f lon=%.8f\n",
                        FOM_sqrcub,socetset_map,ogRef_lat,ogRef_lon);
                writeToScript(isis_script,command);
             }
        
            // Now run maptemplate to define the projection of the standard
            // cube depending on the planet
            // This is also where any change in pos_lon_dir will happen if
            // needed (ie, +E from +W for Mars)
            //
            // For Mars:
            //    cubes with average lat range within +/- 65 deg lat,
            //         project to equi:0.0,180.0,ocentric
            //    cubes with average lat range greater than +/- 65 deg lat,
            //         project to pola:+/-90,0.0
            //
            // Note: you must enter the resolution to the DEG value used above.
            //       If you let lev2tolev2 default to the scale on the labels,
            //       it will use the km/pix value instead of pix/deg, resulting
            //       in a different deg/pix resolution than intended.

            float lat_cutoff;
            if (strstr(uc_ellipsoid,"MARS") != NULL)  // this is Mars
                lat_cutoff = 65.0;
            else
                lat_cutoff = 999.0;
         
            avglat = minlat + (maxlat-minlat)/2.0;

            if (avglat < -lat_cutoff) {
                 // We are at the South Pole
                 sprintf(command,"maptemplate map=%s projection=polarstereographic clon=0.0 clat=-90.0 targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=%s londom=%d resopt=ppd resolution=%.8f\n",standard_map,isisTargName,eqradius,polradius,lattype,ocentricPosLonDir,londom,1.0/ySpacingDG);
               }
               else if (avglat > lat_cutoff) {
                 // We are at the North Pole
                 sprintf(command,"maptemplate map=%s projection=polarstereographic clon=0.0 clat=90.0 targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=%s londom=%d resopt=ppd resolution=%.8f\n",standard_map,isisTargName,eqradius,polradius,lattype,ocentricPosLonDir,londom,1.0/ySpacingDG);
               }
               else {
                 // We are in the equatorial belt
                 sprintf(command,"maptemplate map=%s projection=equirectangular clon=%.1f clat=%.8f targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=%s londom=%d resopt=ppd resolution=%.8f\n",standard_map,clon,clat,isisTargName,eqradius,polradius,lattype,ocentricPosLonDir,londom,1.0/ySpacingDG);
               }

               writeToScript(isis_script,command);

            // Now reproject the DEM/FOM/ORTHO cubes
            // For the FOM cube, use nearest neighbor interpolation to
            // preserve original FOM values as much as possible

            sprintf(command,"map2map from=%s map=%s to=%s pixres=map interp=bilinear\n",
                   sqrcub,standard_map,outcub_tiled);
            writeToScript(isis_script,command);

            if (strstr(productType,"DEM")) {
               sprintf(command,"map2map from=%s map=%s to=%s pixres=map interp=nearestneighbor\n",
                      FOM_sqrcub,standard_map,FOM_outcub_tiled);
               writeToScript(isis_script,command);
            }
         }
         else { // THIS IS A SPHEROID
            // The input cubes are already "square", so replace the
            // sqrcub and FOM_sqrcub file names with the input file names
            strcpy(sqrcub,input_cub);
            strcpy(FOM_sqrcub,input_fomcub);
             
            // For Moon and other spheres:
            //    keep cubes as equi, with original clat, but change clon as specified.
            //    (changing clon does not resample pixels.)
            //
           // Note: you must enter the resolution to the DEG value used above.
            //       If you let lev2tolev2 default to the scale on the labels,
            //       it will use the km/pix value instead of pix/deg, resulting
            //       in a different deg/pix resolution than intended.
             

            strcpy(lattype,"planetocentric");  //For Spheroids, ographic==ocentric latitudes, so call it ocentric in ISIS
		 
            sprintf(command,"maptemplate map=%s projection=%s clon=%.8f clat=%.8f targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=%s londom=%d rngopt=user minlat=%.8f maxlat=%.8f minlon=%.8f maxlon=%.8f resopt=ppd resolution=%.8f\n",standard_map,projection,clon,clat,isisTargName,eqradius,polradius,lattype,ocentricPosLonDir,londom,minlat,maxlat,minlon,maxlon,1.0/ySpacingDG);
            writeToScript(isis_script,command);

            sprintf(command,"maplab from=%s map=%s sample=0.5 line=0.5 coordinates=latlon lat=%.8f lon=%.8f\n",
                    sqrcub,standard_map,ocRef_lat,ocRef_lon);
            writeToScript(isis_script,command);
                
            if (strstr(productType,"DEM")) {
                sprintf(command,"maplab from=%s map=%s sample=0.5 line=0.5 coordinates=latlon lat=%.8f lon=%.8f\n",
                        FOM_sqrcub,standard_map,ocRef_lat,ocRef_lon);
                writeToScript(isis_script,command);
            }
            // The input cubes are the "outcub_tiled" cubes, so replace the
            // outcub_tiled and FOM_outcub_tiled file names with the input file names
            strcpy(outcub_tiled,input_cub);
            strcpy(FOM_outcub_tiled,input_fomcub);

        }

         // For portability, convert tiled cube to BSQ format.
         // NOTE: it is faster run map2map on tiled cubes, and then
         //       run cubeatt to do the reformat (especially for
         //       polar projected Hirise cubes!)
         sprintf(command,"cubeatt from=%s to=%s+BandSequential+Lsb+Attached\n",
                 outcub_tiled,outcub);
         writeToScript(isis_script,command);

         if (strstr(productType,"DEM")) {
            sprintf(command,"cubeatt from=%s to=%s+BandSequential+Lsb+Attached\n",
                    FOM_outcub_tiled,FOM_outcub);
            writeToScript(isis_script,command);
         }

         // Delete remaining temporary files
         if (strstr(productType,"DEM"))
            sprintf(command,"/bin/rm -f %s %s %s %s %s %s %s %s\n",input_cub,input_fomcub,sqrcub,FOM_sqrcub,outcub_tiled,FOM_outcub_tiled,socetset_map,standard_map);
         else
            sprintf(command,"/bin/rm -f %s %s %s %s %s\n",sqrcub,input_cub,outcub_tiled,socetset_map,standard_map);
         writeToScript(isis_script,command);

         /////////////////////////////////////////////////////////
         // If a lower-resolution layout file is need, generate it
         /////////////////////////////////////////////////////////
         
         if (layout_flag[0]=='y' || layout_flag[0]=='Y') {
           // Output description of this section to script file
           sprintf(command,"#############################################################");
           writeToScript(isis_script,command);
           sprintf(command,"## Generate reduced resolution ISIS cubes for ARCMAP layouts");
           writeToScript(isis_script,command);
           sprintf(command,"#############################################################\n");
           writeToScript(isis_script,command);
         
           sprintf(command,"reduce from=%s to=%s sscale=5.0 lscale=5.0 validper=10 vper_replace=nearest\n",outcub,layout_cub);
           writeToScript(isis_script,command);
         }

      break;

   case 6: // Grid (Map-Projected) Coordinates
       
      // NOTE: In ISIS3, we only need to the (x,y) map projected
      // coordinate (in meters) of a line, sample location in the DEM/FOM/ORTHO for
      // geo-referencing.  We will not attempt to walk the DEM/FOM/ORTHO for the
      // lat/lon range of valid pixels since a lat/lon range is not
      // required for further processing in ISIS3 as it was for ISIS2
      // (e.g. for map2map)
      //
      // For our purposes, use the ulcenter_Xlon, ulcenter_Ylat coordinate as
      // our georeferencing point.
      //
      // Also, stay in the +East, 180 londom system.  I.E., stay in the coordinate
      // system used by Socet Set, and don't worry about ISIS conventions for now.

      // Based on the projection, output description of this section to script file
      // or output a report that the projection is not supported in this program
   
      // Get Projection
      ret = parse_label(prj, "PROJECTION_TYPE", projection);
      
      if(strcmp(projection,"POLAR_STEREOGRAPHIC_PROJECTION")==0) {
         sprintf(command,"#####################################################################");
         writeToScript(isis_script,command);   
         sprintf(command,"## Generate Polar Stereographic Cube");
         writeToScript(isis_script,command);
         sprintf(command,"## These cubes are in the native Socet Set format and not resampled.");
         writeToScript(isis_script,command);
         sprintf(command,"#####################################################################\n");
         writeToScript(isis_script,command);
      }
      else if(strcmp(projection,"SINUSOIDAL_PROJECTION")==0) {
         sprintf(command,"#####################################################################");
         writeToScript(isis_script,command);   
         sprintf(command,"## Generate Sinusoidal Cube");
         writeToScript(isis_script,command);
         sprintf(command,"## These cubes are in the native Socet Set format and not resampled.");
         writeToScript(isis_script,command);
         sprintf(command,"#####################################################################\n");
         writeToScript(isis_script,command);
      }
      else {
         printf ("/n/n%s PROJECTION NOT SUPPORTED!!\n",projection);
         exit(1);
      }

      // Now get the specifics on the map projection by parsing the project file...

      if(strcmp(projection,"POLAR_STEREOGRAPHIC_PROJECTION")==0) {
    
         // Get N or S for Central Long
         ret = parse_label(prj, "POLAR_ASPECT", polarAspect);

         // convert polarAspect to Center Lat
         if (strcmp("N",polarAspect) == 0)
           clat = 90.0;
         if (strcmp("S",polarAspect) == 0)
           clat = -90.0; 
    
         // get Center Long
         ret = parse_label(prj, "CENTER_LONGITUDE", value);
         east180_clon = atof(value);
    
         // Get projection scale bounds
         ret = parse_label(prj, "CENTRAL_SCALE_FACTOR", value);
         prj_scale = atof(value);
      }

      if(strcmp(projection,"SINUSOIDAL_PROJECTION")==0) {
    
         // get Center Long and convert from radians to degrees
         ret = parse_label(prj, "CENTRAL_MERIDIAN", value);
         east180_clon = atof(value) * 57.295779513082320876798154814105;
    
      }
      
       // Make sure clon is in +/- 180 system
       if (east180_clon > 180.0) east180_clon = east180_clon - 360.0;

       londom = 180;  // SOCET Set is in +/- 180 system, so stay in it
      
       strcpy(lattype,"planetocentric");  // For spheres, ographic=ocentric
                                          // For ellipsoids, so long as
                                          // clat = +/- 90, we don't have the
                                          // scaling radius issue, and the
                                          // projection equations are ocentric.
       
//************************************************************************
// Generate the ISIS cube in polar stereographic coordinates
//************************************************************************

       // Run maptemplate to set up the map file

       if(strcmp(projection,"POLAR_STEREOGRAPHIC_PROJECTION")==0)
          sprintf(command,"maptemplate map=%s projection=polarstereographic clon=%.8f clat=%.8f targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=positiveEast londom=180 resopt=mpp resolution=%.8f\n",standard_map,east180_clon,clat,isisTargName,eqradius,polradius,lattype,x_realspacing);
       else // Generate sinusoidal map file
          sprintf(command,"maptemplate map=%s projection=sinusoidal clon=%.8f targopt=user targetname=%s eqradius=%.3f polradius=%.3f lattype=%s londir=positiveEast londom=180 resopt=mpp resolution=%.8f\n",standard_map,east180_clon,isisTargName,eqradius,polradius,lattype,x_realspacing);   

       writeToScript(isis_script,command);

       // Run maplab to add the mapping keywords to the DEM/FOM/ORTHO cube labels
       sprintf(command,"maplab from=%s map=%s sample=1.0 line=1.0 x=%.8f y=%.8f\n",input_cub,standard_map,ulcenter_Xlon,ulcenter_Ylat);
       writeToScript(isis_script,command);

       if (strstr(productType,"DEM")) {
          sprintf(command,"maplab from=%s map=%s sample=1.0 line=1.0 x=%.8f y=%.8f\n",input_fomcub,standard_map,ulcenter_Xlon,ulcenter_Ylat);
          writeToScript(isis_script,command);
       }

       // For portability, convert tiled cubes to BSQ format.
       // NOTE: it is faster run ISIS programs on tiled cubes, and then
       //       run cubeatt to do the reformat
       sprintf(command,"cubeatt from=%s to=%s+BandSequential+Lsb+Attached\n",
                 input_cub,outcub);
       writeToScript(isis_script,command);
       
       if (strstr(productType,"DEM")) {
          sprintf(command,"cubeatt from=%s to=%s+BandSequential+Lsb+Attached\n",
                  input_fomcub,FOM_outcub);
          writeToScript(isis_script,command);
       }

       // Delete remaining temporary files
       if (strstr(productType,"DEM"))
          sprintf(command,"/bin/rm -f %s %s %s %s\n",input_cub,input_fomcub,socetset_map,standard_map);
       else
          sprintf(command,"/bin/rm -f %s %s %s\n",input_cub,socetset_map,standard_map);
       writeToScript(isis_script,command);

       /////////////////////////////////////////////////////////
       // If a lower-resolution layout file is need, generate it
       /////////////////////////////////////////////////////////
       
       if (layout_flag[0]=='y' || layout_flag[0]=='Y') {
         // Output description of this section to script file
         sprintf(command,"#############################################################");
         writeToScript(isis_script,command);
         sprintf(command,"## Generate reduced resolution ISIS cubes for ARCMAP layouts");
         writeToScript(isis_script,command);
         sprintf(command,"#############################################################\n");
         writeToScript(isis_script,command);
       
         sprintf(command,"reduce from=%s to=%s sscale=5.0 lscale=5.0 validper=10 vper_replace=nearest\n",outcub,layout_cub);
         writeToScript(isis_script,command);
       }
    
////SOME UTM NOTES KEPT FOR REFERENECE:
////UTM then clon=(ZONE-31)*6+3, SCFA=.9996 for Earth
//if (strcmp(projection,"UTM_PROJECTION") == 0) {
//// Get Zone
//ret = parse_label(prj, "ZONE", value);
//zone = atoi(value);
//clon = (zone - 31) * 6 + 3;
//prj_scale = 0.9996;
//       }
         
   } // end switch (coord_sys)
   
   return(0);

} // END of generate_ss2isis_script



/**************  parse_label.c  ********************
*                                                  *
*  This routine reads a Socet Set project file and *
*                                                  *
*  by Trent Hare       for USGS, flagstaff         *
*                                                  *
****************************************************/
int parse_label(char *file, char *keyword, char *value)
{

//declaration 
   int scan_value;
   char id[FILELEN];

   FILE *fp;


   fp = fopen(file,"r");
   if (fp == NULL)
     {
      printf("\ncan't open the input file %s!\n",file);
      exit(1);
   }   

//Search for sent parameter

  scan_value = fscanf (fp,"%s",id);
  while(scan_value != EOF) {
    if (strcmp(id,keyword) == 0) {
      scan_value = fscanf (fp,"%s",id);
      strcpy(value,id);
      //printf("%s\n",value);
      return(1);
    }
    scan_value = fscanf (fp,"%s",id);
  }
  fclose(fp);
  strcpy(value,"-9999");
  return(0);
}

//**************  getTargetInfo  ******************************************

// ADD COMMENTS
// by Elpitha H-Kraus for USGS, flagstaff
//
// HIST   12 Mar 2008   EHK - Original Version
//                                        (Combination of get* in dem2isis2)
//            19 Mar 2014,  EHK - Changed postive longitude direction from West
//                                         to East for all bodies except TITAN.
//                                          NOTE: this means IAU postitive longitiude
//                                          standards are not longer being met, but the
//                                          ARCMAP defaults
//            10 Apr 2014,  EHK - Removed isisTargDef
//*************************************************************************
int getTargetInfo (char *uc_ellipsoid, char *isisTargName,
                   char *ographicPosLonDir, char *ocentricPosLonDir)
{

  // NOTE: This subroutine supports some "extra" figures from what is in the
  // in the current <ss_install>\internal_dbs\GEODETIC files.  These are
  // legacy figures from past Socet Set GEODETIC  configurations, and kept
  // in case an old project is restored, and data products exported.

  if(strstr(uc_ellipsoid,"ADRASTEA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Adrastea2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Adrastea");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"AMALTHEA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Amalthea2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Amalthea");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
  }
  else if(strstr(uc_ellipsoid,"ANANKE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Ananke2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Ananke");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"ARIEL") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Ariel2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Ariel");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"ATLAS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Atlas2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Atlas");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
  }
  else if(strstr(uc_ellipsoid,"BELINDA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Belinda2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Belinda");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"BIANCA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Bianca2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Bianca");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"CALLISTO") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Callisto2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Callisto");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
  }
  else if(strstr(uc_ellipsoid,"CALYPSO") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Calypso2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Calypso");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"CARME") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Carme2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Carme");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"CHARON") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Charon2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Charon");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"CORDELIA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Cordelia2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Cordelia");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"CRESSIDA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Cressida2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Cressida");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"DEIMOS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Deimos2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Deimos");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"DESDEMONA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Desdemona2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Desdemona");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"DESPINA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Despina2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Despina");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"DIONE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Dione2006
	  // Past Socet Set Ellipsoid name(s): DIONE,
    strcpy(isisTargName, "Dione");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"ELARA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Elara2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Elara");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"ENCELADUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Enceladus2006
	  // Past Socet Set Ellipsoid name(s): ENCELADUS,  ENCELADUS2009
    strcpy(isisTargName, "Enceladus");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"EPIMETHEUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Epimetheus2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Epimetheus");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"EROS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Eros
	  // Past Socet Set Ellipsoid name(s): EROS
    strcpy(isisTargName, "Eros");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"EUROPA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Europa2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Europa");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"GALATEA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Galatea2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Galatea");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"GANYMEDE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Ganymede2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Ganymede");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"HELENE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Helene2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Helene");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"HIMALIA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Himalia2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Himalia");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"HYPERION") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Hyperion2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Hyperion");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"IAPETUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Iapetus2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Iapetus");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"IO") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Io2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Io");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"JANUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Janus2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Janus");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"JULIET") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Juliet2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Juliet");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"JUPITER") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Jupiter2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Jupiter");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"LARISSA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Larissa2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Larissa");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"LEDA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Leda2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Leda");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"LYSITHEA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Lysithea2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Lysithea");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"MARS") != NULL && strstr(uc_ellipsoid,"MARS1991") == NULL) {
	  // Current SOCET Set Ellipsoid name(s):  Mars2000, Mars2000_MOLA_SPHERE
	  // Past Socet Set Ellipsoid name(s):  MARS1991 
    strcpy(isisTargName, "Mars");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"MARS1991") != NULL) {
    strcpy(isisTargName, "Mars");
          // MARS1991 is a legacy ellipsoid, meant to out Mars products in ographic lat, pos west lon
	  // so maintain positive west lon direction
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
  }
  else if(strstr(uc_ellipsoid,"MERCURY") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Mercury_MSGR, Mercury2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Mercury");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"METIS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Metis2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Metis");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"MIMAS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Mimas2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Mimas");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"MIRANDA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Miranda2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Miranda");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"MOON") != NULL && strstr(uc_ellipsoid,"MOONW") == NULL) {
	  // Current SOCET Set Ellipsoid name(s): Moon2000
	  // Past Socet Set Ellipsoid name(s): n/a MOON, MOONW
    strcpy(isisTargName, "Moon");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
    else if(strstr(uc_ellipsoid,"MOONW") != NULL) {
	  // MOONW is a legacy ellipsoid, meant to treat the MOON as positive West longitude,
	  // so maintain positive west lon direction
    strcpy(isisTargName, "Moon");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
  }
  else if(strstr(uc_ellipsoid,"NAIAD") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Naiad2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Naiad");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"NEPTUNE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Neptune2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Neptune");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"NEREID") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Nereid2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Nereid");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"OBERON") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Oberon2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Oberon");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
    else if(strstr(uc_ellipsoid,"OPHELIA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Ophelia2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Ophelia");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PAN") != NULL && strstr(uc_ellipsoid,"PANDORA") == NULL) {
	  // Current SOCET Set Ellipsoid name(s): Pan2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Pan");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PANDORA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Pandora2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Pandora");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PASIPHAE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Pasiphae2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Pasiphae");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PHOBOS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Phobos2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Phobos");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PHOEBE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Phoebe2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Phoebe");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PLUTO") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Pluto2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Pluto");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PORTIA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Portia2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Portia");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PROMETHEUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Prometheus2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Prometheus");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PROTEUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Proteus2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Proteus");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"PUCK") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Puck2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Puck");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"RHEA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Rhea2006
	  // Past Socet Set Ellipsoid name(s): RHEA, RHEA2009
    strcpy(isisTargName, "Rhea");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"ROSALIND") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Rosalind2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Rosalind");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"SATURN") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Saturn2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Saturn");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"SINOPE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Sinope2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Sinope");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"TELESTO") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Telesto2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Telesto");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"TETHYS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Tethys2006
	  // Past Socet Set Ellipsoid name(s): TETHYS
    strcpy(isisTargName, "Tethys");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"THALASSA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Thalassa2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Thalassa");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"THEBE") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Thebe2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Thebe");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"TITAN") != NULL && strstr(uc_ellipsoid,"TITANIA") == NULL) {
	  // Current SOCET Set Ellipsoid name(s): Titan2000
	  // Past Socet Set Ellipsoid name(s): TITAN
    strcpy(isisTargName, "Titan");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
  }
  else if(strstr(uc_ellipsoid,"TITAN2000") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): 
	  // Past Socet Set Ellipsoid name(s): n/a  // TO DO
    strcpy(isisTargName, "Titan");
    strcpy(ographicPosLonDir, "positiveWest");
    strcpy(ocentricPosLonDir, "positiveWest");
    //Set *PosLonDir to positiveEast for compatibility with ISIS3 and ARCMAP ????
    //strcpy(ographicPosLonDir, "positiveEast");
    //strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"TITANIA") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Titania2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Titania");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"TRITON") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Triton2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Triton");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"UMBRIEL") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Umbriel2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Umbriel");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"URANUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Uranus2006
	  // Past Socet Set Ellipsoid name(s): n/a
    strcpy(isisTargName, "Uranus");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"VENUS") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Venus1985, Venus2000
	  // Past Socet Set Ellipsoid name(s): VENUS, VENUSMED  TO DO
    strcpy(isisTargName, "Venus");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"WILD2") != NULL) {
	  // Current SOCET Set Ellipsoid name(s): Wild2
	  // Past Socet Set Ellipsoid name(s): WILD2
    strcpy(isisTargName, "Wild2");
    //PosLonDir undocumented in ISIS, so set to positiveEast for now for compatibility with ISIS3 and ARCMAP
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else if(strstr(uc_ellipsoid,"WGS_84") != NULL) {
    strcpy(isisTargName, "Earth");
    strcpy(ographicPosLonDir, "positiveEast");
    strcpy(ocentricPosLonDir, "positiveEast");
  }
  else
    return(1);

  return(0);

} // End of getTargetInfo

int writeToScript(char *isis_script, char *command)

{
  //output isis command to print.prt

  FILE *fp;

  // Open script file
  fp = fopen(isis_script,"a");
  if (fp == NULL) {
     printf("\ncan't open isis script file %s!\n",isis_script);
     return(1);
  }

  fprintf (fp,"%s\n",command);

  fclose (fp);
  return (0);

} // End of writeToScript

