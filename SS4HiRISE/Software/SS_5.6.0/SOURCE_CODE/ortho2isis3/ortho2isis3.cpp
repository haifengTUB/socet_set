////////////////////////////////////////////////////////////////////////////////
//
//_Title ORTHO2ISIS3 outputs a SOCET orthoimage as a raw 8-bit image plus ISIS3 script 
//                                                  
//_Desc  This is a SOCET Set program that uses SOCET DEV_KIT routines to access
//       a SOCET orthoimage and output a raw 8-bit file and an associated ISIS3
//       processing script.  The raw file and script are to be transferred to
//       an ISIS machine for port to ISIS3.
//
//       Orthoimages in Geographic Coordinates or Polar Stereographic (Grid)
//       coordinates are currently supported.
//
//       Input parameters are:
//
//              SS_project
//              socet_orthoimage.sup
//              isis_ortho.cub
//
//       Output files are:
//
//              ./isis_ortho.raw
//              ./isis_ortho2isis3.sh
//
//       This isis_ortho2isis3.sh script will generate up to three output
//       files:
//              SS_isis.cub (For Geographic projects only)
//              isis.cub
//              isis_layout.cub (optional)
//
//
//       For Orthoimages in the Geographic Coordinate System:
//       -------------------------------------------------------------
//       Orthoimages in Geographic Coordinates in SS are in essentially
//       equirectangular map projection.  However, SS calculates scaling
//       radii such that the resulting pixels are square w.r.t. meter
//       measurements (and non-square w.r.t. degrees.)  In ISIS, a
//       different scaling radius is calculated for the equi-rectangular
//       map projection, and the resulting pixels are square in degree-space
//       but not meter-space.  For bodies treated as a spheroid, SS and ISIS
//       scaling radii are the same, and products from either package are
//       compatable.  However, for ellipsoidal planets such as Mars, the
//       different approaches to scaling radii become evident. Because
//       of these differences, we output two versions of the orthoimage:
//       SS_isis.cub and isis.cub.
//
//       SS_isis.cub maintains the original SS orthoimage geometry with respect
//       to scaled radii values, and is retained in case someone needs to use
//       the non-resampled image data.  Map projection information is added to 
//       the ISIS labels - with the resolution based on the Y spacing of the
//       SOCET orthoimage in dg/px. For ellipsoidal bodies, a
//       SS2ISIS_IMPORT_NOTES keyword is added to the SS_isis.cub labels
//       stating the scaling is not compatible with ISIS.  Furthermore, for
//       Mars, these orthoimages are in the ographic latitude, +West longitude
//       system.
//
//       isis.cub is our standard ISIS distribution cube.  For
//       ellipsoidal bodies, pixel are scaled to be compatible with ISIS
//       radii scaling (so as to have 'square pixels' in ISIS) and map
//       projection labels.  Furthermore, for Mars, orthoimages above
//       +/- 65 degrees latitude are in polar projection.  Between +/-65 degrees
//       latitude, they are in Equirectangular projection with clat=0, clon=180.
//       All standard Mars orthoimages are in the +E 360 degree longitude
//       system, with ocentric lats.
//
//       Finally, for orthoimages in geographic coordinates, the lat/lon range
//       is well known in SOCET and this information is added to the ISIS
//       labels.
//
//       To summarize:
//       Two ISIS cubes will be output when mapping in SS Geographic coordinates,
//       with the following conventions:
//            1) The first will preserve the native SS pixels, and
//                  - have SS_ prefix added to the output file name
//                  - be in equirectangular map projection with clat, clon
//                    as set by user during SS project creation
//                  - have ographic lats
//                  - follow IAU conventions for ographic lats
//                    (ie, MARS products will have +W lons, ographic lats)
//                  - For ellipsoidal bodies, a SS2ISIS_IMPORT_NOTES keyword
//                    warning will be added to the ISIS cube labels that reads
//                    "PIXEL SCALE NOT ISIS COMPATIBLE"
//            2) The second will follow the format of a 'standard' ISIS cube
//               for distribution, and be resampled to
//                  - follow IAU conventions for ocentric lats and positive
//                    longitude direction.
//                  - have ISIS compatible pixel scales
//               Furthermore:
//                  orthoimages between +/- 65 degrees lat will be in:
//                     - equirectangular map projection with clat=0, clon=180
//                  orthoimages from +/- 65 degrees lat to the poles will be in
//                     - polar stereographic map projection.
//                  Note: the center/average lat of the cube is used to
//                  determine if cube falls in 'polar region'.
//
//
//       For orthoimages in Polar Stereographic (Grid) Coordinates:
//       -------------------------------------------------------------
//       SOCET Set orthoimages produced in a Polar Steregraphic Grid coordinate
//       system will not have the same scaling radii issues that we have with
//       geographic coordinates so long as the clat and clon of the projection
//       is at the poles (i.e., +/-90 degrees.)  In this case, we need only
//       produce the standard isis.cub, using the polar stereographic
//       projected coordinate of the upperleft orthoimages post to geo-reference
//       the orthoimage in ISIS.
//
//       Also note that the lat/lon range of polar stereographic orthoimages
//       will *not* be placed on the ISIS cube labels.  The lat/lon range of
//       orthoimages in grid coordinates is not readily known out of SOCET SET
//       (especially for map-projections without straight meridians and/or
//       parallels, such as polar stereographic).  However, unlike ISIS2, ISIS3
//       does not require a lat/lon range on the labels for further processing,
//       so no attempts are made to add a lat/lon range.
//
//
//_Hist Aug 19 2008 EHK  Orig Version - port ortho2isis2 to ISIS3
//      Sep 05 2008 EHK  Output isis commands to a script file
//      Oct 06 2008 EHK  Some HiRISE images processed through hi4socet.pl
//                       prior to June 4, 2008 have HRS pixels.  When running
//                       map2map in isis3, these HRS pixels become NULLs.
//                       To avoid this, added 'stretch' to the processing steps
//                       to map LRS,LIS and HRS,HIS to 1 and 254 respectively.
//      Oct 7 2008  EHK  Split conversion of SS image into a raw image into
//                       two sections to avoid memory allocation problems
//      Oct 14 2008 EHK  map2map is mapping pixels with DNs of 1 to Nulls in
//                       some areas, so changed stretch to mapp LRS and LIS
//                       pixels to 2.0 until the bug us worked out in isis.
//                       Also, to be on the safe side, mapped hrs=his=253 for
//                       now.
//      Nov 03 2008 EHK  Added option to generate layout_cube for ARCMAP
//                       and made modifications to use SS libraries for
//                       building filenames, etc., to match ortho2isis3 interface
//
//                       2) Added targets subdirectory to /farm/geo/.../isisdata
//                          and added new target def files
//                          mars_mola_DLR_sphere.def.1 and
//                          mars_2000_east_ographic.def.1 to getTargDef routine,
//                          and 'replaced' MARS2000 with MARS2000W.
//      Nov 10 2008 EHK  Split conversion of SS into a raw image into 4 sections
//                       (previously I split it into 2 sections, but that was
//                       not small enough for Windows.)
//                       Also, changed interpolation method for map2map to
//                       bilinear, rather than the default of cubic convolution.
//                       (Cubic convolution will map valid pixels to nulls in
//                       areas of significant DN change...this is a feature of
//                       the polynomial fit for cubic convolution.)
//      Jan 30 2009 EHK  Broke conversion of SS into a raw image into 10
//                       sections to insure there are no memory problems...
//                       ...especially when other processes are running.
//      Feb  4 2009 EHK  Missed updating a hardcoded value used when we
//                       broke the image into four sections.  To avoid this
//                       in any future 'break-ups', changed hard-coded
//                       values to num_sec and sec variables.
//      Jun 15 2009 EHK  Added export of polar stereographic orthoimages from
//                       SS for LRO and updated documentation
//      Aug 20 2009 EHK  Added WGS_84 (Earth) as a supported datum
//      Jan 20 2010 EHK  Updated ISIS commands for version ISIS3.1.21
//                       and added setisis command to ISIS script compatible
//                       with ISIS commands.
//      Feb 21 2010 EHK  Handle MOON2000 and TITAN2000 as a Socet Set target names
//                       to conform with SSv55 geodetic files
//      May 26 2010 EHK  Corrected bug that put Moon2000 in positive West system
//                       when exporting to ISIS.
//      Sep 20 2010 EHK  Reorganized code to call subroutines in common with
//                       ortho2isis3, and avoid rerunning raw2isis,etc, multiple times.
//                       Also changed output scripts not to delete raw binary files.
//      Oct 17 2011 EHK  Added SS_ prefix to output *.raw files to avoid confusion with
//                        the 'standard' ISIS cubes having the same core name, but different
//                         number of lines and samples
//
//_End
//
////////////////////////////////////////////////////////////////////////////////

//System includes
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

//SOCET SET
#include <project/get_proj.h>
#include <sens/SensorModel.h>
#include <sens/smplugins/basic_plugin/OrthoSensorModel.h>
#include <img/img.h>
#include <img/img_main.h>
#include <util/string_dpw.h>
#include <util/init_socet_app.h>
#include <ground_point.h>

#define FILELEN 512

// prototypes
extern int parse_label(char *file, char *keyword, char *value);
extern int getTargetInfo (char *ellipsoid, char *isisTargName, char *isisTargDef,
                           char *ographicPosLonDir, char *ocentricPosLonDir);
extern int writeToScript(char *isis_script, char *command);
extern int generate_ss2isis_script (char *isis_script, char *prj,
            char *productType, char *byteOrder, char *outcub_name,
            char *layout_flag, int lines, int samples, double x_realspacing,
            double y_realspacing, double ulcenter_Xlon, double ulcenter_Ylat);

void main(int argc,char *argv[])
{

/*        declaration            */

   // Input variables
   char file[FILELEN];           //holds ortho-image file name
   char ortho[FILELEN];      //temp hold for ortho support file name
   char orthoName[FILELEN];  //temp hold for ortho support name

   int in_img;
   char square_pixel_flag;

   // SOCET declarations
   int iband,bands,lines,samples, tile_x, tile_y;
   int tenth_lines, rest_lines, num_sec, sec;
   SensorModel *is;               //Sensormodel struct for ortho image 
   OrthoSensorModel *ortho_sup;    //OrthoSensormodel struct for ortho image
   double Ylatref; // Y (meters) or lat (radians) at center of ll pixel 
   double Xlonref; // X (meters) or lon (radians) at center of ll pixel
   double x_realspacing;
   double y_realspacing;
   double ulcenter_Xlon, ulcenter_Ylat;
   FILE *out_img;

   // Project File Variables
   img_proj_struct project;        //SS project structure
   img_proj_struct_ptr  proj_ptr;  //Pointer to SS project file
   char prj[FILELEN];              //SS project with full path and extension
   char  projectName[FILELEN];      // SS project w/o path or .prj ext
   char  projectFile[FILELEN];      // SS project w/.prj ext
   int coord_sys, xy_units, z_units;
   double radius, ecc;
   char ellipsoid[FILELEN];
   char projection[FILELEN];
   char polarAspect[FILELEN];
   double prj_scale;

   // ISIS variables
   char rawORTHO[FILELEN];
   char isis_script[FILELEN];
   char socetset_map[FILELEN];
   char standard_map[FILELEN];
   char SS_tempcub[FILELEN];
   char strcub[FILELEN];
   char tempcub[FILELEN];
   char sqrcub[FILELEN];
   char outcub_tiled[FILELEN];
   char outcub[FILELEN];
   char outcub_name[FILELEN];
   char SS_outcub[FILELEN];
   char note[256];
   char layout_flag[1];
   char layout_cub[FILELEN];

   // temporary files and names
   int   prjReadErr;            // Error flag for reading project file

   // Misc declarations
   int ii=0,scan_value;
   unsigned long int i;
   int ret;
   char value[FILELEN];
   char command[512];
   double rad2deg = 180.0/M_PI;  //convert deg to radians and back
   char productType[4];   //productType=DEM,FOM or ORT
   char byteOrder[4];     //needed for DEMs, but here as place holder

  /////////////////////////////////////////////////////////////////////////////
  // Check number of command line args and issue help if needed
  // Otherwise initiate the socet set application
  /////////////////////////////////////////////////////////////////////////////

   if (argc<4) {
    cerr << "\nRun ortho2isis3 as follows:\n";
    cerr << "start_socet -single ortho2isis3 <project> <ortho> <isis>.cub <layout_flag>\n";
    cerr << "\nwhere:\n";
    cerr << "project = SOCET SET project name to export Orthoimage from\n";
    cerr << "          (path and extension is not required)\n";
    cerr << "ortho = SOCET SET support file of orthoimage image to export\n";
    cerr << "        (path and extension is not required)\n";
    cerr << "isis.cub = desired name of standard isis cube\n";
    cerr << "          (path not allowed since the *.raw and *.sh files)\n";
    cerr << "          (are to be copied to an ISIS machine)\n";
    cerr << "layout_flag = flag to generate lower resolution standard cube for\n";
    cerr << "              use in ARCMAP layouts.  Enter y or n, default=n\n";
    exit(1);
   }

   // Top level SOCET SET initialization  routine.
   // Should be  called  before  any  PCI services.
   init_socet_app ( argv[0], argc, argv);

  /////////////////////////////////////////////////////////////////////////////
  //Get input arguments
  /////////////////////////////////////////////////////////////////////////////

   strcpy(prj,argv[1]);
   strcpy(ortho,argv[2]);
   strcpy(outcub,argv[3]);
   if (argc==5)
      strcpy(layout_flag,argv[4]);
   else
      strcpy(layout_flag,"n");

  /////////////////////////////////////////////////////////////////////////////
  // Populate the project structure - with error checking
  // and get variations of project file representation:
  //     a) name only
  //     b) name + extension
  //     c) fullpath + name + extension
  /////////////////////////////////////////////////////////////////////////////

  // Make sure project name contains no path, but has .prj extension
  strcpy(projectName,ReturnFileName(prj));
  StripFileExt(projectName);
  strcpy(projectFile,concat(projectName, ".prj"));

  // Load project structure
  prjReadErr = project.read(projectFile);
  switch (prjReadErr) {
    case 0:
      setCurrentProj(project);
      break;
    case -1:
      cerr << "Failed to open project file " << projectFile << endl;
      exit (1);
      break;
    case -2:
      cerr << "Project file read error: unknow line " << endl;
      exit (1);
      break;
    default:
      cerr << "Project file read error on line #" << prjReadErr << endl;
      exit (1);
      break;
  }

  // Get project with full path and extension
  strcpy(prj, concat(project.project_data_path,".prj"));

  /////////////////////////////////////////////////////////////////////////////
  // Get variations of ORTHO file representation:
  //      a) fullpath + name + ".sup"
  //      b) fullpath + name
  // and make sure it exists
  /////////////////////////////////////////////////////////////////////////////

  // Make sure project name contains no path, but has .prj extension
  strcpy(orthoName,ReturnFileName(ortho));
  StripFileExt(orthoName);
  build_file_name(ortho, project.project_data_path, orthoName, ".sup");

  if (!file_exists(ortho)) {
    cerr << "\northo " << orthoName << " does not exist!\n";
    cerr << "(NOTE: looking for " << ortho << ")\n";
    exit(1);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Generate raw binary output file name
  /////////////////////////////////////////////////////////////////////////////

   // strip extension from outcub to get the base filename
   strcpy(outcub_name,ReturnFileName(outcub));
   StripFileExt(outcub_name);

   strcpy(rawORTHO,"SS_");
   strcat(rawORTHO,outcub_name);
   strcat(rawORTHO,".raw");
   if (file_exists(rawORTHO))
      file_remove(rawORTHO);
   
  /////////////////////////////////////////////////////////////////////////////
  // Generate output isis script name
  /////////////////////////////////////////////////////////////////////////////

   strcpy(isis_script,outcub_name);
   strcat(isis_script,"_ortho2isis3.sh");
   if (file_exists(isis_script))
     file_remove(isis_script);

   /////////////////////////////////////////////////////////////////////////////
   //output ortho2isis3 command that was issued to isis_script file
   /////////////////////////////////////////////////////////////////////////////
   sprintf(command,"#!/bin/csh\n");
   writeToScript(isis_script,command);

   if (argc==4) sprintf(command,"## start_socet -single ortho2isis3 %s %s %s\n",argv[1],argv[2],argv[3]);
   if (argc==5) sprintf(command,"## start_socet -single ortho2isis3  %s %s %s %s\n",argv[1],argv[2],argv[3],argv[4]);
   writeToScript(isis_script,command);

  /////////////////////////////////////////////////////////////////////////////
  // Set the project pointer and read support file
  /////////////////////////////////////////////////////////////////////////////
   
   // Set pointer to project structure.
   proj_ptr = getCurrentProjStruct();

   // Read support file
   is = read_support_file(ortho);
   if ( is == NULL ) {
      printf("Error opening or reading input support %s\n",ortho);
      exit (-1);
   }

   // Open ortho image file and check status
   strcpy(file, is->image_file_name[0]);
   if ((in_img = img_openfile(file,0)) < 0) {
      printf("Failed to open SOCET ortho image %s\n",file);
      free ( proj_ptr );
      exit(-1);
   }

/*************** Read Support File *******/

   bands = img_query_num_bands(in_img);
   img_query_dimension(in_img,&lines,&samples);
   img_query_tile_size(in_img,&tile_x,&tile_y);
/*   depth = img_query_depth(in_img); */

   /////////////////////////////////////////////////////////////////////////////
   // output ORTHO as a raw binary file
   /////////////////////////////////////////////////////////////////////////////

   out_img=fopen(rawORTHO,"wb");
   if (out_img  == NULL)
   {
      fprintf(stderr,"Failed to open %s\n",rawORTHO);
      fclose(out_img);
      exit(-1);
   }

   num_sec = 9;  // ZERO_BASED
   tenth_lines = lines / (num_sec+1);
   rest_lines = lines - num_sec * tenth_lines;

   cout << "Converting ortho image to a raw file...\n";

   unsigned char *buf1 = new unsigned char [tenth_lines*samples];
   unsigned char *save_buf1 = buf1;
   
   for (sec=0; sec <num_sec; sec++)
   {
      buf1 = save_buf1;
      for (iband=0; iband < bands; iband++)
      {
        img_load_buffer(in_img,sec*tenth_lines,0,tenth_lines,samples,0,
                        buf1,samples,(unsigned char *)"\0");
        for (i=0 ; i < tenth_lines*samples ; i++)
           fwrite(buf1++,sizeof(unsigned char),1,out_img);
      }
      cout << "...Conversion " << 10*(sec+1) << "% Done\n";
   }

   delete [] buf1;

   unsigned char *buf2 = new unsigned char [rest_lines*samples];
   for (iband=0; iband < bands; iband++)
   {
     img_load_buffer(in_img,num_sec*tenth_lines,0,rest_lines,samples,0,
                     buf2,samples,(unsigned char *)"\0");
     for (i=0 ; i < rest_lines*samples ; i++)
        fwrite(buf2++,sizeof(unsigned char),1,out_img);
   }
   delete [] buf2;

   img_closefile(in_img);
   fclose(out_img);

   cout << "...Conversion 100% Done\n";
   
   /////////////////////////////////////////////////////////////////////////////
   // Get ortho unique parameters from support file
   /////////////////////////////////////////////////////////////////////////////
   ortho_sup = (OrthoSensorModel*) is;
   Ylatref = ortho_sup->lat_ref_pt;
   Xlonref = ortho_sup->lon_ref_pt;
   y_realspacing = ortho_sup->interline_dist;  // Vert GSD
   x_realspacing = ortho_sup->interpixel_dist;  // Horiz GSD
   
   // An ORTHO support file stores a (lat_ref, lon_ref) coordinated,
   // which is the center of the lower-left pixel (units are either radians
   // or meters, depending on the map projection.)
   // From the (Ylatref, Xlonref) coordinate, and number of lines/samps
   // of the image, compute the georeference coordinate at the center of
   // the upper left pixel for ISIS
   ulcenter_Xlon = Xlonref;
   ulcenter_Ylat = Ylatref + (y_realspacing*(lines-1));  

   /////////////////////////////////////////////////////////////////////////////
   // Generate the ortho2isis3 script
   /////////////////////////////////////////////////////////////////////////////

   strcpy (productType,"ORT");
   strcpy (byteOrder,"   ");
    
   generate_ss2isis_script (isis_script,
                            prj,
                            productType,
                            byteOrder,
                            outcub_name,
                            layout_flag,
                            lines,
                            samples,
                            x_realspacing,
                            y_realspacing,
                            ulcenter_Xlon,
                            ulcenter_Ylat);
                            
} // END MAIN

