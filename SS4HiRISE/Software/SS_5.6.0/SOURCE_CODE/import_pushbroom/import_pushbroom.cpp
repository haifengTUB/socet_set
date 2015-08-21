///////////////////////////////////////////////////////////////////////////////
//
//_Title IMPORT_PUSHBROOM import linescanner/pushbroom image to SOCET Set v5.4.x
//
//_Desc  This program imports an 8-bit linescanner/pushbroom image - in raw
//       format - to SOCET Set, first as a framing camera with unknown
//       position/orientation, and then converts the framing camera support
//       file to that of a generic pushbroom sensor.  Parameter values for
//       image import and the generic pushbroom support file are supplied in
//       the required *_keywords.lis input file.  Note that *_keywords.lis
//       was created by ISIS3 program create_pushbroom_keywords, on a
//       platform running ISIS3.
//
//_Exec
//       start_socet -single import_pushbroom <project> <fullpath/image>.raw <fullpath/keywords>.lis
//
//_Hist  Oct 21 2008 Elpitha H-Kraus, USGS, Flagstaff Original Version
//       Nov 18 2008 EHK - removed some debug statements
//       Feb 18 2011 EHK - added step of deleting existing *.img* files prior
//                         so import.  When re=importing an image, the minified
//                         images are not updated, so we must delete them first
//                         to insure they correspond to the full res image
//       Feb 08 2012 EHK - modified to import 16-bit tiffs also
//
///////////////////////////////////////////////////////////////////////////////

#include <system_includes.h>
#include <project/proj.h>
#include <util/string_dpw.h>
#include <util/init_socet_app.h>
#include <key/handle_key.h>

/* Define internal functions */
int parse_keywords(char *file, char keyword[], char *value);
void create_pushbroom_sup (char *frmSup, char *keywordFile, double gp_origin_z, char *pushSup);

#define FILELEN 512
#define LINELENGTH 200

int
main(int argc, char *argv[])
{

  // Input arguments
  char     inputImg[FILELEN];        // Input linescanner/pushbroom *raw*
                                   // image processed in ISIS and converted
                                   // to 8-bit.  Images are expected to
                                   // have padding such that the boresight
                                   // is at the center sample

  char     keywordFile[FILELEN];   // Input file of import parameters and
                                   // pushbroom sensor keyword values
                                   // associated with the raw input image.
                                   // This keywordFile was created by
                                   // get_pushbroom_keywords on a
                                   // machine/platform running ISIS

  char  projectName[FILELEN];      // SS project w/o path or .prj ext
  char  projectFile[FILELEN];      // SS project w/.prj ext

  // Output Support files
  char     frmSup[FILELEN];   // temporary framing camera support file
  char     pushSup[FILELEN];  // generic pushbroom support file

  // settingsfile keywords
  int      nl;                 // Number of image lines
  int      ns;                 // Number of image samples
  double   sizex, sizey;       // Image dimensions in x/y, mm
  char     camPath[FILELEN];  // Path to internal databases CAM directory
  char     camFile[FILELEN];  // Camera file needed from frame import

  // project file keywords
  img_proj_struct project;  //SS project structure

  // Environment variables
  char  DBDIR_path[FILELEN];   // decoded path to SS internal_dbs directory
  int   unix_os, windows_os;         // flags indicating platform we are on

  // temporary files and names
  int   prjReadErr;            // Error flag for reading project file
  int   fileExistErr;          // Error flag checking for input files
  char  inputImg_path[FILELEN];  // Path to inputimage
  char  inputImg_type[4];  // flag indicating a raw 8-bit or 16-bit tiff
  char  img_name[FILELEN];     // HiRISE image name (w/o extensions)
  char  socet_img[FILELEN];    // Default SOCET Set image name
  char  img_dir[FILELEN];      // SOCET image directory with full path
  char  settingsfile[FILELEN]; // SOCET settings file for batch import
  FILE  *setfp;                // file pointer to settings file

  // System Call variables
  char  msg[FILELEN];
  int   scan_value;
  char  value[200];
  int   stat;
  int   count;

  /////////////////////////////////////////////////////////////////////////////
  // Check number of command line args and issue help if needed
  // Otherwise initiate the socet set application
  /////////////////////////////////////////////////////////////////////////////

  if (argc != 4) {
    //cerr << "\nRun " << argv[0] << " as follows:\n";
    //cerr << "start_socet -single " << argv[0] << " project fullpath/image.raw fullpath/pushbroom_keywords.lis\n";
    cerr << "\nRun import_pushbroom as follows:\n";
    cerr << "start_socet -single import_pushbroom <project> fullpath/<image>.raw|tif fullpath/<pushbroom_keywords.lis>\n";
    cerr << "\nwhere:\n";
    cerr << "project = SOCET SET project name to import images under\n";
    cerr << "          (path and extension is not required)\n";
    cerr << "fullpath/image.raw|tif = linescanner/pushbroom image in *raw* 8-bit format OR 16-bit tiff format\n";
    cerr << "                     (Path to raw|tif image is required)\n";
    cerr << "fullpath/pushbroom_keywords.lis = output file of get_pushbroom_keywords\n";
    cerr << "                                  associated with input image.\n";
    cerr << "                                  (Path to keywords file is required)\n";
    return -1;
  }

  // Top level SOCET SET initialization  routine.
  // Should be  called  before  any  PCI services.
  init_socet_app( argv[0], argc, argv);

  /////////////////////////////////////////////////////////////////////////////
  //Get input arguments
  /////////////////////////////////////////////////////////////////////////////

  strcpy (projectFile,argv[1]);
  strcpy (inputImg,argv[2]);
  strcpy (keywordFile,argv[3]);

  /////////////////////////////////////////////////////////////////////////////
  // Populate the project structure - with error checking
  /////////////////////////////////////////////////////////////////////////////

  // Make sure project name contains no path, but has .prj extension
  strcpy(projectName,ReturnFileName(projectFile));
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

  /////////////////////////////////////////////////////////////////////////////
  // Check for existence of input files 
  /////////////////////////////////////////////////////////////////////////////

  fileExistErr = 0;

  if (!file_exists(inputImg)) {
    cerr << "\nInput image " << inputImg << " does not exist!\n";
    fileExistErr = 1;
  }

  if (!file_exists(keywordFile)) {
    cerr << "\nInput keyword file " << keywordFile << " does not exist!\n";
    fileExistErr = 1;
  }

  if (fileExistErr)
    exit(1);

  /////////////////////////////////////////////////////////////////////////////
  // Decode DBDIR environment variable while alse determining platform
  // (Windows or Unix)
  /////////////////////////////////////////////////////////////////////////////

  unix_os = 0;
  windows_os = 0;
  str_decode_env_path ("$DBDIR",DBDIR_path); // This works for both unix and
                                             // windows?  At the prompt, the
                                             // syntax for environment variables
                                             // in windows is %DBDIR%

  if(strstr(DBDIR_path,"DBDIR")==0)  // the decode worked properly
    if (strstr(DBDIR_path,":\\")!=0) // look for :\ in path for windows
      windows_os = 1;
    else
      unix_os=1;
  else {
   cerr << "Unable to decode DBDIR environment variable...is it missing?" << endl;
   exit (1);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Determine if input image is 8-bit raw or 16-bit tif
  /////////////////////////////////////////////////////////////////////////////
  
  if (strstr(inputImg,"tif")!=0) 
    strcpy(inputImg_type,"tif");
  else
    strcpy( inputImg_type,"raw");

  /////////////////////////////////////////////////////////////////////////////
  // Generate SOCET SET and temporary filenames
  /////////////////////////////////////////////////////////////////////////////

  // extract HiRISE image name from input inputImg name
  ReturnPathAndFile(inputImg,inputImg_path,img_name);
  StripFileExts (img_name);

  // if local directory is defaulted to for location of inputImg, make sure
  // the syntax is correct

  if (strlen(inputImg_path)==0)
    if (unix_os)
      strcpy(inputImg_path,"./");
    else
      strcpy(inputImg_path,".\\");

  // Create a temporary SOCET Set framing import settings file name
  build_file_name(settingsfile, inputImg_path, "xxtempframe", ".set");

  // Create HiRISE support file name, with SOCET Set data path
  build_file_name(pushSup, project.project_data_path, img_name, ".sup");

  // Create temporary framing camera support file name, with SOCET Set data path
  build_file_name(frmSup, project.project_data_path, img_name, ".sup_frame");

  // Build path+file to the default.cam file for framing camera import
  add_path_if_none (DBDIR_path, "CAM", camPath);
  build_file_name (camFile, camPath,"default",".cam");
    
  /////////////////////////////////////////////////////////////////////////////
  // If pushSup already exists, we are re-importing the images, so delete all
  // image files first by running usgs_delete_image
  /////////////////////////////////////////////////////////////////////////////

  if (file_exists(pushSup)) {
    sprintf(msg,"start_socet -single usgs_delete_image %s\n",pushSup);
    printf("%s\n",msg);
    scan_value = system(msg);
  }
  
  /////////////////////////////////////////////////////////////////////////////
  // Retrieve or generate needed framing camera import settingsfile parameters
  /////////////////////////////////////////////////////////////////////////////

  // Get image size from Keyword list (PVL) file
  stat = parse_keywords(keywordFile, "TOTAL_LINES", value);
  nl = atoi(value);

  stat = parse_keywords(keywordFile, "TOTAL_SAMPLES", value);
  ns = atoi(value);

  // Set images size to nl and ns for frame import
  // (We are not scaling to mm since this is a temporary import)
  sizex = ns;
  sizey = nl;

  // Generate the default SOCET Set *.img|tif file name as a 
  // keyword value in the settings file
  if (strcmp(inputImg_type,"tif")==0) 
    strcpy (socet_img,concat(img_name,".tif"));
  else
    strcpy (socet_img,concat(img_name,".img"));

  /////////////////////////////////////////////////////////////////////////////
  // Populate the temporary framing camera settings file
  /////////////////////////////////////////////////////////////////////////////

  setfp = fopen (settingsfile,"w");
  if (setfp == NULL) {
    printf ("Unable to open frame import settings file %s\n",settingsfile);
    exit (1);
  }

  fprintf (setfp,"setting_file                  1.1\n");
  fprintf (setfp,"multi_frame.project                 %s\n",projectFile);
  fprintf (setfp,"multi_frame.cam_calib_filename      %s\n",camFile);
  fprintf (setfp,"multi_frame.create_files            IMAGE_AND_SUPPORT\n");
  fprintf (setfp,"multi_frame.atmos_ref               0\n");
  if (strcmp(inputImg_type,"tif")==0)
    fprintf (setfp,"multi_frame.auto_min                NO\n");
  else
   fprintf (setfp,"multi_frame.auto_min                YES\n");
  fprintf (setfp,"multi_frame.digital_cam             YES\n");
  fprintf (setfp,"multi_frame.input_image_filename    %s\n",inputImg);
  if (strcmp(inputImg_type,"tif")==0)
    fprintf (setfp,"multi_frame.output_format           img_type_tiff_tiled\n");
  else
    fprintf (setfp,"multi_frame.output_format           img_type_vitec\n");
  fprintf (setfp,"multi_frame.output_name             %s\n",socet_img);
  fprintf (setfp,"multi_frame.output_location         %s\n",project.project_image_location);
  fprintf (setfp,"multi_frame.cam_loc_ang_sys         OPK\n");
  fprintf (setfp,"multi_frame.cam_loc_ang_units       UNIT_DEGREES\n");
  fprintf (setfp,"multi_frame.cam_loc_xy_units        UNIT_DEGREES\n");
  fprintf (setfp,"multi_frame.img_size_lines          %d\n",nl);
  fprintf (setfp,"multi_frame.img_size_samps          %d\n",ns);
  fprintf (setfp,"multi_frame.sizex                   %lf\n",sizex);
  fprintf (setfp,"multi_frame.sizey                   %lf\n",sizey);
  fprintf (setfp,"multi_frame.orientation             1\n");

  fclose (setfp);

  /////////////////////////////////////////////////////////////////////////////
  // Import image to Socet Set as a framing camera
  /////////////////////////////////////////////////////////////////////////////

  sprintf(msg,"start_socet -single multi_frame -a frame -batch -s %s\n",
          settingsfile);
  printf("%s\n",msg);
  scan_value = system(msg);
  
  /////////////////////////////////////////////////////////////////////////////
  // If this is a tif image, minify it
  // (Autominification during import results in a single tif image with internal
  // pyramids - an undesirable feature for gdal....)
  /////////////////////////////////////////////////////////////////////////////

   if(strcmp(inputImg_type,"tif")==0) {
    //file_remove(settingsfile);  // delete xxtempframe.set
   
    // Create a temporary SOCET Set minification settings file name
    build_file_name(settingsfile, inputImg_path, "xxtempminify", ".set");
  
    setfp = fopen (settingsfile,"w");
    if (setfp == NULL) {
      printf ("Unable to open frame import settings file %s\n",settingsfile);
      exit (1);
    }

    fprintf (setfp,"setting_file                     1.1\n");
    fprintf (setfp,"min.project                     %s\n",projectFile);
    fprintf (setfp,"min.resampling_method      BILINEAR\n",pushSup);
    fprintf (setfp,"min.output_location          %s\n",project.project_image_location);
    fprintf (setfp,"min.location_method         SINGLE_LOCATION\n",pushSup);
    fprintf (setfp,"min.create_single_file        NO\n",pushSup);
    fprintf (setfp,"min.single_file_format        img_type_tiff_tiled\n",pushSup);
    fprintf (setfp,"min.input_sup                  %s\n",pushSup);
    fprintf (setfp,"min.condor_selected          0\n",pushSup);
    fprintf (setfp,"min.condor_machine          Any Machine\n",pushSup);
    
    fclose (setfp);
    
    sprintf(msg,"start_socet -single minifier -batch -s %s\n",
            settingsfile);
    printf("%s\n",msg);
    scan_value = system(msg);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Rename the framing camera support file to *.sup_frame (command depends on
  // operating system)
  //
  // NOTE: I tried using the rename command for DOS, but had syntax errors
  //       when supplying paths to the file.  Goggling rename, I found a page
  //       that said the move command was available before rename, so I
  //       tried move, and it worked....
  /////////////////////////////////////////////////////////////////////////////

  if (unix_os)
    sprintf(msg,"/bin/mv %s %s\n",pushSup,frmSup);
  else
    sprintf(msg,"move %s %s\n",pushSup,frmSup);

  scan_value = system(msg);

  /////////////////////////////////////////////////////////////////////////////
  // Create the generic pushbroom support file by replacing frameing camera
  // keywords with pushbroom keywords in the keywordFile
  /////////////////////////////////////////////////////////////////////////////

  create_pushbroom_sup (frmSup, keywordFile, project.gp_origin.z, pushSup);

  /////////////////////////////////////////////////////////////////////////////
  // Delete temporary files
  /////////////////////////////////////////////////////////////////////////////

  //file_remove(settingsfile);
  file_remove(frmSup);

  cout << "\nImport complete\n";

} // end of import_pushbroom

/**************  parse_keywords   ******************
*                                                  *
*  This routine reads a list of SS support file    *
*  keywords and grabs the desired value            *
*                                                  *
*  EHK  USGS  Oct 24, 2008                         *
*                                                  *
****************************************************/
int
parse_keywords(char *file, char keyword[], char *value)

{

/***************      declaration           **********/
   int i,ii=0,index=0,scan_value;
   char ch,id[512],key[20];

   FILE *fp;


   fp = fopen(file,"r");
   if (fp == NULL)
     {
      printf("\ncan't open the input file %s!\n",file);
      exit(1);
   }

/*************** Read Header file write label *******/

  fgets(id,512,fp);
  while (strstr(id,keyword) == NULL)
    fgets(id,512,fp);

  sscanf (id,"%s %s",key,value);

  fclose(fp);
  return(-1);

} // end of parse_keywords

void
create_pushbroom_sup (
   char *frmSup,       // Input framing camera support file
   char *keywordFile,  // Input pushbroom keywords list
   double gp_origin_z, //project's ground point origin height, m
   char *pushSup)      // Output generic pushbroom support file
{

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINELENGTH 200

///////////////////////////////////////////////////////////////////////////////
//
//_Title CREATE_PUSHBROOM_SUP generates a SS generic pushbroom support file
//
//_Desc	This subroutine takes as input a framing camera support file and a
//       corresponding pushbroom keyword file, and merges the content of the
//       two files to create a SS Generic Pushbroom support file
//
//_Hist	Oct 23 2008 Elpitha H. Kraus, USGS, Flagstaff Original Version
//
//_End
//
///////////////////////////////////////////////////////////////////////////////


  FILE	*frmfp;      // File pointer to frmSup file
  FILE	*keyfp;      // File pointer to keywordFile
  FILE	*pushfp;     // File pointer to pushSup file

  char	base_line[LINELENGTH]; // Line of data in the base class
                               // sensor model section of the SOCET
                               // SET framing camera support file
  char  push_line[LINELENGTH]; // Line of data in the pushbroom
                               // keyword file

  char  updated_line[LINELENGTH]; // push_line with updated information

  int           count;
  int           len;
  int           ret;
  int		i;

  /////////////////////////////////////////////////////////////////////////////
  // Open input and output files
  /////////////////////////////////////////////////////////////////////////////

  frmfp = fopen (frmSup,"r");
  if (frmfp == NULL) {
    printf ("Unable to open input framing camera support file: %s\n",frmSup);
    exit (1);
  }

  keyfp = fopen (keywordFile,"r");
  if (keyfp == NULL) {
    printf ("Unable to open input pushbroom keyword file: %s\n",keywordFile);
    fclose(frmfp);
    exit (1);
  }

  pushfp = fopen (pushSup,"w");
  if (pushfp == NULL) {
    printf ("Unable to open output generic pushbroom support file: %s\n",pushSup);
    fclose(frmfp);
    fclose(keyfp);
    exit (1);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Now output the Generic Pushbroom support file....
  //
  // First, copy the Base Class SOCET SET Sensor Model keywords and values from 
  // the frame support file to the generic pushbroom support file.  Note
  // that we must update the RECTIFICATION TERMS, GROUND_ZERO and LOAD_PT.
  // In addition, make sure IMAGE_MOTION is 0.
  /////////////////////////////////////////////////////////////////////////////

  // Get base class keywords up to RECTIFICATION_TERMS from frame sup file
  fgets(base_line,LINELENGTH,frmfp);
  while (strstr(base_line,"RECTIFICATION_TERMS") == NULL) {
    fputs (base_line,pushfp);
    fgets(base_line,LINELENGTH,frmfp);
  }

  // Skip Retification terms in framing support file
  for (i=0; i<2; i++) 
    fgets(base_line,LINELENGTH,frmfp);

  // Get Retification terms from keyword file */
  for (i=0; i<=2; i++) { 
    fgets(push_line,LINELENGTH,keyfp);
    fputs (push_line,pushfp);
  }

  // Get base class keywords up to GROUND_ZERO from frame sup file
  fgets(base_line,LINELENGTH,frmfp);
  while (strstr(base_line,"GROUND_ZERO") == NULL) {
    fputs (base_line,pushfp);
    fgets(base_line,LINELENGTH,frmfp);
  }

  // Get GROUND_ZERO from keyword file
  fgets(push_line,LINELENGTH,keyfp);
  len = strlen(push_line);
  while (push_line[len] != ' ')
     len--;
  push_line[len] = '\0';
  sprintf(updated_line,"%s %.14le\n",push_line,gp_origin_z);
  fputs (updated_line,pushfp);
     
  // Get base class keywords up to LOAD_PT from frame sup file
  fgets(base_line,LINELENGTH,frmfp);
  while (strstr(base_line,"LOAD_PT") == NULL) {
    fputs (base_line,pushfp);
    fgets(base_line,LINELENGTH,frmfp);
  }

  // Get LOAD_PT from keyword file
  fgets(push_line,LINELENGTH,keyfp);
  len = strlen(push_line);
  while (push_line[len] != ' ')
     len--;
  push_line[len] = '\0';
  sprintf(updated_line,"%s %.14le\n",push_line,gp_origin_z);
  fputs (updated_line,pushfp);


  // Get base class keywords up to COORD_SYSTEM from frame sup file
  fgets(base_line,LINELENGTH,frmfp);
  while (strstr(base_line,"COORD_SYSTEM") == NULL) {
    fputs (base_line,pushfp);
    fgets(base_line,LINELENGTH,frmfp);
  }

  // Get COORD_SYSTEM from keyword file
  fgets(push_line,LINELENGTH,keyfp);
  fputs (push_line,pushfp);

  // Get base class keywords up to IMAGE_MOTION from frame sup file
  fgets(base_line,LINELENGTH,frmfp);
  while (strstr(base_line,"IMAGE_MOTION") == NULL) {
    fputs (base_line,pushfp);
    fgets(base_line,LINELENGTH,frmfp);
  }

  // Get IMAGE_MOTION from keyword file
  fgets(push_line,LINELENGTH,keyfp);
  fputs (push_line,pushfp);

  // Get remaining base class keywords (up to SENSOR_TYPE) from
  // frame sup file
  fgets(base_line,LINELENGTH,frmfp);
  while (strstr(base_line,"SENSOR_TYPE") == NULL) {
    fputs (base_line,pushfp);
    fgets(base_line,LINELENGTH,frmfp);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Output the SS generic pushbroom specific keywords to the support file
  /////////////////////////////////////////////////////////////////////////////

  while (fgets(push_line,LINELENGTH,keyfp) != NULL)
    fputs (push_line,pushfp);

  /////////////////////////////////////////////////////////////////////////////
  // Copy the "planet parameters" stored in the last section of the framing
  // camera support file to the pushbroom support file
  /////////////////////////////////////////////////////////////////////////////

  while (strstr(base_line,"ELLIPSOID") == NULL)
    fgets(base_line,LINELENGTH,frmfp);
  fputs(base_line,pushfp);
  for(i=1; i<=7; i++) {
    fgets(base_line,LINELENGTH,frmfp);
    fputs(base_line,pushfp);
  }

  /////////////////////////////////////////////////////////////////////////////
  // All done so cleanup and return
  /////////////////////////////////////////////////////////////////////////////

  fclose (frmfp);
  fclose (keyfp);
  fclose (pushfp);

  return;
}
