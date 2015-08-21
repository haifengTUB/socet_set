////////////////////////////////////////////////////////////////////////////////
//
//_Title CALCORTHOBDRY reports the LL and UR coord. range to match a DEM range
//
//_Desc  This is a SOCET Set devkit program.  It calculates the lat/lon or
//       X/Y boundaries (as appropriate) of the DEM used to orthorectify images.
//       This boundary is then used in orthophotogeneration to get 1:1
//       pixel correspondence between the orthoimage and the DEM. 
//
//       Input parameters are:
//
//              SS_project
//              socet_dem
//
//       Output files are:
//
//              <project_data_dir>/calcOrthoBdry.log
//
//_Hist	May 18 2007 Elpitha H. Kraus, USGS, Flagstaff Original Version
//      Jun 18 2007 EHK - corrected formatted output error when reporting
//                        x/y boundary in meters
//      Jul 03 2007 EHK - changed float to double precision on *realspacing
//                        variables to avoid round off errors
//      Nov 03 2008 EHK - Modified to use dthHeader library so as not to
//                        require reading in the entire DEM and run into
//                        memory limitations on Windows
//                        Also removed include of iostream.h (it is not
//                        in Windows, and not necessary under Solaris)
//      Nov 04 2008 EHK - Changed log file from print.prt to calcOrthoBdry.log
//      Oct 14 2010 EHK - Changed log file from calcOrthoBdry.log to calcOrthoBdry_<dem>.log
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
#include <dtmAccess/DtmHeader.h>
#include <dtmAccess/dtm_edit_util.h>
#include <dtmAccess/fom_defs.h>
#include <project/proj.h>
#include <util/string_dpw.h>
#include <util/init_socet_app.h>
#include <ground_point.h>
#include <image_point.h>
#include <coord/coord_system.h>
#include <coord/lat_long_xy.h>

// Set up ISIS NULL values
#define FILELEN 512
#define NO_ERRS 0
#define PARINV_ERR -1

// prototypes
int stripp(char instr[], char outstr[], int position);
int parse_label(char *file, char *keyword, char *value);
int writeToLog(char *msg, char *logfile);
int writeReport(char *msg, char *logfile);

void main(int argc,char *argv[]) 
{
   // DECLARATIONS:

   // Input variables
   char dem[FILELEN];
   char demName[FILELEN];
   char fname[FILELEN];
   char projectFile[FILELEN];
   char projectName[FILELEN];
   char logfile[FILELEN];
   char logfileName[FILELEN];

   // DEM Header Variables
   DtmHeader* di_header;
   double x_realspacing, y_realspacing; //real spacing of a DEM, in project coordinates
   ground_point_struct ll_corner; 
   ground_point_struct ur_corner; 

   // Project File Variables
   img_proj_struct project;
   img_proj_struct_ptr  proj_ptr;
   char prj[FILELEN];
   int   prjReadErr;            // Error flag for reading project file
   int coord_sys;

   // Misc declarations
   int i,ii=0,scan_value;
   int ret;
   char value[FILELEN];
   char msg[256];
   ground_point_struct ortho_ll, ortho_ur;
   char DMS[20];

  /////////////////////////////////////////////////////////////////////////////
  // Check number of command line args and issue help if needed
  // Otherwise initiate the socet set application
  /////////////////////////////////////////////////////////////////////////////

   if (argc!=3) {
    //cerr << "\nRun " << argv[0] << " as follows:\n\n";
    //cerr << "start_socet -single " << argv[0] << " <project> <DEM>\n";
    cerr << "\nRun calcOrthoBdry as follows:\n\n";
    cerr << "start_socet -single calcOrthoBdry <project> <DEM>\n";
    cerr << "\nwhere:\n";
    cerr << "project = SOCET SET project name\n";
    cerr << "          (path and extension is not required)\n";
    cerr << "DEM = SOCET SET DEM to be used for Orthophoto Generation\n";
    cerr << "      (path and extension is not required)\n";
    cerr << "\ncalcOrthoBdry will output:\n";
    cerr << "           Upperleft and Lower Right coordinates to be used in Orthophoto Generation.\n";
    cerr << "           Coordinates will be written to the screen and appended to logfile:\n";
    cerr << "           <project_data_path>/calcOrthoBdry_<DEM>.log.\n";
    exit(1);
   }

   // Top level SOCET SET initialization  routine.
   // Should be  called  before  any  PCI services.
   init_socet_app ( argv[0], argc, argv);

  /////////////////////////////////////////////////////////////////////////////
  //Get input arguments
  /////////////////////////////////////////////////////////////////////////////

   strcpy(prj,argv[1]);
   strcpy(dem,argv[2]);

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
  // Get variations of DEM file representation:
  //      a) fullpath + name + ".dth"
  //      b) fullpath + name
  // and make sure it exists
  /////////////////////////////////////////////////////////////////////////////

  // Make sure project name contains no path, but has .prj extension
  strcpy(demName,ReturnFileName(dem));
  StripFileExt(demName);
  build_file_name(dem, project.project_data_path, demName, ".dth");
  strcpy(fname, dem);
  StripFileExt(fname);

  if (!file_exists(dem)) {
    cerr << "\ndem " << demName << " does not exist!\n";
    cerr << "(NOTE: looking for " << dem << ")\n";
    exit(1);
  }
  // Get logfile name
  strcpy(logfileName,"calcOrthoBdry_");
  strcat(logfileName,demName);
  build_file_name(logfile, project.project_data_path, logfileName, ".log");

   //output calcOrthoBdry command that was issued to calcOrthoBdry.log file
   sprintf(msg,"start_socet -single %s %s %s",
      argv[0],argv[1],argv[2]);
   writeToLog(msg, logfile);
  
  /////////////////////////////////////////////////////////////////////////////
  // Set the DEM header and load DEM structure
  /////////////////////////////////////////////////////////////////////////////
   proj_ptr = getCurrentProjStruct();
   di_header = new DtmHeader(proj_ptr);
   di_header->load(fname);

   // valid DEM?
   if(di_header->dtmFormat() != DTM_GRID) {
      cerr << "ERROR: input DEM is not in GRID format\n";
      free(di_header);
      free(proj_ptr);
      exit(-1);
   }

   // We have a DEM (ie, di.header->dtmFormat() == DTM_GRID)
   // Get needed parameters from SOCET DEM header that pertain 
   // to both geographic and map projected projects

   x_realspacing = di_header->xRealSpacing();
   y_realspacing = di_header->yRealSpacing();
   ll_corner = di_header->llCorner(); 
   ur_corner = di_header->urCorner(); 

   // Parse Project file for project's coordinate system

   ret = parse_label(prj, "COORD_SYS", value);
   coord_sys = atoi(value);

   // Get boundary of DEM and expand it by half the x/y spacing
   // as appropriate to match DEM range of a corresponding ISIS
   // cube (this works for radians or meters)
   // ------------------------------------------------------
   // The LOWER_LEFT_XYZ and UPPER_RIGHT_XYZ coordinates stored in
   // the .dth file are the coordinates of the upper right and lower
   // left posts of the DEM, which is equivalent to the real-world
   // coordinate at the center of the upper right and lower left
   // pixels in an ISIS cube.  Old-ISIS (version 1) expects the
   // lat/lon or X/Y range to cover from "edge-to-edge" of the image
   // space, so offset the LOWER_LEFT_XYZ and UPPER_RIGHT_XYZ by half
   // the  x/y spacing to get an edge-to-edge range

   ortho_ll.y = (ll_corner.y - y_realspacing/2.0 ); 
   ortho_ll.x = (ll_corner.x - x_realspacing/2.0 ); 
   ortho_ur.y = (ur_corner.y + y_realspacing/2.0 ); 
   ortho_ur.x = (ur_corner.x + x_realspacing/2.0 ); 

   //If this is a geographic project, convert radians to DMS and report
   //corners.  If this project in a map projection (meters), just report
   //corners.

   if (coord_sys == 1) { //Geographic coords...convert to DMS and report
     convert_radians_to_deg_min_sec_string(ortho_ll.x,1,0,5,DMS);
     sprintf(msg,"\nLL: Lon %s",DMS);
     writeReport(msg, logfile);

     convert_radians_to_deg_min_sec_string(ortho_ll.y,0,0,5,DMS);
     sprintf(msg,"LL: Lat %s",DMS);
     writeReport(msg, logfile);

     convert_radians_to_deg_min_sec_string(ortho_ur.x,1,0,5,DMS);
     sprintf(msg,"\nUR: Lon %s",DMS);
     writeReport(msg, logfile);

     convert_radians_to_deg_min_sec_string(ortho_ur.y,0,0,5,DMS);
     sprintf(msg,"UR: Lat %s",DMS);
     writeReport(msg, logfile);

   }
   else {
     sprintf(msg,"\nLL: X %.5lf",ortho_ll.x); 
     writeReport(msg, logfile);
     sprintf(msg,"LL: Y %.5lf",ortho_ll.y); 
     writeReport(msg, logfile);
     sprintf(msg,"\nUR: X %.5lf",ortho_ur.x); 
     writeReport(msg, logfile);
     sprintf(msg,"UR: Y %.5lf",ortho_ur.y); 
     writeReport(msg, logfile);
   } 

} // END MAIN


//*****************************************************************************
//*****************************************************************************

int writeToLog(char *msg, char *logfile)

{
  //output isis command to log file

  FILE *fp;

  // Open log file
  fp = fopen(logfile,"a");
  if (fp == NULL) {
     printf("\ncan't open %s file!\n", logfile);
     return(1);
  }

  fprintf (fp,"\nCALCORTHOBDRY command issued:\n");
  fprintf (fp,"-----------------------------\n");
  fprintf (fp,"%s\n",msg);

  fclose (fp);
  return (0);

} // End of writeToLog

int writeReport(char *msg, char *logfile)

{
  //output report to logfile and screen

  FILE *fp;

  // Open log file
  fp = fopen(logfile,"a");
  if (fp == NULL) {
     printf("\ncan't open %s file!\n",logfile);
     return(1);
  }

  fprintf (fp,"%s\n",msg);

  fclose (fp);

  printf("%s\n",msg);

  return (0);

} // End of writeReport




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

int stripp(char instr[], char outstr[], int position)
/*************************************************************************
*_Title stripp Trim .xxx off file name
*_Nobinding
*_Args  Type  Name             I/O  Description
*_Parm  char  instr[];          I   String to trim blanks from
*_Parm  char  outstr[];         O   Buffer to receive trimmed string
*_Parm  int   position;         I   how many non-char to skip
*_Parm  int   strippara;        O   Returns number of characters copied
*_Desc  strippara strips non-character like - , ( from string
*       and returns the result.
*_Keys  STRING
*       May 22 1997 Trent Hare - Changed to isistypes
*_END
*************************************************************************/
{
    int i=1;             /* Number characters copied */
    int ii=1;             /* Number characters copied */
    int j=1;             /* Number characters copied */

/****************************************************************
  Find the first character from start of string
*****************************************************************/
    ii=0;
    for (j=0; j<position; j++) {
      i=0;
      while((instr[ii] != '.'))
      {
          outstr[i] = instr[ii];
          i++;
          ii++;
      }
      ii++;
    }
    outstr[i] = '\0';
    return (ii);
}

