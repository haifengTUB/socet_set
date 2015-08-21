////////////////////////////////////////////////////////////////////////////////
//
//_Title DEM2ISIS3 outputs a SOCET DEM as a raw 32-bit image plus ISIS3 script.
//
//_Desc  This is a SOCET Set program that uses SOCET DEV_KIT routines to access
//       a SOCET DEM and output a raw 32-bit DEM file, a raw 8-bit FOM file,
//       and an associated ISIS3 processing script.  The raw files and script
//       are to be transferred to an ISIS machine for port to ISIS3.
//
//       DEMs in Geographic Coordinates or Polar Stereographic (Grid) coordinates
//       are currently supported.
//
//       Note that TINs are not supported in ISIS, so this program only works
//       for DEMs.
//
//       Input parameters are:
//
//              SS_project
//              socet_dem.dth
//              isis_dem.cub
//              layout_flag
//
//
//       Output files are:
//
//              ./isis_dem.raw
//              ./FOM_isis_dem.raw (FOM filename auto created by either:
//                                  1) replace DEM with FOM in <isis_dem>, if "DEM" string exists, or
//                                  2) adding FOM prefix t<isis_dem> if "DEM" string does not exist
//              ./isis_dem2isis3.sh
//
//       This isis_dem2isis3.sh script will generate up to four output
//       files:
//              SS_isis_dem.cub (For Geographic projects of ellipsoids only)
//              isis_dem.cub
//              FOM_isis_dem.cub
//              isis_dem_layout.cub (optional)
//
//
//       For DEMs in the Geographic Coordinate System:
//       -------------------------------------------------------------
//       DEMs in Geographic Coordinates in SS are in essentially
//       equirectangular map projection.  However, SS calculates scaling
//       radii such that the resulting pixels are square w.r.t. meter
//       measurements (and non-square w.r.t. degrees.)  In ISIS, a
//       different scaling radius is calculated for the equi-rectangular
//       map projection, and the resulting pixels are square in degree-space
//       but not meter-space.  For bodies treated as a spheroid, SS and ISIS
//       scaling radii are the same, and products from either package are
//       compatable.  However, for ellipsoidal planets such as Mars, the
//       different approaches to scaling radii become evident. Because
//       of these differences, we output two versions of the DEM:
//       SS_isis.cub and isis.cub.
//
//       SS_isis.cub maintains the original SS DEM geometry with respect to
//       scaled radii values, and is retained in case someone needs to use the
//       non-resampled DEM data.  Map projection information is added to the
//       ISIS labels - with the resolution based on the Y spacing of the SOCET
//       DEM in dg/px. For ellipsoidal bodies, a SS2ISIS_IMPORT_NOTES keyword
//       is added to the SS_isis.cub labels stating the scaling is not
//       compatible with ISIS.  Furthermore, for Mars, these DEMs are in
//       the ographic latitude, +West longitude system.
//
//       isis_dem.cub is our standard ISIS distribution cube.  For
//       ellipsoidal bodies, pixels are scaled to be compatible with ISIS
//       radii scaling (so as to have 'square pixels' in ISIS) and map
//       projection labels.  Furthermore, for Mars, DEMs above +/- 65 degrees
//       latitude are in polar projection.  Between -65 to + 65 degrees latitude,
//       they are in Equirectangular projection with clat=0, clon=180.  All
//       standard Mars DEMs are in the +E 360 degree longitude system, with
//       ocentric lats.
//
//	 Finally, for DEMs in geographic coordinates, the lat/lon range is well
//       known in SOCET and this information is added to the ISIS labels.
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
//                  DEMs between +/- 65 degrees lat will be in:
//                     - equirectangular map projection with clat=0, clon=180
//                  DEMs from +/- 65 degrees lat to the poles will be in
//                     - polar stereographic map projection.
//                  Note: the center/average lat of the cube is used to
//                  determine if cube falls in 'polar region'.
//
//
//       For DEMs in Polar Stereographic (Grid) Coordinates:
//       -------------------------------------------------------------
//       SOCET Set DEMs produced in a Polar Steregraphic Grid coordinate
//       system will not have the same scaling radii issues that we have with
//       geographic coordinates so long as the clat and clon of the projection
//       is at the poles (i.e., +/-90 degrees.)  In this case, we need only
//       produce the standard isis.cub, using the polar stereographic
//       projected coordinate of the upperleft DEM post to geo-reference the
//       DEM in ISIS.
//
//       Also note that the lat/lon range of polar stereographic DEMs will *not*
//       be placed on the ISIS cube labels.  The lat/lon range of DEMs in
//       grid coordinates is not readily known out of SOCET SET (especially
//       for map-projections without straight meridians and parallels, such as
//       polar stereographic).  However, unlike ISIS2, ISIS3 does not require
//       a lat/lon range on the labels for further processing, so no
//       attempts are made to add a lat/lon range.
//
//_Hist Sep 26 2008 EHK  Orig Version - port dem2isis2 to ISIS3
//                       Also output isis commands to a script file
//      Nov 03 2008 EHK  Windows apparently has round off error that is
//                       mapping what we want as 32-bit NULL pixels as LRS
//                       pixels.  The solution for now is to stretch the
//                       DEM cubes so that LRS become NULL.
//      Nov 10 2008 EHK  Changed interpolation method for map2map to
//                       bilinear, rather than the default of cubic convolution.
//                       (Cubic convolution will map valid pixels to nulls in
//                       areas of significant DN change...this is a feature of
//                       the polynomial fit for cubic convolution.)
//      Nov 12 2008 EHK  Cleaned up deletion of temporary files
//      Jun 11 2009 EHK  Added export of polar stereographic DEMs from SS for
//                       LRO
//      Jun 15 2009 EHK  Updated documentation.
//      Aug 20 2009 EHK  Added WGS_84 (Earth) as a supported datum
//      Jan 20 2010 EHK  Updated ISIS commands for version ISIS3.1.21
//                       and added setisis command to ISIS script compatible
//                       with ISIS commands.
//      Feb 21 2010 EHK  Handle MOON2000 and TITAN2000 as a Socet Set target names
//                       to conform with SSv55 geodetic files
//      May 26 2010 EHK  Corrected bug that put Moon2000 in positive West system
//                       when exporting to ISIS.
//      Jun  3 2010 EHK  Added export of FOM file as a standard product and
//                       reorganized creation of raw files to ISIS cubes
//                       to avoid rerunning raw2isis, etc, multiple times.
//      Sep 20 2010 EHK  Reorganized code to call subroutines in common with
//                       ortho2isis3.  Also changed output scripts not to delete
//                       raw binary files.
//      Sep 28 2010 EHK  When creating FOM filename, added search for DTM
//                       along with DEM (to be replaced by FOM)
//      Oct 17 2011 EHK  Added SS_ prefix to output *.raw files to avoid confusion with
//                        the 'standard' ISIS cubes having the same core name, but different
//                         number of lines and samples
//      Jul 08 2011 EHK chnaged dem2isis3 to dem2isis3.exe in help for start_socet command
//                         (the .exe extension will insure program will run in a command prompt),
//                         Changed variable name outcub_name to outcubName
//
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

// For the record, these are the ISIS NULL values
// I tried to use under windows, but what should have
// resulted as NULL pixels were LRS (round-off???)
//#define NULL3 -0.3402822655089E+39       (value from dem2isis2)
#define NULL3 -3.4028226550889044521e+38 // (modified value Trent uses in gdal) 
//Set ISIS3 NULL value
//const int INULL4 = 0xFF7FFFFB;
//const float NULL4 =(*((const float *) &INULL4));

// prototypes
extern int parse_label(char *file, char *keyword, char *value);
extern int getTargetInfo (char *ellipsoid, char *isisTargName, char *isisTargDef,
	        char *ographicPosLonDir, char *ocentricPosLonDir);
extern int writeToScript(char *isis_script, char *command);
extern int generate_ss2isis_script (char *isis_script, char *prj,
            char *productType, char *byteOrder, char *outcubName,
            char *layout_flag, int lines, int samples, double x_realspacing,
            double y_realspacing, double ulcenter_Xlon, double ulcenter_Ylat);
            
void main(int argc, char *argv[]) {
	// DECLARATIONS:

	// Input variables
	char dem[FILELEN];
	char demName[FILELEN];
	char fname[FILELEN];
	char outcub[FILELEN];
	char outcubName[FILELEN];
	char layout_flag[1];

	// DEM Header Variables
	unsigned char dem_loaded;
	DtmGrid* di;
	DtmHeader* di_header;
	int ncols, nrows;
	double x_realspacing, y_realspacing; //real spacing of a DEM, in project coordinates
	ground_point_struct ll_corner;
	ground_point_struct ur_corner;
	double ulcenter_Xlon, ulcenter_Ylat;

	// Project File Variables
	img_proj_struct project;        //SS project structure
	img_proj_struct_ptr proj_ptr;  //Pointer to SS project file
	char prj[FILELEN];         //SS project with full path and extension
	char projectName[FILELEN];      // SS project w/o path or .prj ext
	char projectFile[FILELEN];      // SS project w/.prj ext

	// ISIS variables
	char isis_script[FILELEN];
    float null;

	// Environment variables
	char DBDIR_path[FILELEN];   // decoded path to SS internal_dbs directory
	int unix_os, windows_os;         // flags indicating platform we are on

	// temporary files and names
	int prjReadErr;            // Error flag for reading project file
	int demReadErr;            // Error flag for reading project file
	int fileExistErr;          // Error flag checking for input files

	// Misc declarations
    char rawDEM[FILELEN];
    char rawFOM[FILELEN];
    char FOM_outcubName[FILELEN];
	int i, ii = 0, scan_value;
	int ret;
	char value[FILELEN];
	int x, y;
	char command[512];
	float demz[1];
	float z[1];
	char fom[1];
	double rad2deg = 180.0 / M_PI;  //convert deg to radians and back
	ground_point_struct ul_corner;
	ground_point_struct gp_ul, gp_lr;
	char productType[4];   //productType=DEM,FOM or ORT
	char byteOrder[4];     //needed for DEMs only
	FILE *ofp_DEM;
	FILE *ofp_FOM;

	//Set what NULL to use
	null = (float) NULL3;

	/////////////////////////////////////////////////////////////////////////////
	// Check number of command line args and issue help if needed
	// Otherwise initiate the socet set application
	/////////////////////////////////////////////////////////////////////////////

	if (argc < 4) {
		cerr << "\nRun dem2isis3 as follows:\n";
		cerr << "start_socet -single dem2isis3.exe <project> <socet_dem> <isis.cub> <layout_flag>\n";
		cerr << "\nwhere:\n";
		cerr << "project = SOCET SET project name to export DEM from\n";
		cerr << "          (path and extension is not required)\n";
		cerr << "socet_dem = SOCET SET dem to export\n";
		cerr << "          (path and extension is not required)\n";
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

	strcpy(prj, argv[1]);
	strcpy(dem, argv[2]);
	strcpy(outcub, argv[3]);
	if (argc == 5)
		strcpy(layout_flag, argv[4]);
	else
		strcpy(layout_flag, "n");

	/////////////////////////////////////////////////////////////////////////////
	// Populate the project structure - with error checking
	// and get variations of project file representation:
	//     a) name only
	//     b) name + extension
	//     c) fullpath + name + extension
	/////////////////////////////////////////////////////////////////////////////

	// Make sure project name contains no path, but has .prj extension
	strcpy(projectName, ReturnFileName(prj));
	StripFileExt(projectName);
	strcpy(projectFile, concat(projectName, ".prj"));

	// Load project structure
	prjReadErr = project.read(projectFile);
	switch (prjReadErr) {
		case 0:
		setCurrentProj(project);
		break;
		case - 1:
		cerr << "Failed to open project file " << projectFile << endl;
		exit (1);
		break;
		case - 2:
		cerr << "Project file read error: unknow line " << endl;
		exit (1);
		break;
		default:
		cerr << "Project file read error on line #" << prjReadErr << endl;
		exit (1);
		break;
	}

	// Get project with full path and extension
	strcpy(prj, concat(project.project_data_path, ".prj"));

	/////////////////////////////////////////////////////////////////////////////
	// Get variations of DEM file representation:
	//      a) fullpath + name + ".dth"
	//      b) fullpath + name
	// and make sure it exists
	/////////////////////////////////////////////////////////////////////////////

	// Make sure project name contains no path, but has .prj extension
	strcpy(demName, ReturnFileName(dem));
	StripFileExt(demName);
	build_file_name(dem, project.project_data_path, demName, ".dth");
	strcpy(fname, dem);
	StripFileExt(fname);

	if (!file_exists(dem)) {
		cerr << "\ndem " << demName << " does not exist!\n";
		cerr << "(NOTE: looking for " << dem << ")\n";
		exit(1);
	}

	/////////////////////////////////////////////////////////////////////////////
	// Generate raw binary output file names
	/////////////////////////////////////////////////////////////////////////////

	// strip extension from outcub to get the base filename
	strcpy(outcubName, ReturnFileName(outcub));
	StripFileExt(outcubName);

	strcpy(rawDEM, "SS_");
	strcat(rawDEM, outcubName);
	strcat(rawDEM, ".raw");
	if (file_exists(rawDEM))
		file_remove(rawDEM);

        // Form output raw FOM file name...add the SS_ prefix and *.raw exension after
        // the FOM prefix is added
        strcpy(FOM_outcubName, outcubName);
        upper_case(FOM_outcubName);
        char* pos=NULL;
        if (pos=strstr(FOM_outcubName,"DEM")) {
           strcpy (FOM_outcubName,outcubName); // Reset original case of filename
           strncpy (pos,"FOM",3);  // Replace "DEM" with FOM, if DEM is in filename
        }
        else if (pos=strstr(FOM_outcubName,"DTM")) {
           strcpy (FOM_outcubName,outcubName); // Reset original case of filename
           strncpy (pos,"FOM",3);  // Replace "DTM" with FOM, if DTM is in filename
        }
        else {
           strcpy (FOM_outcubName,"FOM_"); // else, add FOM_ prefix
           strcat (FOM_outcubName,outcubName);
        }
    
        strcpy(rawFOM, "SS_");  // add SS_  prefix
        strcat(rawFOM,FOM_outcubName);
        strcat(rawFOM, ".raw");  // add *.raw file exension

        if (file_exists(rawFOM))
                 file_remove(rawFOM);

	/////////////////////////////////////////////////////////////////////////////
	// Generate output isis script name
	/////////////////////////////////////////////////////////////////////////////

	strcpy(isis_script, outcubName);
	strcat(isis_script, "_dem2isis3.sh");
	if (file_exists(isis_script))
		file_remove(isis_script);

	/////////////////////////////////////////////////////////////////////////////
	//output dem2isis3 command that was issued to isis_script file
	/////////////////////////////////////////////////////////////////////////////
	sprintf(command, "#!/bin/csh\n");
	writeToScript(isis_script, command);

	if (argc == 4)
		sprintf(command, "## start_socet -single dem2isis3 %s %s %s\n", argv[1], argv[2], argv[3]);
	if (argc == 5)
		sprintf(command, "## start_socet -single dem2isis3 %s %s %s %s\n", argv[1], argv[2], argv[3], argv[4]);
	writeToScript(isis_script, command);

	/////////////////////////////////////////////////////////////////////////////
	// Decode DBDIR environment variable to determine platform
	// (Windows or Unix)
	/////////////////////////////////////////////////////////////////////////////

	unix_os = 0;
	windows_os = 0;
	str_decode_env_path ("$DBDIR", DBDIR_path); // This works for both unix and
	// windows?  At the prompt, the
	// syntax for environment variables
	// in windows is %DBDIR%

	if (strstr(DBDIR_path, "DBDIR") == 0)  // the decode worked properly
		if (strstr(DBDIR_path, ":\\") != 0) // look for :\ in path for windows
			windows_os = 1;
		else
			unix_os = 1;
	else {
		cerr << "Unable to decode DBDIR environment variable...is it missing?" << endl;
		exit (1);
	}

	/////////////////////////////////////////////////////////////////////////////
	// Set the DEM header and load DEM structure
	/////////////////////////////////////////////////////////////////////////////
	proj_ptr = getCurrentProjStruct();
	di_header = new DtmHeader(proj_ptr);
	di_header->load(fname);

	// valid DEM?
	di = new DtmGrid(di_header);
	if (demReadErr = di->openDtm(fname, FALSE, O_RDONLY, FALSE, FALSE)) {
		cerr << "DEM READ ERROR #" << demReadErr << " reading DEM file.\n";
		free(di_header);
		free(proj_ptr);
		exit( -1);
	}

	if (di_header->dtmFormat() != DTM_GRID) {
		cerr << "ERROR: input DEM is not in GRID format\n";
		free(di_header);
		free(proj_ptr);
		exit( -1);
	}
    
	/////////////////////////////////////////////////////////////////////////////
	// We have a DEM (ie, di.header->dtmFormat() == DTM_GRID)
	// Get needed parameters from SOCET DEM header that pertain
	// to both geographic and map projected projects
	/////////////////////////////////////////////////////////////////////////////
	ncols = di_header->numXPosts();
	nrows = di_header->numYPosts();
	x_realspacing = di_header->xRealSpacing();
	y_realspacing = di_header->yRealSpacing();
	ll_corner = di_header->llCorner();
	ur_corner = di_header->urCorner();

	// The LOWER_LEFT_XYZ and UPPER_RIGHT_XYZ coordinates stored in
	// the .dth file are the coordinates of the upper right and lower left
	// posts of the DEM, which is equivalent to the center of the upper right
	// and lower left pixels in an ISIS cube (units are either radians or meters,
	// depending on the map projection.)
	// From the ll and ur coordinates, compute the georeference
	// coordinate at center of the upper left pixel for ISIS
	ulcenter_Xlon = ll_corner.x;
	ulcenter_Ylat = ur_corner.y;
    
	/////////////////////////////////////////////////////////////////////////////
	// Generate a temporary raw binary DEM file
    /////////////////////////////////////////////////////////////////////////////

	ofp_DEM = fopen(rawDEM, "wb");
	if (ofp_DEM == NULL) {
		printf("\ncan't open the output raw DEM file: %s!\n", rawDEM);
		exit(1);
	}

	ofp_FOM = fopen(rawFOM, "wb");
	if (ofp_FOM == NULL) {
		printf("\ncan't open the output raw FOM file: %s!\n", rawFOM);
		exit(1);
	}

	cout << "Converting DEM and FOM to raw files...\n";
	int quarter_rows = nrows / 4;
	int half_rows = nrows / 2;
	int threequarter_rows = 3 * quarter_rows;

	//Read in data backwards (flips on mid horizontal line)
	int index_y;
	demz[0] = null;

	float *elev_buf = new float [ncols];
	char *fom_buf = new char [ncols];

	for (index_y = nrows - 1; index_y >= 0; index_y--) {

		di->getElevationBlock(0, index_y, ncols - 1, index_y, elev_buf);
		di->getFomBlock(0, index_y, ncols - 1, index_y, fom_buf);

		for (i = 0; i < ncols; i++) {

			fom[0] = *(fom_buf + i);
			fwrite(&fom[0], sizeof(char), 1, ofp_FOM);

			if (*(fom_buf + i) < 2 )
				fwrite(&demz[0], sizeof(float), 1, ofp_DEM);
			else {
				z[0] = *(elev_buf + i);
				fwrite(&z[0], sizeof(float), 1, ofp_DEM);
			}
		}
		if (index_y == threequarter_rows )
			cout << "...Conversion 25% Done\n";
		if (index_y == half_rows )
			cout << "...Conversion 50% Done\n";
		if (index_y == quarter_rows )
			cout << "...Conversion 75% Done\n";
	}

	cout << "...Conversion 100% Done\n";

	free(elev_buf);
	free(fom_buf);
	fclose(ofp_DEM);
	fclose(ofp_FOM);

	/////////////////////////////////////////////////////////////////////////////
	// Set the byte order for the DEM based on the Socet Set platform
	/////////////////////////////////////////////////////////////////////////////
	if (unix_os == 1)
		strcpy(byteOrder,"msb");
	else
		strcpy(byteOrder,"lsb");

	/////////////////////////////////////////////////////////////////////////////
	// Generate the dem2isis3 script
	/////////////////////////////////////////////////////////////////////////////

    strcpy (productType,"DEM");
    
	generate_ss2isis_script (isis_script,
	        prj,
	        productType,
	        byteOrder,
	        outcubName,
	        layout_flag,
	        nrows,
	        ncols,
	        x_realspacing,
	        y_realspacing,
	        ulcenter_Xlon,
	        ulcenter_Ylat);
} // END MAIN

