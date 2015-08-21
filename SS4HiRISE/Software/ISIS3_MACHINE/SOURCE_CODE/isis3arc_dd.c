#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**************  isis3arc_dd.c *********************
*                                                  *
*  This routine reads an ISIS3 cub label and       *
*  creates a generic 32 compatible Arc ASCII file  *
*  Registration in geographic, cellsize = deg/pix  *
   Only used for Simple Cylindrical Projections
*                                                  *
*  by Trent Hare       for USGS, flagstaff         *
*                                                  *
* Nov 2008, rewrite from isis2arc_dd.c             *
****************************************************/

// Set up ISIS NULL values
#define NULL1 0
#define NULL2 -32768
#define NULL3 -0.3402822655089E+39 /*0xFF7FFFFB*/


/* routines */
int parse_label(char *file, char *keyword, char *value);
int strstrip(char instr[], char outstr[], int position);
int stripp(char instr[], char outstr[], int position);
void swap_float_4(float *tnf4)              /* 4 byte floating point numbers */
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));  
}

main(int argc,char *argv[])
{

 int i=0, j=0, ii;
 int record_bytes;
 int skipbytes;
 int itype;
 int nsamples, nlines, nbands;
 int positiveWest = 1; /* FALSE */
 char value[60];
 char value_strip[60];
 char type[10];
 char file[256];
 char str[256];
 char str2[256];
 char byteorder[20];
 float xulcenter,yulcenter, yllcenter,cellsize, usernodata, nodata;
 float xulcorner,yulcorner, yllcorner;
 float xi, yi, x, y, yul, z;
 unsigned char z8[1];
 short int z16[1];
 float z32[1];
 FILE *ifp,*outfp,*VRTfp;

 if ((argc<3) || (argc>4))
    {
     printf ("\nUSAGE: %s infile.cub outfile.asc {NULL}\n",argv[0]);
     printf ("NULL value optional. It should be less than the minumum valid pixel value and integer. i.e. -99999\n");
     exit(1);
    } 

/***********   Grab Format ************/
 i = parse_label(argv[1], "Format", value);
 if (i == 0) {
   printf("\nFormat not found. Abort\n\n");
   exit(1);
 }
 
 if (strcmp(value,"Tile") == 0) { 
    printf("\nFormat of '%s' is not supported. Please only send BandSequential (BSQ) ISIS3 files (e.g. in any program use TO=output.cub+BandSequential). Exiting\n\n", value); 
    exit(1);
 }

 /***********   Grab samples ************/
 i = parse_label(argv[1], "Samples", value);
 if (i == 0) {
   printf("\nSamples not found. Abort\n\n");
   exit(1);
 }
 nsamples = atoi(value);
 
 
/***********   Grab lines ************/
 i = parse_label(argv[1], "Lines", value);
 if (i == 0) {
   printf("\nLines not found. Abort\n\n");
   exit(1);
 }
 nlines = atoi(value);

/***********   Grab bands ************/
 i = parse_label(argv[1], "Bands", value);
 if (i == 0) {
   printf("\nBands not found. Abort\n\n");
   exit(1);
 }
 nbands = atoi(value);

/***********   Grab StartByte data start point **********/
 i = parse_label(argv[1], "StartByte", value);
 if (i == 0) {
   printf("\nStartByte not found. Abort\n\n");
   exit(1);
 }
 skipbytes = (atoi(value) - 1);

/***********   Grab Type ************/
 i = parse_label(argv[1], "Type", value);
 if (i == 0) {
   printf("\nType not found. Abort\n\n");
   exit(1);
 }
 
 if (strcmp(value,"UnsignedByte") == 0) { 
    itype = 8;
    nodata = NULL1;
    strcpy(type, "Byte");
    if (argc > 3)
     usernodata = (float) atoi(argv[3]);
 } else if (strcmp(value,"SignedWord") == 0) {
    itype = 16;
    nodata = NULL2;
    strcpy(type, "Int16");
    if (argc > 3)
     usernodata = (float) atoi(argv[3]);
 } else if (strcmp(value,"Real") == 0) {
    itype = 32;
    nodata = NULL3;
    strcpy(type, "Float32");
    if (argc > 3)
	 usernodata = (float) atof(argv[3]);
 } else {
    printf("\nType of '%s' is not supported in ISIS3. Exiting\n\n", value); 
    exit(1);
 }
 
/***********   Grab ByteOrder ************/
 i = parse_label(argv[1], "ByteOrder", value);
 if (i == 0) {
   printf("\nByteOrder not found. Abort\n\n");
   exit(1);
 }
 if (strcmp(value,"Lsb") == 0) 
   strcpy(byteorder,"LSBFIRST");
 else
   strcpy(byteorder,"MSBFIRST");

/***********   Grab Scale ************/
 i = parse_label(argv[1], "Scale", value);
 if (i == 0) {
   printf("\nIt appears this cub has no projection. Using defaults:\n");
   printf("\nMAP_RESOLUTION not found. Setting to 1.0\n");
   cellsize = 1.0;
 } else {
   cellsize = 1 / atof(value);
 }

/***********   Grab LongitudeDirection  ************/
 i = parse_label(argv[1], "LongitudeDirection", value);
 if (i == 0) {
   printf("LongitudeDirection not found.\n");
 } else {
    if (strcmp(value,"WEST") == 0) 
       positiveWest = 0;
 }

/***********   Grab  MinimumLatitude ************/
 i = parse_label(argv[1], "MinimumLatitude", value);
 if (i == 0) {
   printf(" MinimumLatitude not found. Setting LL_Y to %f\n",nlines + 0.5);
   yllcenter = nlines + 0.5;
   yllcorner = nlines;
 } else {
   yllcorner = atof(value);
   yllcenter = ((yllcorner) + (cellsize/2));
 }

 /***********   Grab  MaximumLatitude ************/
 i = parse_label(argv[1], "MaximumLatitude", value);
 if (i == 0) {
   printf(" MaximumLatitude not found. Setting UL_Y to 0.0\n");
   yulcorner = 0.0;
 } else {
   yulcorner = atof(value);
 } 
  
/***********   Grab MinimumLongitude ************/
 i = parse_label(argv[1], "MinimumLongitude", value);
 if (i == 0) {
   printf("MinimumLongitude not found. Setting UL_X to 0.5\n\n");
   xulcorner = 0.0;
   xulcenter = 0.5;
 } else {
   xulcorner = atof(value);
   xulcenter = ((xulcorner) + (cellsize/2));
 }

if (positiveWest == 0)
  if (xulcorner > 180) 
     xulcorner =  360 - xulcorner;
  else
    xulcorner = xulcorner * -1;
else 
  if (xulcorner > 180) 
     xulcorner = xulcorner - 360;

  /* open input cub file */
  ifp=fopen(argv[1],"rb");
	if(!ifp)
	{
	   printf("Could not open '%s'.\n",argv[1]);
	   exit(1);
	}
        
  /** skip header ISIS header */
  for (i=1;i<=skipbytes;i++)
	fread(z8,sizeof(char),1,ifp);

  printf("\nLines %d /Samples %d /Bands %d /Type %d\n",nlines,nsamples,nbands,itype);
  if (nbands > 1) {
    printf("Multiple Bands. Program will add _b# for each band...\n");
    printf("Creating a .VRT file for multi-band gdal_translate conversion\n\n");

    /*gdal vrt OUTPUT */
    strcpy(file,argv[2]);
    i = stripp(file, file, 1);
    strcat(file,".vrt");
    VRTfp=fopen(file,"w");
    if(!VRTfp)
    {
        printf("Could not open '%s'.\n",file);
        exit(1);
    }
    fprintf(VRTfp,"<VRTDataset rasterXSize=\"%d\" rasterYSize=\"%d\">\n",nsamples,nlines);
    fprintf(VRTfp,"  <GeoTransform>%.8f,  %.8f,  0.0, %.8f,  0.0, %.8f</GeoTransform>\n",
                  xulcorner,cellsize,yulcorner,cellsize * -1);
  }

  if (argc > 3)
    printf("User set No Data to: %.f\n", usernodata);
     
  for (ii=1; ii <= nbands; ii++) {
	 strcpy(file,argv[2]);
    if (nbands > 1) {
        i = stripp(file, str2, 2);
        i = stripp(file, file, 1);
        if (ii < 10)
            sprintf(str,"_b0%1d",ii);
        else
            sprintf(str,"_b%2d",ii);
        strcat(file,str);
        strcat(file,".");
        strcat(file,str2);
        printf("Current band: %d Name: %s\n",ii,file);
        
        /*gdal vrt OUTPUT */
        fprintf(VRTfp,"  <VRTRasterBand dataType=\"%s\" band=\"1\">\n",type);
        if (argc > 3)
	   fprintf(VRTfp,"    <NoDataValue>%.f</NoDataValue>\n",usernodata);
	else
	   fprintf(VRTfp,"    <NoDataValue>%E</NoDataValue>\n",nodata);
        
        fprintf(VRTfp,"    <SimpleSource>\n");
        fprintf(VRTfp,"      <SourceFilename>%s</SourceFilename>\n",file);
        fprintf(VRTfp,"      <SourceBand>1</SourceBand>\n",file);
        fprintf(VRTfp,"      <SrcRect xOff=\"0\" yOff=\"0\" xSize=\"%d\" ySize=\"%d\"/>\n",nsamples,nlines);
        fprintf(VRTfp,"      <DstRect xOff=\"0\" yOff=\"0\" xSize=\"%d\" ySize=\"%d\"/>\n",nsamples,nlines);
        fprintf(VRTfp,"    </SimpleSource>\n");
        fprintf(VRTfp,"  </VRTRasterBand>\n",type);        
    }

     outfp=fopen(file,"w");
	 if(!outfp)
	 {
	   printf("Could not open '%s'.\n",file);
	   exit(1);
	 }

	 /** write out header **/
	 fprintf(outfp,"nrows %d\n",nlines);
	 fprintf(outfp,"ncols %d\n",nsamples);
	 fprintf(outfp,"xllcorner %.8f\n",xulcorner);
	 fprintf(outfp,"yllcorner %.8f\n",yllcorner);
	 fprintf(outfp,"cellsize %.8f\n",cellsize);
	 /* fprintf(outfp,"byteorder %s\n",byteorder); */
	 if (argc > 3)
	   fprintf(outfp,"NODATA_value %.f\n",usernodata);
	 else
	   fprintf(outfp,"NODATA_value %E\n", nodata);

	  if (itype == 8) {
		/**** Loop through 8 bit binary DEM ******/
		y = yllcenter;
		for (xi=1; xi <= nlines; xi++) {
		  x=xulcenter;
		  for (yi=1; yi <= nsamples; yi++) {
			fread(z8,sizeof(unsigned char),1,ifp);
            if (argc > 3) {
		  	  z = atoi(z8[0]);
			  if (z > nodata)
  			    fprintf(outfp,"%d ", z8[0]);
			  else
 			    fprintf(outfp,"%d ", usernodata);
			} else {
			  fprintf(outfp,"%d ", z8[0]);
			}
			x = x + cellsize;
			j++;
		  }
		  fprintf(outfp,"\n");
		  y = y - cellsize;
		}
	  } else if (itype == 16) {
		/**** Loop through 16 bit binary DEM ******/
		y = yllcenter;
		for (xi=1; xi <= nlines; xi++) {
		  x=xulcenter;
		  for (yi=1; yi <= nsamples; yi++) {
			fread(z16,sizeof(short int),1,ifp);
            if (argc > 3) {
		  	  z = atoi(z16[0]);
			  if (z > nodata)
  			    fprintf(outfp,"%d ", z16[0]);
			  else
 			    fprintf(outfp,"%d ", usernodata);
			} else {
			  fprintf(outfp,"%d ", z16[0]);
			}
			x = x + cellsize;
			j++;
		  }
		  fprintf(outfp,"\n");
		  y = y - cellsize;
		}
	  } else { 
		/**** Loop through 32 bit binary Image ******/
		y = yllcenter;
		for (xi=1; xi <= nlines; xi++) {
		  x=xulcenter;
		  for (yi=1; yi <= nsamples; yi++) {
			fread(z32,sizeof(float),1,ifp);
 			  sprintf(str,"%.8f",z32[0]);
			  /*swap_float_4(z32);*/                      
			  /* adding 10 because you cannot compare floats */
			  z = (float) atof(str) + 10;
			  if (z > nodata)
  			    fprintf(outfp,"%.8f ", z32[0]);
			  else {
                if (argc > 3) 
				  fprintf(outfp,"%.f ", usernodata);
			    else 
				  fprintf(outfp,"%E ", z32[0]);
			  }
			x = x + cellsize;
			j++;
		  }
		  fprintf(outfp,"\n");
		  y = y - cellsize;
		}
	  } /* if (itype == 8) */
      fclose(outfp);
  }

  ending:
  fclose(ifp);
  if (nbands > 1) {
    fprintf(VRTfp,"</VRTDataset>\n");
    fclose(VRTfp);
  }
  printf("Complete. %d pixels written in %d files(bands).\n\n",j,ii - 1);
    
}


/**************  parse_label.c    ******************
*                                                  *
*  This routine reads an ISIS cub label            *
*  by Trent Hare       for USGS, flagstaff         *
*                                                  *
****************************************************/
int parse_label(char *file, char *keyword, char *value)
{

/***************      declaration           **********/

   int ii=0,index=0,scan_value;
   char id[512];

   FILE *fp;


   fp = fopen(file,"r");
   if (fp == NULL)
     {
      printf("\ncan't open the input file %s!\n",file);
      exit(1);
   }   


/*************** Read Header file write label *******/

  fscanf (fp,"%s",id);
  while(strcmp(id,"End") != 0) {
    if (strcmp(id,keyword) == 0) {
      fscanf (fp,"%s",id);
      fscanf (fp,"%s",id);
      strcpy(value,id);
      /* printf("%s\n",value); */
      return(1);
    }
    fscanf (fp,"%s",id);
  }
  fclose(fp);
  return(0);
}


int strstrip(char instr[], char outstr[], int position)
/*************************************************************************
*_Title strstrip Trim leading and trailing blanks from string
*_Nobinding
*_Args  Type   Name             I/O  Description
*_Parm  char   instr[];          I   String to trim blanks from
*_Parm  char   outstr[];         O   Buffer to receive trimmed string
*_Parm  int    position;         I   how many non-char to skip
*_Parm  int   strstrip;        O   Returns number of characters copied
*_Desc  strstrip strips non-character like - , ( from string
*       and returns the result.
*_Keys  STRING
*       May 22 1997 Trent Hare - Changed to isistypes
*_END
*************************************************************************/
{
    char *ip;          /* Pointer to first non-blank character */
    char *op;          /* Pointer to output string */
    int i=1;             /* Number characters copied */
    int j=1;             /* Number characters copied */

/****************************************************************
  Find the first character from start of string
*****************************************************************/
    ip = instr;
    for (j=0; j<position; j++) {
      op = outstr;
      while (((*ip<'0')&&(*ip>'9')) || ((*ip<'a')&&(*ip>'z')) ||
             ((*ip<'A')&&(*ip>'Z'))) ip++;
      ip++;
      while (((*ip>='0')&&(*ip<='9')) || ((*ip>='a')&&(*ip<='z')) ||
             ((*ip>='A')&&(*ip<='Z'))) { *op++ = *ip++; i++; }
    }
    *op = '\0';
    return (i);
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
