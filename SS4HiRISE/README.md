-------------------
README_JULY2015.txt
-------------------
Disclaimer: Any use of trade, product, or firm names in this document is for 
            descriptive purposes only and does not imply endorsement by the
            U.S. Government.

Table of Contents:
   Section 1: Overview
   Section 2: Distribution Files Content
   Section 3: ISIS Machine Set-Up
   Section 4: SOCET SET Workstation Set-Up

To easily download all these files, you grab this github repository from [here] (https://github.com/USGS-Astrogeology/socet_set). Use the "Download Zip" button on right-side of page.

====================
SECTION 1: OVERIVIEW
====================

   This distribution site contains data files, in-house software, and tutorials
   used by the USGS Astrogeology Science Center for stereo processing of HiRISE
   images using ISIS3 and the WINDOWS version of SOCET SETv5.6.0 (SOCET SET
   is a registered trademark of BAE).

   The site is divided into the three zip areas for download:

    1) Data\ 
          Data.zip contains global MOLA gridded dataset as ISIS3 cubes, and
          procedures to download the MOLA Track (PEDR) data.  This is a one-
          time download -- independent of software or tutorial versions.

    2) Software\
          SOFTWARE_JULY2015.zip contains our in-house software, PERL scripts
          and the SOCET SET plugin sensor needed for HiRISE stereoprocessing
          using SOCET SET v5.6.0 under WINDOWS, and ISIS3.4.9. 

          To follow our procedures, you will also need to download FWTools 
          software (instructions below), and install the GNU FORTRAN and C
          compilers on your ISIS processing machine.

    3) Tutorials\
          PDFs and TUTORIALS.zip contains the tutorial 
          presentations describing our stereo processing procedures using 
          HiRISE images, ISIS3 and SOCET SET v5.6.0 

======================================
Section 2: DISTRIBUTION FILES CONTENTS
======================================

   The content of each zip file is listed below.  SOFTWARE_JULY2015.zip is
   organized by software package.  Pleae note that the USGS Astrogeology
   in-house SOCET SET programs have been compiled for you, with associated
   binaries placed in the bin and smplugins folders of the zip file.  The
   SOCET SET source code for our non-proprietary software is provided for
   reference in folder SS_5.6.0/SOURCE_CODE.  The source code provided for
   the ISIS machine will require compilation at your facility using the
   GNU FORTRAN or C compilers.  See the documents in TUTORIALS.zip
   for details on software use and procedures for stereo processing.

   --------
   Data
   --------
      GLOBAL_MOLA_DEMS/ --> (data used by hidata4socet.pl)
         mola_128ppd_north_simp_88lat.isis3.cub.gz
         mola_128ppd_south_simp_88lat.isis3.cub.gz

      GLOBAL_MOLA_PEDR/ --> (data used by hidata4socet.pl)
         PEDR2TAB.PRM
         wget_mola.script
         Example_mola_files.txt

   ---------------------
   Software
   ---------------------

      ISIS_MACHINE/
         PERL_SCRIPTS/
            hi_isis2geotiff.pl
            hi4socet.pl
            hidata4socet.pl
            hinoproj.pl--> (Run by hi4socet.pl)
            isis3gdal_jp2.pl
            isis3world.pl
            pedrTAB2SHP_og.pl --> (Run by hidata4socet.pl)

         SOURCE_CODE/
            isis3arc_dd.c  --> a C-program that converts an ISIS3 DEM into an
                               ascii ARC GRID, in decimal degrees.
                               (Run by hidata4socet.pl)

            pedr2tab.PCLINUX.f --> a FORTRAN program that generates an ASCII table
                                of PEDR data.  Note that this is a revised
                                version of the pedr2tab.f program originally written
                                by the MOLA team modified to compile under LINUX or
                                WINDOWS.  (Run by hidata4socet.pl)

      SS_5.6.0/
         bin/
            dem2isis3.exe
            import_frame.exe
            import_pushbroom.exe
            ortho2isis3.exe
            USGS_calcOrthoBdry.exe
            USGS_CMD.exe
            USGS_CMD.exe.BAT
            usgs_delete_image.exe
            USGS_Export_DEM_to_ISIS.exe
            USGS_Export_Orthos_to_ISIS.exe
            USGS_FOM_Masking.exe
            USGS_GPF_Leveling.exe
            USGS_Import_Frame_Shell.exe
            USGS_Import_PushBroom_Shell.exe
            USGS_Orientation_Backup.exe

         internal_dbs/     
            DTM_STRATEGY/  
               ngate_HIRISE.strategy  --> NGATE strategy file for HiRISE images
               adapt.strat.onepassAfterNGATE --> AATE strategy file used to smooth 
                                                 NGATE results
               filterpass.strat  --> AATE strategy file used with the Difference of
                                     Gaussian method described in
                                     SOCET_SET_for_HiRISE_June2009.ppt

            GEODETIC/ --> database files updated with planetary figures
               datum_name_mapping.dat
               ellipsiod.dat
               geodetic.dat
               geodetic.doc

         lib/
            smplugins/ --> USGS Astro Sensor Models
               libSS_USGSAstroLineScanner_pluginSM.dll

         SOURCE_CODE/
            calcOrthoBdry/ -->Report extent of a DEM for orthoimage generation
               calcOrthoBdry.cpp
               makefile_calcOrthoBdry.win

            dem2isis3/ -->Exports SS DEM of Mars to ISIS3 cubes
               dem2isis3.cpp
               makefile_dem2isis3.win

            export_subs/ --> subroutines shared by dem2isis3 and ortho2isis3
               export_subroutines.cpp

            import_pushbroom/ --> Imports a linescanner/pushbroom image 
                                   Currently supports cameras: HiRISE, MRO CTX,
                                   LRO-NAC and LRO-WAC 
               import_pushbroom.cpp
               makefile_import_pushbroom.win

            ortho2isis3/ -->Exports SS Orthoimages of Mars to ISIS3 cubes
               makefile_ortho2isis3.win
               ortho2isis3.cpp

   ----------------------
   Tutorials
   ----------------------
      HiRISE_StereoProcessing_Tutorial_Aug_2015.pdf
      SOCETSET_for_HiRISE_July_2015_Intro.pdf
      SOCETSET_for_HiRISE_July_2015_Training.pdf
      
===============================
SECTION 3: ISIS MACHINE SET-UP
===============================

   1) INSTALL ISIS3
   ----------------
   Instructions for installing ISIS3 are found at:
   http://isis.astrogeology.usgs.gov/documents/InstallGuide/index.html

   2) DOWNLOAD AND INSTALL FWTools
   ------------------------------
      FWTools is need for programs ogr2ogr and gdal.  ogr2ogr is run by 
      PERL scrpt pedrTAB2SHP_og.pl to convert MOLA tracks from ocentric to 
      ographic latitude coordinates.  gdal is run by PERL script 
      isis3gdal_jp2.pl to convert ISIS3 images to JPEG 2000.

      Download and install ftools/ogr2ogr from:
           Source: http://trac.osgeo.org/gdal/wiki/DownloadSource 
           Instructions: http://trac.osgeo.org/gdal/wiki/BuildingOnUnix 

      Binaries for linux: 
           http://fwtools.maptools.org/  (after uncompressing type "./install" 
           and then place the "bin_safe" directory in your $path).

      Binaries for MAC:
           Get GDAL Framework from KyngChaos
           http://www.kyngchaos.com/

   2) INSTALL GLOBAL MOLA DEMs
   ----------------------------
      Copy the contents of DATA/GLOBAL_MOLA_DEMS to a local directory, and 
      Gunzip each ISIS cube.  When gunzip'ed, each cube will be 991 MB.

   3) DOWNLOAD MOLA PEDR FILES
   ----------------------------
      The MOLA PEDR files contain the MOLA Track data.  Copy 
      DATA/GLOBAL_MOLA_PEDR to a local directory with at least 23 GB of space.  
      Then download all pedr files into that directory:  Either source the 
      provided script (DATA/GLOBAL_MOLA_PEDR/wget_mola.script), or download 
      the files from 
      http://pds-geosciences.wustl.edu/missions/mgs/pedr.html

      After files are downloaded, create a list file of all PEDR files, with 
      full path. Name the list file mola_files.txt and place it in your local 
      directory containing the PEDR files.  (For an example listing, see 
      GLOBAL_MOLA_PEDR/Example_mola_files.txt)

      NOTE: DATA/GLOBAL_MOLA_PEDR contains a PEDR parameters file named 
      PEDR2TAB.PRM.  **Do not modify or replace this file** -  it is a master file
      used by PERL script hidata4socet.pl, and differs from the one provided
      by the PDS-GEOSCIENCES web site.

   4) INSTALL ISIS PROCESSING SOURCE CODE
   --------------------------------------
      Place contents of ISIS3_MACHINE/SOURCE_CODE in a local directory, and add 
      the location of that directory to your $path

      Using the GNU Fortran compiler, compile pedr2tab.f as follows:
         gfortran -o pedr2tab.PCLINUX pedr2tab.PCLINUX.f

      Using the GNU C compiler, compile isis3arc_dd.c as follows:
         gcc -o isis3arc_dd isis3arc_dd.c

   6) INSTALL PERL SCRIPTS
   -----------------------
      Place contents of ISIS3_MACHINE/PERL_SCRIPTS in a local directory, and add 
      the location of that directory to your $path

      Change your working directory into the local directory the perl
      scripts are stored, and chmod them so they are exectuable, as follows:
         chmod uog+rwx *.pl

   7) PERL SCRIPT MODIFICATIONS
   ----------------------------
      Our perl scripts contain paths to non-ISIS programs, data files and
      other perl scripts.  For portability, all paths are treated as variables
      in the scripts, and defined in a block labeled "Location-dependent paths,
      fix for your system:".  (If this block does not exist in a perl
      script, then no modifications are needed.)  For the following perl scripts,
      the following modifications are required:
   
      hi4socet.pl:
      ------------
         modify line 21 to replace /usgs/dev/contrib/bin/ with the
         path to your local location of hinoproj.pl.

      hidata4socet.pl:
      ---------------
         modify line 23 to replace /archive/projects/SOCET_SET/DATABASES/MOLA/
         with the path to your local location of the files installed from
         DATA/GLOBAL_MOLA_DEMS in step 2, above.

         modify line 24 to replace /archive/projects/SOCET_SET/DATABASES/MOLA/
         with the path to your local location of file DATA/GLOBAL_MOLA_PEDR/PEDR2TAB.PRM
         installed in step 3, above.
 
         modify line 25 to replace /archive/projects/mola/ with the path to
         your local location of file mola_files.txt, and the PDER files created
         in step 3, above

         modify line 26 to replace /home/thare/bin/linux with the path to your
         local location of program pedrtab.PCLINUX complied in step 4 (above.)

         modify line 27 to replace /usgs/cdev/contrib/bin/ with the path to your
         local location perl script pedrTAB2SHP_og.pl installed from
         ISIS3/PERL_SCRIPTS in step 6, above.

         modify line 28 to replace /home/thare/bin/ with the path to your
         local location of program isis3arc_dd compiled in step 4, above.

=======================================
SECTION 4: SOCET SET WORKSTATION SET-UP
=======================================

   1) INSTALL SOCET SET v5.6.0
   ----------------------------

   2) MAKE A BACKUP OF adapt.strat
   -------------------------------
      Make a copy of your <ss_install_dir>/internal_dbs/DTM_STRATEGY/adapt.strat
      file and name it adapt.strat.original. (Keep this backup copy in your
      <ss_install_dir>/internal_dbs/DTM_STRATEGY/adapt.strat folder.)
  
   3) ADD USGS Astrogeology BINARY FILES
   -------------------------------------
      Copy the contents of SOFTWARE_JULY2015/SS_5.6.0/bin into your 
      <socet_set_install>/bin folder.

   4) ADD USGSASTROLineScanner SENSOR MODEL
   --------------------------------------------
      Copy the contents of SOFTWARE_JULY2015/SS_5.6.0/lib/smplugins into your
      <socet_set_install>/lib/smplugins folder.

   5) UPDATE CONFIG FILES
   -------------------------
      Replace <socet_set_install_dir>/internal_dbs/CONFIG files with the 
      contents of SOFTWARE_JULY2015/SS_5.6.0/internal_dbs/CONFIG.

   6) ADD CUSTOMIZED DTM STRATEGY FILES
   ------------------------------------
      Copy the contents of SOFTWARE_JULY2015/SS_5.6.0/internal_dbs/DTM_STRATEGY/
      to your <socet_set_install_dir>/internal_dbs/DTM_STRATEGY folder.

   7) UPDATE GEODETIC FILES
   -------------------------
      Replace <socet_set_install_dir>/internal_dbs/GEODETIC files with the 
      contents of SOFTWARE_JULY2015/SS_5.6.0/internal_dbs/GEODETIC.

