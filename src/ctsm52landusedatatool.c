/* Generates Global CTSM52 Surface Data from LUH2 format time series and MODIS and EarthStat current day reference data */
/* Author Peter Lawrence - Terrestrial Sciences Section - National Center for Atmospheric Research */
/* Email:  lawrence@ucar.edu */ 
/* Web:    https://www.cgd.ucar.edu/staff/lawrence */
/* GitHub: https://github.com/lawrencepj1 */
/* Phone:  +1 303-4971727 (work) +1 303-9566932 (mobile) */

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define MAXCTSMPIX 1440
#define MAXCTSMLIN 720

#define CTSMPIXSIZE 0.25
#define CTSMLLX -180.0
#define CTSMLLY -90.0

#define PI 4.0*atan(1.0)
#define EarthCir 40075.017

#define MAXPFT 15
#define MAXCFTRAW 32
#define MAXCFT 64

#define firsttreepft 1
#define lasttreepft  8

#define RANK_natpft 1
#define RANK_cft 1
#define RANK_EDGEN 0
#define RANK_EDGEE 0
#define RANK_EDGES 0
#define RANK_EDGEW 0
#define RANK_LAT 1
#define RANK_LATIXY 2
#define RANK_LON 1
#define RANK_LONGXY 2
#define RANK_LANDMASK 2
#define RANK_LANDFRAC 2
#define RANK_AREA 2
#define RANK_PCT_GLACIER 2
#define RANK_PCT_LAKE 2
#define RANK_PCT_WETLAND 2
#define RANK_PCT_URBAN 2
#define RANK_PCT_NATVEG 2
#define RANK_PCT_CROP 2
#define RANK_PCT_NAT_PFT 3
#define RANK_PCT_CFT 3
#define RANK_FERTNITRO_CFT 3
#define RANK_HARVEST_VH1 2
#define RANK_HARVEST_VH2 2
#define RANK_HARVEST_SH1 2
#define RANK_HARVEST_SH2 2
#define RANK_HARVEST_SH3 2
#define RANK_GRAZING 2
#define RANK_UNREPRESENTED_PFT_LULCC 3
#define RANK_UNREPRESENTED_CFT_LULCC 3

long MAXOUTPIX = MAXCTSMPIX;
long MAXOUTLIN = MAXCTSMLIN;
long OUTLONOFFSET = 0;
long OUTLATOFFSET = 0;
float OUTPIXSIZE = CTSMPIXSIZE;
float OUTLLX = CTSMLLX;
float OUTLLY = CTSMLLY;

long OUTDATASIZE = sizeof(float) * MAXCTSMPIX * MAXCTSMLIN;
long OUTDBLDATASIZE = sizeof(double) * MAXCTSMPIX * MAXCTSMLIN;

/* Namelist Variables */

char authorname[1024];
char regionfilename[1024];
char outputdir[1024];
char outputseries[1024];
char timestamp[1024];
int refyear;
int startyear;
int endyear;
int ctsmcurrentsurfstartyear;
int ctsmcurrentsurfendyear;
int ctsmcurrentsurfreadyear = -99999;
char ctsmcurrentsurfdb[1024];
int ctsmLUHforeststartyear;
int ctsmLUHforestendyear;
int ctsmLUHforestreadyear = -99999;
char ctsmLUHforestdb[1024];
int ctsmLUHpasturestartyear;
int ctsmLUHpastureendyear;
int ctsmLUHpasturereadyear = -99999;
char ctsmLUHpasturedb[1024];
int ctsmLUHotherstartyear;
int ctsmLUHotherendyear;
int ctsmLUHotherreadyear = -99999;
char ctsmLUHotherdb[1024];
int ctsmLUHc3annstartyear;
int ctsmLUHc3annendyear;
int ctsmLUHc3annreadyear = -99999;
char ctsmLUHc3anndb[1024];
int ctsmLUHc4annstartyear;
int ctsmLUHc4annendyear;
int ctsmLUHc4annreadyear = -99999;
char ctsmLUHc4anndb[1024];
int ctsmLUHc3perstartyear;
int ctsmLUHc3perendyear;
int ctsmLUHc3perreadyear = -99999;
char ctsmLUHc3perdb[1024];
int ctsmLUHc4perstartyear;
int ctsmLUHc4perendyear;
int ctsmLUHc4perreadyear = -99999;
char ctsmLUHc4perdb[1024];
int ctsmLUHc3nfxstartyear;
int ctsmLUHc3nfxendyear;
int ctsmLUHc3nfxreadyear = -99999;
char ctsmLUHc3nfxdb[1024];
int refstatesstartyear;
int refstatesendyear;
int refstatesreadyear = -99999;
char refstatesdb[1024];
int luhstatesstartyear;
int luhstatesendyear;
int luhcurrentstatesreadyear = -99999;
int luhprevstatesreadyear = -99999;
char luhstatesdb[1024];
int luhmanagementstartyear;
int luhmanagementendyear;
int luhwoodharvestreadyear = -99999;
int luhcropmanagementreadyear = -99999;
char luhmanagementdb[1024];
int luhtransitionsstartyear;
int luhtransitionsendyear;
int luhsecdfunrepreadyear = -99999;
int luhsecdnunrepreadyear = -99999;
char luhtransitionsdb[1024];
char pftparamfile[1024];
char cftrawparamfile[1024];
char cftparamfile[1024];
int flipLUHgrids;
int includeOcean;

char PFTluhtype[MAXPFT][256];
char CFTRAWluhtype[MAXCFTRAW][256];
char CFTluhtype[MAXCFT][256];

float *tempGrid;
float *tempflipGrid;
float *tempextrapGrid;
float *tempoutGrid;
float *secdfCROPINGrid;
float *secdfOTHERINGrid;
float *secdfOTHEROUTGrid;
float *secdnCROPINGrid;
float *secdnOTHERINGrid;
float *secdnOTHEROUTGrid;

int *innatpft;
int *incft;
float inEDGEN;
float inEDGEE;
float inEDGES;
float inEDGEW;
float *inLAT;
float *inLATIXY;
float *inLON;
float *inLONGXY;

float *inLANDMASKGrid;
float *inLANDFRACGrid;
float *inAREAGrid;
float *inPCTGLACIERGrid;
float *inPCTLAKEGrid;
float *inPCTWETLANDGrid;
float *inPCTURBANGrid;
float *inPCTNATVEGGrid;
float *inPCTCROPGrid;
float *inCURRENTPCTPFTGrid[MAXPFT];
float *inCURRENTPCTCFTGrid[MAXCFT];

float *inFORESTPCTPFTGrid[MAXPFT];
float *inPASTUREPCTPFTGrid[MAXPFT];
float *inOTHERPCTPFTGrid[MAXPFT];

float *inC3ANNPCTCFTGrid[MAXCFTRAW];
float *inC4ANNPCTCFTGrid[MAXCFTRAW];
float *inC3PERPCTCFTGrid[MAXCFTRAW];
float *inC4PERPCTCFTGrid[MAXCFTRAW];
float *inC3NFXPCTCFTGrid[MAXCFTRAW];

float *inBASEPRIMFGrid;
float *inBASEPRIMNGrid;
float *inBASESECDFGrid;
float *inBASESECDNGrid;
float *inBASEPASTRGrid;
float *inBASERANGEGrid;
float *inBASEC3ANNGrid;
float *inBASEC4ANNGrid;
float *inBASEC3PERGrid;
float *inBASEC4PERGrid;
float *inBASEC3NFXGrid;
float *inBASEURBANGrid;

float *inCURRPRIMFGrid;
float *inCURRPRIMNGrid;
float *inCURRSECDFGrid;
float *inCURRSECDNGrid;
float *inCURRPASTRGrid;
float *inCURRRANGEGrid;
float *inCURRCropGrid;
float *inCURRC3ANNGrid;
float *inCURRC4ANNGrid;
float *inCURRC3PERGrid;
float *inCURRC4PERGrid;
float *inCURRC3NFXGrid;
float *inCURRURBANGrid;

float *inPREVDELTASECDFGrid;
float *inPREVDELTASECDNGrid;
float *inPREVDELTAPASTRGrid;
float *inPREVDELTARANGEGrid;
float *inPREVDELTAC3ANNGrid;
float *inPREVDELTAC4ANNGrid;
float *inPREVDELTAC3PERGrid;
float *inPREVDELTAC4PERGrid;
float *inPREVDELTAC3NFXGrid;

float *inHARVESTVH1Grid;
float *inHARVESTVH2Grid;
float *inHARVESTSH1Grid;
float *inHARVESTSH2Grid;
float *inHARVESTSH3Grid;

float *inBIOHVH1Grid;
float *inBIOHVH2Grid;
float *inBIOHSH1Grid;
float *inBIOHSH2Grid;
float *inBIOHSH3Grid;

float *inUNREPSECDFGrid;
float *inUNREPSECDNGrid;
float *inUNREPPASTRGrid;
float *inUNREPRANGEGrid;
float *inUNREPC3ANNGrid;
float *inUNREPC4ANNGrid;
float *inUNREPC3PERGrid;
float *inUNREPC4PERGrid;
float *inUNREPC3NFXGrid;

float *inFERTC3ANNGrid;
float *inFERTC4ANNGrid;
float *inFERTC3PERGrid;
float *inFERTC4PERGrid;
float *inFERTC3NFXGrid;

float *inIRRIGC3ANNGrid;
float *inIRRIGC4ANNGrid;
float *inIRRIGC3PERGrid;
float *inIRRIGC4PERGrid;
float *inIRRIGC3NFXGrid;

float *inBASEFORESTTOTALGrid;
float *inBASENONFORESTTOTALGrid;
float *inBASECROPTOTALGrid;
float *inBASEURBANTOTALGrid;
float *inBASEMISSINGGrid;
float *inBASEOTHERGrid;
float *inBASENATVEGGrid;

float *inCURRFORESTTOTALGrid;
float *inCURRNONFORESTTOTALGrid;
float *inCURRCROPTOTALGrid;
float *inPREVCROPTOTALGrid;
float *inCURRURBANTOTALGrid;
float *inCURRMISSINGGrid;
float *inCURROTHERGrid;
float *inCURRNATVEGGrid;

float *inUNREPFORESTGrid;
float *inUNREPOTHERGrid;

float *outLANDMASKGrid;
float *outPCTURBANGrid;
float *outPCTNATVEGGrid;
float *outPCTCROPGrid;
float *outPCTPFTGrid[MAXPFT];
float *outPCTCFTGrid[MAXCFT];

float *outFERTNITROGrid[MAXCFT];

float *outUNREPPFTGrid[MAXPFT];
float *outUNREPCFTGrid[MAXCFT];

float *outHARVESTVH1Grid;
float *outHARVESTVH2Grid;
float *outHARVESTSH1Grid;
float *outHARVESTSH2Grid;
float *outHARVESTSH3Grid;

float *outBIOHVH1Grid;
float *outBIOHVH2Grid;
float *outBIOHSH1Grid;
float *outBIOHSH2Grid;
float *outBIOHSH3Grid;

float *outRBIOHVH1Grid;
float *outRBIOHVH2Grid;
float *outRBIOHSH1Grid;
float *outRBIOHSH2Grid;
float *outRBIOHSH3Grid;

float *outRBIOHTreePFTAreaGrid;
float *outRBIOHTreePFTWeightedAreaGrid;
float *outRBIOHPFTAreaGrid;
float *outRBIOHORIGTOTALGrid;
float *outRBIOHCURRTOTALGrid;

double *outLANDFRACdblGrid;
double *outAREAdblGrid;
double *outPCTGLACIERdblGrid;
double *outPCTLAKEdblGrid;
double *outPCTWETLANDdblGrid;

double *outPCTURBANdblGrid;
double *outPCTNATVEGdblGrid;
double *outPCTCROPdblGrid;
double *outPCTPFTdblGrid[MAXPFT];
double *outPCTCFTdblGrid[MAXCFT];

double *outFERTNITROdblGrid[MAXCFT];

double *outUNREPPFTdblGrid[MAXPFT];
double *outUNREPCFTdblGrid[MAXCFT];

double *outRBIOHVH1dblGrid;
double *outRBIOHVH2dblGrid;
double *outRBIOHSH1dblGrid;
double *outRBIOHSH2dblGrid;
double *outRBIOHSH3dblGrid;

/* Out Surface Data NetCDF variables */
int  stat;  /* return status */
int  ncid;  /* netCDF id */

/* dimension ids */
int natpft_dim;
int cft_dim;
int lon_dim;
int lat_dim;
int nchar_dim;

/* dimension lengths */
size_t natpft_len = 15;
size_t cft_len = 64;
size_t lon_len = MAXCTSMPIX;
size_t lat_len = MAXCTSMLIN;
size_t nchar_len = 128;

/* variable ids */
int natpft_id;
int cft_id;
int EDGEN_id;
int EDGEE_id;
int EDGES_id;
int EDGEW_id;
int LAT_id;
int LATIXY_id;
int LON_id;
int LONGXY_id;
int LANDMASK_id;
int LANDFRAC_id;
int AREA_id;
int PCT_GLACIER_id;
int PCT_LAKE_id;
int PCT_WETLAND_id;
int PCT_URBAN_id;
int PCT_NATVEG_id;
int PCT_CROP_id;
int PCT_NAT_PFT_id;
int PCT_CFT_id;
int FERTNITRO_CFT_id;
int HARVEST_VH1_id;
int HARVEST_VH2_id;
int HARVEST_SH1_id;
int HARVEST_SH2_id;
int HARVEST_SH3_id;
int GRAZING_id;
int UNREPRESENTED_PFT_LULCC_id;
int UNREPRESENTED_CFT_LULCC_id;

/* variable shapes */
int natpft_dims[RANK_natpft];
int cft_dims[RANK_cft];
int LAT_dims[RANK_LAT];
int LATIXY_dims[RANK_LATIXY];
int LON_dims[RANK_LON];
int LONGXY_dims[RANK_LONGXY];
int LANDMASK_dims[RANK_LANDMASK];
int LANDFRAC_dims[RANK_LANDFRAC];
int AREA_dims[RANK_AREA];
int PCT_GLACIER_dims[RANK_PCT_GLACIER];
int PCT_LAKE_dims[RANK_PCT_LAKE];
int PCT_WETLAND_dims[RANK_PCT_WETLAND];
int PCT_URBAN_dims[RANK_PCT_URBAN];
int PCT_NATVEG_dims[RANK_PCT_NATVEG];
int PCT_CROP_dims[RANK_PCT_CROP];
int PCT_NAT_PFT_dims[RANK_PCT_NAT_PFT];
int PCT_CFT_dims[RANK_PCT_CFT];
int FERTNITRO_CFT_dims[RANK_FERTNITRO_CFT];
int HARVEST_VH1_dims[RANK_HARVEST_VH1];
int HARVEST_VH2_dims[RANK_HARVEST_VH2];
int HARVEST_SH1_dims[RANK_HARVEST_SH1];
int HARVEST_SH2_dims[RANK_HARVEST_SH2];
int HARVEST_SH3_dims[RANK_HARVEST_SH3];
int GRAZING_dims[RANK_GRAZING];
int UNREPRESENTED_PFT_LULCC_dims[RANK_UNREPRESENTED_PFT_LULCC];
int UNREPRESENTED_CFT_LULCC_dims[RANK_UNREPRESENTED_CFT_LULCC];

int readnamelist(char *namelist) {

  FILE *namelistfile;
  char templine[1024];
  char fieldname[256];
  char *token;
  int tokencount;

  printf("Reading Namelist: %s\n",namelist);
  namelistfile = fopen(namelist,"r");
  
  sprintf(authorname,"");
  fgets(templine,sizeof(templine),namelistfile);
  token = strtok(templine," ");
  tokencount = 1;
  while( token != NULL ) {
     if (tokencount == 2) {
         sprintf(authorname,"%s",token);
     }
     if (tokencount > 2 && strcmp(token,"\n") != 0) {
        sprintf(authorname,"%s %s",authorname,token);
     }
     token = strtok(NULL," ");
     tokencount++;
  }
  fscanf(namelistfile,"%s %s",fieldname,regionfilename);
  fscanf(namelistfile,"%s %s",fieldname,outputdir);
  fscanf(namelistfile,"%s %s",fieldname,outputseries);
  fscanf(namelistfile,"%s %s",fieldname,timestamp);
  fscanf(namelistfile,"%s %d",fieldname,&refyear);
  fscanf(namelistfile,"%s %d",fieldname,&startyear);
  fscanf(namelistfile,"%s %d",fieldname,&endyear);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmcurrentsurfstartyear,&ctsmcurrentsurfendyear,ctsmcurrentsurfdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHforeststartyear,&ctsmLUHforestendyear,ctsmLUHforestdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHpasturestartyear,&ctsmLUHpastureendyear,ctsmLUHpasturedb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHotherstartyear,&ctsmLUHotherendyear,ctsmLUHotherdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHc3annstartyear,&ctsmLUHc3annendyear,ctsmLUHc3anndb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHc4annstartyear,&ctsmLUHc4annendyear,ctsmLUHc4anndb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHc3perstartyear,&ctsmLUHc3perendyear,ctsmLUHc3perdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHc4perstartyear,&ctsmLUHc4perendyear,ctsmLUHc4perdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&ctsmLUHc3nfxstartyear,&ctsmLUHc3nfxendyear,ctsmLUHc3nfxdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&refstatesstartyear,&refstatesendyear,refstatesdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&luhstatesstartyear,&luhstatesendyear,luhstatesdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&luhmanagementstartyear,&luhmanagementendyear,luhmanagementdb);
  fscanf(namelistfile,"%s %d %d %s",fieldname,&luhtransitionsstartyear,&luhtransitionsendyear,luhtransitionsdb);
  fscanf(namelistfile,"%s %s",fieldname,pftparamfile);
  fscanf(namelistfile,"%s %s",fieldname,cftrawparamfile);
  fscanf(namelistfile,"%s %s",fieldname,cftparamfile);  
  fscanf(namelistfile,"%s %d",fieldname,&flipLUHgrids);
  fscanf(namelistfile,"%s %d",fieldname,&includeOcean);

  return 0;

}

int setregionoptions() {

  FILE *pftparamfile;
  float lllon, lllat, urlon, urlat, pixsize;
  
  printf("Reading %s\n",regionfilename);
  pftparamfile = fopen(regionfilename,"r");

  fscanf(pftparamfile,"%f",&lllon);
  fscanf(pftparamfile,"%f",&lllat);
  fscanf(pftparamfile,"%f",&urlon);
  fscanf(pftparamfile,"%f",&urlat);
  fscanf(pftparamfile,"%f",&pixsize);
  
  OUTPIXSIZE = pixsize;
  MAXOUTPIX = (long) (urlon - lllon) / OUTPIXSIZE;
  MAXOUTLIN = (long) (urlat - lllat) / OUTPIXSIZE;

  lon_len = MAXOUTPIX;
  lat_len = MAXOUTLIN;
  
  OUTLLX = lllon;
  OUTLLY = lllat;
  
  OUTLATOFFSET = (long) (90.0 - urlat) / OUTPIXSIZE;
  OUTLONOFFSET = (long) (lllon + 180.0) / OUTPIXSIZE;

  OUTDATASIZE = MAXOUTPIX * MAXOUTLIN * sizeof(float);
  OUTDBLDATASIZE = MAXOUTPIX * MAXOUTLIN * sizeof(double);
  
  return 0;

}

int createallgrids() {

  int pftid, cftid;

  tempGrid = (float *) malloc(OUTDATASIZE);
  tempoutGrid = (float *) malloc(OUTDATASIZE);
  tempflipGrid = (float *) malloc(OUTDATASIZE);
  tempextrapGrid = (float *) malloc(OUTDATASIZE);
  secdfCROPINGrid = (float *) malloc(OUTDATASIZE);
  secdfOTHERINGrid = (float *) malloc(OUTDATASIZE);
  secdfOTHEROUTGrid = (float *) malloc(OUTDATASIZE);
  secdnCROPINGrid = (float *) malloc(OUTDATASIZE);
  secdnOTHERINGrid = (float *) malloc(OUTDATASIZE);
  secdnOTHEROUTGrid = (float *) malloc(OUTDATASIZE);

  innatpft = (int *) malloc(MAXPFT * sizeof(int));
  incft = (int *) malloc(MAXCFT * sizeof(int));
  inLAT = (float *) malloc(MAXOUTLIN * sizeof(float));
  inLATIXY = (float *) malloc(OUTDATASIZE);
  inLON = (float *) malloc(MAXOUTPIX * sizeof(float));
  inLONGXY = (float *) malloc(OUTDATASIZE);

  inLANDMASKGrid = (float *) malloc(OUTDATASIZE);
  inLANDFRACGrid = (float *) malloc(OUTDATASIZE);
  inAREAGrid = (float *) malloc(OUTDATASIZE);
  inPCTGLACIERGrid = (float *) malloc(OUTDATASIZE);
  inPCTLAKEGrid = (float *) malloc(OUTDATASIZE);
  inPCTWETLANDGrid = (float *) malloc(OUTDATASIZE);
  inPCTURBANGrid = (float *) malloc(OUTDATASIZE);
  inPCTNATVEGGrid = (float *) malloc(OUTDATASIZE);
  inPCTCROPGrid = (float *) malloc(OUTDATASIZE);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      inCURRENTPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {  
      inCURRENTPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      inFORESTPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
      inPASTUREPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
      inOTHERPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      inC3ANNPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC4ANNPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC3PERPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC4PERPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC3NFXPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }

  inBASEPRIMFGrid = (float *) malloc(OUTDATASIZE);
  inBASEPRIMNGrid = (float *) malloc(OUTDATASIZE);
  inBASESECDFGrid = (float *) malloc(OUTDATASIZE);
  inBASESECDNGrid = (float *) malloc(OUTDATASIZE);
  inBASEPASTRGrid = (float *) malloc(OUTDATASIZE);
  inBASERANGEGrid = (float *) malloc(OUTDATASIZE);
  inBASEC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inBASEC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inBASEC3PERGrid = (float *) malloc(OUTDATASIZE);
  inBASEC4PERGrid = (float *) malloc(OUTDATASIZE);
  inBASEC3NFXGrid = (float *) malloc(OUTDATASIZE);
  inBASEURBANGrid = (float *) malloc(OUTDATASIZE);

  inCURRPRIMFGrid = (float *) malloc(OUTDATASIZE);
  inCURRPRIMNGrid = (float *) malloc(OUTDATASIZE);
  inCURRSECDFGrid = (float *) malloc(OUTDATASIZE);
  inCURRSECDNGrid = (float *) malloc(OUTDATASIZE);
  inCURRPASTRGrid = (float *) malloc(OUTDATASIZE);
  inCURRRANGEGrid = (float *) malloc(OUTDATASIZE);
  inCURRCropGrid = (float *) malloc(OUTDATASIZE);
  inCURRC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inCURRC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inCURRC3PERGrid = (float *) malloc(OUTDATASIZE);
  inCURRC4PERGrid = (float *) malloc(OUTDATASIZE);
  inCURRC3NFXGrid = (float *) malloc(OUTDATASIZE);
  inCURRURBANGrid = (float *) malloc(OUTDATASIZE);

  inPREVDELTASECDFGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTASECDNGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTAPASTRGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTARANGEGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTAC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTAC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTAC3PERGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTAC4PERGrid = (float *) malloc(OUTDATASIZE);
  inPREVDELTAC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inHARVESTVH1Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTVH2Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTSH1Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTSH2Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTSH3Grid = (float *) malloc(OUTDATASIZE);

  inBIOHVH1Grid = (float *) malloc(OUTDATASIZE);
  inBIOHVH2Grid = (float *) malloc(OUTDATASIZE);
  inBIOHSH1Grid = (float *) malloc(OUTDATASIZE);
  inBIOHSH2Grid = (float *) malloc(OUTDATASIZE);
  inBIOHSH3Grid = (float *) malloc(OUTDATASIZE);

  inUNREPSECDFGrid = (float *) malloc(OUTDATASIZE);
  inUNREPSECDNGrid = (float *) malloc(OUTDATASIZE);
  inUNREPPASTRGrid = (float *) malloc(OUTDATASIZE);
  inUNREPRANGEGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC3PERGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC4PERGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inFERTC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inFERTC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inFERTC3PERGrid = (float *) malloc(OUTDATASIZE);
  inFERTC4PERGrid = (float *) malloc(OUTDATASIZE);
  inFERTC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inIRRIGC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC3PERGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC4PERGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inBASEFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASENONFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASECROPTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASEURBANTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASEMISSINGGrid = (float *) malloc(OUTDATASIZE);
  inBASEOTHERGrid = (float *) malloc(OUTDATASIZE);
  inBASENATVEGGrid = (float *) malloc(OUTDATASIZE);

  inCURRFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRNONFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRCROPTOTALGrid = (float *) malloc(OUTDATASIZE);
  inPREVCROPTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRURBANTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRMISSINGGrid = (float *) malloc(OUTDATASIZE);
  inCURROTHERGrid = (float *) malloc(OUTDATASIZE);
  inCURRNATVEGGrid = (float *) malloc(OUTDATASIZE);

  inUNREPFORESTGrid = (float *) malloc(OUTDATASIZE);
  inUNREPOTHERGrid = (float *) malloc(OUTDATASIZE);

  outLANDMASKGrid = (float *) malloc(OUTDATASIZE);
  outPCTURBANGrid = (float *) malloc(OUTDATASIZE);
  outPCTNATVEGGrid = (float *) malloc(OUTDATASIZE);
  outPCTCROPGrid = (float *) malloc(OUTDATASIZE);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outUNREPPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outUNREPCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outFERTNITROGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  outHARVESTVH1Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTVH2Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTSH1Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTSH2Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTSH3Grid = (float *) malloc(OUTDATASIZE);

  outBIOHVH1Grid = (float *) malloc(OUTDATASIZE);
  outBIOHVH2Grid = (float *) malloc(OUTDATASIZE);
  outBIOHSH1Grid = (float *) malloc(OUTDATASIZE);
  outBIOHSH2Grid = (float *) malloc(OUTDATASIZE);
  outBIOHSH3Grid = (float *) malloc(OUTDATASIZE);

  outRBIOHVH1Grid = (float *) malloc(OUTDATASIZE);
  outRBIOHVH2Grid = (float *) malloc(OUTDATASIZE);
  outRBIOHSH1Grid = (float *) malloc(OUTDATASIZE);
  outRBIOHSH2Grid = (float *) malloc(OUTDATASIZE);
  outRBIOHSH3Grid = (float *) malloc(OUTDATASIZE);

  outRBIOHTreePFTAreaGrid = (float *) malloc(OUTDATASIZE);
  outRBIOHTreePFTWeightedAreaGrid = (float *) malloc(OUTDATASIZE);
  outRBIOHPFTAreaGrid = (float *) malloc(OUTDATASIZE);
  outRBIOHORIGTOTALGrid = (float *) malloc(OUTDATASIZE);
  outRBIOHCURRTOTALGrid = (float *) malloc(OUTDATASIZE);

  outLANDFRACdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outAREAdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTGLACIERdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTLAKEdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTWETLANDdblGrid = (double *) malloc(OUTDBLDATASIZE);

  outPCTURBANdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTNATVEGdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTCROPdblGrid = (double *) malloc(OUTDBLDATASIZE);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outPCTPFTdblGrid[pftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outPCTCFTdblGrid[cftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outUNREPPFTdblGrid[pftid] = (double *) malloc(OUTDBLDATASIZE);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outUNREPCFTdblGrid[cftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outFERTNITROdblGrid[cftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  outRBIOHVH1dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outRBIOHVH2dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outRBIOHSH1dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outRBIOHSH2dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outRBIOHSH3dblGrid = (double *) malloc(OUTDBLDATASIZE);

  return 0;

}


int
readpftparamfile() {

  FILE *pftparaminfile;
  int inpft, inpftid;
  char inPFTluhtype[256], inPFTname[256];
  
  printf("Reading %s\n",pftparamfile);
  pftparaminfile = fopen(pftparamfile,"r");

  for (inpft = 0; inpft < MAXPFT; inpft++) {
      fscanf(pftparaminfile,"%d%s%s",&inpftid,inPFTluhtype,inPFTname);
      sprintf(PFTluhtype[inpft],"%s",inPFTluhtype);
  }  
  
  return 0;

}

int
readcftrawparamfile() {

  FILE *cftrawparaminfile;
  int incftraw, incftrawid;
  char inCFTRAWluhtype[256], inCFTRAWname[256];
  
  printf("Reading %s\n",cftrawparamfile);
  cftrawparaminfile = fopen(cftrawparamfile,"r");

  for (incftraw = 0; incftraw < MAXCFTRAW; incftraw++) {
      fscanf(cftrawparaminfile,"%d%s%s",&incftrawid,inCFTRAWluhtype,inCFTRAWname);
      sprintf(CFTRAWluhtype[incftraw],"%s",inCFTRAWluhtype);
  }  
  
  return 0;

}

int
readcftparamfile() {

  FILE *cftparaminfile;
  int incft, incftid;
  char inCFTluhtype[256], inCFTname[256];
  
  printf("Reading %s\n",cftparamfile);
  cftparaminfile = fopen(cftparamfile,"r");

  for (incft = 0; incft < MAXCFT; incft++) {
      fscanf(cftparaminfile,"%d%s%s",&incftid,inCFTluhtype,inCFTname);
      sprintf(CFTluhtype[incft],"%s",inCFTluhtype);
  }  
  
  return 0;

}


int initializeGrids() {

  long ctsmlin, ctsmpix;
  long pftid, cftid;
  float ctsmlat, ctsmlatdistance, ctsmlondistance, ctsmarea;
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          outLANDMASKGrid[MAXOUTPIX * MAXOUTLIN] = 0.0;
          outPCTURBANGrid[MAXOUTPIX * MAXOUTLIN] = 0.0;
          outPCTNATVEGGrid[MAXOUTPIX * MAXOUTLIN] = 0.0;
          outPCTCROPGrid[MAXOUTPIX * MAXOUTLIN] = 0.0;
          for (pftid = 0; pftid < MAXPFT; pftid++) {
              outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outUNREPPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }

          outHARVESTVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outHARVESTVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outHARVESTSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outHARVESTSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outHARVESTSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;

          outBIOHVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outBIOHVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outBIOHSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outBIOHSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          outBIOHSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;

          for (cftid = 0; cftid < MAXCFT; cftid++) {
              outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outFERTNITROGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outUNREPCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
      }
  }

  return 0;
                
}

void
check_err(const int stat, const int line, const char *file) {

    if (stat != NC_NOERR) {
        (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
        exit(1);
    }
}

int
openncinputfile(char *netcdffilename) {

    printf("Opening NetCDF File: %s\n",netcdffilename); 
    stat = nc_open(netcdffilename, NC_NOWRITE, &ncid);
    check_err(stat,__LINE__,__FILE__);

    return 0;

}

int
openncoutputfile(char *netcdffilename) {

    printf("Opening NetCDF File: %s\n",netcdffilename); 
    stat = nc_open(netcdffilename, NC_WRITE, &ncid);
    check_err(stat,__LINE__,__FILE__);

    return 0;

}

int
createncoutputfile(char *netcdffilename) {

    printf("Creating NetCDF File: %s\n",netcdffilename); 

    /* enter define mode */
    stat = nc_create(netcdffilename, NC_CLOBBER|NC_CDF5, &ncid);
    check_err(stat,__LINE__,__FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "natpft", natpft_len, &natpft_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "cft", cft_len, &cft_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "lon", lon_len, &lon_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "lat", lat_len, &lat_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "nchar", nchar_len, &nchar_dim);
    check_err(stat,__LINE__,__FILE__);

    /* define variables */

    natpft_dims[0] = natpft_dim;
    stat = nc_def_var(ncid, "natpft", NC_INT, RANK_natpft, natpft_dims, &natpft_id);
    check_err(stat,__LINE__,__FILE__);

    cft_dims[0] = cft_dim;
    stat = nc_def_var(ncid, "cft", NC_INT, RANK_cft, cft_dims, &cft_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGEN", NC_FLOAT, RANK_EDGEN, 0, &EDGEN_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGEE", NC_FLOAT, RANK_EDGEE, 0, &EDGEE_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGES", NC_FLOAT, RANK_EDGES, 0, &EDGES_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGEW", NC_FLOAT, RANK_EDGEW, 0, &EDGEW_id);
    check_err(stat,__LINE__,__FILE__);

    LAT_dims[0] = lat_dim;
    stat = nc_def_var(ncid, "LAT", NC_FLOAT, RANK_LAT, LAT_dims, &LAT_id);
    check_err(stat,__LINE__,__FILE__);

    LATIXY_dims[0] = lat_dim;
    LATIXY_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LATIXY", NC_FLOAT, RANK_LATIXY, LATIXY_dims, &LATIXY_id);
    check_err(stat,__LINE__,__FILE__);

    LON_dims[0] = lon_dim;
    stat = nc_def_var(ncid, "LON", NC_FLOAT, RANK_LON, LON_dims, &LON_id);
    check_err(stat,__LINE__,__FILE__);

    LONGXY_dims[0] = lat_dim;
    LONGXY_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LONGXY", NC_FLOAT, RANK_LONGXY, LONGXY_dims, &LONGXY_id);
    check_err(stat,__LINE__,__FILE__);

    LANDMASK_dims[0] = lat_dim;
    LANDMASK_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LANDMASK", NC_FLOAT, RANK_LANDMASK, LANDMASK_dims, &LANDMASK_id);
    check_err(stat,__LINE__,__FILE__);

    LANDFRAC_dims[0] = lat_dim;
    LANDFRAC_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LANDFRAC", NC_DOUBLE, RANK_LANDFRAC, LANDFRAC_dims, &LANDFRAC_id);
    check_err(stat,__LINE__,__FILE__);

    AREA_dims[0] = lat_dim;
    AREA_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "AREA", NC_DOUBLE, RANK_AREA, AREA_dims, &AREA_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_GLACIER_dims[0] = lat_dim;
    PCT_GLACIER_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_GLACIER", NC_DOUBLE, RANK_PCT_GLACIER, PCT_GLACIER_dims, &PCT_GLACIER_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_LAKE_dims[0] = lat_dim;
    PCT_LAKE_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_LAKE", NC_DOUBLE, RANK_PCT_LAKE, PCT_LAKE_dims, &PCT_LAKE_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_WETLAND_dims[0] = lat_dim;
    PCT_WETLAND_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_WETLAND", NC_DOUBLE, RANK_PCT_WETLAND, PCT_WETLAND_dims, &PCT_WETLAND_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_URBAN_dims[0] = lat_dim;
    PCT_URBAN_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_URBAN", NC_DOUBLE, RANK_PCT_URBAN, PCT_URBAN_dims, &PCT_URBAN_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_NATVEG_dims[0] = lat_dim;
    PCT_NATVEG_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_NATVEG", NC_DOUBLE, RANK_PCT_NATVEG, PCT_NATVEG_dims, &PCT_NATVEG_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_CROP_dims[0] = lat_dim;
    PCT_CROP_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_CROP", NC_DOUBLE, RANK_PCT_CROP, PCT_CROP_dims, &PCT_CROP_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_NAT_PFT_dims[0] = natpft_dim;
    PCT_NAT_PFT_dims[1] = lat_dim;
    PCT_NAT_PFT_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "PCT_NAT_PFT", NC_DOUBLE, RANK_PCT_NAT_PFT, PCT_NAT_PFT_dims, &PCT_NAT_PFT_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_CFT_dims[0] = cft_dim;
    PCT_CFT_dims[1] = lat_dim;
    PCT_CFT_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "PCT_CFT", NC_DOUBLE, RANK_PCT_CFT, PCT_CFT_dims, &PCT_CFT_id);
    check_err(stat,__LINE__,__FILE__);

    FERTNITRO_CFT_dims[0] = cft_dim;
    FERTNITRO_CFT_dims[1] = lat_dim;
    FERTNITRO_CFT_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "FERTNITRO_CFT", NC_DOUBLE, RANK_FERTNITRO_CFT, FERTNITRO_CFT_dims, &FERTNITRO_CFT_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_VH1_dims[0] = lat_dim;
    HARVEST_VH1_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_VH1", NC_DOUBLE, RANK_HARVEST_VH1, HARVEST_VH1_dims, &HARVEST_VH1_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_VH2_dims[0] = lat_dim;
    HARVEST_VH2_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_VH2", NC_DOUBLE, RANK_HARVEST_VH2, HARVEST_VH2_dims, &HARVEST_VH2_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_SH1_dims[0] = lat_dim;
    HARVEST_SH1_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_SH1", NC_DOUBLE, RANK_HARVEST_SH1, HARVEST_SH1_dims, &HARVEST_SH1_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_SH2_dims[0] = lat_dim;
    HARVEST_SH2_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_SH2", NC_DOUBLE, RANK_HARVEST_SH2, HARVEST_SH2_dims, &HARVEST_SH2_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_SH3_dims[0] = lat_dim;
    HARVEST_SH3_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_SH3", NC_DOUBLE, RANK_HARVEST_SH3, HARVEST_SH3_dims, &HARVEST_SH3_id);
    check_err(stat,__LINE__,__FILE__);

    GRAZING_dims[0] = lat_dim;
    GRAZING_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "GRAZING", NC_DOUBLE, RANK_GRAZING, GRAZING_dims, &GRAZING_id);
    check_err(stat,__LINE__,__FILE__);

    UNREPRESENTED_PFT_LULCC_dims[0] = natpft_dim;
    UNREPRESENTED_PFT_LULCC_dims[1] = lat_dim;
    UNREPRESENTED_PFT_LULCC_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "UNREPRESENTED_PFT_LULCC", NC_DOUBLE, RANK_UNREPRESENTED_PFT_LULCC, UNREPRESENTED_PFT_LULCC_dims, &UNREPRESENTED_PFT_LULCC_id);
    check_err(stat,__LINE__,__FILE__);

    UNREPRESENTED_CFT_LULCC_dims[0] = cft_dim;
    UNREPRESENTED_CFT_LULCC_dims[1] = lat_dim;
    UNREPRESENTED_CFT_LULCC_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "UNREPRESENTED_CFT_LULCC", NC_DOUBLE, RANK_UNREPRESENTED_CFT_LULCC, UNREPRESENTED_CFT_LULCC_dims, &UNREPRESENTED_CFT_LULCC_id);
    check_err(stat,__LINE__,__FILE__);

    /* assign global attributes */

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 8, "NCAR-CSM");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "Author", strlen(authorname), authorname);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "History_Log", strlen(timestamp), timestamp);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "Region", strlen(regionfilename), regionfilename);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMCurrentDB", strlen(ctsmcurrentsurfdb), ctsmcurrentsurfdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMForestDB", strlen(ctsmLUHforestdb), ctsmLUHforestdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMPastureDB", strlen(ctsmLUHpasturedb), ctsmLUHpasturedb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMOtherDB", strlen(ctsmLUHotherdb), ctsmLUHotherdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMC3AnnDB", strlen(ctsmLUHc3anndb), ctsmLUHc3anndb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMC4AnnDB", strlen(ctsmLUHc4anndb), ctsmLUHc4anndb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMC3PerDB", strlen(ctsmLUHc3perdb), ctsmLUHc3perdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMC4PerDB", strlen(ctsmLUHc4perdb), ctsmLUHc4perdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "CTSMC3NfxDB", strlen(ctsmLUHc3nfxdb), ctsmLUHc3nfxdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "LUH2StatesDB",strlen(luhstatesdb),luhstatesdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "LUH2ManagementDB",strlen(luhmanagementdb),luhmanagementdb);
    check_err(stat,__LINE__,__FILE__);
    }
    
    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "LUH2TransitionsDB",strlen(luhtransitionsdb),luhtransitionsdb);
    check_err(stat,__LINE__,__FILE__);
    }
    

    /* assign per-variable attributes */

    {
    stat = nc_put_att_text(ncid, natpft_id, "long_name", 23, "indices of natural PFTs");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, natpft_id, "units", 5, "index");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, cft_id, "long_name", 15, "indices of CFTs");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, cft_id, "units", 5, "index");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEN_id, "long_name", 29, "northern edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEN_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEE_id, "long_name", 28, "eastern edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEE_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGES_id, "long_name", 29, "southern edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGES_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEW_id, "long_name", 28, "western edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEW_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LAT_id, "long_name", 3, "lat");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LAT_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const float mksrf_file__FillValue_att[1] = {((float)9.96921e+36)} ;
    stat = nc_put_att_float(ncid, LATIXY_id, "_FillValue", NC_FLOAT, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LATIXY_id, "long_name", 11, "latitude-2d");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LATIXY_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LON_id, "long_name", 3, "lon");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LON_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const float mksrf_file__FillValue_att[1] = {((float)9.96921e+36)} ;
    stat = nc_put_att_float(ncid, LONGXY_id, "_FillValue", NC_FLOAT, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LONGXY_id, "long_name", 12, "longitude-2d");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LONGXY_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDMASK_id, "long_name", 9, "land mask");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDMASK_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDFRAC_id, "long_name", 25, "land fraction of gridcell");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDFRAC_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, AREA_id, "long_name", 16, "area of gridcell");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, AREA_id, "units", 4, "km^2");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_GLACIER_id, "long_name", 30, "total percent glacier landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_GLACIER_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_GLACIER_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_LAKE_id, "long_name", 27, "total percent lake landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_LAKE_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_LAKE_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_WETLAND_id, "long_name", 30, "total percent wetland landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_WETLAND_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_WETLAND_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_URBAN_id, "long_name", 28, "total percent urban landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_URBAN_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_URBAN_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NATVEG_id, "long_name", 41, "total percent natural vegetation landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NATVEG_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_NATVEG_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CROP_id, "long_name", 27, "total percent crop landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CROP_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_CROP_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NAT_PFT_id, "long_name", 73, "percent plant functional type on the natural veg landunit (% of landunit)");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NAT_PFT_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_NAT_PFT_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CFT_id, "long_name", 65, "percent crop functional type on the crop landunit (% of landunit)");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CFT_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_CFT_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, FERTNITRO_CFT_id, "long_name", 33, "nitrogen fertilizer for each crop");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, FERTNITRO_CFT_id, "units", 8, "gN/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, FERTNITRO_CFT_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH1_id, "long_name", 27, "harvest from primary forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH1_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_VH1_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH2_id, "long_name", 31, "harvest from primary non-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH2_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_VH2_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH1_id, "long_name", 36, "harvest from secondary mature-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH1_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_SH1_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH2_id, "long_name", 35, "harvest from secondary young-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH2_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_SH2_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH3_id, "long_name", 33, "harvest from secondary non-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH3_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_SH3_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, GRAZING_id, "long_name", 25, "grazing of herbacous pfts");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, GRAZING_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, GRAZING_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_PFT_LULCC_id, "long_name", 41, "unrepresented PFT gross LULCC transitions");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_PFT_LULCC_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, UNREPRESENTED_PFT_LULCC_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_CFT_LULCC_id, "long_name", 42, "unrepresented crop gross LULCC transitions");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_CFT_LULCC_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, UNREPRESENTED_CFT_LULCC_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }


    /* leave define mode */
    stat = nc_enddef (ncid);
    check_err(stat,__LINE__,__FILE__);

    /* assign variable data */

    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);
    return 0;
}

int
closencfile() {

    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);

    return 0;

}

int readnc0dfield(char *FieldName, float *targetvalue) {

    int varid;
        
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_get_var_float(ncid, varid, targetvalue);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int readnc1dfield(char *FieldName, float *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_get_var_float(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int readnc1dintfield(char *FieldName, int *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_get_var_int(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int readnc2dfield(char *FieldName, float *targetgrid, int flipgrid) {

    int varid;
    long ctsmlin, ctsmpix, fliplin;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    if (flipgrid == 0) {
        stat =  nc_get_var_float(ncid, varid, targetgrid);
        check_err(stat,__LINE__,__FILE__);
    }
    else {
        stat =  nc_get_var_float(ncid, varid, tempflipGrid);
        check_err(stat,__LINE__,__FILE__);
        for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
            fliplin = MAXOUTLIN - ctsmlin - 1;
            for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
                targetgrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempflipGrid[fliplin * MAXOUTPIX + ctsmpix];
            }
	}
    }
    
    return 0;

}


int readnc3dfield(char *FieldName, int index1d, float *targetgrid, int flipgrid) {

    int varid;
    long ctsmlin, ctsmpix, fliplin;
    size_t start[3], count[3];
    
    count[0] = 1;
    count[1] = MAXOUTLIN;
    count[2] = MAXOUTPIX;
    start[0] = index1d;
    start[1] = 0;
    start[2] = 0;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    if (flipgrid == 0) {
        stat =  nc_get_vara_float(ncid, varid, start, count, targetgrid);
        check_err(stat,__LINE__,__FILE__);
    }
    else {
        stat =  nc_get_vara_float(ncid, varid, start, count, tempflipGrid);
        check_err(stat,__LINE__,__FILE__);
        for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
            fliplin = MAXOUTLIN - ctsmlin - 1;
            for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
                targetgrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempflipGrid[fliplin * MAXOUTPIX + ctsmpix];
            }
	}
    }
        
    return 0;
    
}


int readnc4dfield(char *FieldName, int index1d, int index2d, float *targetgrid, int flipgrid) {

    int varid;
    long ctsmlin, ctsmpix, fliplin;
    size_t start[4], count[4];
    
    count[0] = 1;
    count[1] = 1;
    count[2] = MAXOUTLIN;
    count[3] = MAXOUTPIX;
    start[0] = index1d;
    start[1] = index2d;
    start[2] = 0;
    start[3] = 0;
       
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    if (flipgrid == 0) {
        stat =  nc_get_vara_float(ncid, varid, start, count, targetgrid);
        check_err(stat,__LINE__,__FILE__);
    }
    else {
        stat =  nc_get_vara_float(ncid, varid, start, count, tempflipGrid);
        check_err(stat,__LINE__,__FILE__);
        for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
            fliplin = MAXOUTLIN - ctsmlin - 1;
            for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
                targetgrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempflipGrid[fliplin * MAXOUTPIX + ctsmpix];
            }
	}
    }
        
    return 0;
    
}

int writenc0dfield(char *FieldName, float *targetvalue) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_float(ncid, varid, targetvalue);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc1dfield(char *FieldName, float *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_float(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc1dintfield(char *FieldName, int *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_int(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc2dfield(char *FieldName, float *targetgrid) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_float(ncid, varid, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc3dfield(char *FieldName, int index1d, float *targetgrid) {

    int varid;
    size_t start[3], count[3];
    
    count[0] = 1;
    count[1] = MAXOUTLIN;
    count[2] = MAXOUTPIX;
    start[0] = index1d;
    start[1] = 0;
    start[2] = 0;
        
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_vara_float(ncid, varid, start, count, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;
    
}

int writenc2ddblfield(char *FieldName, double *targetgrid) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_double(ncid, varid, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc3ddblfield(char *FieldName, int index3d, double *targetgrid) {

    int varid;
    size_t start[3], count[3];
    
    count[0] = 1;
    count[1] = MAXOUTLIN;
    count[2] = MAXOUTPIX;
    start[0] = index3d;
    start[1] = 0;
    start[2] = 0;
        
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_vara_double(ncid, varid, start, count, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;
    
}

int readctsmcurrentGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmcurrentsurfreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmcurrentsurfstartyear && ctsmcurrentsurfreadyear == ctsmcurrentsurfstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmcurrentsurfendyear && ctsmcurrentsurfreadyear == ctsmcurrentsurfendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmcurrentsurfstartyear) {
      yearindex = 0;
      ctsmcurrentsurfreadyear = ctsmcurrentsurfstartyear;
  }
  else {
      if (currentyear >= ctsmcurrentsurfendyear) {
          yearindex = ctsmcurrentsurfendyear - ctsmcurrentsurfstartyear;
	  ctsmcurrentsurfreadyear = ctsmcurrentsurfendyear;
      }
      else {
          yearindex = currentyear - ctsmcurrentsurfstartyear;
	  ctsmcurrentsurfreadyear = currentyear;
      }
  }

  openncinputfile(ctsmcurrentsurfdb);
  
  readnc1dintfield("natpft",innatpft);
  readnc1dintfield("cft",incft);
  readnc0dfield("EDGEN",&inEDGEN);
  readnc0dfield("EDGEE",&inEDGEE);
  readnc0dfield("EDGES",&inEDGES);
  readnc0dfield("EDGEW",&inEDGEW);
  readnc1dfield("LAT",inLAT);
  readnc2dfield("LATIXY",inLATIXY,0);
  readnc1dfield("LON",inLON);
  readnc2dfield("LONGXY",inLONGXY,0);
  readnc2dfield("LANDMASK",inLANDMASKGrid,0);
  readnc2dfield("LANDFRAC",inLANDFRACGrid,0);
  readnc2dfield("AREA",inAREAGrid,0);
  readnc2dfield("PCT_GLACIER",inPCTGLACIERGrid,0);
  readnc2dfield("PCT_LAKE",inPCTLAKEGrid,0);
  readnc2dfield("PCT_WETLAND",inPCTWETLANDGrid,0);
  readnc2dfield("PCT_URBAN",inPCTURBANGrid,0);
  readnc3dfield("PCT_NATVEG",yearindex,inPCTNATVEGGrid,0);
  readnc3dfield("PCT_CROP",yearindex,inPCTCROPGrid,0);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc4dfield("PCT_NAT_PFT",pftid,yearindex,inCURRENTPCTPFTGrid[pftid],0);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      readnc4dfield("PCT_CFT",cftid,yearindex,inCURRENTPCTCFTGrid[cftid],0);
  }
  
  closencfile();

  return 0;
  
}


int readctsmLUHforestGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHforestreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHforeststartyear && ctsmLUHforestreadyear == ctsmLUHforeststartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHforestendyear && ctsmLUHforestreadyear == ctsmLUHforestendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHforeststartyear) {
      yearindex = 0;
      ctsmLUHforestreadyear = ctsmLUHforeststartyear;
  }
  else {
      if (currentyear >= ctsmLUHforestendyear) {
          yearindex = ctsmLUHforestendyear - ctsmLUHforeststartyear;
	  ctsmLUHforestreadyear = ctsmLUHforestendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHforeststartyear;
	  ctsmLUHforestreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHforestdb);  

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc4dfield("PCT_NAT_PFT",pftid,yearindex,inFORESTPCTPFTGrid[pftid],0);
  }
  
  closencfile();
  
  return 0;
  
}
  
  
int readctsmLUHpastureGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHpasturereadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHpasturestartyear && ctsmLUHpasturereadyear == ctsmLUHpasturestartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHpastureendyear && ctsmLUHpasturereadyear == ctsmLUHpastureendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHpasturestartyear) {
      yearindex = 0;
      ctsmLUHpasturereadyear = ctsmLUHpasturestartyear;
  }
  else {
      if (currentyear >= ctsmLUHpastureendyear) {
          yearindex = ctsmLUHpastureendyear - ctsmLUHpasturestartyear;
	  ctsmLUHpasturereadyear = ctsmLUHpastureendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHpasturestartyear;
	  ctsmLUHpasturereadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHpasturedb);  

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc4dfield("PCT_NAT_PFT",pftid,yearindex,inPASTUREPCTPFTGrid[pftid],0);
  }
  
  closencfile();
  
  return 0;
  
}


int readctsmLUHotherGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHotherreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHotherstartyear && ctsmLUHotherreadyear == ctsmLUHotherstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHotherendyear && ctsmLUHotherreadyear == ctsmLUHotherendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHotherstartyear) {
      yearindex = 0;
      ctsmLUHotherreadyear = ctsmLUHotherstartyear;
  }
  else {
      if (currentyear >= ctsmLUHotherendyear) {
          yearindex = ctsmLUHotherendyear - ctsmLUHotherstartyear;
	  ctsmLUHotherreadyear = ctsmLUHotherendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHotherstartyear;
	  ctsmLUHotherreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHotherdb);  

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc4dfield("PCT_NAT_PFT",pftid,yearindex,inOTHERPCTPFTGrid[pftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readctsmLUHc3annGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHc3annreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc3annstartyear && ctsmLUHc3annreadyear == ctsmLUHc3annstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHc3annendyear && ctsmLUHc3annreadyear == ctsmLUHc3annendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc3annstartyear) {
      yearindex = 0;
      ctsmLUHc3annreadyear = ctsmLUHc3annstartyear;
  }
  else {
      if (currentyear >= ctsmLUHc3annendyear) {
          yearindex = ctsmLUHc3annendyear - ctsmLUHc3annstartyear;
	  ctsmLUHc3annreadyear = ctsmLUHc3annendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHc3annstartyear;
	  ctsmLUHc3annreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHc3anndb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc4dfield("PCT_CFT",cftid,yearindex,inC3ANNPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readctsmLUHc4annGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHc4annreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc4annstartyear && ctsmLUHc4annreadyear == ctsmLUHc4annstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHc4annendyear && ctsmLUHc4annreadyear == ctsmLUHc4annendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc4annstartyear) {
      yearindex = 0;
      ctsmLUHc4annreadyear = ctsmLUHc4annstartyear;
  }
  else {
      if (currentyear >= ctsmLUHc4annendyear) {
          yearindex = ctsmLUHc4annendyear - ctsmLUHc4annstartyear;
	  ctsmLUHc4annreadyear = ctsmLUHc4annendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHc4annstartyear;
	  ctsmLUHc4annreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHc4anndb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc4dfield("PCT_CFT",cftid,yearindex,inC4ANNPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readctsmLUHc3perGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHc3perreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc3perstartyear && ctsmLUHc3perreadyear == ctsmLUHc3perstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHc3perendyear && ctsmLUHc3perreadyear == ctsmLUHc3perendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc3perstartyear) {
      yearindex = 0;
      ctsmLUHc3perreadyear = ctsmLUHc3perstartyear;
  }
  else {
      if (currentyear >= ctsmLUHc3perendyear) {
          yearindex = ctsmLUHc3perendyear - ctsmLUHc3perstartyear;
	  ctsmLUHc3perreadyear = ctsmLUHc3perendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHc3perstartyear;
	  ctsmLUHc3perreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHc3perdb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc4dfield("PCT_CFT",cftid,yearindex,inC3PERPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readctsmLUHc4perGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHc4perreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc4perstartyear && ctsmLUHc4perreadyear == ctsmLUHc4perstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHc4perendyear && ctsmLUHc4perreadyear == ctsmLUHc4perendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc4perstartyear) {
      yearindex = 0;
      ctsmLUHc4perreadyear = ctsmLUHc4perstartyear;
  }
  else {
      if (currentyear >= ctsmLUHc4perendyear) {
          yearindex = ctsmLUHc4perendyear - ctsmLUHc4perstartyear;
	  ctsmLUHc4perreadyear = ctsmLUHc4perendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHc4perstartyear;
	  ctsmLUHc4perreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHc4perdb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc4dfield("PCT_CFT",cftid,yearindex,inC4PERPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readctsmLUHc3nfxGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == ctsmLUHc3nfxreadyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc3nfxstartyear && ctsmLUHc3nfxreadyear == ctsmLUHc3nfxstartyear) {
      return 0;
  }
  
  if (currentyear >= ctsmLUHc3nfxendyear && ctsmLUHc3nfxreadyear == ctsmLUHc3nfxendyear) {
      return 0;
  }
  
  if (currentyear <= ctsmLUHc3nfxstartyear) {
      yearindex = 0;
      ctsmLUHc3nfxreadyear = ctsmLUHc3nfxstartyear;
  }
  else {
      if (currentyear >= ctsmLUHc3nfxendyear) {
          yearindex = ctsmLUHc3nfxendyear - ctsmLUHc3nfxstartyear;
	  ctsmLUHc3nfxreadyear = ctsmLUHc3nfxendyear;
      }
      else {
          yearindex = currentyear - ctsmLUHc3nfxstartyear;
	  ctsmLUHc3nfxreadyear = currentyear;
      }
  }
  
  openncinputfile(ctsmLUHc3nfxdb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc4dfield("PCT_CFT",cftid,yearindex,inC3NFXPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}


int readLUHbasestateGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == refstatesreadyear) {
      return 0;
  }
  
  if (currentyear <= refstatesstartyear && refstatesreadyear == refstatesstartyear) {
      return 0;
  }
  
  if (currentyear >= refstatesendyear && refstatesreadyear == refstatesendyear) {
      return 0;
  }
  
  if (currentyear <= refstatesstartyear) {
      yearindex = 0;
      refstatesreadyear = refstatesstartyear;
  }
  else {
      if (currentyear >= refstatesendyear) {
          yearindex = refstatesendyear - refstatesstartyear;
	  refstatesreadyear = refstatesendyear;
      }
      else {
          yearindex = currentyear - refstatesstartyear;
	  refstatesreadyear = currentyear;
      }
  }
  
  openncinputfile(refstatesdb); 

  readnc3dfield("primf",yearindex,inBASEPRIMFGrid,flipLUHgrids);
  readnc3dfield("primn",yearindex,inBASEPRIMNGrid,flipLUHgrids);
  readnc3dfield("secdf",yearindex,inBASESECDFGrid,flipLUHgrids);
  readnc3dfield("secdn",yearindex,inBASESECDNGrid,flipLUHgrids);
  readnc3dfield("pastr",yearindex,inBASEPASTRGrid,flipLUHgrids);
  readnc3dfield("range",yearindex,inBASERANGEGrid,flipLUHgrids);
  readnc3dfield("c3ann",yearindex,inBASEC3ANNGrid,flipLUHgrids);
  readnc3dfield("c4ann",yearindex,inBASEC4ANNGrid,flipLUHgrids);
  readnc3dfield("c3per",yearindex,inBASEC3PERGrid,flipLUHgrids);
  readnc3dfield("c4per",yearindex,inBASEC4PERGrid,flipLUHgrids);
  readnc3dfield("c3nfx",yearindex,inBASEC3NFXGrid,flipLUHgrids);
  readnc3dfield("urban",yearindex,inBASEURBANGrid,flipLUHgrids);
  
  closencfile();

  return 0;
  
}


int readLUHcurrstateGrids(int currentyear) {

  int pftid, cftid, yearindex;
  
  if (currentyear == luhcurrentstatesreadyear) {
      return 0;
  }
  
  if (currentyear <= luhstatesstartyear && luhcurrentstatesreadyear == luhstatesstartyear) {
      return 0;
  }
  
  if (currentyear >= luhstatesendyear && luhcurrentstatesreadyear == luhstatesendyear) {
      return 0;
  }
  
  if (currentyear <= luhstatesstartyear) {
      yearindex = 0;
      luhcurrentstatesreadyear = luhstatesstartyear;
  }
  else {
      if (currentyear >= luhstatesendyear) {
          yearindex = luhstatesendyear - luhstatesstartyear;
	  luhcurrentstatesreadyear = luhstatesendyear;
      }
      else {
          yearindex = currentyear - luhstatesstartyear;
	  luhcurrentstatesreadyear = currentyear;
      }
  }
  
  openncinputfile(luhstatesdb); 

  readnc3dfield("primf",yearindex,inCURRPRIMFGrid,flipLUHgrids);
  readnc3dfield("primn",yearindex,inCURRPRIMNGrid,flipLUHgrids);
  readnc3dfield("secdf",yearindex,inCURRSECDFGrid,flipLUHgrids);
  readnc3dfield("secdn",yearindex,inCURRSECDNGrid,flipLUHgrids);
  readnc3dfield("pastr",yearindex,inCURRPASTRGrid,flipLUHgrids);
  readnc3dfield("range",yearindex,inCURRRANGEGrid,flipLUHgrids);
  readnc3dfield("c3ann",yearindex,inCURRC3ANNGrid,flipLUHgrids);
  readnc3dfield("c4ann",yearindex,inCURRC4ANNGrid,flipLUHgrids);
  readnc3dfield("c3per",yearindex,inCURRC3PERGrid,flipLUHgrids);
  readnc3dfield("c4per",yearindex,inCURRC4PERGrid,flipLUHgrids);
  readnc3dfield("c3nfx",yearindex,inCURRC3NFXGrid,flipLUHgrids);
  readnc3dfield("urban",yearindex,inCURRURBANGrid,flipLUHgrids);
  
  closencfile();

  return 0;
  
}

int readLUHprevdeltastateGrids(int prevyear) {

  int pftid, cftid, curryear, yearindex1, yearindex2;
  long ctsmlin, ctsmpix;
  
  if (prevyear == luhprevstatesreadyear) {
      return 0;
  }
  
  if (prevyear <= luhstatesstartyear && luhprevstatesreadyear == luhstatesstartyear) {
      return 0;
  }
  
  if (prevyear >= luhstatesendyear && luhprevstatesreadyear == luhstatesendyear) {
      return 0;
  }
  
  if (prevyear <= luhstatesstartyear) {
      yearindex1 = 0;
      luhprevstatesreadyear = luhstatesstartyear;
  }
  else {
      if (prevyear >= luhstatesendyear) {
          yearindex1 = luhstatesendyear - luhstatesstartyear - 1;
	  luhprevstatesreadyear = luhstatesendyear;
      }
      else {
          yearindex1 = prevyear - luhstatesstartyear;
	  luhprevstatesreadyear = prevyear;
      }
  }
  
  curryear = prevyear + 1;
  if (curryear < luhstatesstartyear) {
      yearindex2 = 0;
  }
  else {
      if (curryear >= luhstatesendyear) {
          yearindex2 = luhstatesendyear - luhstatesstartyear;
      }
      else {
          yearindex2 = yearindex1 + 1;
      }
  }
    
  openncinputfile(luhstatesdb); 

  readnc3dfield("secdf",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("secdf",yearindex2,inPREVDELTASECDFGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTASECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTASECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdn",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("secdn",yearindex2,inPREVDELTASECDNGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTASECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTASECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("pastr",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("pastr",yearindex2,inPREVDELTAPASTRGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTAPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTAPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("range",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("range",yearindex2,inPREVDELTARANGEGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTARANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTARANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3ann",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("c3ann",yearindex2,inPREVDELTAC3ANNGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTAC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTAC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c4ann",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("c4ann",yearindex2,inPREVDELTAC4ANNGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTAC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTAC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3per",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("c3per",yearindex2,inPREVDELTAC3PERGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTAC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTAC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c4per",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("c4per",yearindex2,inPREVDELTAC4PERGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTAC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTAC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3nfx",yearindex1,tempGrid,flipLUHgrids);
  readnc3dfield("c3nfx",yearindex2,inPREVDELTAC3NFXGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              inPREVDELTAC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVDELTAC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix] - tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }
  
  closencfile();

  return 0;
  
}


int readLUHwoodharvestGrids(int prevyear) {

  int pftid, cftid, yearindex;
  
  if (prevyear == luhwoodharvestreadyear) {
      return 0;
  }
  
  if (prevyear <= luhtransitionsstartyear && luhwoodharvestreadyear == luhtransitionsstartyear) {
      return 0;
  }
  
  if (prevyear >= luhtransitionsendyear && luhwoodharvestreadyear == luhtransitionsendyear) {
      return 0;
  }
  
  if (prevyear <= luhtransitionsstartyear) {
      yearindex = 0;
      luhwoodharvestreadyear = luhtransitionsstartyear;
  }
  else {
      if (prevyear >= luhtransitionsendyear) {
          yearindex = luhtransitionsendyear - luhtransitionsstartyear;
	  luhwoodharvestreadyear = luhtransitionsendyear;
      }
      else {
          yearindex = prevyear - luhtransitionsstartyear;
	  luhwoodharvestreadyear = prevyear;
      }
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("primf_harv",yearindex,inHARVESTVH1Grid,flipLUHgrids);
  readnc3dfield("primn_harv",yearindex,inHARVESTVH2Grid,flipLUHgrids);
  readnc3dfield("secmf_harv",yearindex,inHARVESTSH1Grid,flipLUHgrids);
  readnc3dfield("secyf_harv",yearindex,inHARVESTSH2Grid,flipLUHgrids);
  readnc3dfield("secnf_harv",yearindex,inHARVESTSH3Grid,flipLUHgrids);
  readnc3dfield("primf_bioh",yearindex,inBIOHVH1Grid,flipLUHgrids);
  readnc3dfield("primn_bioh",yearindex,inBIOHVH2Grid,flipLUHgrids);
  readnc3dfield("secmf_bioh",yearindex,inBIOHSH1Grid,flipLUHgrids);
  readnc3dfield("secyf_bioh",yearindex,inBIOHSH2Grid,flipLUHgrids);
  readnc3dfield("secnf_bioh",yearindex,inBIOHSH3Grid,flipLUHgrids);
  
  closencfile();
  
  return 0;
  
}


int readUNREPSECDFGrids(int prevyear) {

  int yearindex;
  long ctsmlin, ctsmpix;
  float cropstatechange, secdfstatechange, secdfotherchange, secdfresidualchange; 
  float secdfcropinval, secdfotherinval, secdfotheroutval, latval;
  float unreploss;
  
  if (prevyear == luhsecdfunrepreadyear) {
      return 0;
  }
  
  if (prevyear <= luhtransitionsstartyear && luhsecdfunrepreadyear == luhtransitionsstartyear) {
      return 0;
  }
  
  if (prevyear >= luhtransitionsendyear && luhsecdfunrepreadyear == luhtransitionsendyear) {
      return 0;
  }
  
  if (prevyear <= luhtransitionsstartyear) {
      yearindex = 0;
      luhsecdfunrepreadyear = luhtransitionsstartyear;
  }
  else {
      if (prevyear >= luhtransitionsendyear) {
          yearindex = luhtransitionsendyear - luhtransitionsstartyear;
	  luhsecdfunrepreadyear = luhtransitionsendyear;
      }
      else {
          yearindex = prevyear - luhtransitionsstartyear;
	  luhsecdfunrepreadyear = prevyear;
      }
  }
  
  openncinputfile(luhtransitionsdb); 

  readnc3dfield("c3ann_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
	  else {
	      secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }

  readnc3dfield("c4ann_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3per_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c4per_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3nfx_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("primf_harv",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
	  else {
              secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }

  readnc3dfield("secdn_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("pastr_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("range_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("urban_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdf_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
	  else {
              secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }

  readnc3dfield("secdf_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdf_to_range",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdf_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }


  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
	  cropstatechange = inPREVDELTAC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (cropstatechange > 0.0) {
	      cropstatechange = 0.0;
	  }
          secdfstatechange = inPREVDELTASECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          secdfcropinval = secdfCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  secdfotherinval = secdfOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  secdfotheroutval = secdfOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  secdfotherchange = secdfotherinval - secdfotheroutval;
	  secdfresidualchange = secdfotherchange - secdfstatechange;
	  if (secdfresidualchange > 0.0) {
	      secdfresidualchange = 0.0;
	  }
	  unreploss = secdfcropinval + secdfresidualchange + cropstatechange;
	  if (unreploss < 0.0) {
	      unreploss = 0.0;
	  }

          latval = inLATIXY[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (latval > 30.0 || latval < -30.0) {
	      unreploss = 0.0;
	  }

          inUNREPSECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix] = unreploss;
      }
  }
  
  closencfile();

  return 0;
  
}


int readUNREPSECDNGrids(int prevyear) {

  int yearindex;
  long ctsmlin, ctsmpix;
  float cropstatechange, secdnstatechange, secdnotherchange, secdnresidualchange; 
  float secdncropinval, secdnotherinval, secdnotheroutval, latval;
  float unreploss;
  
  if (prevyear == luhsecdnunrepreadyear) {
      return 0;
  }
  
  if (prevyear <= luhtransitionsstartyear && luhsecdnunrepreadyear == luhtransitionsstartyear) {
      return 0;
  }
  
  if (prevyear >= luhtransitionsendyear && luhsecdnunrepreadyear == luhtransitionsendyear) {
      return 0;
  }
  
  if (prevyear <= luhtransitionsstartyear) {
      yearindex = 0;
      luhsecdnunrepreadyear = luhtransitionsstartyear;
  }
  else {
      if (prevyear >= luhtransitionsendyear) {
          yearindex = luhtransitionsendyear - luhtransitionsstartyear;
	  luhsecdnunrepreadyear = luhtransitionsendyear;
      }
      else {
          yearindex = prevyear - luhtransitionsstartyear;
	  luhsecdnunrepreadyear = prevyear;
      }
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("c3ann_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
	  else {
	      secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }

  readnc3dfield("c4ann_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3per_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c4per_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("c3nfx_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("primn_harv",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
	  else {
              secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }

  readnc3dfield("secdf_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("pastr_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("range_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("urban_to_secdn",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdn_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
	  else {
              secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }

  readnc3dfield("secdn_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdn_to_range",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }

  readnc3dfield("secdn_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (tempGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= 1.0) {
              secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + tempGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  }
      }
  }


  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
	  cropstatechange = inPREVDELTAC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  cropstatechange += inPREVDELTAC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (cropstatechange > 0.0) {
	      cropstatechange = 0.0;
	  }
          secdnstatechange = inPREVDELTASECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          secdncropinval = secdnCROPINGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  secdnotherinval = secdnOTHERINGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  secdnotheroutval = secdnOTHEROUTGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  secdnotherchange = secdnotherinval - secdnotheroutval;
	  secdnresidualchange = secdnotherchange - secdnstatechange;
	  if (secdnresidualchange > 0.0) {
	      secdnresidualchange = 0.0;
	  }
	  unreploss = secdncropinval + secdnresidualchange + cropstatechange;
	  if (unreploss < 0.0) {
	      unreploss = 0.0;
	  }

          latval = inLATIXY[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (latval > 30.0 || latval < -30.0) {
	      unreploss = 0.0;
	  }

          inUNREPSECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readLUHcropmanagementGrids(int currentyear) {

  int yearindex;
  long ctsmlin, ctsmpix;
  
  if (currentyear == luhcropmanagementreadyear) {
      return 0;
  }
  
  if (currentyear <= luhmanagementstartyear && luhcropmanagementreadyear == luhmanagementstartyear) {
      return 0;
  }
  
  if (currentyear >= luhmanagementendyear && luhcropmanagementreadyear == luhmanagementendyear) {
      return 0;
  }
  
  if (currentyear <= luhmanagementstartyear) {
      yearindex = 0;
      luhcropmanagementreadyear = luhmanagementstartyear;
  }
  else {
      if (currentyear >= luhmanagementendyear) {
          yearindex = luhmanagementendyear - luhmanagementstartyear;
	  luhcropmanagementreadyear = luhmanagementendyear;
      }
      else {
          yearindex = currentyear - luhmanagementstartyear;
	  luhcropmanagementreadyear = currentyear;
      }
  }
  
  openncinputfile(luhmanagementdb); 
  
  readnc3dfield("irrig_c3ann",yearindex,inIRRIGC3ANNGrid,flipLUHgrids);
  readnc3dfield("irrig_c4ann",yearindex,inIRRIGC4ANNGrid,flipLUHgrids);
  readnc3dfield("irrig_c3per",yearindex,inIRRIGC3PERGrid,flipLUHgrids);
  readnc3dfield("irrig_c4per",yearindex,inIRRIGC4PERGrid,flipLUHgrids);
  readnc3dfield("irrig_c3nfx",yearindex,inIRRIGC3NFXGrid,flipLUHgrids);

  readnc3dfield("fertl_c3ann",yearindex,inFERTC3ANNGrid,flipLUHgrids);
  readnc3dfield("fertl_c4ann",yearindex,inFERTC4ANNGrid,flipLUHgrids);
  readnc3dfield("fertl_c3per",yearindex,inFERTC3PERGrid,flipLUHgrids);
  readnc3dfield("fertl_c4per",yearindex,inFERTC4PERGrid,flipLUHgrids);
  readnc3dfield("fertl_c3nfx",yearindex,inFERTC3NFXGrid,flipLUHgrids);

  closencfile();

  return 0;
  
}


int extrapindvidualfertGrid(float *cropgrid, float *fertgrid) {

  long ctsmlin, ctsmpix;
  long searchlin, searchpix;
  long searchboxinside, searchboxoutside;
  float allcropfraction, cropfraction, searchcrop, searchfert, searchcropsum, searchfertsum;
  int searchlinokay, searchpixokay;
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {

          tempextrapGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  
	  allcropfraction = inCURRC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
      
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1 && allcropfraction >= 0.0 && allcropfraction <= 1.0) {
	      cropfraction = cropgrid[ctsmlin * MAXOUTPIX + ctsmpix];
              if (cropfraction > 0.0 && cropfraction <= 1.0) {
	          tempextrapGrid[ctsmlin * MAXOUTPIX + ctsmpix] = fertgrid[ctsmlin * MAXOUTPIX + ctsmpix];
              }
	      else {
	          searchcropsum = 0.0;
		  searchfertsum = 0.0;
		  searchboxinside = 0;
		  searchboxoutside = 2;
		  while (searchcropsum == 0.0) {
                      for (searchlin = ctsmlin - searchboxoutside; searchlin <= ctsmlin + searchboxoutside; searchlin++) {
                          searchlinokay = 1;
	                  if (searchlin < 0 || searchlin >= MAXOUTLIN) {
                              searchlinokay = 0;
                          }
                          if (searchlin >= ctsmlin - searchboxinside && searchlin <= ctsmlin + searchboxinside) {
                              searchlinokay = 0;
                          }
                          if (searchlinokay == 1) {
                              for (searchpix = ctsmpix - searchboxoutside; searchpix <= ctsmpix + searchboxoutside; searchpix++) {
                                  searchpixokay = 1;
	                          if (searchpix < 0 || searchpix >= MAXOUTPIX) {
                                      searchpixokay = 0;
                                  }
                                  if (searchpix >= ctsmpix - searchboxinside && searchpix <= ctsmpix + searchboxinside) {
                                      searchpixokay = 0;
                                  }
                                  if (searchpixokay) {
                                      searchcrop = cropgrid[searchlin * MAXOUTPIX + searchpix];
                                      if (searchcrop > 0.0 && searchcrop <= 1.0) {
                                          searchfert = fertgrid[searchlin * MAXOUTPIX + searchpix];
                                          if (searchfert >= 0.0 && searchfert < 10000.0) {
                                              searchcropsum += searchcrop;
                                              searchfertsum += searchcrop * searchfert;
                                          }
                                      }
                                  }
                              }
                          }
                      }
		      searchboxinside = searchboxoutside;
		      searchboxoutside = searchboxoutside * 2; 
                      if (searchboxoutside > 16) {
		          searchfertsum = fertgrid[ctsmlin * MAXOUTPIX + ctsmpix];
		          searchcropsum = 1.0;
                      }
		  }
                  tempextrapGrid[ctsmlin * MAXOUTPIX + ctsmpix] = searchfertsum / searchcropsum;		      
              }
          }
      }
  }
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          fertgrid[ctsmlin * MAXOUTPIX + ctsmpix] = tempextrapGrid[ctsmlin * MAXOUTPIX + ctsmpix];
      }
  }

  return 0;

}


int extrapfertGrids() {

  long ctsmlin, ctsmpix;
  float allcropfraction;

  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {

          allcropfraction = inCURRC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          allcropfraction += inCURRC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          allcropfraction += inCURRC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          allcropfraction += inCURRC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          allcropfraction += inCURRC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  
	  if (allcropfraction >= 0.0 && allcropfraction <= 1.0) {
              inCURRCropGrid[ctsmlin * MAXOUTPIX + ctsmpix] = allcropfraction;
	  }
	  else {
	      inCURRCropGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
      }
  }
  
  extrapindvidualfertGrid(inCURRC3ANNGrid,inFERTC3ANNGrid);
  extrapindvidualfertGrid(inCURRC4ANNGrid,inFERTC4ANNGrid);
  extrapindvidualfertGrid(inCURRC3PERGrid,inFERTC3PERGrid);
  extrapindvidualfertGrid(inCURRC4PERGrid,inFERTC4PERGrid);
  extrapindvidualfertGrid(inCURRC3NFXGrid,inFERTC3NFXGrid);
  
  return 0;
  
}



int generateLUHcollectionGrids() {

  long ctsmlin, ctsmpix;
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
      
          inBASEFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inBASEPRIMFGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASESECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inBASEFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASEFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inBASENONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inBASEPRIMNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASESECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inBASENONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASENONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inBASECROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inBASEC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASEC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASEC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASEC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASEC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inBASECROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASECROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inBASEURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inBASEURBANGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inBASEURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASEURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
	  if (inBASEPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASEPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
	  if (inBASERANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASERANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inBASEOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inBASEPRIMNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASESECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASERANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inBASEOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inBASEOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inBASEMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0 - inBASEFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inBASENONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inBASEPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inBASERANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inBASECROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          if (inBASEMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] < 0.0) {
              inBASEMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
          if (inBASEMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 1.0) {
              inBASEMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
          }
          inBASENATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inBASEFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASEPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inBASEOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix];

          inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRPRIMFGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRSECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inCURRNONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRPRIMNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRSECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inCURRNONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inCURRNONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inCURRCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inCURRCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inCURRCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inPREVCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inPREVDELTAC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inPREVDELTAC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inPREVDELTAC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inPREVDELTAC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inPREVDELTAC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inPREVCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inPREVCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inCURRURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRURBANGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inCURRURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inCURRURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inCURRMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0 - inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inCURRNONFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inCURRPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inCURRRANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix] - inCURRCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          if (inCURRMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] < 0.0) {
              inCURRMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
          if (inCURRMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 1.0) {
              inCURRMISSINGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
          }
          inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRPRIMNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRSECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRRANGEGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	  if (inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 10.0) {
	      inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
	  }
          inCURRNATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURRPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix];

          inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inUNREPSECDFGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          if (inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] < 0.0) {
              inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
          if (inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 1.0) {
              inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
          }
	  
          inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inUNREPSECDNGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          if (inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] < 0.0) {
              inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
          if (inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] > 1.0) {
              inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
          }
	  
	  if (inPREVCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] < (inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix])) {
	      if ((inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix]) > 0.0) {
	          inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] * inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] / (inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix]);
	          inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] = inPREVCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] * inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] / (inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] + inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix]);
              }
          }

      }
  }
            
  return 0;
  
}


int generatectsmURBANGrids() {

  long ctsmlin, ctsmpix;
  float pctglacierval, pctlakeval, pctwetlandval, pctavail;
  float pcturbanval;

  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1.0) {
              pctglacierval = inPCTGLACIERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              pctlakeval = inPCTLAKEGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              pctwetlandval = inPCTWETLANDGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	      pctavail = 100.0 - pctglacierval + pctlakeval + pctwetlandval;
	      pcturbanval = inCURRURBANTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] * 100.0;
	      if (pcturbanval > pctavail) {
	          pcturbanval = pctavail;
              }
	      outPCTURBANGrid[ctsmlin * MAXOUTPIX + ctsmpix] = pcturbanval;
	  }
      }
  }

  return 0;

}


int generatectsmPFTGrids() {

  long ctsmlin, ctsmpix;
  int pftid;
  float pctnatvegval, pctnatvegbase;
  float foresttotalbaseval, foresttotalfracval, foresttotalfracdelta, foresttotalcurrentval;
  float pasturebaseval, pasturefracval, pasturecurrentval, pasturefracdelta;
  float otherbaseval, otherfracval, othercurrentval, otherfracdelta, missingbaseval;
  float currentpctforestpft, deltapctforestpft, unreppctforestpft;
  float currentpctpasturepft, deltapctpasturepft;
  float currentpctotherpft, deltapctotherpft;
  float currentpctmissingpft, unreppctotherpft;
  float forestunrepfrac, otherunrepfrac;
  float newpctpft, unrepfrac, newpctpfttotal;
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1) {
              pctnatvegval = inCURRNATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] * 100.0;
              forestunrepfrac = 0.0;
              otherunrepfrac = 0.0;
              if (pctnatvegval > 0.0) {
                  outPCTNATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = pctnatvegval;
                  pctnatvegbase = inBASENATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] * 100.0;
                  if (pctnatvegbase > 0.0) {
                      foresttotalbaseval = inBASEFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] / pctnatvegbase * 100.0;
                      foresttotalfracval = inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] / pctnatvegval * 100.0;
                      foresttotalfracdelta = foresttotalfracval - foresttotalbaseval;
                      if (foresttotalfracdelta >= 0.0) {
                          foresttotalcurrentval = foresttotalbaseval;
                      }
                      else {
                          foresttotalcurrentval = foresttotalbaseval + foresttotalfracdelta;
                          foresttotalfracdelta = 0.0;
                      }
		      if (inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] >= 0.01) {
		          if (inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix]) {
                              forestunrepfrac = inUNREPFORESTGrid[ctsmlin * MAXOUTPIX + ctsmpix] / inCURRFORESTTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix];
			  }
			  else {
			      forestunrepfrac = 0.0;
			  }
		      }
                      pasturebaseval = inBASEPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] / pctnatvegbase * 100.0;
                      pasturefracval = inCURRPASTRGrid[ctsmlin * MAXOUTPIX + ctsmpix] / pctnatvegval * 100.0;
                      pasturecurrentval = 0.0;
                      pasturefracdelta = pasturefracval;
                      otherbaseval = inBASEOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] / pctnatvegbase * 100.0;
                      otherfracval = inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] / pctnatvegval * 100.0;
                      otherfracdelta = otherfracval - otherbaseval;
		      missingbaseval = 0.0;

                      if (otherfracdelta >= 0.0) {
                          othercurrentval = otherbaseval;
                      }
                      else {
                          othercurrentval = otherbaseval + otherfracdelta;
                          otherfracdelta = 0.0;
                      }
		      if (inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] >= 0.01) {
		          if (inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] <= inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix]) {
                              otherunrepfrac = inUNREPOTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix] / inCURROTHERGrid[ctsmlin * MAXOUTPIX + ctsmpix];
			  }
			  else {
			      otherunrepfrac = 0.0;
			  }
		      }
                  }
                  else {
                      foresttotalcurrentval = 0.0;
                      foresttotalfracdelta = 0.0;
                      pasturecurrentval = 0.0;
                      pasturefracdelta = 0.0;
                      othercurrentval = 0.0;
                      otherfracdelta = 0.0;
		      missingbaseval = 1.0;
		      
                  }
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      currentpctforestpft = foresttotalcurrentval * inCURRENTPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      deltapctforestpft = foresttotalfracdelta * inFORESTPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      unreppctforestpft = forestunrepfrac * (currentpctforestpft + deltapctforestpft);
                      currentpctpasturepft = pasturecurrentval * inCURRENTPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      deltapctpasturepft = pasturefracdelta * inPASTUREPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      currentpctotherpft = othercurrentval * inCURRENTPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      deltapctotherpft = otherfracdelta * inOTHERPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
		      currentpctmissingpft = missingbaseval * inCURRENTPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
		      unreppctotherpft = otherunrepfrac * (currentpctotherpft + deltapctotherpft);
                      newpctpft = currentpctforestpft + deltapctforestpft + currentpctpasturepft + deltapctpasturepft + currentpctotherpft + deltapctotherpft + currentpctmissingpft;
		      if (pftid > 0 && newpctpft > 0.0) {
		          unrepfrac = (unreppctforestpft + unreppctotherpft) / newpctpft; 
			  if (unrepfrac < 0.001) {
			      unrepfrac = 0.0;
			  }
			  if (unrepfrac > 0.25) { 
			      unrepfrac = 0.25;
			  }
		      }
		      else {
		          unrepfrac = 0.0;
		      }
                      outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = newpctpft;
                      outUNREPPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = unrepfrac;
                  }
                  newpctpfttotal = 0.0;
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      newpctpfttotal = newpctpfttotal + outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                  }
                  if (newpctpfttotal > 0.0) {
                      for (pftid = 0; pftid < MAXPFT; pftid++) {
                          newpctpft = outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                          if (newpctpft > 0.0) {
                              newpctpft = newpctpft / newpctpfttotal * 100.0;
                              outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = newpctpft;
                          }
                          else {
                              outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                          }
                      }
                  }
/*		  else {
		      printf("No newpctpfttotal %f at %d %d\n",newpctpfttotal,ctsmlin,ctsmpix);
		  } */
              }
              else {
                  outPCTNATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = inCURRENTPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      outUNREPPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
              }
          }
      }
  }

  return 0;
  
}


int generatectsmCFTGrids() {

  long ctsmlin, ctsmpix;
  int cftid, rawcftid, rainfedcftid, irrigcftid;
  float pctcropval, c3annunrepval, c4annunrepval, c3perunrepval, c4perunrepval, c3nfxunrepval;
  float newpctrainfedcft, newpctirrigcft, newunreprainfedval, newunrepirrigval;
  float newpctcroptotal, newpctcft;

  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1) {
              pctcropval = inCURRCROPTOTALGrid[ctsmlin * MAXOUTPIX + ctsmpix] * 100.0;
              if (pctcropval > 0.0 && pctcropval <= 100.0) {
                  outPCTCROPGrid[ctsmlin * MAXOUTPIX + ctsmpix] = pctcropval;
                  c3annunrepval = 0.0; /* inUNREPC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix]; */
                  c4annunrepval = 0.0; /* inUNREPC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix]; */
                  c3perunrepval = 0.0; /* inUNREPC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix]; */
                  c4perunrepval = 0.0; /* inUNREPC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix]; */
                  c3nfxunrepval = 0.0; /* inUNREPC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix]; */
                  for (rawcftid = 0; rawcftid < MAXCFTRAW; rawcftid++) {
                      rainfedcftid = 2 * (rawcftid);
                      irrigcftid = 2 * (rawcftid) + 1;
                      newpctrainfedcft = inCURRC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (1.0 - inIRRIGC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC3ANNPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newpctirrigcft = inCURRC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (inIRRIGC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC3ANNPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newunreprainfedval = c3annunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c3annunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunreprainfedval;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctirrigcft;
                          outUNREPCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunrepirrigval;
                      }
                      newpctrainfedcft = inCURRC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (1.0 - inIRRIGC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC4ANNPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newpctirrigcft = inCURRC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (inIRRIGC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC4ANNPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newunreprainfedval = c4annunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c4annunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunreprainfedval;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunrepirrigval;
                      }
                      newpctrainfedcft = inCURRC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (1.0 - inIRRIGC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC3PERPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newpctirrigcft = inCURRC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (inIRRIGC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC3PERPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newunreprainfedval = c3perunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c3perunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunreprainfedval;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunrepirrigval;
                      }
                      newpctrainfedcft = inCURRC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (1.0 - inIRRIGC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC4PERPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newpctirrigcft = inCURRC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (inIRRIGC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC4PERPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newunreprainfedval = c4perunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c4perunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunreprainfedval;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunrepirrigval;
                      }
                      newpctrainfedcft = inCURRC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (1.0 - inIRRIGC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC3NFXPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newpctirrigcft = inCURRC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix] * (inIRRIGC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix]) * inC3NFXPCTCFTGrid[rawcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                      newunreprainfedval = c3nfxunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c3nfxunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunreprainfedval;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = outUNREPCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] + newunrepirrigval;
                      }
                  }
                  newpctcroptotal = 0.0;
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      newpctcroptotal = newpctcroptotal + outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix];
                  }
                  if (newpctcroptotal > 0.0) {
                      for (cftid = 0; cftid < MAXCFT; cftid++) {
                          newpctcft = outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix];
                          if (newpctcft > 0.0) {
                              newpctcft = newpctcft / newpctcroptotal * 100.0;
                              outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = newpctcft;
                          }
                          else {
                              outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                          }
                      }
                  }
              }
              else {
                  outPCTCROPGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  outPCTCFTGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
                  outUNREPCFTGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  for (cftid = 1; cftid < MAXCFT; cftid++) {
                      outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                      outUNREPCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
              }
          }
      }
  }

  return 0;
  
}


int generatectsmwoodharvestGrids() {

  long ctsmlin, ctsmpix;
  int pftid;
  float TreePFTArea, TreeFrac, TreeScale, PFTArea;
  float newharvestvh1, newharvestvh2, newharvestsh1, newharvestsh2, newharvestsh3;
  float newbiohvh1, newbiohvh2, newbiohsh1, newbiohsh2, newbiohsh3;

  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1.0) {
              TreePFTArea = 0.0;
              TreeFrac = 0.0;
              PFTArea = inAREAGrid[ctsmlin * MAXOUTPIX + ctsmpix] * inLANDFRACGrid[ctsmlin * MAXOUTPIX + ctsmpix] * outPCTNATVEGGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 100.0 * 1.0e6;              
              for (pftid = firsttreepft; pftid <= lasttreepft; pftid++) {
                  TreePFTArea = TreePFTArea + PFTArea * outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] / 100.0;
                  TreeFrac = TreeFrac + outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] / 100.0;
              }
              TreeScale = 1.0;
              if (TreePFTArea > 1.0e6) {
                  newharvestvh1 = inHARVESTVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] * TreeScale;
                  if (newharvestvh1 < 0.0 || newharvestvh1 > 9.0e4) {
                      newharvestvh1 = 0.0;
                  }
                  if (newharvestvh1 > 0.98) {
                      newharvestvh1 = 0.98;
                  }
                  outHARVESTVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newharvestvh1;
                  newbiohvh1 = inBIOHVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohvh1 < 0.0) {
                      newbiohvh1 = 0.0;
                  }
                  if (newbiohvh1 > 10000.0) {
                      newbiohvh1 = 10000.0;
                  }
                  outBIOHVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newbiohvh1;
                  newharvestvh2 = inHARVESTVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] * TreeScale;
                  if (newharvestvh2 < 0.0 || newharvestvh2 > 9.0e4) {
                      newharvestvh2 = 0.0;
                  }
                  if (newharvestvh2 > 0.98) {
                      newharvestvh2 = 0.98;
                  }
                  outHARVESTVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newharvestvh2;
                  newbiohvh2 = inBIOHVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohvh2 < 0.0) {
                      newbiohvh2 = 0.0;
                  }
                  if (newbiohvh2 > 10000.0) {
                      newbiohvh2 = 10000.0;
                  }
                  outBIOHVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newbiohvh2;
                  newharvestsh1 = inHARVESTSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] * TreeScale;
                  if (newharvestsh1 < 0.0 || newharvestsh1 > 9.0e4) {
                      newharvestsh1 = 0.0;
                  }
                  if (newharvestsh1 > 0.98) {
                      newharvestsh1 = 0.98;
                  }
                  outHARVESTSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newharvestsh1;
                  newbiohsh1 = inBIOHSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohsh1 < 0.0) {
                      newbiohsh1 = 0.0;
                  }
                  if (newbiohsh1 > 10000.0) {
                      newbiohsh1 = 10000.0;
                  }
                  outBIOHSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newbiohsh1;
                  newharvestsh2 = inHARVESTSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] * TreeScale;
                  if (newharvestsh2 < 0.0 || newharvestsh2 > 9.0e4) {
                      newharvestsh2 = 0.0;
                  }
                  if (newharvestsh2 > 0.98) {
                      newharvestsh2 = 0.98;
                  }
                  outHARVESTSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newharvestsh2;
                  newbiohsh2 = inBIOHSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohsh2 < 0.0) {
                      newbiohsh2 = 0.0;
                  }
                  if (newbiohsh2 > 10000.0) {
                      newbiohsh2 = 10000.0;
                  }
                  outBIOHSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newbiohsh2;
                  newharvestsh3 = inHARVESTSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] * TreeScale;
                  if (newharvestsh3 < 0.0 || newharvestsh3 > 9.0e4) {
                      newharvestsh3 = 0.0;
                  }
                  if (newharvestsh3 > 0.98) {
                      newharvestsh3 = 0.98;
                  }
                  outHARVESTSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newharvestsh3;
                  newbiohsh3 = inBIOHSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohsh3 < 0.0) {
                      newbiohsh3 = 0.0;
                  }
                  if (newbiohsh3 > 10000.0) {
                      newbiohsh3 = 10000.0;
                  }
                  outBIOHSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] = newbiohsh3;
              }
          }
      }
  }

  return 0;
  
}


float calcmaxtotalbioh(float TreePFTArea, float TreePFTWeightedArea, float PFTArea) {

  float maxbioh, maxtotalbioh;
  
  if (PFTArea > 0.0) {
      maxbioh = 10.0 + 790.0 * TreePFTWeightedArea / PFTArea;
  }
  else {
      maxbioh = 0.0;
  }
  
  maxtotalbioh = TreePFTArea * maxbioh;
  
  return maxtotalbioh;
  
}


int generatectsmbiohdirectGrids() {

  long ctsmlin, ctsmpix;  

  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          outRBIOHVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = outBIOHVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix];
          outRBIOHVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = outBIOHVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix];
          outRBIOHSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix] = outBIOHSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix];
          outRBIOHSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix] = outBIOHSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix];
          outRBIOHSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix] = outBIOHSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix];
      }
  }

  return 0;

}


int generatectsmfertGrids() {

  long ctsmlin, ctsmpix;
  int cftid, rawcftid, rainfedcftid, irrigcftid;
  float pctcropval, pctrainfedcft, pctirrigcft;
  float fertamount;

  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1) {
              for (rawcftid = 0; rawcftid < MAXCFTRAW; rawcftid++) {
                  rainfedcftid = 2 * (rawcftid);
                  irrigcftid = 2 * (rawcftid) + 1;
	          fertamount = 0.0;
	          if (strcmp(CFTRAWluhtype[rawcftid],"C3ANN") == 0) {
		      fertamount = inFERTC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 10.0;
		  }
	          if (strcmp(CFTRAWluhtype[rawcftid],"C4ANN") == 0) {
		      fertamount = inFERTC4ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 10.0;
		  }
	          if (strcmp(CFTRAWluhtype[rawcftid],"C3PER") == 0) {
		      fertamount = inFERTC3PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 10.0;
		  }
	          if (strcmp(CFTRAWluhtype[rawcftid],"C4PER") == 0) {
		      fertamount = inFERTC4PERGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 10.0;
		  }
	          if (strcmp(CFTRAWluhtype[rawcftid],"C3NFX") == 0) {
		      fertamount = inFERTC3NFXGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 10.0;
		  }
	          if (strcmp(CFTRAWluhtype[rawcftid],"EXCLD") == 0) {
		      fertamount = inFERTC3ANNGrid[ctsmlin * MAXOUTPIX + ctsmpix] / 10.0;
		  }
		  pctrainfedcft = outPCTCFTGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                  if (pctrainfedcft >= 0.0) {
                      outFERTNITROGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = fertamount;
                  }
		  else {
                      outFERTNITROGrid[rainfedcftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
		  pctirrigcft = outPCTCFTGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix];
                  if (pctirrigcft >= 0.0) {
                      outFERTNITROGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = fertamount;
                  }
		  else {
                      outFERTNITROGrid[irrigcftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
              }
          }
      }
  }

  return 0;
  
}



double truncCTSMValues(double CTSMValueIn, double CTSMMaxValue, double CTSMThreshold) {

  double CTSMValueOut = 0.0;
  
  if (CTSMValueIn >= CTSMThreshold) {
      CTSMValueOut = ((double) ((int) (CTSMValueIn / CTSMThreshold))) * CTSMThreshold;
  }
  
  if (CTSMValueOut > CTSMMaxValue) {
      CTSMValueOut = CTSMMaxValue;
  }
  
  return CTSMValueOut;
    
}


int generatedblGrids() {

  long ctsmlin, ctsmpix;
  int pftid, cftid;
  double LandFRAC, AvailPCT, OtherPCT, AllPFTs, AllCFTs, CropPCT, NatVegPCT, tempdblPCT;
  double LargestPCTPFT, RescaledPCTPFT, LargestPCTCFT, RescaledPCTCFT;
  int LargestPFT, LargestCFT;
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          outAREAdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(inAREAGrid[ctsmlin * MAXOUTPIX + ctsmpix],1000000.0,0.001);
          if (inLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 1.0) {
              outLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
              LandFRAC = truncCTSMValues(inLANDFRACGrid[ctsmlin * MAXOUTPIX + ctsmpix],1.0,0.0001);
              outLANDFRACdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC;
              outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(inPCTGLACIERGrid[ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
              outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(inPCTLAKEGrid[ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
              outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(inPCTWETLANDGrid[ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
              outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outPCTURBANGrid[ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
	      OtherPCT = outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] + outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] + outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] + outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              AvailPCT = truncCTSMValues(100.0 - OtherPCT,100.0,0.01);
	      if (AvailPCT < 0.0) {
	          AvailPCT = 0.0;
              }
              CropPCT = truncCTSMValues(outPCTCROPGrid[ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
	      if (CropPCT > AvailPCT) {
	          CropPCT = AvailPCT;
	      }
	      NatVegPCT = AvailPCT - CropPCT;
              outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = CropPCT;
              outPCTNATVEGdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = NatVegPCT;
              AllPFTs = 0.0;
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  AllPFTs += truncCTSMValues(outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
              }
              AllCFTs = 0.0;
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  AllCFTs += truncCTSMValues(outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
              }
              if (AvailPCT == 0.0) {
                  outPCTPFTdblGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
                  for (pftid = 1; pftid < MAXPFT; pftid++) {
                      outPCTPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
                  outPCTCFTdblGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
                  for (cftid = 1; cftid < MAXCFT; cftid++) {
                      outPCTCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outFERTNITROdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      outUNREPPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outUNREPCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  }
                  outRBIOHVH1dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  outRBIOHVH2dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  outRBIOHSH1dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  outRBIOHSH2dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                  outRBIOHSH3dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              else {
                  if (AllPFTs == 0.0) {
                      outPCTPFTdblGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
                      for (pftid = 1; pftid < MAXPFT; pftid++) {
                          outPCTPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                      }
                  }
                  else {
	              LargestPCTPFT = 0.0;
	              LargestPFT = 0;
	              RescaledPCTPFT = 0.0;
                      for (pftid = 0; pftid < MAXPFT; pftid++) {
                          tempdblPCT = truncCTSMValues(outPCTPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
			  tempdblPCT = truncCTSMValues(tempdblPCT * 100.0 / AllPFTs,100.0,0.01);
                          outPCTPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = tempdblPCT;
		          RescaledPCTPFT += tempdblPCT;
		          if (tempdblPCT > LargestPCTPFT) {
		              LargestPCTPFT = tempdblPCT;
		              LargestPFT = pftid;
                          }
                      }
		      outPCTPFTdblGrid[LargestPFT][ctsmlin * MAXOUTPIX + ctsmpix] += 100.0 - RescaledPCTPFT;
                  }
                  if (AllCFTs == 0.0) {
                      outPCTCFTdblGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
                      for (cftid = 1; cftid < MAXCFT; cftid++) {
                          outPCTCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
                      }
                  }
                  else {
	              LargestPCTCFT = 0.0;
	              LargestCFT = 0;
	              RescaledPCTCFT = 0.0;
                      for (cftid = 0; cftid < MAXCFT; cftid++) {
                          tempdblPCT = truncCTSMValues(outPCTCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix],100.0,0.01);
			  tempdblPCT = truncCTSMValues(tempdblPCT * 100.0 / AllCFTs,100.0,0.01);
                          outPCTCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = tempdblPCT;
		          RescaledPCTCFT += tempdblPCT;
		          if (tempdblPCT > LargestPCTCFT) {
		              LargestPCTCFT = tempdblPCT;
		              LargestCFT = cftid;
                          }
                      }
		      outPCTCFTdblGrid[LargestCFT][ctsmlin * MAXOUTPIX + ctsmpix] += 100.0 - RescaledPCTCFT;
		  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outFERTNITROdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outFERTNITROGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix],100000.0,0.01);
                  }
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      outUNREPPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outUNREPPFTGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix],1.0,0.0001);
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outUNREPCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outUNREPCFTGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix],1.0,0.0001);
                  }
                  outRBIOHVH1dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outRBIOHVH1Grid[ctsmlin * MAXOUTPIX + ctsmpix],100000.0,0.01);
                  outRBIOHVH2dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outRBIOHVH2Grid[ctsmlin * MAXOUTPIX + ctsmpix],100000.0,0.01);
                  outRBIOHSH1dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outRBIOHSH1Grid[ctsmlin * MAXOUTPIX + ctsmpix],100000.0,0.01);
                  outRBIOHSH2dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outRBIOHSH2Grid[ctsmlin * MAXOUTPIX + ctsmpix],100000.0,0.01);
                  outRBIOHSH3dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = truncCTSMValues(outRBIOHSH3Grid[ctsmlin * MAXOUTPIX + ctsmpix],100000.0,0.01);
              }        
              outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC * outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC * outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC * outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC * outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC * outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
              outPCTNATVEGdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = LandFRAC * outPCTNATVEGdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
          }
	  else { 
              outLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outLANDFRACdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
          if (outLANDFRACdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 0.0) {
              outLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTNATVEGdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  outPCTPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outPCTCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outFERTNITROdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  outUNREPPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outUNREPCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              outRBIOHVH1dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outRBIOHVH2dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outRBIOHSH1dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outRBIOHSH2dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outRBIOHSH3dblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
          }
      }
  }
  
  return 0;

}


int swapoceanGrids() {

  double scalelandunits;
  long ctsmlin, ctsmpix;
  int pftid, cftid;
  
  for (ctsmlin = 0; ctsmlin < MAXOUTLIN; ctsmlin++) {
      for (ctsmpix = 0; ctsmpix < MAXOUTPIX; ctsmpix++) {
          if (outLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] == 0.0) {
              outLANDMASKGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
	      outLANDFRACdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
	      outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
	      outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              outPCTPFTdblGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
              for (pftid = 1; pftid < MAXPFT; pftid++) {
                  outPCTPFTdblGrid[pftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              outPCTCFTdblGrid[0][ctsmlin * MAXOUTPIX + ctsmpix] = 100.0;
              for (cftid = 1; cftid < MAXCFT; cftid++) {
                  outPCTCFTdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outFERTNITROdblGrid[cftid][ctsmlin * MAXOUTPIX + ctsmpix] = 0.0;
              }	      
	  }
	  else {
	      scalelandunits = outLANDFRACdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	      outLANDFRACdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = 1.0;
	      outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = scalelandunits * outPCTGLACIERdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	      outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = scalelandunits * outPCTLAKEdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	      outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = scalelandunits * outPCTWETLANDdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] + (1.0 - scalelandunits) * 100.0;
	      outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = scalelandunits * outPCTURBANdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	      outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = scalelandunits * outPCTCROPdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	      outPCTNATVEGdblGrid[ctsmlin * MAXOUTPIX + ctsmpix] = scalelandunits * outPCTNATVEGdblGrid[ctsmlin * MAXOUTPIX + ctsmpix];
	 }
      }
  }
	      
  return 0;
  
}


int writegrids(int currentyear) {

  char outncfilename[1024];
  long ctsmlin, ctsmpix;
  int pftid, cftid;
  char pftidstr[256];
  char cftidstr[256];

  sprintf(outncfilename,"%s/%s_%d.%s.nc",outputdir,outputseries,currentyear,timestamp);
  createncoutputfile(outncfilename);
  openncoutputfile(outncfilename);
  
  writenc1dintfield("natpft",innatpft);
  writenc1dintfield("cft",incft);
  writenc0dfield("EDGEN",&inEDGEN);
  writenc0dfield("EDGEE",&inEDGEE);
  writenc0dfield("EDGES",&inEDGES);
  writenc0dfield("EDGEW",&inEDGEW);
  writenc1dfield("LAT",inLAT);
  writenc2dfield("LATIXY",inLATIXY);
  writenc1dfield("LON",inLON);
  writenc2dfield("LONGXY",inLONGXY);
  writenc2dfield("LANDMASK",outLANDMASKGrid);
  writenc2ddblfield("LANDFRAC",outLANDFRACdblGrid);
  writenc2ddblfield("AREA",outAREAdblGrid);
  writenc2ddblfield("PCT_GLACIER",outPCTGLACIERdblGrid);
  writenc2ddblfield("PCT_LAKE",outPCTLAKEdblGrid);
  writenc2ddblfield("PCT_WETLAND",outPCTWETLANDdblGrid);
  writenc2ddblfield("PCT_URBAN",outPCTURBANdblGrid);
  writenc2ddblfield("PCT_NATVEG",outPCTNATVEGdblGrid);
  writenc2ddblfield("PCT_CROP",outPCTCROPdblGrid);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      writenc3ddblfield("PCT_NAT_PFT",pftid,outPCTPFTdblGrid[pftid]);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      writenc3ddblfield("PCT_CFT",cftid,outPCTCFTdblGrid[cftid]);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {
      writenc3ddblfield("FERTNITRO_CFT",cftid,outFERTNITROdblGrid[cftid]);
  }

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      writenc3ddblfield("UNREPRESENTED_PFT_LULCC",pftid,outUNREPPFTdblGrid[pftid]);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      writenc3ddblfield("UNREPRESENTED_CFT_LULCC",cftid,outUNREPCFTdblGrid[cftid]);
  }

  writenc2ddblfield("HARVEST_VH1",outRBIOHVH1dblGrid);
  writenc2ddblfield("HARVEST_VH2",outRBIOHVH2dblGrid);
  writenc2ddblfield("HARVEST_SH1",outRBIOHSH1dblGrid);
  writenc2ddblfield("HARVEST_SH2",outRBIOHSH2dblGrid);
  writenc2ddblfield("HARVEST_SH3",outRBIOHSH3dblGrid);

  closencfile();
  
  return 0;

}


main(long narg, char **argv) {

  int yearnumber;
    
  if(narg != 2){
        printf("Usage ctsm5landdatatool namelistfile\n");
        return 0;
  }
  
  readnamelist(argv[1]);
  setregionoptions();
  readpftparamfile();
  readcftrawparamfile();
  readcftparamfile();

  createallgrids();

  for (yearnumber = startyear; yearnumber <= endyear; yearnumber++) {
  
      initializeGrids();
      
      readctsmcurrentGrids(yearnumber);
      readctsmLUHforestGrids(yearnumber);
      readctsmLUHpastureGrids(yearnumber);
      readctsmLUHotherGrids(yearnumber);
      readctsmLUHc3annGrids(yearnumber);
      readctsmLUHc4annGrids(yearnumber);
      readctsmLUHc3perGrids(yearnumber);
      readctsmLUHc4perGrids(yearnumber);
      readctsmLUHc3nfxGrids(yearnumber);

      readLUHbasestateGrids(refyear);
  
      readLUHcurrstateGrids(yearnumber);
      readLUHprevdeltastateGrids(yearnumber-1);
  
      readLUHwoodharvestGrids(yearnumber-1);
  
      readUNREPSECDFGrids(yearnumber-1);
      readUNREPSECDNGrids(yearnumber-1);

      readLUHcropmanagementGrids(yearnumber);
      extrapfertGrids();

      generateLUHcollectionGrids();
      generatectsmURBANGrids();
      generatectsmPFTGrids();
      generatectsmCFTGrids();
      generatectsmwoodharvestGrids();
      generatectsmbiohdirectGrids();
      generatectsmfertGrids();
            
      generatedblGrids();
      
      if (includeOcean != 1) {
          swapoceanGrids();
      }
      
      writegrids(yearnumber);

  }
  
  return 1;
  
}
