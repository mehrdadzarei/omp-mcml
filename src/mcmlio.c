/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Input/output of data.
 ****/

#include "mcml.h"

/***********************************************************
 *	Structure used to check against duplicated file names.
 ****/
struct NameList {
  char name[STRLEN];
  struct NameList * next;
};

typedef struct NameList NameNode;
typedef NameNode * NameLink;

/***********************************************************
 *	Get a filename and open it for reading, retry until 
 *	the file can be opened.  '.' terminates the program.
 *      
 *	If Fname != NULL, try Fname first.
 ****/
FILE *GetFile(char *Fname)
{
  FILE * file=NULL;
  Boolean firsttime=1;
  
  do {
    if(firsttime && Fname[0]!='\0') { 
	  /* use the filename from command line */
      firsttime = 0;
    }
    else {
      printf("Input filename(or . to exit):");
      scanf("%s", Fname);
      firsttime = 0;
    }

    if(strlen(Fname) == 1 && Fname[0] == '.') 
      exit(1);			/* exit if no filename entered. */
    
    file = fopen(Fname, "r");
  }  while(file == NULL);
  
  return(file);
}

/***********************************************************
 *	Kill the ith char (counting from 0), push the following 
 *	chars forward by one.
 ****/
void KillChar(size_t i, char * Str)
{
  size_t sl = strlen(Str);
  
  for(;i<sl;i++) Str[i] = Str[i+1];
}

/***********************************************************
 *	Eliminate the chars in a string which are not printing 
 *	chars or spaces.
 *
 *	Spaces include ' ', '\f', '\t' etc.
 *
 *	Return 1 if no nonprinting chars found, otherwise 
 *	return 0.
 ****/
Boolean CheckChar(char * Str)
{
  Boolean found = 0;	/* found bad char. */
  size_t sl = strlen(Str);
  size_t i=0;
  
  while(i<sl) 
    if (Str[i]<0 || Str[i]>255)
      nrerror("Non-ASCII file\n");
    else if(isprint(Str[i]) || isspace(Str[i])) 
      i++;
    else {
      found = 1;
      KillChar(i, Str);
      sl--;
    }
  
  return(found);	
}

/***********************************************************
 *	Return 1 if this line is a comment line in which the 
 *	first non-space character is "#".
 *
 *	Also return 1 if this line is space line.
 ****/
Boolean CommentLine(char *Buf)
{
  size_t spn, cspn;
  
  spn = strspn(Buf, " \t");	
  /* length spanned by space or tab chars. */

  cspn = strcspn(Buf, "#\n");
  /* length before the 1st # or return. */

  if(spn == cspn) 	/* comment line or space line. */
	return(1);
  else				/* the line has data. */	 
	return(0);		
}

/***********************************************************
 *	Skip space or comment lines and return a data line only.
 ****/
char * FindDataLine(FILE *File_Ptr)
{
  char buf[STRLEN];
  
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  return(buf);
}

/***********************************************************
 *	Skip file version, then read number of runs.
 ****/
short ReadNumRuns(FILE* File_Ptr)
{
  char buf[STRLEN];
  short n=0;
  
  FindDataLine(File_Ptr); /* skip file version. */
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') nrerror("Reading number of runs\n");
  sscanf(buf, "%hd",&n);
  return(n);
}
  
/***********************************************************
 *	Read the file name and the file format.
 *
 *	The file format can be either A for ASCII or B for
 *	binary.
 ****/
void ReadFnameFormat(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in file name and format. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading file name and format.\n");
  sscanf(buf, "%s %c", In_Ptr->out_fname, &(In_Ptr->out_fformat) );
  if(toupper(In_Ptr->out_fformat) != 'B') 
	In_Ptr->out_fformat = 'A';
}

/***********************************************************
 *	Read the number of photons.
 ****/
void ReadNumPhotons(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of photons. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading number of photons.\n");
  sscanf(buf, "%ld", &In_Ptr->num_photons);
  if(In_Ptr->num_photons<=0) 
	nrerror("Nonpositive number of photons.\n");
}

void ReadNumThreads(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of photons. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading number of threads.\n");
  sscanf(buf, "%hd", &In_Ptr->num_threads);
  if(In_Ptr->num_threads<=0) 
	nrerror("Nonpositive number of threads.\n");
}

/***********************************************************
 *	Read the seed for random generation.
 ****/
void ReadSeed(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of photons. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading seed for random generation.\n");
  sscanf(buf, "%ld", &In_Ptr->seed);
  if(In_Ptr->seed<0) 
	nrerror("Nonpositive seed for random generation.\n");
}

/***********************************************************
 *	Read the members dz and dr.
 ****/
void ReadDzDrDxDy(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in dz, dr. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') nrerror("Reading dx, dy, dz, dr.\n");
  sscanf(buf, "%lf%lf%lf%lf", &In_Ptr->dx, &In_Ptr->dy, &In_Ptr->dz, &In_Ptr->dr);
  if(In_Ptr->dx<=0) nrerror("Nonpositive dx.\n");
  if(In_Ptr->dy<=0) nrerror("Nonpositive dy.\n");
  if(In_Ptr->dz<=0) nrerror("Nonpositive dz.\n");
  if(In_Ptr->dr<=0) nrerror("Nonpositive dr.\n");
}

/***********************************************************
 *	Read the members nz, nr, na.
 ****/
void ReadxyzNa(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];
  double x, y, z;
  
  /** read in number of dz, dr, da. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading number of dx, dy, dz, da's.\n");
  sscanf(buf, "%lf%lf%lf%hd", &x, &y, &z, &In_Ptr->na);
  In_Ptr->nx =  (int)(x/In_Ptr->dx);
  In_Ptr->ny =  (int)(y/In_Ptr->dy);
  In_Ptr->nz =  (int)(z/In_Ptr->dz);
  In_Ptr->nr =  (int)(sqrt((x*x+y*y))/In_Ptr->dr);
  if(In_Ptr->nx<=0) 
      nrerror("Nonpositive number of dx's.\n");
  if(In_Ptr->ny<=0) 
      nrerror("Nonpositive number of dy's.\n");
  if(In_Ptr->nz<=0) 
      nrerror("Nonpositive number of dz's.\n");
  if(In_Ptr->nr<=0) 
	nrerror("Nonpositive number of dr's.\n");
  if(In_Ptr->na<=0) 
	nrerror("Nonpositive number of da's.\n");
  In_Ptr->da = 0.5*PI/In_Ptr->na;
}

void ReadSourcePosition(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in srcx, srcy, srcz. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') nrerror("Reading source position.\n");
  sscanf(buf, "%lf%lf%lf", &In_Ptr->srcx, &In_Ptr->srcy, &In_Ptr->srcz);  
 }

void ReadSourceDirection(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in source direction. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') nrerror("Reading source direction.\n");
  sscanf(buf, "%lf%lf%lf", &In_Ptr->dirx, &In_Ptr->diry, &In_Ptr->dirz);
  if(In_Ptr->dirx<-1 || In_Ptr->dirx>1) nrerror("source direction is beetwin -1 and 1.\n");  
  if(In_Ptr->diry<-1 || In_Ptr->diry>1) nrerror("source direction is beetwin -1 and 1.\n");  
  if(In_Ptr->dirz<-1 || In_Ptr->dirz>1) nrerror("source direction is beetwin -1 and 1.\n");  
 }
 
void ReadSourceType(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in source type. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') nrerror("Reading source type.\n");
  sscanf(buf, "%d%lf%lf", &In_Ptr->srctyp, &In_Ptr->srcpar1, &In_Ptr->srcpar2);
  if(In_Ptr->srctyp<0) nrerror("Nonpositive source type.\n");
  if(In_Ptr->srcpar1<0.0) nrerror("Nonpositive source parameters.\n");
  if(In_Ptr->srcpar2<0.0) nrerror("Nonpositive source parameters.\n");
 }

void ReadNumDetec(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of detectors and raduis. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading number of detectors and raduis.\n");
  sscanf(buf, "%hd%lf%lf", &In_Ptr->num_detectors, &In_Ptr->rad_detectors, &In_Ptr->NA_detectors);
  if(In_Ptr->num_detectors<0) 
	nrerror("Nonpositive number of detectors.\n");
  if(In_Ptr->rad_detectors<=0) 
	nrerror("Nonpositive raduis of detectors.\n");
  if(In_Ptr->NA_detectors<0) 
	nrerror("Nonpositive Numerical Aperture of detectors.\n");
}

/***********************************************************
 *	Read the position of detectors.
 *
 ****/
Boolean ReadOneDetector(FILE *File_Ptr, LayerStruct * Detector_Ptr, InputStruct *In_Ptr, short * index_detector){
  char buf[STRLEN], msg[STRLEN];  

  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') return(1);	/* error. */

  sscanf(buf, "%lf%lf", &Detector_Ptr->xD, &Detector_Ptr->yD);
  if(Detector_Ptr->xD<0 || Detector_Ptr->yD<0)
    return(1);			/* error. */
  
  return(0);
}

/***********************************************************
 *	Read the position of detectors at a time.
 ****/
void ReadDetectors(FILE *File_Ptr, short Num_Detectors, LayerStruct ** Detectors_PP, InputStruct *In_Ptr){
  char msg[STRLEN];
  short i;

  /* Allocate an array for the detectors parameters. */
  *Detectors_PP = (LayerStruct *)
    malloc((unsigned) (Num_Detectors)*sizeof(LayerStruct));
  if (!(*Detectors_PP))
    nrerror("allocation failure in ReadDetectors()");

  for(i=0; i<Num_Detectors; i++)  
    if(ReadOneDetector(File_Ptr, &((*Detectors_PP)[i]), In_Ptr, &i)){
      sprintf(msg, "Error reading %hd of %hd detectors\n", i, In_Ptr->num_layers);
	  nrerror(msg);
    }
}

/***********************************************************
 *	Read the number of layers.
 ****/
void ReadNumLayers(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of layers. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading number of layers.\n");
  sscanf(buf, "%hd", &In_Ptr->num_layers);
  if(In_Ptr->num_layers<=0) 
	nrerror("Nonpositive number of layers.\n");
}

/***********************************************************
 *	Read the refractive index n of the ambient.
 ****/
void ReadAmbient(FILE *File_Ptr, LayerStruct * Layer_Ptr, char *side)
{
  char buf[STRLEN], msg[STRLEN];
  double n;

  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') {
    sprintf(msg, "Rading n of %s ambient.\n", side);
	nrerror(msg);
  }

  sscanf(buf, "%lf", &n );
  if(n<=0) nrerror("Wrong n.\n");
  Layer_Ptr->n = n;
}

/***********************************************************
 *	Read the parameters of one layer.
 *
 *	Return 1 if error detected.
 *	Return 0 otherwise.
 *
 *	*Z_Ptr is the z coordinate of the current layer, which
 *	is used to convert thickness of layer to z coordinates
 *	of the two boundaries of the layer.
 ****/
Boolean ReadOneLayer(FILE *File_Ptr, LayerStruct * Layer_Ptr, double *Z_Ptr)
{
  char buf[STRLEN], msg[STRLEN];
  double d, n, mua, mus, g;	/* d is thickness. */

  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') return(1);	/* error. */

  sscanf(buf, "%lf%lf%lf%lf%lf", &n, &mua, &mus, &g, &d);
  if(d<0 || n<=0 || mua<0 || mus<0 || g<0 || g>1) 
    return(1);			/* error. */
    
  Layer_Ptr->n	= n;
  Layer_Ptr->mua = mua;	
  Layer_Ptr->mus = mus;	
  Layer_Ptr->g   = g;
  Layer_Ptr->z0	= *Z_Ptr;
  *Z_Ptr += d;
  Layer_Ptr->z1	= *Z_Ptr;

  return(0);
}

/***********************************************************
 *	Read the parameters of one layer at a time.
 ****/
void ReadLayerSpecs(FILE *File_Ptr, short Num_Layers, LayerStruct ** Layerspecs_PP)
{
  char msg[STRLEN];
  short i=0;
  double z = 0.0;	/* z coordinate of the current layer. */
  
  /* Allocate an array for the layer parameters. */
  /* layer 0 and layer Num_Layers + 1 are for ambient. */
  *Layerspecs_PP = (LayerStruct *)
    malloc((unsigned) (Num_Layers+2)*sizeof(LayerStruct));
  if (!(*Layerspecs_PP)) 
    nrerror("allocation failure in ReadLayerSpecs()");
  
  ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "top"); 
  for(i=1; i<=Num_Layers; i++)  
    if(ReadOneLayer(File_Ptr, &((*Layerspecs_PP)[i]), &z)) {
      sprintf(msg, "Error reading %hd of %hd layers\n", i, Num_Layers);
	  nrerror(msg);
    }
  ReadAmbient(File_Ptr, &((*Layerspecs_PP)[i]), "bottom"); 
}

/***********************************************************
 *	Compute the critical angles for total internal
 *	reflection according to the relative refractive index
 *	of the layer.
 *	All layers are processed.
 ****/
void CriticalAngle( short Num_Layers, LayerStruct ** Layerspecs_PP)
{
  short i=0;
  double n1, n2;
  
  for(i=1; i<=Num_Layers; i++)  {
    n1 = (*Layerspecs_PP)[i].n;
    n2 = (*Layerspecs_PP)[i-1].n;
    (*Layerspecs_PP)[i].cos_crit0 = n1>n2 ? 
		sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
    
    n2 = (*Layerspecs_PP)[i+1].n;
    (*Layerspecs_PP)[i].cos_crit1 = n1>n2 ? 
		sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
  }
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void ReadRsp(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRsp);
  if(In_Ptr->isRsp<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void ReadRd_ra(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRd_ra);
  if(In_Ptr->isRd_ra<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void ReadRd_lr(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRd_lr);
  if(In_Ptr->isRd_lr<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void ReadRd_r(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRd_r);
  if(In_Ptr->isRd_r<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void ReadRd_d(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRd_d);
  if(In_Ptr->isRd_d<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void  ReadRd_a(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRd_a);
  if(In_Ptr->isRd_a<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state of diffuse reflectance.
 ****/
void ReadRd(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of diffuse reflectance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of diffuse reflectance.\n");
  sscanf(buf, "%ld", &In_Ptr->isRd);
  if(In_Ptr->isRd<0) 
	nrerror("Nonpositive state of diffuse reflectance.\n");
}

/***********************************************************
 *	Read the state absorption.
 ****/
void ReadA_xyz(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state absorption. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state absorption.\n");
  sscanf(buf, "%ld", &In_Ptr->isA_xyz);
  if(In_Ptr->isA_xyz<0) 
	nrerror("Nonpositive state absorption.\n");
}

/***********************************************************
 *	Read the state absorption.
 ****/
void ReadA_rz(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state absorption. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state absorption.\n");
  sscanf(buf, "%ld", &In_Ptr->isA_rz);
  if(In_Ptr->isA_rz<0) 
	nrerror("Nonpositive state absorption.\n");
}

/***********************************************************
 *	Read the state absorption.
 ****/
void ReadA_z(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state absorption. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state absorption.\n");
  sscanf(buf, "%ld", &In_Ptr->isA_z);
  if(In_Ptr->isA_z<0) 
	nrerror("Nonpositive state absorption.\n");
}

/***********************************************************
 *	Read the state absorption.
 ****/
void ReadA_l(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state absorption. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state absorption.\n");
  sscanf(buf, "%ld", &In_Ptr->isA_l);
  if(In_Ptr->isA_l<0) 
	nrerror("Nonpositive state absorption.\n");
}

/***********************************************************
 *	Read the state of total absorption.
 ****/
void ReadA(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of total absorption. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of total absorption.\n");
  sscanf(buf, "%ld", &In_Ptr->isA);
  if(In_Ptr->isA<0) 
	nrerror("Nonpositive state of total absorption.\n");
}

/***********************************************************
 *	Read the state transmittance.
 ****/
void ReadTt_ra(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state transmittance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state transmittance.\n");
  sscanf(buf, "%ld", &In_Ptr->isTt_ra);
  if(In_Ptr->isTt_ra<0) 
	nrerror("Nonpositive state transmittance.\n");
}

/***********************************************************
 *	Read the state transmittance.
 ****/
void ReadTt_r(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of photons. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state transmittance.\n");
  sscanf(buf, "%ld", &In_Ptr->isTt_r);
  if(In_Ptr->isTt_r<0) 
	nrerror("Nonpositive state transmittance.\n");
}

/***********************************************************
 *	Read the state transmittance.
 ****/
void ReadTt_d(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of photons. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state transmittance.\n");
  sscanf(buf, "%ld", &In_Ptr->isTt_d);
  if(In_Ptr->isTt_d<0) 
	nrerror("Nonpositive state transmittance.\n");
}

/***********************************************************
 *	Read the state transmittance.
 ****/
void ReadTt_a(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state transmittance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state transmittance.\n");
  sscanf(buf, "%ld", &In_Ptr->isTt_a);
  if(In_Ptr->isTt_a<0) 
	nrerror("Nonpositive state transmittance.\n");
}

/***********************************************************
 *	Read the state total transmittance.
 ****/
void ReadTt(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state total transmittance. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state total transmittance.\n");
  sscanf(buf, "%ld", &In_Ptr->isTt);
  if(In_Ptr->isTt<0) 
	nrerror("Nonpositive stat total transmittance.\n");
}

/***********************************************************
 *	Read the state of tracing.
 ****/
void ReadBannana(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in state of tracing. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading state of tracing.\n");
  sscanf(buf, "%ld", &In_Ptr->isBannana);
  if(In_Ptr->isBannana<0) 
	nrerror("Nonpositive state of tracing.\n");
}

/***********************************************************
 *	Read the state of PPL.
 ****/
void ReadPPL(FILE *File_Ptr, InputStruct *In_Ptr)
{
  char buf[STRLEN];

  /** read in number of photons. **/
  buf[0] = '\0';
  do {	/* skip space or comment lines. */    
    if(fgets(buf, 255, File_Ptr) == NULL)  {
      printf("Incomplete data\n");
      buf[0]='\0';
      break;
    }
    else        
      CheckChar(buf);
  } while(CommentLine(buf));
  // strcpy(buf, FindDataLine(File_Ptr));
  if(buf[0]=='\0') 
	nrerror("Reading is PPL.\n");
  sscanf(buf, "%d", &In_Ptr->isPPL);
  if(In_Ptr->isPPL<0) 
	nrerror("Nonpositive for PPL.\n");
}

/***********************************************************
 *	Read in the input parameters for one run.
 ****/
void ReadParm(FILE* File_Ptr, InputStruct * In_Ptr)
{
  char buf[STRLEN];
  
  In_Ptr->Wth = WEIGHT;

  ReadFnameFormat(File_Ptr, In_Ptr);
  ReadNumPhotons(File_Ptr, In_Ptr);
  ReadNumThreads(File_Ptr, In_Ptr);
  ReadSeed(File_Ptr, In_Ptr);
  ReadDzDrDxDy(File_Ptr, In_Ptr);
  ReadxyzNa(File_Ptr, In_Ptr);
  ReadSourcePosition(File_Ptr, In_Ptr);
  ReadSourceDirection(File_Ptr, In_Ptr);
  ReadSourceType(File_Ptr, In_Ptr);
  ReadNumDetec(File_Ptr, In_Ptr);
  ReadDetectors(File_Ptr, In_Ptr->num_detectors, &In_Ptr->detectors, In_Ptr);
  ReadNumLayers(File_Ptr, In_Ptr);
  ReadLayerSpecs(File_Ptr, In_Ptr->num_layers, &In_Ptr->layerspecs);
  CriticalAngle(In_Ptr->num_layers, &In_Ptr->layerspecs);
  /* Read state of outputs */
  ReadRsp(File_Ptr, In_Ptr);
  ReadRd_ra(File_Ptr, In_Ptr);
  ReadRd_lr(File_Ptr, In_Ptr);
  ReadRd_r(File_Ptr, In_Ptr);
  ReadRd_d(File_Ptr, In_Ptr);
  ReadRd_a(File_Ptr, In_Ptr);
  ReadRd(File_Ptr, In_Ptr);
  ReadA_xyz(File_Ptr, In_Ptr);
  ReadA_rz(File_Ptr, In_Ptr);
  ReadA_z(File_Ptr, In_Ptr);
  ReadA_l(File_Ptr, In_Ptr);
  ReadA(File_Ptr, In_Ptr);
  ReadTt_ra(File_Ptr, In_Ptr);
  ReadTt_r(File_Ptr, In_Ptr);
  ReadTt_d(File_Ptr, In_Ptr);
  ReadTt_a(File_Ptr, In_Ptr);
  ReadTt(File_Ptr, In_Ptr);
  ReadBannana(File_Ptr, In_Ptr);
  ReadPPL(File_Ptr, In_Ptr);
}

/***********************************************************
 *	Return 1, if the name in the name list.
 *	Return 0, otherwise.
 ****/
Boolean NameInList(char *Name, NameLink List)
{
  while (List != NULL) {
	if(strcmp(Name, List->name) == 0) 
	  return(1);
	List = List->next;
  };
  return(0);
}

/***********************************************************
 *	Add the name to the name list.
 ****/
void AddNameToList(char *Name, NameLink * List_Ptr)
{
  NameLink list = *List_Ptr;

  if(list == NULL) {	/* first node. */
	*List_Ptr = list = (NameLink)malloc(sizeof(NameNode));
	strcpy(list->name, Name);
	list->next = NULL;
  }
  else {				/* subsequent nodes. */
	/* Move to the last node. */
	while(list->next != NULL)
	  list = list->next;

	/* Append a node to the list. */
	list->next = (NameLink)malloc(sizeof(NameNode));
	list = list->next;
	strcpy(list->name, Name);
	list->next = NULL;
  }
}

/***********************************************************
 *	Check against duplicated file names.
 *
 *	A linked list is set up to store the file names used
 *	in this input data file.
 ****/
Boolean FnameTaken(char *fname, NameLink * List_Ptr)
{
  if(NameInList(fname, *List_Ptr))
	return(1);
  else {
	AddNameToList(fname, List_Ptr);
	return(0);
  }
}

/***********************************************************
 *	Free each node in the file name list.
 ****/
void FreeFnameList(NameLink List)
{
  NameLink next;

  while(List != NULL) {
	next = List->next;
	free(List);
	List = next;
  }
}

/***********************************************************
 *	Check the input parameters for each run.
 ****/
void CheckParm(FILE* File_Ptr, InputStruct * In_Ptr)
{
  short i_run;
  short num_runs;	/* number of independent runs. */
  NameLink head = NULL;
  Boolean name_taken;/* output files share the same file name.*/
  char msg[STRLEN];
  
  num_runs = ReadNumRuns(File_Ptr);
  for(i_run=1; i_run<=num_runs; i_run++)  {
    printf("Checking input data for run %hd\n", i_run);
    ReadParm(File_Ptr, In_Ptr);

	name_taken = FnameTaken(In_Ptr->out_fname, &head);
	if(name_taken) 
	  sprintf(msg, "file name %s duplicated.\n", In_Ptr->out_fname);

    free(In_Ptr->layerspecs);
	if(name_taken) nrerror(msg);
  }
  FreeFnameList(head);
  rewind(File_Ptr);
}

/***********************************************************
 *	Allocate the arrays in OutStruct for one run, and 
 *	array elements are automatically initialized to zeros.
 ****/
void InitOutputData(InputStruct In_Parm, OutStruct * Out_Ptr, Boolean m)
{
  short nx = In_Parm.nx;
  short ny = In_Parm.ny;
  short nz = In_Parm.nz;
  short nr = In_Parm.nr;
  short na = In_Parm.na;
  short nd = In_Parm.num_detectors;
  short nl = In_Parm.num_layers;	/* remember to use nl+2 because of 2 for ambient. */
  
  if(nz<=0 || nr<=0 || na<=0 || nl<=0 || nx<=0 || ny<=0 || nd<=0) 
    nrerror("Wrong grid parameters.\n");
  
  /* Init pure numbers. */
  if(m){    
    if(In_Parm.isRd)
        Out_Ptr->Rd  = 0.0;
    if(In_Parm.isA)
        Out_Ptr->A   = 0.0;
    if(In_Parm.isTt)
        Out_Ptr->Tt  = 0.0;
    Out_Ptr->nthread = 0;
  }
  Out_Ptr->Rsp = 0.0;
  
  /* Allocate the arrays and the matrices. */
  if(In_Parm.isRd_ra || In_Parm.isRd_r || In_Parm.isRd_a || In_Parm.isRd)
    Out_Ptr->Rd_ra = AllocMatrix(0,nr-1,0,na-1);
  if(In_Parm.isRd_lr)
    Out_Ptr->Rd_lr  = AllocMatrix(1,nl,0,nr-1);
  if(m){
    if(In_Parm.isRd_r)
        Out_Ptr->Rd_r  = AllocVector(0,nr-1);
    if(In_Parm.isRd_a)
        Out_Ptr->Rd_a  = AllocVector(0,na-1);
  }
  if(In_Parm.isRd_d)
    Out_Ptr->Rd_d  = AllocVector(0,nd-1);
  
  if(In_Parm.isA_xyz)
    Out_Ptr->A_xyz  = AllocCubic(-nx,nx-1,-ny,ny-1,0,nz-1);
  if(In_Parm.isA_rz || In_Parm.isA_z || In_Parm.isA_l || In_Parm.isA)
    Out_Ptr->A_rz  = AllocMatrix(0,nr-1,0,nz-1);
  if(m){
    if(In_Parm.isA_z)
        Out_Ptr->A_z   = AllocVector(0,nz-1);
    if(In_Parm.isA_l)
        Out_Ptr->A_l   = AllocVector(0,nl+1);
  }
        
  if(In_Parm.isTt_ra || In_Parm.isTt_r || In_Parm.isTt_a || In_Parm.isTt)
    Out_Ptr->Tt_ra = AllocMatrix(0,nr-1,0,na-1);
  if(m){
    if(In_Parm.isTt_r)
        Out_Ptr->Tt_r  = AllocVector(0,nr-1);
    if(In_Parm.isTt_a)
        Out_Ptr->Tt_a  = AllocVector(0,na-1);
  }
  if(In_Parm.isTt_d)
    Out_Ptr->Tt_d  = AllocVector(0,nd-1);
  
  if(In_Parm.isBannana){
    Out_Ptr->Bannana  = AllocCubic(0,nd-1,-nx,nx-1,0,nz-1);
    if(!m)
        Out_Ptr->Bannana_inner  = AllocMatrix(-nx,nx-1,0,nz-1);
  }
  
  if(In_Parm.isPPL){
    Out_Ptr->Li  = AllocMatrix(0, nd, 1,nl); 
    if(!m)
        Out_Ptr->li  = AllocVector(1,nl);
    Out_Ptr->num_Ps  = (long *) calloc((unsigned) (nd), sizeof(long));
  }
}

/***********************************************************
 *	Undo what InitOutputData did.
 *  i.e. free the data allocations.
 ****/
void FreeData(InputStruct In_Parm, OutStruct * Out_Ptr, short i)
{
  short nx = In_Parm.nx;
  short ny = In_Parm.ny;
  short nz = In_Parm.nz;
  short nr = In_Parm.nr;
  short na = In_Parm.na;
  short nd = In_Parm.num_detectors;
  short nl = In_Parm.num_layers;
  /* remember to use nl+2 because of 2 for ambient. */
  
  if(i)
    free(In_Parm.layerspecs);
  
  if(In_Parm.isRd_ra || In_Parm.isRd_r || In_Parm.isRd_a || In_Parm.isRd)
    FreeMatrix(Out_Ptr->Rd_ra, 0,nr-1,0,na-1);
  if(In_Parm.isRd_lr)
    FreeMatrix(Out_Ptr->Rd_lr, 1,nl,0,nr-1);
  if(i){
    if(In_Parm.isRd_r)
        FreeVector(Out_Ptr->Rd_r, 0,nr-1);
    if(In_Parm.isRd_a)
        FreeVector(Out_Ptr->Rd_a, 0,na-1);
  }
  if(In_Parm.isRd_d)
    FreeVector(Out_Ptr->Rd_d, 0,nd-1);
  
  if(In_Parm.isA_xyz)
    FreeCubic(Out_Ptr->A_xyz, -nx, nx-1, -ny, ny-1, 0,nz-1);
  if(In_Parm.isA_rz || In_Parm.isA_z || In_Parm.isA_l || In_Parm.isA)
    FreeMatrix(Out_Ptr->A_rz, 0, nr-1, 0,nz-1);
  if(i){
    if(In_Parm.isA_z)
        FreeVector(Out_Ptr->A_z, 0, nz-1);
    if(In_Parm.isA_l)
        FreeVector(Out_Ptr->A_l, 0,nl+1);
  }
  if(In_Parm.isTt_ra || In_Parm.isTt_r || In_Parm.isTt_a || In_Parm.isTt)
    FreeMatrix(Out_Ptr->Tt_ra, 0,nr-1,0,na-1);
  if(i){
    if(In_Parm.isTt_r)
        FreeVector(Out_Ptr->Tt_r, 0,nr-1);
    if(In_Parm.isTt_a)
        FreeVector(Out_Ptr->Tt_a, 0,na-1);
  }
  if(In_Parm.isTt_d)
      FreeVector(Out_Ptr->Tt_d, 0,nd-1);
  if(In_Parm.isPPL){
    free(Out_Ptr->num_Ps);
    if(!i)
        FreeVector(Out_Ptr->li, 1,nl);    
    FreeMatrix(Out_Ptr->Li, 0, nd, 1,nl);
  }
  if(In_Parm.isBannana){
    if(!i)
        FreeMatrix(Out_Ptr->Bannana_inner, -nx,nx-1,0,nz-1);
    FreeCubic(Out_Ptr->Bannana, 0, nd-1, -nx, nx-1, 0, nz-1);
  }
}

/***********************************************************
 *	Get 1D array elements by summing the 2D array elements.
 ****/
void Sum2DRd(InputStruct In_Parm, OutStruct * Out_Ptr, OutStruct * Out_Ptr_Inner)
{
  short nr = In_Parm.nr;
  short na = In_Parm.na;
  short ir,ia, i;
  short num_layers = In_Parm.num_layers;
  double sum;
  
  Out_Ptr->Rsp += Out_Ptr_Inner->Rsp;
  
  if(In_Parm.isRd_ra || In_Parm.isRd_r || In_Parm.isRd_a || In_Parm.isRd)
      for(ir=0; ir<nr; ir++)
          for(ia=0; ia<na; ia++)
              Out_Ptr->Rd_ra[ir][ia] += Out_Ptr_Inner->Rd_ra[ir][ia];
  if(In_Parm.isRd_lr)
      for(i=1; i<=num_layers; i++)
        for(ir=0; ir<nr; ir++)
            Out_Ptr->Rd_lr[i][ir] += Out_Ptr_Inner->Rd_lr[i][ir];
  if(In_Parm.isRd_d)
    for(i=0;i<In_Parm.num_detectors;i++)
        Out_Ptr->Rd_d[i] += Out_Ptr_Inner->Rd_d[i];
  if(In_Parm.isRd_r)
    for(ir=0; ir<nr; ir++){
        sum = 0.0;
        for(ia=0; ia<na; ia++) sum += Out_Ptr_Inner->Rd_ra[ir][ia];
        Out_Ptr->Rd_r[ir] += sum;
    }
  if(In_Parm.isRd_a)
    for(ia=0; ia<na; ia++)  {
        sum = 0.0;
        for(ir=0; ir<nr; ir++) sum += Out_Ptr_Inner->Rd_ra[ir][ia];
        Out_Ptr->Rd_a[ia] += sum;
    }
  if(In_Parm.isRd){
    sum = 0.0;
    for(ir=0; ir<nr; ir++)
        for(ia=0; ia<na; ia++)
            sum += Out_Ptr_Inner->Rd_ra[ir][ia];
    Out_Ptr->Rd += sum;
  }
}

/***********************************************************
 *	Return the index to the layer according to the index
 *	to the grid line system in z direction (Iz).
 *
 *	Use the center of box.
 ****/
short IzToLayer(short Iz, InputStruct In_Parm)
{
  short i=1;	/* index to layer. */
  short num_layers = In_Parm.num_layers;
  double dz = In_Parm.dz;
  
  while( (Iz+0.5)*dz >= In_Parm.layerspecs[i].z1 
	&& i<num_layers) i++;
  
  return(i);
}

/***********************************************************
 *	Get 1D array elements by summing the 2D array elements.
 ****/
void Sum2DA(InputStruct In_Parm, OutStruct * Out_Ptr, OutStruct * Out_Ptr_Inner)
{
  short nz = In_Parm.nz;
  short nr = In_Parm.nr;
  short nx = In_Parm.nx;
  short ny = In_Parm.ny;
  short iz,ix,iy,ir,nthrds;
  double sum;
  
  if(In_Parm.isA_xyz)
      for(iz=0;iz<nz;iz++)
          for(iy=-ny;iy<ny;iy++)
              for(ix=-nx;ix<nx;ix++)
                  Out_Ptr->A_xyz[ix][iy][iz] += Out_Ptr_Inner->A_xyz[ix][iy][iz];
  if(In_Parm.isA_rz || In_Parm.isA_z || In_Parm.isA_l || In_Parm.isA)
      for(iz=0; iz<nz; iz++)
        for(ir=0; ir<nr; ir++)
            Out_Ptr->A_rz[ir][iz]+= Out_Ptr_Inner->A_rz[ir][iz];
  if(In_Parm.isA_z)
    for(iz=0; iz<nz; iz++){
        sum = 0.0;
        for(ir=0; ir<nr; ir++) sum += Out_Ptr_Inner->A_rz[ir][iz];
        Out_Ptr->A_z[iz] += sum;
    }
  if(In_Parm.isA_l || In_Parm.isA){
    double sumaz;
    sum = 0.0;
    for(iz=0; iz<nz; iz++){
        sumaz = 0.0;
        for(ir=0; ir<nr; ir++) sumaz += Out_Ptr_Inner->A_rz[ir][iz];        
        sum += sumaz;
        if(In_Parm.isA_l)
            Out_Ptr->A_l[IzToLayer(iz, In_Parm)] += sumaz;
    }
    if(In_Parm.isA)
        Out_Ptr->A += sum;
  }
}

/***********************************************************
 *	Get 1D array elements by summing the 2D array elements.
 ****/
void Sum2DTt(InputStruct In_Parm, OutStruct * Out_Ptr, OutStruct * Out_Ptr_Inner)
{
  short nr = In_Parm.nr;
  short na = In_Parm.na;
  short ir,ia;
  double sum;
  
  if(In_Parm.isTt_ra || In_Parm.isTt_r || In_Parm.isTt_a || In_Parm.isTt)
      for(ir=0; ir<nr; ir++)
          for(ia=0; ia<na; ia++)
              Out_Ptr->Tt_ra[ir][ia] += Out_Ptr_Inner->Tt_ra[ir][ia];
  if(In_Parm.isTt_d)
    for(short i=0;i<In_Parm.num_detectors;i++)
        Out_Ptr->Tt_d[i] += Out_Ptr_Inner->Tt_d[i];
  if(In_Parm.isTt_r)
    for(ir=0; ir<nr; ir++){
        sum = 0.0;
        for(ia=0; ia<na; ia++) sum += Out_Ptr_Inner->Tt_ra[ir][ia];
        Out_Ptr->Tt_r[ir] += sum;
    }
  if(In_Parm.isTt_a)
    for(ia=0; ia<na; ia++){
        sum = 0.0;
        for(ir=0; ir<nr; ir++) sum += Out_Ptr_Inner->Tt_ra[ir][ia];
        Out_Ptr->Tt_a[ia] += sum;
    }
  if(In_Parm.isTt){
    sum = 0.0;
    for(ir=0; ir<nr; ir++)
        for(ia=0; ia<na; ia++)
            sum += Out_Ptr_Inner->Tt_ra[ir][ia];
    Out_Ptr->Tt += sum;
  }
}

void SumLi(InputStruct In_Parm, OutStruct * Out_Ptr, OutStruct * Out_Ptr_Inner)
{ 
    short index_layer, i;
    short num_layers = In_Parm.num_layers;
    short num_det = In_Parm.num_detectors;

	for(i = 0; i < num_det; i++){
        for(index_layer = 1; index_layer<= num_layers; index_layer++)
            Out_Ptr->Li[i][index_layer] += Out_Ptr_Inner->Li[i][index_layer];
        Out_Ptr->num_Ps[i] += Out_Ptr_Inner->num_Ps[i];
	}
}

void SumBannana(InputStruct In_Parm, OutStruct *Out_Ptr, OutStruct *Out_Ptr_Inner){
  short nx = In_Parm.nx;
  short nz = In_Parm.nz;
  short nd = In_Parm.num_detectors;
  short ix, iz, i;

  for(iz = 0; iz < nz; iz++)
     for(ix = -nx; ix < nx; ix++)
        for( i = 0; i < nd; i++)
           Out_Ptr->Bannana[i][ix][iz] += Out_Ptr_Inner->Bannana[i][ix][iz];
}

/***********************************************************
 *	Scale Rd and Tt properly.
 *
 *	"a" stands for angle alpha.
 ****
 *	Scale Rd(r,a) and Tt(r,a) by
 *      (area perpendicular to photon direction)
 *		x(solid angle)x(No. of photons).
 *	or
 *		[2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
 *	or
 *		[2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
 ****
 *	Scale Rd(r) and Tt(r) by
 *		(area on the surface)x(No. of photons).
 ****
 *	Scale Rd(a) and Tt(a) by
 *		(solid angle)x(No. of photons).
 ****/
void ScaleRdTt(InputStruct In_Parm, OutStruct *	Out_Ptr)
{
  short nr = In_Parm.nr;
  short na = In_Parm.na;
  double dr = In_Parm.dr;
  double da = In_Parm.da;
  short nl = In_Parm.num_layers;
  short ir,ia;
  double scale1, scale2;
  
  Out_Ptr->Rsp /= In_Parm.num_photons;
  if(In_Parm.isTt_ra || In_Parm.isRd_ra){
    scale1 = 4.0*PI*PI*dr*sin(da/2)*dr*In_Parm.num_photons;
	/* The factor (ir+0.5)*sin(2a) to be added. */

    for(ir=0; ir<nr; ir++)  
        for(ia=0; ia<na; ia++) {
            scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);
            if(In_Parm.isRd_ra)
                Out_Ptr->Rd_ra[ir][ia] *= scale2;
            if(In_Parm.isTt_ra)
                Out_Ptr->Tt_ra[ir][ia] *= scale2;
        }
  }
  if(In_Parm.isTt_r || In_Parm.isRd_r || In_Parm.isRd_lr){
    scale1 = 2.0*PI*dr*dr*In_Parm.num_photons;  
	/* area is 2*PI*[(ir+0.5)*dr]*dr. ir+0.5 to be added. */

    for(ir=0; ir<nr; ir++) {
        scale2 = 1.0/((ir+0.5)*scale1);
        if(In_Parm.isRd_r)
            Out_Ptr->Rd_r[ir] *= scale2;
        if(In_Parm.isTt_r)
            Out_Ptr->Tt_r[ir] *= scale2;
        if(In_Parm.isRd_lr)
            for(short i=1; i<=nl; i++)
                Out_Ptr->Rd_lr[i][ir] *= scale2;
    }
  }
  if(In_Parm.isTt_a || In_Parm.isRd_a){
    scale1  = 2.0*PI*da*In_Parm.num_photons;
    /* solid angle is 2*PI*sin(a)*da. sin(a) to be added. */

    for(ia=0; ia<na; ia++) {
        scale2 = 1.0/(sin((ia+0.5)*da)*scale1);
        if(In_Parm.isRd_a)
            Out_Ptr->Rd_a[ia] *= scale2;
        if(In_Parm.isTt_a)
            Out_Ptr->Tt_a[ia] *= scale2;
    }
  }
  if(In_Parm.isTt_d || In_Parm.isRd_d){
      scale1  = 1.0/(double)(PI*In_Parm.rad_detectors*In_Parm.rad_detectors*In_Parm.num_photons);
      for(short i=0;i<In_Parm.num_detectors;i++){
          if(In_Parm.isTt_d)
              Out_Ptr->Tt_d[i] *= scale1;
          if(In_Parm.isRd_d)
              Out_Ptr->Rd_d[i] *= scale1;
      }
  }
  if(In_Parm.isTt || In_Parm.isRd){
    scale2 = 1.0/(double)In_Parm.num_photons;
    if(In_Parm.isRd)
        Out_Ptr->Rd *= scale2;
    if(In_Parm.isTt)
        Out_Ptr->Tt *= scale2;
  }
}

/***********************************************************
 *	Scale absorption arrays properly.
 ****/
void ScaleA(InputStruct In_Parm, OutStruct * Out_Ptr)
{
  short nx = In_Parm.nx;
  short ny = In_Parm.ny;
  short nz = In_Parm.nz;
  short nr = In_Parm.nr;
  double dx = In_Parm.dx;
  double dy = In_Parm.dy;
  double dz = In_Parm.dz;
  double dr = In_Parm.dr;
  short nl = In_Parm.num_layers;
  short ix,iy,iz,ir;
  short il;
  double scale1;
  
  if(In_Parm.isA_xyz){
      scale1 = 1.0/(dx*dy*dz*In_Parm.num_photons);
      for(ix=-nx; ix<nx; ix++) 
        for(iy=-ny; iy<ny; iy++) 
            for(iz=0; iz<nz; iz++)
                Out_Ptr->A_xyz[ix][iy][iz] *= scale1;
  }
  if(In_Parm.isA_rz){
    /* Scale A_rz. */
    scale1 = 2.0*PI*dr*dr*dz*In_Parm.num_photons;	
	/* volume is 2*pi*(ir+0.5)*dr*dr*dz.*/ 
	/* ir+0.5 to be added. */
    for(iz=0; iz<nz; iz++) 
        for(ir=0; ir<nr; ir++) 
            Out_Ptr->A_rz[ir][iz] /= (ir+0.5)*scale1;
  }
  if(In_Parm.isA_z){
    /* Scale A_z. */
    scale1 = 1.0/(dz*In_Parm.num_photons);
    for(iz=0; iz<nz; iz++) 
        Out_Ptr->A_z[iz] *= scale1;
  }
  if(In_Parm.isA_l){
    /* Scale A_l. Avoid int/int. */
    scale1 = 1.0/(double)In_Parm.num_photons;	
    for(il=0; il<=nl+1; il++)
        Out_Ptr->A_l[il] *= scale1;
  }
  if(In_Parm.isA){
    scale1 = 1.0/(double)In_Parm.num_photons;	
    Out_Ptr->A *=scale1;
  }
}

void ScaleLi(InputStruct In_Parm, OutStruct *	Out_Ptr)
{
  short nl = In_Parm.num_layers;
  short nd = In_Parm.num_detectors;
  short il, id;

  for(id = 0; id < nd; id++)
     for(il = 1; il <= nl; il++)
        Out_Ptr->Li[id][il] /= (double)Out_Ptr->num_Ps[id];
}

void ScaleBannana(InputStruct In_Parm, OutStruct * Out_Ptr){
  short nz = In_Parm.nz;
  short nx = In_Parm.nx;
  short nd = In_Parm.num_detectors;
  short ix, iz, i;
  double scale;

  scale = 1.0/100.0;

  for(iz = 0; iz < nz; iz++)
     for(ix = -nx; ix < nx; ix++)
        for( i = 0; i < nd; i++)
           Out_Ptr->Bannana[i][ix][iz] *= scale;
}

/***********************************************************
 *	Sum and scale results of current run.
 ****/
void SumScaleResult(InputStruct In_Parm, OutStruct * Out_Ptr, OutStruct * Out_Ptr_Inner)
{
  short nthrds=1;
  /* Get 1D & 0D results. */
  Out_Ptr->nthread++;  
  Sum2DRd(In_Parm, Out_Ptr, Out_Ptr_Inner);
  Sum2DA(In_Parm,  Out_Ptr, Out_Ptr_Inner);
  Sum2DTt(In_Parm, Out_Ptr, Out_Ptr_Inner);
  if(In_Parm.isPPL)
    SumLi(In_Parm, Out_Ptr, Out_Ptr_Inner);
  if(In_Parm.isBannana)
    SumBannana(In_Parm, Out_Ptr, Out_Ptr_Inner);
  FreeData(In_Parm, Out_Ptr_Inner, 0);
#ifdef _OPENMP
  nthrds = omp_get_num_threads();
#endif  
  if(Out_Ptr->nthread == nthrds){
    ScaleRdTt(In_Parm, Out_Ptr);
    ScaleA(In_Parm, Out_Ptr);
    if(In_Parm.isPPL)
        ScaleLi(In_Parm, Out_Ptr);
    if(In_Parm.isBannana)
        ScaleBannana(In_Parm, Out_Ptr);
  }
}

/***********************************************************
 *	Write the input parameters to the file.
 ****/
void WriteInParm(FILE *file, InputStruct In_Parm)
{
  short i;
  
  fprintf(file, "InParm \t\t\t# Input parameters. cm is used.\n");  
  fprintf(file, "%ld \t\t\t# No. of photons\n", In_Parm.num_photons);  
  fprintf(file, "%ld \t\t\t# RNG seed\n", In_Parm.seed);  
  fprintf(file, "%G\t%G\t%G\t%G# dx, dy, dz, dr [cm]\n", In_Parm.dx,In_Parm.dy,In_Parm.dz,In_Parm.dr);
  fprintf(file, "%.3f\t%.3f\t%.3f\t# size of x, y and z.\n\n",
          (In_Parm.nx+1)*In_Parm.dx, (In_Parm.ny+1)*In_Parm.dy, (In_Parm.nz+1)*In_Parm.dz);
  
  fprintf(file, "%.3f\t%.3f\t%.3f\t# source position\n", In_Parm.srcx, In_Parm.srcy, In_Parm.srcz);
  fprintf(file, "%.3f\t%.3f\t%.3f\t# source direction\n", In_Parm.dirx, In_Parm.diry, In_Parm.dirz);
  switch(In_Parm.srctyp){
      case 0: fprintf(file, "source type is pencil\n\n"); break;
      case 1: fprintf(file, "source type is disk distribution\n\n"); break;
      case 2: fprintf(file, "source type is gaussian\n\n"); break;
      case 3: fprintf(file, "source type is fiber\n\n"); break;
  }  
  
  fprintf(file, "%hd\t%.3f\t%.3f\t\t# Number, raduis and NA of detectors\n\n", In_Parm.num_detectors, In_Parm.rad_detectors, In_Parm.NA_detectors);
  
  fprintf(file, "%hd\t\t\t\t\t# Number of layers\n", In_Parm.num_layers);  
  fprintf(file, "#n\tmua\tmus\tg\td\t# One line for each layer\n"); 
  fprintf(file, "%G\t\t\t\t\t# n for medium above\n", In_Parm.layerspecs[0].n); 
  for(i=1; i<=In_Parm.num_layers; i++){
    LayerStruct s;
    s = In_Parm.layerspecs[i];
    fprintf(file, "%G\t%G\t%G\t%G\t%G\t# layer %hd\n",
	        s.n, s.mua, s.mus, s.g, s.z1-s.z0, i);
  }
  fprintf(file, "%G\t\t\t\t\t# n for medium below\n\n", 
	In_Parm.layerspecs[i].n); 
}

/***********************************************************
 *	Write reflectance, absorption, transmission. 
 ****/
void WriteRAT(FILE * file, OutStruct Out_Parm, InputStruct In_Parm)
{
  fprintf(file, "RAT #Reflectance, absorption, transmission. \n");
	/* flag. */
  if(In_Parm.isRsp)
    fprintf(file, "%-14.6G \t#Specular reflectance [-]\n", Out_Parm.Rsp);
  if(In_Parm.isRd)
    fprintf(file, "%-14.6G \t#Diffuse reflectance [-]\n", Out_Parm.Rd);
  if(In_Parm.isA)
    fprintf(file, "%-14.6G \t#Absorbed fraction [-]\n", Out_Parm.A);
  if(In_Parm.isTt)
    fprintf(file, "%-14.6G \t#Transmittance [-]\n", Out_Parm.Tt);
  
  fprintf(file, "\n");
}

/***********************************************************
 *	Write absorption as a function of layer. 
 ****/
void WriteA_layer(FILE * file, short Num_Layers, OutStruct Out_Parm)
{
  short i;
  
  fprintf(file, "A_l #Absorption as a function of layer. [-]\n layer_index\tAbsorption\n");
	/* flag. */

  for(i=1; i<=Num_Layers; i++)
    fprintf(file, "%hd\t%12.4G\n", i, Out_Parm.A_l[i]);
  fprintf(file, "\n");
}

/***********************************************************
 *	5 numbers each line.
 ****/
void WriteRd_ra(FILE * file, short Nr, short Na, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ir, ia;
  
  fprintf(file, 
	  "\n\n%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	  "# position   Rd[r][angle]. [1/(cm2sr)].",
	  "# Rd[0][0], [0][1],..[0][na-1]",
	  "# Rd[1][0], [1][1],..[1][na-1]",
	  "# ...",
	  "# Rd[nr-1][0], [nr-1][1],..[nr-1][na-1]",
	  "r    angle   Rd_ra\n\n");
  
  for(ir=0;ir<Nr;ir++)
    for(ia=0;ia<Na;ia++)
      fprintf(file, "%.3f\t%.3f\t%12.4E\n", In_Parm.dr*ir, In_Parm.da*ia, Out_Parm.Rd_ra[ir][ia]);    
}

void WriteRd_lr(FILE * file, short Nl, short Nr, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ir, il;
  
  fprintf(file, 
	  "%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	  "# Rd[nl][x]. [1/cm].",
	  "# Rd[0][0], [0][1],..[0][nr-1]",
	  "# Rd[1][0], [1][1],..[1][nr-1]",
	  "# ...",
	  "# Rd[nl-1][0], [nl-1][1],..[nl-1][nr-1]",
	  "Rd_lr\n");
  
  for(il=1;il<=Nl;il++){
    fprintf(file, "layer %hd\n",il);
    for(ir=0;ir<Nr;ir++)
      fprintf(file, "%.3f\t%12.4E\n", In_Parm.dr*ir, Out_Parm.Rd_lr[il][ir]);
  }
}

/***********************************************************
 *	1 number each line.
 ****/
void WriteRd_r(FILE * file, short Nr, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ir;
  
  fprintf(file, "Rd_r #Rd[0], [1],..Rd[nr-1]. [1/cm2]\n");	/* flag. */
  
  for(ir=0;ir<Nr;ir++) {
    fprintf(file, "%.3f\t%12.4E\n", In_Parm.dr*ir, Out_Parm.Rd_r[ir]);
  }  
}

/***********************************************************
 *	1 number each line.
 ****/
void WriteRd_a(FILE * file, short Na, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ia;
  
  fprintf(file, "\n\nRd_a #Rd[0], [1],..Rd[na-1]. [sr-1]\n\n");	/* flag. */
  
  for(ia=0;ia<Na;ia++) {
    fprintf(file, "%.3f\t%12.4E\n", In_Parm.da*ia, Out_Parm.Rd_a[ia]);
  }
}

/***********************************************************
 *	5 numbers each line.
 ****/
void WriteTt_ra(FILE * file, short Nr, short Na, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ir, ia;
  
  fprintf(file, 
	  "\n\n%s\n%s\n%s\n%s\n%s\n%s\n",	/* flag. */
	  "# position   Tt[r][angle]. [1/(cm2sr)].",
	  "# Tt[0][0], [0][1],..[0][na-1]",
	  "# Tt[1][0], [1][1],..[1][na-1]",
	  "# ...",
	  "# Tt[nr-1][0], [nr-1][1],..[nr-1][na-1]",
	  "r    angle   Tt_ra\n");
  
  for(ir=0;ir<Nr;ir++)
    for(ia=0;ia<Na;ia++)
      fprintf(file, "%.3f\t%.3f\t%12.4E\n ", In_Parm.dr*ir, In_Parm.da*ia, Out_Parm.Tt_ra[ir][ia]);
}

void WriteA_xyz(FILE *file, short Nx, short Ny, short Nz, OutStruct Out_Parm, InputStruct In_Parm)
{
  short iz, ix, iy;
  
  fprintf(file, 
	  "\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n", /* flag. */
	  "# position   A[x][y][z]. [1/cm3]",
	  "# A[0][0][0], [1][0][0],..[nx-1][0][0]",
	  "# A[0][1][0], [1][1][0],..[nx-1][1][0]",
	  "# ...",
	  "# A[0][ny-1][0], [1][ny-1][0],..[nx-1][ny-1][0]",
      "# ...",
      "# A[0][ny-1][nz-1], [1][ny-1][nz-1],..[nx-1][ny-1][nz-1]",
	  "x    y   z   A_xyz\n");
      
  for(iz=0;iz<Nz;iz++)
    for(iy=-Ny;iy<Ny;iy++)
        for(ix=-Nx;ix<Nx;ix++)  {
            fprintf(file, "%.3f\t%.3f\t%.3f\t%12.4E\n", In_Parm.dx*ix, In_Parm.dy*iy , In_Parm.dz*iz, Out_Parm.A_xyz[ix][iy][iz]);
        }
  
}

/***********************************************************
 *	5 numbers each line.
 ****/
void WriteA_rz(FILE * file, short Nr, short Nz, OutStruct Out_Parm, InputStruct In_Parm)
{
  short iz, ir;
  
  fprintf(file, 
	  "%s\n%s\n%s\n%s\n%s\n%s\n", /* flag. */
	  "# position   A[r][z]. [1/cm3]",
	  "# A[0][0], [1][0],..[nr-1][0]",
	  "# A[0][1], [1][1],..[nr-1][1]",
	  "# ...",
	  "# A[0][nz-1], [1][nz-1],..[nr-1][nz-1]",
	  "r    z   A_rz\n");
  
  for(iz=0;iz<Nz;iz++)
    for(ir=0;ir<Nr;ir++){
      fprintf(file, "%.3f\t%.3f\t%12.4E\n", In_Parm.dr*ir, In_Parm.dz*iz, Out_Parm.A_rz[ir][iz]);
    } 
}

/***********************************************************
 *	1 number each line.
 ****/
void WriteA_z(FILE * file, short Nz, OutStruct Out_Parm, InputStruct In_Parm)
{
  short iz;
  
  fprintf(file, "A_z #A[0], [1],..A[nz-1]. [1/cm]\n\n");	/* flag. */
  
  for(iz=0;iz<Nz;iz++) {
    fprintf(file, "%.3f\t%12.4E\n", iz*In_Parm.dz, Out_Parm.A_z[iz]);
  }
  
  fprintf(file, "\n");
}

/***********************************************************
 *	1 number each line.
 ****/
void WriteTt_r(FILE * file, short Nr, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ir;
  
  fprintf(file, "Tt_r #Tt[0], [1],..Tt[nr-1]. [1/cm2]\n\n"); /* flag. */
  
  for(ir=0;ir<Nr;ir++) {
    fprintf(file, "%.3f\t%12.4E\n", In_Parm.dr*ir, Out_Parm.Tt_r[ir]);
  }
}

/***********************************************************
 *	1 number each line.
 ****/
void WriteTt_a(FILE * file, short Na, OutStruct Out_Parm, InputStruct In_Parm)
{
  short ia;
  
  fprintf(file, 
	"\n\nTt_a #Tt[0], [1],..Tt[na-1]. [sr-1]\n"); /* flag. */
  
  for(ia=0;ia<Na;ia++) {
    fprintf(file, "%.3f\t%12.4E\n", In_Parm.da*ia, Out_Parm.Tt_a[ia]);
  }
}

void WriteBannana(FILE * file, InputStruct In_Parm, OutStruct Out_Parm){
  short ix, iz, id;
  short nz = In_Parm.nz;
  short nx = In_Parm.nx;
  short nd = In_Parm.num_detectors;

  for(id = 0; id < nd; id++){
    for(iz = 0; iz < nz; iz++){
     for(ix = -nx; ix < nx; ix++)
        fprintf(file, "%12.4E\t",Out_Parm.Bannana[id][ix][iz]);
     fprintf(file,"\n");
    }
    fprintf(file,"\n\n");
  }
}

/***********************************************************
 ****/
void WriteResult(InputStruct In_Parm, OutStruct Out_Parm, char * TimeuserReport, char * ClockReport)
{
  FILE *file, *filea;
  char filewrite[STRLEN];
  
  file = fopen(In_Parm.out_fname, "w");
  if(file == NULL) nrerror("Cannot open file to write.\n");
  
  fprintf(file, "# %s\n", TimeuserReport);
  fprintf(file, "# %s\n", ClockReport);
  
  WriteInParm(file, In_Parm);
  if(In_Parm.isRd || In_Parm.isRsp || In_Parm.isA || In_Parm.isTt)
    WriteRAT(file, Out_Parm, In_Parm);	/* reflectance, absorption, transmittance. */
  if(In_Parm.isPPL){
    fprintf(file,"PPL (cm)\n\n");
    for(short i=0;i<In_Parm.num_detectors;i++)
        for(short j=1;j<=In_Parm.num_layers;j++)
            fprintf(file,"Detector(%hd)\t\tLayer(%hd)\t%.3f\tNo. of Phothons(%d)\n",i+1,j,Out_Parm.Li[i][j],Out_Parm.num_Ps[i]);
  }
  if(In_Parm.isRd_d){
    fprintf(file,"\nRecorded reflected photons in each detector\n\n");
    for(short i=0;i<In_Parm.num_detectors;i++)
        fprintf(file,"Detector(%d)\t\t%12.4E\n",i+1,Out_Parm.Rd_d[i]);
  }
  if(In_Parm.isTt_d){
    fprintf(file,"\nRecorded transmited photons in each detector\n\n");
    for(short i=0;i<In_Parm.num_detectors;i++)
        fprintf(file,"Detector(%d)\t\t%12.4E\n",i+1,Out_Parm.Tt_d[i]);
  }
  fclose(file);
  if(In_Parm.isA_l || In_Parm.isA_z || In_Parm.isA_rz || In_Parm.isA_xyz){
    sprintf(filewrite,"absorption_%s",In_Parm.out_fname);
    if((filea = fopen(filewrite, "w"))== NULL) nrerror("Cannot open file to write.\n");
    if(In_Parm.isA_l)
        WriteA_layer(filea, In_Parm.num_layers, Out_Parm);
    if(In_Parm.isA_z)
        WriteA_z(filea, In_Parm.nz, Out_Parm, In_Parm);
    if(In_Parm.isA_rz)
        WriteA_rz(filea, In_Parm.nr, In_Parm.nz, Out_Parm, In_Parm);
    if(In_Parm.isA_xyz)
        WriteA_xyz(file, In_Parm.nx, In_Parm.ny, In_Parm.nz, Out_Parm, In_Parm);
    fclose(filea);
  }
  if(In_Parm.isRd_r || In_Parm.isRd_a || In_Parm.isRd_ra || In_Parm.isRd_lr){
    sprintf(filewrite,"reflectance_%s",In_Parm.out_fname);
    if((filea = fopen(filewrite, "w"))== NULL) nrerror("Cannot open file to write.\n");
    if(In_Parm.isRd_r)
        WriteRd_r(file, In_Parm.nr, Out_Parm, In_Parm);
    if(In_Parm.isRd_a)
        WriteRd_a(file, In_Parm.na, Out_Parm, In_Parm);
    if(In_Parm.isRd_ra)
        WriteRd_ra(file, In_Parm.nr, In_Parm.na, Out_Parm, In_Parm);
    if(In_Parm.isRd_lr)
        WriteRd_lr(file, In_Parm.num_layers, In_Parm.nr, Out_Parm, In_Parm);
    fclose(filea);
  }
  if(In_Parm.isTt_r || In_Parm.isTt_a || In_Parm.isTt_ra){
      sprintf(filewrite,"transmittance_%s",In_Parm.out_fname);
      if((filea = fopen(filewrite, "w"))== NULL) nrerror("Cannot open file to write.\n");
      if(In_Parm.isTt_r)
        WriteTt_r(file, In_Parm.nr, Out_Parm, In_Parm);
      if(In_Parm.isTt_r)
        WriteTt_a(file, In_Parm.na, Out_Parm, In_Parm);
      if(In_Parm.isTt_r)
        WriteTt_ra(file, In_Parm.nr, In_Parm.na, Out_Parm, In_Parm);
      fclose(filea);
  }
  if(In_Parm.isBannana){
    sprintf(filewrite,"bannana_%s",In_Parm.out_fname);
    if((filea = fopen(filewrite, "w"))== NULL) nrerror("Cannot open file to write.\n");
    if(In_Parm.isBannana)
        WriteBannana(file, In_Parm, Out_Parm);
    fclose(filea);
  }
}