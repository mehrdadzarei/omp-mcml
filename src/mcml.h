/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Monte Carlo simulation of photon distribution in 
 *	multi-layered turbid media in ANSI Standard C.
 ****
 *	Starting Date:	10/1991.
 *	Current Date:	6/1992.
 *
 *	Lihong Wang, Ph. D.
 *	Steven L. Jacques, Ph. D.
 *	Laser Biology Research Laboratory - 17
 *	M.D. Anderson Cancer Center
 *	University of Texas
 *	1515 Holcombe Blvd.
 *	Houston, TX 77030
 *	USA
 *
 *	This program was based on:
 *	(1) The Pascal code written by Marleen Keijzer and 
 *	Steven L. Jacques in this laboratory in 1989, which
 *	deals with multi-layered turbid media.
 *
 *	(2) Algorithm for semi-infinite turbid medium by 
 *	S.A. Prahl, M. Keijzer, S.L. Jacques, A.J. Welch, 
 *	SPIE Institute Series Vol. IS 5 (1989), and by 
 *	A.N. Witt, The Astrophysical journal Supplement
 *	Series 35, 1-6 (1977).
 *	
 *	Major modifications include:
 *		. Conform to ANSI Standard C.
 *		. Removal of limit on number of array elements, 
 *		  because arrays in this program are dynamically 
 *		  allocated. This means that the program can accept 
 *		  any number of layers or gridlines as long as the 
 *		  memory permits.
 *		. Avoiding global variables whenever possible.  This
 *		  program has not used global variables so far.
 *		. Grouping variables logically using structures.
 *		. Top-down design, keep each subroutine clear & 
 *		  short.	
 *		. Reflectance and transmittance are angularly 
 *		  resolved.
 ****
 *	General Naming Conventions:
 *	Preprocessor names: all capital letters, 
 *		e.g. #define PREPROCESSORS
 *	Globals: first letter of each word is capital, no 
 *		underscores, 
 *		e.g. short GlobalVar;
 *	Dummy variables:  first letter of each word is capital,
 *		and words are connected by underscores, 
 *		e.g. void NiceFunction(char Dummy_Var);
 *	Local variables:  all lower cases, words are connected 
 *		by underscores,
 *		e.g. short local_var;
 *	Function names or data types:  same as Globals.
 *
 ****
 *	Dimension of length: cm.
 *
 ****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "posix_randr.c"

#ifdef _OPENMP
  #include <omp.h>
#endif

#define PI 3.1415926
#define WEIGHT 1E-4		/* Critical weight for roulette. */
#define CHANCE 0.1		/* Chance of roulette survival. */
#define STRLEN 256		/* String length. */
#define EPS        1e-6f

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)

/****************** Stuctures *****************************/

/****
 *	Structure used to describe a photon packet.
 ****/
typedef struct {
  double x, y ,z;	/* Cartesian coordinates.[cm] */
  double ux, uy, uz;/* directional cosines of a photon. */
  double w;			/* weight. */
  Boolean dead;		/* 1 if photon is terminated. */
  short layer;		/* index to layer where the photon packet resides. */
  short layer_index; /* layer index for each photon */
  double s;			/* current step size. [cm]. */
  double sleft;		/* step size left. dimensionless [-]. */
} PhotonStruct;

/****
 *	Structure used to describe the geometry and optical
 *	properties of a layer.
 *	z0 and z1 are the z coordinates for the upper boundary
 *	and lower boundary respectively.
 *
 *	cos_crit0 and cos_crit1 are the cosines of the 
 *	critical angle of total internal reflection for the
 *	upper boundary and lower boundary respectively.
 *	They are set to zero if no total internal reflection
 *	exists.
 *	They are used for computation speed.
 ****/
typedef struct {
  double z0, z1;	/* z coordinates of a layer. [cm] */
  double xD;        /* detector position (cm) */
  double yD;        /* detector position (cm) */
  double n;			/* refractive index of a layer. */
  double mua;	    /* absorption coefficient. [1/cm] */
  double mus;	    /* scattering coefficient. [1/cm] */
  double g;		    /* anisotropy. */
  
  double cos_crit0,	cos_crit1;	
} LayerStruct;

/****
 *	Input parameters for each independent run.
 *
 *	z and r are for the cylindrical coordinate system. [cm]
 *	a is for the angle alpha between the photon exiting 
 *	direction and the surface normal. [radian]
 *
 *	The grid line separations in z, r, and alpha
 *	directions are dz, dr, and da respectively.  The numbers 
 *	of grid lines in z, r, and alpha directions are
 *	nz, nr, and na respectively.
 *
 *	The member layerspecs will point to an array of 
 *	structures which store parameters of each layer. 
 *	This array has (number_layers + 2) elements. One
 *	element is for a layer.
 *	The layers 0 and (num_layers + 1) are for top ambient 
 *	medium and the bottom ambient medium respectively.
 ****/
typedef struct {
  char	 out_fname[STRLEN];	/* output file name. */
  char	 out_fformat;		/* output file format. */
							/* 'A' for ASCII, */
							/* 'B' for binary. */
  long	 num_photons; 		/* to be traced. */
  short	 num_threads; 		/* to be traced. */
  long	 seed;       		/* for random generation. */
  double Wth; 				/* play roulette if photon weight < Wth.*/
  
  double dx;				/* x grid separation.[cm] */
  double dy;				/* y grid separation [cm] */
  double dz;				/* z grid separation.[cm] for cylindrical & cartesian*/ 
  double dr;				/* r grid separation.[cm] for cylindrical */
  double da;				/* alpha grid separation. for cylindrical[radian] */
  short nx;			    	/* array range 0..nx-1. */
  short ny;				    /* array range 0..ny-1. */
  short nz;					/* array range 0..nz-1. */
  short nr;					/* array range 0..nr-1. */
  short na;					/* array range 0..na-1. */
  
  double srcx,srcy,srcz;    /* source position*/
  double dirx,diry,dirz;    /* source direction*/
  int srctyp;               /* source type; 0: pencil beam, 1: disc distribution with parameters raduis,
                            2: gaussian beam with parameters raduis, 3: fiber with diameter and NA.*/
  double srcpar1;            /* source parameter*/ 
  double srcpar2;            /* source parameter*/ 
  
  LayerStruct *detectors;   /* position of detectors */
  short	num_detectors;	    /* number of detectors. */
  double rad_detectors;	    /* raduise of detectors. */
  double NA_detectors;	    /* Numerical Aperture of detectors. */
  
  short	num_layers;			/* number of layers. */
  LayerStruct * layerspecs;	/* layer parameters. */	
  
  Boolean isRsp;            /* 1 save Rsp.*/
  Boolean isRd_ra;          /* 1 save Rd_ra.*/
  Boolean isRd_lr;          /* 1 save Rd_lr.*/
  Boolean isRd_r;           /* 1 save Rd_r.*/
  Boolean isRd_d;           /* 1 save Rd_d.*/
  Boolean isRd_a;           /* 1 save Rd_a.*/
  Boolean isRd;             /* 1 save Rd.*/
  Boolean isA_xyz;          /* 1 save A_xyz.*/
  Boolean isA_rz;           /* 1 save A_rz.*/
  Boolean isA_z;            /* 1 save A_z.*/
  Boolean isA_l;            /* 1 save A_l.*/
  Boolean isA;              /* 1 save A.*/
  Boolean isTt_ra;          /* 1 save Tt_ra.*/
  Boolean isTt_r;           /* 1 save Tt_r.*/
  Boolean isTt_d;           /* 1 save Tt_d.*/
  Boolean isTt_a;           /* 1 save Tt_a.*/
  Boolean isTt;             /* 1 save Tt.*/
  Boolean isBannana;        /* 1 save Bannana.*/
  Boolean isPPL;            /* 1 save PPL.*/
} InputStruct;

/****
 *	Structures for scoring physical quantities. 
 *	z and r represent z and r coordinates of the 
 *	cylindrical coordinate system. [cm]
 *	a is the angle alpha between the photon exiting 
 *	direction and the normal to the surfaces. [radian]
 *	See comments of the InputStruct.
 *	See manual for the physcial quantities.
 ****/
typedef struct {
  double    Rsp;	/* specular reflectance. [-] */
  double ** Rd_ra;	/* 2D distribution of diffuse reflectance. [1/(cm2 sr)] */
  double ** Rd_lr;	/* 1D radial distribution of diffuse reflectance for each layer. [1/cm2] */  
  double *  Rd_r;	/* 1D radial distribution of diffuse reflectance. [1/cm2] */
  double *  Rd_d;	/* diffuse reflectance for each detector. [1/cm2] */
  double *  Rd_a;	/* 1D angular distribution of diffuse reflectance. [1/sr] */
  double    Rd;		/* total diffuse reflectance. [-] */
  
  double *** A_xyz;	/* 3D probability density in turbid media over x, y & z. [1/cm3] */
  double ** A_rz;	/* 2D probability density in turbid media over r & z. [1/cm3] */
  double *  A_z;	/* 1D probability density over z. [1/cm] */
  double *  A_l;	/* each layer's absorption probability. [-] */
  double    A;		/* total absorption probability. [-] */
  
  double ** Tt_ra;	/* 2D distribution of total transmittance. [1/(cm2 sr)] */
  double *  Tt_r;	/* 1D radial distribution of transmittance. [1/cm2] */
  double *  Tt_d;	/* transmittance for each detector. [1/cm2] */
  double *  Tt_a;	/* 1D angular distribution of transmittance. [1/sr] */
  double    Tt;		/* total transmittance. [-] */
  
  double **Bannana_inner;/* Bannana shape's array.Dimention for [x,z] */ 
  double ***Bannana;/* Bannana shape's array. One dimention for number of detectors and tow dimention for [x,z] */

  double * li;       /* length of path that photon will move in each layer*/
  double ** Li;     /* length of path that photon will move in each layer for each detector*/
  long * num_Ps ;   /* number of photons which come out for each detector*/

  short nthread;
} OutStruct;

/***********************************************************
 *	Routine prototypes for dynamic memory allocation and 
 *	release of arrays and matrices.
 *	Modified from Numerical Recipes in C.
 ****/
double  *AllocVector(short, short);
double **AllocMatrix(short, short,short, short);
double ***AllocCubic(short, short,short, short,short, short);
void 	 FreeVector(double *, short, short);
void 	 FreeMatrix(double **, short, short, short, short);
void 	 FreeCubic(double ***, short, short, short, short, short, short);
void 	 nrerror(char *);