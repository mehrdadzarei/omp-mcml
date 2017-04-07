/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

#include "mcml.h"

/*	Declare before they are used in main(). */
FILE *GetFile(char *);
short ReadNumRuns(FILE* );
void ReadParm(FILE* , InputStruct * );
void CheckParm(FILE* , InputStruct * );
void InitOutputData(InputStruct, OutStruct *, Boolean m);
void FreeData(InputStruct, OutStruct *, short);
double Rspecular(LayerStruct * );
double LaunchPhoton(/*double, */LayerStruct *, PhotonStruct *, InputStruct *, RandType *, RandType *);
void HopDropSpin(InputStruct  *,PhotonStruct *,OutStruct *, RandType *);
void SumScaleResult(InputStruct, OutStruct *, OutStruct *);
void WriteResult(InputStruct, OutStruct, char *, char *);


/***********************************************************
 *	If F = 0, reset the clock and return 0.
 *
 *	If F = 1, pass the user time to Msg and print Msg on 
 *	screen, return the real time since F=0. 
 *
 *	If F = 2, same as F=1 except no printing.  
 *
 *	Note that clock() and time() return user time and real 
 *	time respectively.
 *	User time is whatever the system allocates to the 
 *	running of the program; 
 *	real time is wall-clock time.  In a time-shared system,
 *	they need not be the same.
 *	
 *	clock() only hold 16 bit integer, which is about 32768 
 *	clock ticks.
 ****/
time_t PunchTime(char F, char *Msguser, char *Msgclock)
{
  static clock_t ut0;	/* user time reference. */
  static time_t  rt0;	/* real time reference. */
  double secsclock, secsuser;
  char s[STRLEN], s1[STRLEN];
  
  if(F==0) {
    ut0 = clock();
    rt0 = time(NULL);
    return(0);
  }
  else if(F==1)  {
    secsclock = (clock() - ut0)/(double)CLOCKS_PER_SEC;
    secsuser = time(NULL) - rt0;
    if (secsclock<0) secsclock=0;	/* clock() can overflow. */
    sprintf(s, "clock: %8.0lf sec = %8.2lf min.  %s\n", 
	    secsclock, secsclock/60.0, Msgclock);
    sprintf(s1, "User time: %8.0lf sec = %8.2lf min.  %s\n", 
	    secsuser, secsuser/60.0, Msguser);

    strcpy(Msguser, s1);
    strcpy(Msgclock, s);
    return(difftime(time(NULL), rt0));
  }
  else if(F==2) return(difftime(time(NULL), rt0));
  else return(0);
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, long Pt)	
{
  time_t now, done_time;
  struct tm *date;
  char s[80];
  
  now = time(NULL);
  date = localtime(&now);
  strftime(s, 80, "%H:%M %x", date);
  printf("Now %s, ", s);
  
  done_time = now + (time_t) (PunchTime(2,"","")*(double)(Pt-P1)/(double)(P1));
  date = localtime(&done_time);
  strftime(s, 80, "%H:%M %x", date);
  printf("End %s\n", s);
}

/***********************************************************
 *	Get the file name of the input data file from the 
 *	argument to the command line.
 ****/
void GetFnameFromArgv(int argc, char * argv[], char * input_filename)
{
  if(argc>=2) {			/* filename in command line */
    strcpy(input_filename, argv[1]);
  }
  else
    input_filename[0] = '\0';
} 
    
/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 ****/
void DoOneRun(InputStruct *In_Ptr)
{
  char timeuser_report[STRLEN];
  char clock_report[STRLEN];
  register long i_photon;	/* index to photon. register for speed.*/
  OutStruct out_parm;		/* distribution of photons.*/  
  long num_photons = In_Ptr->num_photons, photon_rep=10000;
  RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  unsigned int id=0, nthrds=1;
  int seed = In_Ptr->seed;
  Boolean m = 1;        /* 1 master thread*/
  
  InitOutputData(*In_Ptr, &out_parm, m);
  // out_parm.Rsp = Rspecular(In_Ptr->layerspecs);	
  PunchTime(0, "", "");
#ifdef _OPENMP
  omp_set_num_threads(In_Ptr->num_threads);
#endif

#pragma omp parallel private(ran0,ran1,id)
{
  PhotonStruct photon;
  OutStruct out_parm_inner;
  short index_layer;
  short num_layer = In_Ptr->num_layers;

  InitOutputData(*In_Ptr, &out_parm_inner, m=0);
  // rng_init(ran0,ran1,(unsigned int *)&(seed),id);

#ifdef _OPENMP
{
	id=omp_get_thread_num();
    nthrds = omp_get_num_threads();
}
#endif
  
#pragma omp for  
  for(i_photon = 0; i_photon < num_photons; i_photon++){
    if(i_photon==photon_rep && id == 0){
      PredictDoneTime(i_photon*nthrds, num_photons);
      printf("%d photons finished., \n", i_photon*nthrds);      
    }
    if(In_Ptr->isPPL)
        for(index_layer = 1; index_layer <= num_layer; index_layer++)
            out_parm_inner.li[index_layer] = 0;
    out_parm_inner.Rsp += LaunchPhoton(/*out_parm.Rsp, */In_Ptr->layerspecs, &photon, In_Ptr, ran0, ran1);

    do HopDropSpin(In_Ptr, &photon, &out_parm_inner, ran0);
    while (!photon.dead);
  } 

#pragma omp critical
  SumScaleResult(*In_Ptr, &out_parm, &out_parm_inner);
}      
  strcpy(timeuser_report, " Simulation time of this run (user).");
  strcpy(clock_report, " Simulation time of this run (clock).");
  PunchTime(1, timeuser_report, clock_report);
  WriteResult(*In_Ptr, out_parm, timeuser_report, clock_report);      
  FreeData(*In_Ptr, &out_parm, 1);
}

/***********************************************************
 *	The argument to the command line is filename, if any.
 *	Macintosh does not support command line.
 ****/
int main(int argc, char *argv[]) 
{
  char input_filename[STRLEN];
  FILE *input_file_ptr;
  short num_runs;	/* number of independent runs. */
  InputStruct in_parm;

  GetFnameFromArgv(argc, argv, input_filename);
  input_file_ptr = GetFile(input_filename);
  CheckParm(input_file_ptr, &in_parm);	
  num_runs = ReadNumRuns(input_file_ptr);
  
  while(num_runs--)  {
    ReadParm(input_file_ptr, &in_parm);
	DoOneRun(&in_parm);
  }
  
  fclose(input_file_ptr);
  return 0;
}