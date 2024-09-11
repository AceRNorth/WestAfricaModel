#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib> // for exit function
#include <tr1/random>
#include <numeric>
			
/*-----------------------------------------------------Header definitions---------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------*/

#ifndef RANDOMC_H
#define RANDOMC_H

const int nx=133;//grid cells for rain input data
const int ny=49;//grid cells for rain input data
const int TL=20;//juvenile development time (egg to adult)
const int NumGen=6;

using namespace std;
//using std::ifstream;
using namespace std::tr1;

const double PI=3.14159265;
const double TWOPI=6.28318531;
const long long int LONG_MAX=3147483647000;


/*----------------------struct for keeping track of global numbers------------------------*/
struct totals
	{
	
	long long int J[NumGen];
	long long int M[NumGen];
	long long int V[NumGen];
	long long int F[NumGen];
	long long int CentF[NumGen];
	long long int Em[NumGen];
	long long int JTot,MTot,VTot,FTot;
	int CentSqVils;
	double CentSqHum;
	double distW,distD,distR;
	};	
/*----------------------------------------------------------------------------------------*/


/*----------------------struct combining initial condition parameters---------------------*/
struct initials
{
	int NumAdultsWM; int NumAdultsWV; int NumAdultsWF;
	int NumAdultsDM; int NumAdultsDV; int NumAdultsDF;
	int NumAdultsRM; int NumAdultsRV; int NumAdultsRF;
	int NumJW[TL]; int NumJD[TL]; int NumJR[TL];
	vector<pair<int,int>> relpatches;
	//vector<int> reltimes;
	string inputfile;
	string polyconnect;
	string polyedges;
	string rainfile;
	string mortfile;
	double driver_time;
	double driver_start;
	double driver_end;
	double r_time;
	int NumDriver;
	double NumDriverSites;
	double NumDriverD;
	double NumDriverSitesD;
	int NumR;
	int NumRes;
	char dist;
	int recSitesFreq;
	};	
/*----------------------------------------------------------------------------------------*/

/*---------------------------number of settlements in each grid cell----------------------*/
vector<int>SetsPerCell[nx][ny]; 
/*----------------------------------------------------------------------------------------*/
		
/*------------------------state of population in each settlement---------------------------*/
struct Patch
{
	double x;
	double y;
	double Nhum;
	double gam;
	double arab;
	double fun;
	double area;
	double dens;
	int sqx,sqy;
	int numsites;
	long long int J[NumGen][TL];
	long long int JTot;
	long long int MTot;
	long long int M[NumGen];
	long long int V[NumGen];
	//int Fww,Fwd,Fdd,Fwr,Frr,Fdr;
	long long int F[NumGen][NumGen];	
	long long int AesF[NumGen][NumGen];	
	long long int MoveF[NumGen][NumGen];	
	long long int MoveM[NumGen];
	long double comp;
	char arrive;
	char release;
	long double mate_rate;
	vector<int> connecIND;
	vector<double> connecW;
	double TotW;
	vector<int> EDGEconnecIND;
	vector<double> EDGEconnecW;
	double EDGETotW;
	char type;
	char Fixed;
	int CentSq;
	double distOrigin;
	double WaterTemp,WaterPerm;
	double alpha0;
	int clusterID;
	long long int yearTotF;//total biting females produced during year in that site
};
/*----------------------------------------------------------------------------------------*/



/*----------------struct containt simulation timekeeping parameters-------------------------*/
struct Times
{
	int maxT;
	int start;
	int recfreq;
	int recstart;
	int recend;
	int NumRuns;
	int interval;
	int yearnow;
};



/*----------------------------------------------------------------------------------------*/

/*---------------------struct containing model parameters----------------------------------*/
struct Pars
	{
	double muJ,muA,muAF,d,gamma,Fgamma,Mgamma,beta,theta,Frho,Mrho,xi,e,em,ef,LD,dL,muLD,psi,muAES,bias;
	double alpha0,alpha1,alpha2,delta,phi,kappa,al0var,mu,sig;
	double dispa,dispb;
	double EDGEd;
	double omega,s;
	double LarvProbs[TL];
	double meanTL;
	double U;
	double CentRad;
	int kernel;
	int singlepop;
	int NumPat;
	char species;
	int set,index;
	int t_disp1,t_disp2,t_disp3,t_disp4;
	int t_hide1,t_hide2,t_wake1,t_wake2;
	double f[NumGen][NumGen][NumGen];
	int dispStart,dispLength,dispSoFar,dispSoFarEDGE;
	};
/*----------------------------------------------------------------------------------------*/

struct combo
	{
	int num;
	int site;
	};	
	
/*-------------------------------function declarations-------------------------------------*/
		void RunOnceInt(double);
		void record(int);
		void RunMaxT();
		void RunNReps(int);
		void ResetSites();
		void initiate();
		void SitesPopulate(int);
		void OneStep(int);
		void PutDriver();
		void PutDriverSites(int);
		void JuvGetOlder();
		void VirginsMate();
		void AdultsMove(int);
		void Hide();
		void Wake(int);
		void LayEggs();
		void JuvEmerge();
		void AdultsDie();
		void UpdateComp(int);
		void UpdateMate();
		void UpdateConnec();
		void SetFertility();
		long long int random_poisson(double);
		long long int random_binomial(long long int,double);
		double random_normal(double,double);
		int* random_multinom(int,long long int[6]);
		int* random_multinomEqualProb(int,int);
		int* random_multinom_var(int,int,double*,double);
		vector<combo> random_multinom_varB(int,int,double*,double);
		double dist(double,double,double,double);
		double distB(double,double,double,double,double);
		double abso(double);
		long long longmin(long long int, long long int);
		long long longmax(long long int, long long int);

/*----------------------------------------------------------------------------------------*/

/*--------------------------code for random number generation------------------------------*/
// Define 32 bit signed and unsigned integers.
// GieIange these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
   // 16 bit systems use long int for 32 bit integer
   typedef long int           int32;   // 32 bit signed integer
   typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
   // Most other systems use int for 32 bit integer
   typedef int                int32;   // 32 bit signed integer
   typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
   // Microsoft and other compilers under Windows use __int64
   typedef __int64            int64;   // 64 bit signed integer
   typedef unsigned __int64   uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
   // Gnu and other compilers under Linux etc. use long long
   typedef long long          int64;   // 64 bit signed integer
   typedef unsigned long long uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#else
   // 64 bit integers not defined
   // You may include definitions for other platforms here
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(char * ErrorText);     // System-specific error reporting (userintf.cpp)


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else
   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
public:
   CRandomMersenne(uint32 seed) {      // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(uint32 seed);       // Re-seedSetDirectory["~/Dropbox/YDriveBurkina/codesB/DoubleSex/Mosaic"];
   void RandomInitByArray(uint32 seeds[], int length); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32 BRandom();                   // Output random bits
private:
   void Init0(uint32 seed);            // Basic initialization procedure
   uint32 mt[MERS_N];                  // State vector
   int mti;                            // Index into mt
   uint32 LastInterval;                // Last interval length for IRandomX
   uint32 RLimit;                      // Rejection limit used by IRandomX
   enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
   TArch Architecture;                 // Conversion to float depends on architecture
};


class CRandomMother {             // Encapsulate random number generator
public:
   void RandomInit(uint32 seed);       // Initialization
   int IRandom(int min, int max);      // Get integer random number in desired interval
   double Random();                    // Get floating point random number
   uint32 BRandom();                   // Output random bits
   CRandomMother(uint32 seed) {   // Constructor
      RandomInit(seed);}
protected:
   uint32 x[5];                        // History buffer
};

#endif

