#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "freespa.h"
#include "spa.h"

//#define _DEBUG_
void TIC(double *r)
{
	(*r)=(double)clock();
}
double TOC(double *r)
{
	return ((double)clock()-(*r))/CLOCKS_PER_SEC;
}


double AngleBetween(sky_pos p1, sky_pos p2)
{
	double p=fmod(p1.a-p2.a,2*M_PI);
	double f;
	f=cos(p1.z)*cos(p2.z)+sin(p1.z)*sin(p2.z)*cos(p);
	return acos(f);
}



sol_pos SPA_Wrapper(time_t t, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double a_refr, double T)
{
	spa_data spa;	
	struct tm *ut;
	sol_pos P;
	
	
	ut=gmtime(&t);	
	if (!ut)
	{
		P.E=1;
	}
	P.E=0;
	
	spa.year=ut->tm_year+1900;
	spa.month=ut->tm_mon+1;
	spa.day=ut->tm_mday;
	spa.hour=ut->tm_hour;
	spa.minute=ut->tm_min;
	spa.second=(double)ut->tm_sec;
	spa.delta_t=delta_t;
	spa.delta_ut1=delta_ut1;
	spa.timezone=0;
	spa.longitude=180.0*lon/M_PI;
	spa.latitude=180.0*lat/M_PI;
	spa.elevation=e;
	spa.pressure=p;
	spa.temperature=T;
	spa.slope=0;
	spa.azm_rotation=0;
	spa.atmos_refract=180.0*a_refr/M_PI;
	spa.function=SPA_ZA;
	spa_calculate(&spa);
#ifdef _DEBUG_
	#define PPAR(name,par) (printf("%s:\t%.12e\n",name,deg2rad(par)))
	PPAR("x[0]",spa.x0);
	PPAR("x[1]",spa.x1);
	PPAR("x[2]",spa.x2);
	PPAR("x[3]",spa.x3);
	PPAR("x[4]",spa.x4);
	printf("JD.JCE:\t%.12e\n",spa.jce);
	PPAR("dtau",spa.del_tau);
	PPAR("dpsi",spa.del_psi);
	PPAR("deps",spa.del_epsilon);
	PPAR("eps",spa.epsilon);
	PPAR("v",spa.nu);
	PPAR("alpha",spa.alpha);
	PPAR("delta",spa.delta);
	PPAR("H",spa.h);
	PPAR("delta_prime",spa.delta_prime);
	PPAR("H_prime",spa.h_prime);
	PPAR("ele",spa.e0);
	PPAR("dele",spa.del_e);
#endif		
	
	P.sa.z=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.sa.a=fmod(M_PI*spa.azimuth/180, 2*M_PI);
	P.s=P.sa;
	P.s.z+=M_PI*spa.del_e/180;
	return P;
}

#define MIN_EPOCH -125197920000 // year -2000
#define MAX_EPOCH 127090080000 // year +6000
time_t RandEpoch()
{
	/* generate a random epoch between MIN_EPOCH and MAX_EPOCH
	 * Unfortunarely this requires a minimum 38 bit integer
	 * where RAND_MAX is often only 31 bits...
	 */
	 time_t rt;
	 rt=(time_t)rand();
	 rt<<=31;
	 rt|=(time_t)rand();
	 rt%=(MAX_EPOCH-MIN_EPOCH);
	 rt+=MIN_EPOCH;
	 return rt;
}
double RandLon()
{
	return 2*M_PI*((double)rand()/(double)(RAND_MAX));
} 
double RandLat()
{
	return M_PI*((double)rand()/(double)(RAND_MAX))-M_PI/2;
} 
// spa is accurate to 0.0003 degrees, i.e. about 5e-6 radians
// small numerical differences from constant conversions are OK
#define RAD_EPS 2e-7
int SpecificTester(time_t tc, double lat, double lon, int verb)
{
	sol_pos P1, P2;
	double d;
	char* timestr;
	struct tm *ut;
	P1=SPA_NREL(tc, 0, 0, M_PI*lon/180,  M_PI*lat/180, 0, 1010, M_PI*0.5667/180, 20);
	P2=SPA_Wrapper(tc, 0, 0, M_PI*lon/180,  M_PI*lat/180, 0, 1010, M_PI*0.5667/180, 20);
	d=AngleBetween(P1.sa, P2.sa); // note that simply comparing azimuth and zenith has a problem for small zenith angles
								  // (azimuth has no effect for zenith=0) thus we compute the angle between the two 
								  // solar vectors
	if ((fabs(d)>RAD_EPS)||verb)
	{
		timestr=malloc(50*sizeof(char));
		ut=gmtime(&tc);	
		strftime(timestr, 50, "%Y/%m/%d %T %Z",ut);
		printf("%s: %ld %.12e %.12e %.12e\n", timestr, tc, lat, lon,d);
		printf("zenith:  %e\t%e\t%e\n", P1.sa.z, P2.sa.z, P1.sa.z-P2.sa.z);
		printf("azimuth: %e\t%e\t%e\n\n", P1.sa.a, P2.sa.a, P1.sa.a-P2.sa.a);
		free(timestr);
			return 1;
	}
	return 0;
}
int RandomTester()
{
	time_t tc;
	double lat, lon, d;
	sol_pos P1, P2;
	lat=RandLat();
	lon=RandLon();
	tc=RandEpoch();
	return SpecificTester(tc, lat, lon,0);
}
#define LAT 50.902996388
#define LON 6.407165038
#define N 100000000
#define NN 100000
// benchmark routine speed
double Perf(int Nc, sol_pos (*sparoutine)(time_t, double, double, double, double , double, double, double , double))
{
	int i;
	double t;
	time_t tc;
	double lat, lon, s=0;
	sol_pos P;
	TIC(&t);
	lon=RandLon();
	lat=RandLat();// only one random parameter in loop, saves time
	tc=RandEpoch();
	for (i=0;i<Nc;i++)
	{
		//lat=RandLat();// only one random parameter in loop, saves time
		P=sparoutine(tc, 0, 0, M_PI*lon/180,  M_PI*lat/180, 0, 1010, M_PI*0.5667/180, 20);
		s+=P.sa.z; // do not optimize this loop out
	}
	printf("bogus number %e\n",s);
	return TOC(&t);
}
int main()
{
	sol_pos P1, P2;
	spa_data spa;
	double t, e;
	time_t tc;
	char *curtz = getenv("TZ"); // Make a copy of the timezone variable
	char *old=NULL;
	int sum=0,i;
	int NE=0;
	
	if (curtz)
		old=strdup(curtz);
    setenv("TZ", ":/usr/share/zoneinfo/Etc/UTC", 1); // always use UTC
    tzset();
	
	
	if (!testjulian())
		printf("julian date routines OK\n");
	if (!testheliocentricpos())
		printf("heliocentric coordinates OK\n");
		
	SpecificTester(123631357169,-5.410342009552e-01,4.542216060972e+00,1);

	srand((unsigned) time(&tc));
	//srand(0);
	TIC(&t);
	for (i=0;i<N;i++)
	{
		if (RandomTester())
		{
			fprintf(stderr,"Error: spa mismachtch\n");
			NE++;
		}
	}
	printf("tested %d coodinates and times\n", N);
	printf("%d errors\n", NE);
	t=TOC(&t);
	printf("used %f s (%.1f us/test)\n", t, 1e6*t/N);
	//exit(0);
	t=Perf(NN, &SPA_NREL);
	printf("used %f s (%.1f us/call)\n", NN, t, 1e6*t/NN);
	t=Perf(NN, &SPA_Wrapper);
	printf("used %f s (%.1f us/call)\n", NN, t, 1e6*t/NN);
	
    if (old)
    {
		printf("%s\n", old);
		setenv("TZ", old, 1); // Restore old PATH
		free(old); // Don't forget to free!
	}
	else
		unsetenv("TZ");
    tzset();
	return 0;
}
