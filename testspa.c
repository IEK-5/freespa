#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "freespa.h"
#include "spa.h"

void TIC(double *r)
{
	(*r)=(double)clock();
}
double TOC(double *r)
{
	return ((double)clock()-(*r))/CLOCKS_PER_SEC;
}

double AngleBetween(double z1, double a1, double z2, double a2)
{
	double p=fmod(a1-a2,2*M_PI);
	double f;
	f=cos(z1)*cos(z2)+sin(z1)*sin(z2)*cos(p);
	return acos(f);
}

#define ABSOLUTEZERO -273.0 //convert C to K (aparently everyone ignores the .15)
#define AP0 1010.0 // standard sea level air pressure
#define AT0 10.0 // standard sea level air temperature
double Refr(const double coeff[], double p, double T, double h)
{
	//converts true and aparent solar elevation
	return (p/AP0)*((AT0-ABSOLUTEZERO)/(T-ABSOLUTEZERO))*coeff[0]/tan(h+coeff[1]/(h+coeff[2]));
}

double Bennet(double p, double T, double h)
{
	const double BENNET[] = {2.9088820866572158e-04,2.2267533386408395e-03,7.6794487087750510e-02};
	return Refr(BENNET, p, T, h);
}
// simple wrapper around NREL's spa code
sol_pos SPA_Wrapper(struct tm *ut, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T)
{
	spa_data spa;
	sol_pos P;
	
	
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
	spa.atmos_refract=Bennet(p, T, 0)*180/M_PI;
	spa.function=SPA_ZA;
	spa_calculate(&spa);
	
	P.az=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.aa=fmod(M_PI*spa.azimuth/180, 2*M_PI);
	P.a=P.aa;
	P.z=P.az;
	P.z+=M_PI*spa.del_e/180;
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
// spa is allegibly accurate to 0.0003 degrees, i.e. about 5e-6 radians
// small numerical differences from constant conversions are OK
#define RAD_EPS 2e-7
int SpecificTester(time_t tc, double lat, double lon, int verb)
{
	sol_pos P1, P2;
	double d;
	char* timestr;
	struct tm *ut;
	ut=gmtime(&tc);
		
	P1=        SPA(ut, 0, 0, M_PI*lon/180,  M_PI*lat/180, 0, 1010, 10);
	P2=SPA_Wrapper(ut, 0, 0, M_PI*lon/180,  M_PI*lat/180, 0, 1010, 10);
	d=AngleBetween(P1.az, P1.aa, P2.az, P2.aa); // note that simply comparing azimuth and zenith has a problem for small zenith angles
								  // (azimuth has no effect for zenith=0) thus we compute the angle between the two 
								  // solar vectors
	
	if ((fabs(d)>RAD_EPS)||verb)
	{
		timestr=malloc(50*sizeof(char));
		strftime(timestr, 50, "%Y/%m/%d %T %Z",ut);
		printf("%s: %ld %.12e %.12e %.12e\n", timestr, tc, lat, lon,d);
		printf("a zenith:  %e\t%e\t%e\n", P1.az, P2.az, P1.az-P2.az);
		printf("a azimuth: %e\t%e\t%e\n", P1.aa, P2.aa, P1.aa-P2.aa);
		printf("t zenith:  %e\t%e\t%e\n", P1.z, P2.z, P1.z-P2.z);
		printf("t azimuth: %e\t%e\t%e\n\n", P1.a, P2.a, P1.a-P2.a);
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
#define N 10
#define NN 1
// benchmark routine speed
double Perf(int Nc, sol_pos (*sparoutine)(struct tm *, double, double, double, double, double, double , double))
{
	int i;
	double t;
	time_t tc;
	double lat, lon, s=0;
	sol_pos P;
	struct tm *ut;
	TIC(&t);
	lon=RandLon();
	lat=RandLat();// only one random parameter in loop, saves time
	tc=RandEpoch();
	ut=gmtime(&tc);
	for (i=0;i<Nc;i++)
	{
		//lat=RandLat();// only one random parameter in loop, saves time
		P=sparoutine(ut, 0, 0, M_PI*lon/180,  M_PI*lat/180, 0, 1010, 10);
		s+=P.az; // do not optimize this loop out
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
	t=Perf(NN, &SPA);
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
