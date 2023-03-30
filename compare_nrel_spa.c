/*
    freespa
    Copyright (C) 2022  B. E. Pieters,
    IEK-5 Photovoltaik, Forschunszentrum Juelich

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program.  If not, see
    <http://www.gnu.org/licenses/>.
*/
/*****************************************************************
 *  INSTITUT FUER ENERGIE- UND KLIMAFORSCHUNG                    *
 +  IEK-5 PHOTOVOLTAIK                                           *
 *                                                               *
 *        ########                _   _                          *
 *     ##########                |_| |_|                         *
 *    ##########     ##         _ _   _ _     ___ ____ _   _     *
 *   ##########     ####       | | | | | |   |_ _/ ___| | | |    *
 *   #########     #####    _  | | | | | |    | | |   | |_| |    *
 *   #    ###     ######   | |_| | |_| | |___ | | |___|  _  |    *
 *    #          ######     \___/ \___/|_____|___\____|_| |_|    *
 *     ##      #######      F o r s c h u n g s z e n t r u m    *
 *       ##########                                              *
 *                                                               *
 *   http://www.fz-juelich.de/iek/iek-5/DE/Home/home_node.html   *
 *****************************************************************
 *                                                               *
 *    Dr. Bart E. Pieters 2022                                   *
 *                                                               *
 *****************************************************************/
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
#define EARTH_R 6378136.6 //m at the equator [IERS Numerical Standards (IAG 1999)] 

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
sol_pos SPA_Wrapper(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T)
{
	spa_data spa;
	sol_pos P;
	int r;
	
	
	P.E=0;
	
	spa.year=ut->tm_year+1900;
	spa.month=ut->tm_mon+1;
	spa.day=ut->tm_mday;
	spa.hour=ut->tm_hour;
	spa.minute=ut->tm_min;
	spa.second=(double)ut->tm_sec;
	
/* 
 * Morrison, & Stephenson (2004) provide the following equation for Δt
 * 
 *                    2
 *            ⎡y-1820⎤
 * Δt(y) = 32 ⎢──────⎥  - 20
 *            ⎣ 100  ⎦
 * 
 * However, for some reason NREL spa limits Δt to rather small values 
 * (Δt<8000). This limits the validity of NREL spa to approximately the 
 * years between 245 -- 3400. This is still a large range but 
 * considerably smaller than the claimed -2000 -- 6000. 
 * 
 * If you use NREL spa, and need a broader range than 245 -- 3400, you 
 * should remove the 8000 limit (see the "validate_inputs" routine). 
 * 
 * Morrison, L. V., & Stephenson, F. R. (2004). Historical Values of 
 * the Earth’s Clock Error ΔT and the Calculation of Eclipses. 
 * Journal for the History of Astronomy, 35(3), 327–336. 
 * https://doi.org/10.1177/002182860403500305
 * 
 * Ps. freespa does not limit Δt, so knock your self out :)
 */	
	if (delta_t)
		spa.delta_t=*delta_t;
	else
		spa.delta_t=get_delta_t(ut);
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
	if (r)
		fprintf(stderr,"NREL spa returned %d\n", r);
	P.az=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.aa=fmod(M_PI*spa.azimuth/180, 2*M_PI);
	P.a=P.aa;
	P.z=P.az;
	P.z+=M_PI*spa.del_e/180;
	return P;
}

#define MIN_EPOCH -125282592000 // year -2000
#define MAX_EPOCH 127090080000 // year +6000

//#define MIN_EPOCH -54277908808 // year ~250
//#define MAX_EPOCH 45126500400 // year 3400

//#define MIN_EPOCH -12000000000 // year ~250
//#define MAX_EPOCH 45126500400 // year 3400

//#define MIN_EPOCH 946681200 // year 2000
//#define MAX_EPOCH 1640991600 // year 2022
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
	return 2*M_PI*((double)rand()/(double)(RAND_MAX))-M_PI;
} 
double RandLat()
{
	return M_PI*((double)rand()/(double)(RAND_MAX))-M_PI/2;
}  

// there is also elevations below 0
// I ignore that
double RandE()
{
	return 8000*((double)rand()/(double)(RAND_MAX));
} 
double Randp()
{
	/* sea level records:
	 * 870
	 * 1085
	 */
	return 250*((double)rand()/(double)(RAND_MAX))-850;
} 
double RandT()
{
	return 40*((double)rand()/(double)(RAND_MAX))-10;
} 
// spa is allegibly accurate to 0.0003 degrees, i.e. about 5e-6 radians
// small numerical differences from constant conversions are OK
#define RAD_EPS 2e-7
int SpecificTester(time_t tc, double lat, double lon, double E, double p, double T, int verb)
{
	sol_pos P1, P2;
	double d, dt;
	char* timestr;
	struct tm *ut;
	ut=gmjtime(&tc);
	dt=get_delta_t(ut);
	if (dt>=8000) // NREL's spa has this limit, delta t does not ...
		dt=0;
	P1=        SPA(ut, &dt, 0, lon,  lat, E, p, T);
	P2=SPA_Wrapper(ut, &dt , 0, lon,  lat, E, p, T);
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
	double lat, lon, E, T,p;
	lat=RandLat();
	lon=RandLon();
	E=RandE();
	T=RandT();
	p=Randp();
	tc=RandEpoch();
	return SpecificTester(tc, lat, lon,E,p,T, 0);
}
#define LAT 50.902996388
#define LON 6.407165038
#define N 10000 // number of coordinates to test
#define NN 10000 // number of performance test iterations
#define NNN 10000 // number of solar times to test
// benchmark routine speed
double Perf(int Nc, sol_pos (*sparoutine)(struct tm *, double *, double, double, double, double, double , double))
{
	int i;
	double t, dt;
	time_t tc;
	double lat, lon, s=0;
	sol_pos P;
	struct tm *ut;
	TIC(&t);
	lon=RandLon();
	lat=RandLat();// only one random parameter in loop, saves time
	tc=RandEpoch();
	ut=gmjtime(&tc);
	dt=0;
	for (i=0;i<Nc;i++)
	{
		P=sparoutine(ut, &dt, 0, lon,  lat, 0, 1010, 10);
		s+=P.az; // do not optimize this loop out
	}
	printf("bogus number %e\n",s);
	return TOC(&t);
}


/* test time routines */

// simple wrapper around NREL's spa code
int SPA_SunTimes(struct tm ut, double *delta_t, double delta_ut1, 
			double lon, double lat, double e, double p, double T, 
			struct tm *sunrise, struct tm *transit, struct tm *sunset)
{
	spa_data spa;
	time_t t, t0, tt;	
	
	spa.year=ut.tm_year+1900;
	spa.month=ut.tm_mon+1;
	spa.day=ut.tm_mday;
	spa.hour=ut.tm_hour;
	spa.minute=ut.tm_min;
	spa.second=(double)ut.tm_sec;
	if (delta_t)
		spa.delta_t=*delta_t;
	else
		spa.delta_t=get_delta_t(&ut);
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
	spa.function=SPA_ZA_RTS;
	spa_calculate(&spa);
	if ((spa.sunrise<0)||(spa.sunset<0)||(spa.suntransit<0))
	{
		if ((spa.latitude<67.9)&&(spa.latitude>-67.9))
			return 10;
		if (spa.zenith<90)
			return 1;		
		return -1;
	}
	
	ut.tm_hour=0;
	ut.tm_min=0;
	ut.tm_sec=0;
	t0=mkgmjtime(&ut);
	
	tt=t0+(time_t)(3600*spa.suntransit);
	transit=gmjtime_r(&tt, transit);
	
	t=t0+(time_t)(3600*spa.sunrise);
	if (t>tt)
		t-=86400;
	sunrise=gmjtime_r(&t, sunrise);
	t=t0+(time_t)(3600*spa.sunset);
	if (t<tt)
		t+=86400;
	sunset=gmjtime_r(&t, sunset);
	// TODO match freespa return value
	return 0;
}
#define SUN_RADIUS 4.6542695162932789e-03 // in radians

// this value is OK for freespa but counts most NREL spa results as wrong
// #define XEPS deg2rad(0.05) // in radians

// this value should be OK for both
#define XEPS deg2rad(1.0) // in radians
// note that this value is the result of just trying what is still OK
// it is in fact much larger than 

#define SRF 1	//sun rise freespa
#define SSF 2	//sun set freespa
#define STF 4	//sun transit freespa
#define SRN 8	//sun rise NREL spa
#define SSN 16	//sun set NREL spa
#define STN 32	//sun transit NREL spa

int TimeTester(time_t tc, double lat, double lon, double E, int verb)
/* tests whether sunset/sunrise/transit times are accurate
 * Considers errors in radians, not seconds
 * The errors are computed using the own solar position routines
 * return value integer. bits 1&2: number of errors in freespa
 *                       bits 3&4: number of errors in NREL spa
 * 
 * Note:  NREL's routines ignore elevation, freespa does not
 * Note2: freespa does a more expensive optimization than NREL spa
 *        and is generally more accurate.
 */
{
	int R=0, r1, r2;
	double dt;
	sol_pos P;
	double E1, E2;
	struct tm ut;
	struct tm *p;
	struct tm sunset_nrel, sunrise_nrel, transit_nrel;
	struct tm sunset_free, sunrise_free, transit_free;
	double dipg;
	
	// freespa takes the geometric dip into account
	
	dipg=acos(EARTH_R/(EARTH_R+E));
	
	p=gmjtime_r(&tc, &ut);
	dt=get_delta_t(p);
	//if (dt>=8000) // NREL's spa has this limit, delta t does not ...
	//	dt=0;
	
	r1=SunTimes(ut, &dt, 0, lon, lat, E, 1010.0, 10.0, &sunrise_free, &transit_free, &sunset_free);
	r2=SPA_SunTimes(ut, &dt, 0, lon, lat, E, 1010.0, 10.0, &sunrise_nrel, &transit_nrel, &sunset_nrel);
	
	// sunrise & set	
	if (r1==0)
	{
		P=SPA(&sunrise_free, &dt, 0, lon,  lat, E, 1010, 10);
		E1=P.az-M_PI/2-SUN_RADIUS-dipg;
		if (fabs(E1)>XEPS)
			R|=SRF;
			
		P=SPA(&sunset_free, &dt, 0, lon,  lat, E, 1010, 10);
		E1=P.az-M_PI/2-SUN_RADIUS-dipg;
		if (fabs(E1)>XEPS)
			R|=SSF;
	}
	if (r2==0)
	{
		P=SPA_Wrapper(&sunrise_nrel, &dt, 0, lon,  lat, E, 1010, 10);
		E2=P.az-M_PI/2-SUN_RADIUS;
		if (fabs(E2)>XEPS)
			R|=SRN;
		P=SPA_Wrapper(&sunset_nrel, &dt, 0, lon,  lat, E, 1010, 10);
		E2=P.az-M_PI/2-SUN_RADIUS;
		if (fabs(E2)>XEPS)
			R|=SSN;		
	}
	
	//transit
	P=SPA(&transit_free, &dt, 0, lon,  lat, E, 1010, 10);
	E1=atan(sin(P.aa)*fabs(tan(P.az)));
	P=SPA_Wrapper(&transit_nrel, &dt, 0, lon,  lat, E, 1010, 10);
	E2=atan(sin(P.aa)*fabs(tan(P.az)));
	if (fabs(E1)>XEPS)
		R|=STF;
	if (fabs(E2)>XEPS)
		R|=STN;
	
	return R;
}


/* minimum and maximum latitudes
 * Computing sunrise/set times is hard near the poles where the
 * sub may cross the horizon at a very shallow angle (or not at all)
 * For this reason we limit the range a bit
 */

#define MAXLAT 65
#define MINLAT -65

int RandomTimeTester()
{
	time_t tc;
	double lat, lon, E;
	lat=deg2rad(MINLAT+(MAXLAT-MINLAT)*((double)rand()/(double)(RAND_MAX)));
	lon=RandLon();
	E=RandE();
	tc=RandEpoch();
	return TimeTester(tc, lat, lon,E,0);
}
int main()
{
	double t;
	time_t tc;
	char *curtz = getenv("TZ"); // Make a copy of the timezone variable
	char *old=NULL;
	int i, r;
	int NE=0;
	int SREF=0;
	int SSEF=0;
	int STEF=0;
	int SREN=0;
	int SSEN=0;
	int STEN=0;
	struct tm ut={0}, *p;
	tc=1652254722;
	p=gmjtime_r(&tc, &ut);
	
	if (curtz)
		old=strdup(curtz);
    setenv("TZ", ":/usr/share/zoneinfo/Etc/UTC", 1); // always use UTC
    tzset();
	
	
	
	srand((unsigned) time(&tc));
	printf("---------------------------------\n");
	printf("Testing freespa against NREL spa for %d random inputs:\n", N);
	TIC(&t);
	for (i=0;i<N;i++)
	{
		if ((RandomTester()))
		{
			fprintf(stderr,"Error: the spa's do not match!\n");
			NE++;
		}
	}
	printf("tested %d coodinates and times\n", N);
	printf("%d errors\n", NE);
	t=TOC(&t);
	printf("used %f s (%.1f us/test)\n", t, 1e6*t/N);
	printf("---------------------------------\n\n");
	
	printf("---------------------------------\n");
	printf("Benchmarking the spa's\n");
	t=Perf(NN, &SPA);
	printf("freespa  used %f s (%.1f us/call)\n", t, 1e6*t/NN);
	t=Perf(NN, &SPA_Wrapper);
	printf("NREL spa used %f s (%.1f us/call)\n", t, 1e6*t/NN);
	printf("---------------------------------\n\n");
	// TimeTester(1395402298,-0.6832,-2.9471, 1);
	// TimeTester(1508838461,0.2944,-3.0712, 1);
	printf("---------------------------------\n");
	printf("testing sunrise/transit/set times\n");
	TIC(&t);
	
	for (i=0;i<NNN;i++)
	{
		if (r=RandomTimeTester())
		{
			SREF+=((r&SRF)>0);
			SSEF+=((r&SSF)>0);
			STEF+=((r&STF)>0);
			SREN+=((r&SRN)>0);
			SSEN+=((r&SSN)>0);
			STEN+=((r&STN)>0);
		}
	}
	printf("tested %d coodinates and times\n", N);
	printf("freespa:  %10d sun-rise error > %.2e degrees\n", SREF, rad2deg(XEPS));
	printf("freespa:  %10d sun-set  error > %.2e degrees\n", SSEF, rad2deg(XEPS));
	printf("freespa:  %10d transit  error > %.2e degrees\n", STEF, rad2deg(XEPS));
	printf("NREL spa: %10d sun-rise error > %.2e degrees\n", SREN, rad2deg(XEPS));
	printf("NREL spa: %10d sun-set  error > %.2e degrees\n", SSEN, rad2deg(XEPS));
	printf("NREL spa: %10d transit  error > %.2e degrees\n", STEN, rad2deg(XEPS));
	t=TOC(&t);
	printf("---------------------------------\n\n");
	
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
