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
	
	
	P.E=0;
	
	spa.year=ut->tm_year+1900;
	spa.month=ut->tm_mon+1;
	spa.day=ut->tm_mday;
	spa.hour=ut->tm_hour;
	spa.minute=ut->tm_min;
	spa.second=(double)ut->tm_sec;
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
	P.az=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.aa=fmod(M_PI*spa.azimuth/180, 2*M_PI);
	P.a=P.aa;
	P.z=P.az;
	P.z+=M_PI*spa.del_e/180;
	return P;
}
//#define MIN_EPOCH -125197920000 // year -2000
//#define MIN_EPOCH -125282592000 // year -2000
//#define MAX_EPOCH 127090080000 // year +6000
#define MIN_EPOCH 946681200 // year 2000
#define MAX_EPOCH 1640991600 // year 2022
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
// spa is allegibly accurate to 0.0003 degrees, i.e. about 5e-6 radians
// small numerical differences from constant conversions are OK
#define RAD_EPS 2e-7
int SpecificTester(time_t tc, double lat, double lon, int verb)
{
	sol_pos P1, P2;
	double d, dt;
	char* timestr;
	struct tm *ut;
	ut=gmjtime(&tc);
	dt=get_delta_t(ut);
	if (dt>=8000) // NREL's spa has this limit, delta t does not ...
		dt=0;
	P1=        SPA(ut, &dt, 0, lon,  lat, 0, 1010, 10);
	P2=SPA_Wrapper(ut, &dt , 0, lon,  lat, 0, 1010, 10);
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
#define N 10000 // number of coordinates to test
#define NN 10000 // number of performancve test iterations
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
			double lon, double lat, double p, double T, 
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
	spa.elevation=0;
	spa.pressure=p;
	spa.temperature=T;
	spa.slope=0;
	spa.azm_rotation=0;
	spa.atmos_refract=Bennet(p, T, 0)*180/M_PI;
	spa.function=SPA_ZA_RTS;
	spa_calculate(&spa);
	if ((spa.sunrise<0)||(spa.sunset<0)||(spa.suntransit<0))
		return 1;
	
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
	return 0;
}
#define TIMEEPS 60
int TimeTester(time_t tc, double lat, double lon, int verb)
{
	int R=0;
	double dt;
	char buffer [80];
	struct tm ut;
	struct tm *p;
	struct tm sunset_nrel, sunrise_nrel, transit_nrel;
	struct tm sunset_free, sunrise_free, transit_free;
	time_t t1, t2;
	
	p=gmjtime_r(&tc, &ut);
	dt=get_delta_t(p);
	if (dt>=8000) // NREL's spa has this limit, delta t does not ...
		dt=0;
	
	R=SunTimes(ut, &dt, 0, lon, lat, 1010.0, 10.0, &sunrise_free, &transit_free, &sunset_free);
	if (R)
	{
		sol_pos P;
		double E;
		SPA_SunTimes(ut, &dt, 0, lon, lat, 1010.0, 10.0, &sunrise_nrel, &transit_nrel, &sunset_nrel);
		P=SPA(&sunrise_free, &dt, 0, lon,  lat, 0, 1010, 10);
		E=P.az-M_PI/2;
		R=0;
		if (fabs(E)>deg2rad(1))
		{
			printf("Sunrise at %.3g degrees\n", rad2deg(P.az-M_PI/2));
			R=1;
		}		
	}
	return R;
}

int RandomTimeTester()
{
	time_t tc;
	double lat, lon, d;
	sol_pos P1, P2;
	lat=RandLat();
	lon=RandLon();
	tc=RandEpoch();
	return TimeTester(tc, lat, lon,0);
}
int main()
{
	sol_pos P1, P2;
	spa_data spa;
	double t, e, dt, sr, ss, tr, h, m, s;
	time_t tc;
	char *curtz = getenv("TZ"); // Make a copy of the timezone variable
	char *old=NULL;
	int sum=0,i;
	int NE=0;
	struct tm ut={}, *p;
	tc=1652254722;
	p=gmjtime_r(&tc, &ut);
	dt=get_delta_t(p);
	
	if (curtz)
		old=strdup(curtz);
    setenv("TZ", ":/usr/share/zoneinfo/Etc/UTC", 1); // always use UTC
    tzset();
	
	
	
	srand((unsigned) time(&tc));
	printf("Testing freespa against NREL spa for %d random inputs\n", N);
	TIC(&t);
	for (i=0;i<N;i++)
	{
		if (RandomTester())
		{
			fprintf(stderr,"Error: the spa's do not match!\n");
			NE++;
		}
	}
	printf("tested %d coodinates and times\n", N);
	printf("%d errors\n", NE);
	t=TOC(&t);
	printf("used %f s (%.1f us/test)\n", t, 1e6*t/N);
	
	printf("Benchmarking the spa's\n", N);
	t=Perf(NN, &SPA);
	printf("freespa  used %f s (%.1f us/call)\n", NN, t, 1e6*t/NN);
	t=Perf(NN, &SPA_Wrapper);
	printf("NREL spa used %f s (%.1f us/call)\n", NN, t, 1e6*t/NN);
	TimeTester(1395402298,-0.6832,-2.9471, 1);
	TimeTester(1508838461,0.2944,-3.0712, 1);
	
	TIC(&t);
	
	for (i=0;i<N;i++)
	{
		if (RandomTimeTester())
		{
			fprintf(stderr,"Error: the spa's do not match!\n");
			NE++;
		}
	}
	printf("tested %d coodinates and times\n", N);
	printf("%d errors\n", NE);
	t=TOC(&t);
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
