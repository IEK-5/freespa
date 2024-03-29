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
            double lat, double e)
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
	spa.pressure=1010;
	spa.temperature=10;
	spa.slope=0;
	spa.azm_rotation=0;

	spa.atmos_refract=Bennet(1010, 10, 0)*180/M_PI;
	spa.function=SPA_ZA;
	r=spa_calculate(&spa);
	if (r)
		fprintf(stderr,"NREL spa returned %d\n", r);
	
	P.z=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.a=fmod(M_PI*spa.azimuth/180, 2*M_PI);
	P.z+=M_PI*spa.del_e/180;
	return P;
}

// simple wrapper around NREL's spa code
sol_pos SPA_WrapperApp(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
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
	r=spa_calculate(&spa);
	if (r)
		fprintf(stderr,"NREL spa returned %d\n", r);
	
	P.z=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.a=fmod(M_PI*spa.azimuth/180, 2*M_PI);
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
	return 250*((double)rand()/(double)(RAND_MAX))+850;
} 
double RandT()
{
	return 40*((double)rand()/(double)(RAND_MAX))-10;
} 
// spa is allegibly accurate to 0.0003 degrees, i.e. about 5e-6 radians
// small numerical differences from constant conversions are OK
#define RAD_EPS 5e-7
int SpecificTester(time_t tc, double lat, double lon, double E, double p, double T, int verb, double *err)
{
	sol_pos P1, P2, P1a, P2a;
	int ERR=0;
	double d, dt, dip=0;
	char* timestr;
	struct tm *ut;
	ut=gmjtime(&tc);
	dt=get_delta_t(ut);
	if (dt>=8000) // NREL's spa has this limit, delta t does not ...
		dt=0;
	P1=SPA(ut, &dt, 0, lon,  lat, E);
	P1a=ApSolposBennet(P1,&dip,E,p,T); // NREL spa uses Bennet
	
	P2=SPA_Wrapper(ut, &dt , 0, lon,  lat, E);
	P2a=SPA_WrapperApp(ut, &dt , 0, lon,  lat, E, p, T);
	//P2a=AparentSolpos(P2,&dip,E,p,T);
	
	d=AngleBetween(P1.z, P1.a, P2.z, P2.a); // note that simply comparing azimuth and zenith has a problem for small zenith angles
								  // (azimuth has no effect for zenith=0) thus we compute the angle between the two 
								  // solar vectors	
	(*err)=fabs(d);
	if ((fabs(d)>RAD_EPS)||verb)
	{
		timestr=malloc(50*sizeof(char));
		strftime(timestr, 50, "%Y/%m/%d %T %Z",ut);
		printf("%s: %ld %.12e %.12e %.12e %.12e\n", timestr, tc, lat, lon,p,d);
		printf("r zenith:  %e\t%e\t%e\n", P1.z, P2.z, P1.z-P2.z);
		printf("r azimuth: %e\t%e\t%e\n", P1.a, P2.a, P1.a-P2.a);
		free(timestr);
		ERR=1;
	}
	d=AngleBetween(P1a.z, P1a.a, P2a.z, P2a.a); // note that simply comparing azimuth and zenith has a problem for small zenith angles
								  // (azimuth has no effect for zenith=0) thus we compute the angle between the two 
								  // solar vectors
	if (fabs(d)>(*err))	
		(*err)=fabs(d);
		
	if ((fabs(d)>RAD_EPS)||verb)
	{
		timestr=malloc(50*sizeof(char));
		strftime(timestr, 50, "%Y/%m/%d %T %Z",ut);
		printf("%s: %ld %.12e %.12e %.12e\n", timestr, tc, lat, lon,d);
		printf("a zenith:  %e\t%e\t%e\n", P1a.z, P2a.z, P1a.z-P2a.z);
		printf("a azimuth: %e\t%e\t%e\n", P1a.a, P2a.a, P1a.a-P2a.a);
		free(timestr);
		ERR=1;
	}
	return ERR;
}
int RandomTester(double *err)
{
	time_t tc;
	double lat, lon, E, T,p;
	lat=RandLat();
	lon=RandLon();
	E=RandE();
	T=RandT();
	p=Randp();
	tc=RandEpoch();
	return SpecificTester(tc, lat, lon,E,p,T, 0, err);
}
#define LAT 50.902996388
#define LON 6.407165038
#define N 10000 // number of coordinates to test
#define NN 10000 // number of performance test iterations
#define NNN 10000 // number of solar times to test
// benchmark routine speed
double Perf(int Nc, sol_pos (*sparoutine)(struct tm *, double *, double, double, double, double))
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
		P=sparoutine(ut, &dt, 0, lon,  lat, 0);
		s+=P.z; // do not optimize this loop out
	}
	printf("bogus number %e\n",s);
	return TOC(&t);
}


int main()
{
	double t, ee, err;
	time_t tc;
	char *curtz = getenv("TZ"); // Make a copy of the timezone variable
	char *old=NULL;
	int i;
	int NE=0;
	struct tm ut={0}, *p;
	tc=1652254722;
	p=gmjtime_r(&tc, &ut);
	
	if (curtz)
		old=strdup(curtz);
    setenv("TZ", ":/usr/share/zoneinfo/Etc/UTC", 1); // always use UTC
    tzset();
	
	
	
	srand((unsigned) time(&tc));
	printf("--------------------------------------\n");
	printf("Testing freespa against NREL spa for \n%d random inputs:\n", N);
	TIC(&t);
	err=0;
	for (i=0;i<N;i++)
	{
		if ((RandomTester(&ee)))
		{
			fprintf(stderr,"Error: the spa's do not match!\n");
			NE++;
		}
		if (ee>err)
			err=ee;
	}
	printf("%d errors\n", NE);
	printf("Max deviation %.8f°\n", rad2deg(err));
	t=TOC(&t);
	printf("used %f s (%.1f us/test)\n", t, 1e6*t/N);
	printf("--------------------------------------\n\n");
	
	printf("--------------------------------------\n");
	printf("Benchmarking the spa's\n");
	t=Perf(NN, &SPA);
	printf("freespa  used %f s (%.1f us/call)\n", t, 1e6*t/NN);
	t=Perf(NN, &SPA_Wrapper);
	printf("NREL spa used %f s (%.1f us/call)\n", t, 1e6*t/NN);
	printf("--------------------------------------\n\n");
	
	
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
