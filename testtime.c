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

#define deg2rad(a) (M_PI*(a)/180.0)
#define rad2deg(a) (180.0*(a)/M_PI)
#define SUN_RADIUS 4.6542695162932789e-03 // in radians

//#define MIN_EPOCH -125197920000 // year -2000
#define MIN_EPOCH -12219292800 // Oct 15 1582
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
int comparetm(struct tm *p1, struct tm *p2)
{
	int r=0;
	if (p1->tm_year!=p2->tm_year)
		r|=(1);
	if (p1->tm_mon!=p2->tm_mon)
		r|=(1<<1);
	if (p1->tm_mday!=p2->tm_mday)
		r|=(1<<2);
	if (p1->tm_hour!=p2->tm_hour)
		r|=(1<<3);
	if (p1->tm_min!=p2->tm_min)
		r|=(1<<4);
	if (abs(p1->tm_sec-p2->tm_sec)>1) // allow for rounding errors
		r|=(1<<4);
	return r;	
}
#define N 10000
int TestGM()
{
	int i;
	struct tm ut1={0};
	struct tm ut2={0};
	struct tm *p1, *p2;
	char buffer [80];
	time_t tc;
	srand((unsigned) time(&tc));
	tc=0;
	i=0;
	for (i=0;i<N;i++)
	{
		tc=RandEpoch();
		p1=gmtime_r(&tc, &ut1);
		p2=gmjtime_r(&tc, &ut2);
		if (comparetm(p1, p2))
			break;
	}
	if (comparetm(p1, p2))
	{
		printf("%ld\n", tc);
		strftime (buffer,80,"Time UTC: %Y-%m-%d %H:%M:%S",&ut1);
		puts(buffer);
		strftime (buffer,80,"fSpa UTC: %Y-%m-%d %H:%M:%S",&ut2);
		puts(buffer);
		return 1;
	}
	else
		printf("All OK\n");
	return 0;
}
int TestMK()
{
	int i;
	struct tm ut={0};
	struct tm *p;
	char buffer [80];
	time_t tc, tc2;
	srand((unsigned) time(&tc));
	i=0;
	for (i=0;i<N;i++)
	{
		tc=RandEpoch();
		p=gmjtime_r(&tc, &ut);
		tc2=mkgmjtime(p);
		
		if (tc!=tc2)
			break;
	}
	if (tc!=tc2)
	{
		printf("i: %d\n", i);
		printf("in:  %ld\n", tc);
		printf("out: %ld\n", tc2);
		p=gmjtime_r(&tc, &ut);
		strftime (buffer,80,"In  UTC: %Y-%m-%d %H:%M:%S",&ut);
		puts(buffer);
		p=gmjtime_r(&tc2, &ut);
		strftime (buffer,80,"Out UTC: %Y-%m-%d %H:%M:%S",&ut);
		puts(buffer);
		return 1;
	}
	else
		printf("All OK\n");
	return 0;
}

#define ST_EPS deg2rad(2)

int TimeTester(time_t tc, double lat, double lon, double e)
{
	int R=0;
	char buffer [80];
	struct tm ut={0};
	struct tm sunset={0}, sunrise={0}, transit={0};
	
	gmjtime_r(&tc, &ut);
	
	R=SunTimes(ut, NULL, 0, lon, lat, e,1010.0, 10.0, &sunrise, &transit, &sunset);
	//R=SunTimes_(ut, NULL, 0, lon, lat,1010.0, 10.0, &sunrise, &transit, &sunset);
	if (R==0)
	{
		sol_pos P;
		double E;
		P=SPA(&sunrise, NULL, 0, lon,  lat, e, 1010, 10);
		E=P.az-M_PI/2-SUN_RADIUS;
		if (fabs(E)>ST_EPS)
		{
			strftime (buffer,80,"UTC:     %Y-%m-%d %H:%M:%S",&ut);
			puts(buffer);		
			printf("%ld %e %e\n", tc, lon, lat);
			printf("%.3fE %.3fN\n", rad2deg(lon), rad2deg(lat));
			strftime (buffer,80,"sunrise: %Y-%m-%d %H:%M:%S",&sunrise);
			puts(buffer);
			printf("Sunrise at %.3g degrees\n\n", rad2deg(E));
			R|=1;
		}	
		P=SPA(&sunset, NULL, 0, lon,  lat, e, 1010, 10);
		E=P.az-M_PI/2-SUN_RADIUS;
		if (fabs(E)>ST_EPS)
		{
			strftime (buffer,80,"UTC:     %Y-%m-%d %H:%M:%S",&ut);
			puts(buffer);		
			printf("%ld %e %e\n", tc, lon, lat);
			printf("%.3fE %.3fN\n", rad2deg(lon), rad2deg(lat));
			strftime (buffer,80,"sunset: %Y-%m-%d %H:%M:%S",&sunset);
			puts(buffer);	
			printf("Sunset at %.3g degrees\n\n", rad2deg(E));
			R|=2;
		}	
		P=SPA(&transit, NULL, 0, lon,  lat, e, 1010, 10);
		E=fabs(P.aa-M_PI);
		if (E>M_PI/2)
			E=fabs(E-M_PI);
		if ((P.aa>0)&&(P.aa<M_PI))
			E=-E;
		E*=fabs(sin(P.az));
		if (fabs(E)>ST_EPS)
		{
			strftime (buffer,80,"UTC:     %Y-%m-%d %H:%M:%S",&ut);
			puts(buffer);		
			printf("%ld %e %e\n", tc, lon, lat);
			printf("%.3fE %.3fN\n", rad2deg(lon), rad2deg(lat));
			strftime (buffer,80,"transit: %Y-%m-%d %H:%M:%S",&transit);
			puts(buffer);		
			printf("SunTransit at %.3g degrees off\n\n", rad2deg(E));
			R|=4;
		}		
	}
	else 
	{	
		sol_pos P;
		double E;
		R=0;
		P=SPA(&transit, NULL, 0, lon,  lat, e, 1010, 10);
		E=fabs(P.aa-M_PI);
		if (E>M_PI/2)
			E=fabs(E-M_PI);
		if ((P.aa>0)&&(P.aa<M_PI))
			E=-E;
		E*=fabs(sin(P.az));
		if (fabs(E)>ST_EPS)
		{
			strftime (buffer,80,"UTC:     %Y-%m-%d %H:%M:%S",&ut);
			puts(buffer);		
			printf("%ld %e %e\n", tc, lon, lat);
			printf("%.3fE %.3fN\n", rad2deg(lon), rad2deg(lat));
			
			strftime (buffer,80,"transit: %Y-%m-%d %H:%M:%S",&transit);
			puts(buffer);		
			printf("SunTransit at %.3g degrees off\n\n", rad2deg(E));
			R|=4;
		}	
	}
	return R;
}

double RandLon()
{
	return 2*M_PI*((double)rand()/(double)(RAND_MAX))-M_PI;
} 
double RandLat()
{
	return M_PI*((double)rand()/(double)(RAND_MAX))-M_PI/2;
} 
#define NN 1000
int TestSuntimes()
{
	int i, e=0;
	struct tm ut={0};
	double lat, lon;
	time_t tc;
	srand((unsigned) time(&tc));
	tc=0;
	i=0;
	for (i=0;i<NN;i++)
	{
		tc=RandEpoch();
		gmtime_r(&tc, &ut);
		lat=RandLat();
		lon=RandLon();
		if (TimeTester(tc, lat, lon, 0))
			e++;
	}
	if (e==0)
		printf("All OK\n");
	else
		printf("%d errors in %d runs\n", e, NN);
	return 0;
}


int main()
{
	struct tm ut={0};
	ut.tm_year=-2000-1900;
	ut.tm_mday=1;
	printf("The beginning: %ld\n", mkgmjtime(&ut));
	ut.tm_year=1582-1900;
	ut.tm_mon=9;
	ut.tm_mday=15;
	printf("The Trouble: %ld\n", mkgmjtime(&ut));
	
	
	printf("Test gmjtime:\n");
	if (!TestGM())
	{
		printf("Test mkgmjtime:\n");
		TestMK();
	}
		
	TestSuntimes();
	
	exit(0);
}
