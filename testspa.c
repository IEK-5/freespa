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

void TIC(double *r)
{
	(*r)=(double)clock();
}
double TOC(double *r)
{
	return ((double)clock()-(*r))/CLOCKS_PER_SEC;
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
	return 2*M_PI*((double)rand()/(double)(RAND_MAX))-M_PI;
} 
double RandLat()
{
	return M_PI*((double)rand()/(double)(RAND_MAX))-M_PI/2;
} 
 
double RandE()
{
	return 8400*((double)rand()/(double)(RAND_MAX))-400;
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
int LogSPA(char *fn, int N)
{
	sol_pos P, Pa;
	int E=0;
	int i;
	struct tm ut, *p;
	FILE *f;
	time_t tc;
	double lat, lon, e, pp, T;
	
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Error: failed to open %s for writing\n", fn);
		return 1;
	}
	
	srand((unsigned) time(&tc));
	for (i=0;i<N;i++)
	{		
		lat=RandLat();
		lon=RandLon();
		tc=RandEpoch();
		e=RandE();
		pp=Randp();
		T=RandT();
		p=gmtime_r(&tc, &ut);
		if (p)
		{
			P=SPA(p, 0, 0, lon,  lat, e);
			Pa=AparentSolpos(P,NULL, e, pp, T);
			fprintf(f,"%16ld\t%19.16e\t%19.16e\t%19.16e\t%19.16e\t%19.16e\t%19.16e\t%19.16e\t%19.16e\t%19.16e\n", tc, lat,lon,e, pp, T, P.a,P.z,Pa.a, Pa.z);
		}
		else
		{
			E++;
			if (E<10)
			{
				fprintf(stderr,"Warning: failed to convert unix time %ld to broken down UTC\n", tc);
			}
			else if (E==10)
				fprintf(stderr,"More than 10 unix time conversion errors occurred!\n");			
		}
	}
	fclose(f);
	if (E)
		fprintf(stderr,"Warning: %d unix time conversion errors\n", E);
	return 0;
}

double AngleBetween(double z1, double a1, double z2, double a2)
{
	double p=fmod(a1-a2,2*M_PI);
	double f;
	f=cos(z1)*cos(z2)+sin(z1)*sin(z2)*cos(p);
	return acos(f);
}

#define LLEN 255
#define RAD_EPS 2e-7
int TestSPA(char *fn)
{
	sol_pos P, Pa, Pr, Pra;
	int Et=0, Ec=0;
	int i,j=0;
	char line[LLEN];
	struct tm ut, *p;
	FILE *f;
	time_t tc;
	double lat, lon, d, e, pp, T;
	char *timestr;
	timestr=malloc(50*sizeof(char));
	
	if ((f=fopen(fn,"r"))==NULL)
	{
		fprintf(stderr,"Error: failed to open %s for reading\n", fn);
		return 1;
	}
	while (!feof(f))
	{
		fgets(line,LLEN,f);
		if (feof(f))
			break;
		if (line[0]=='#')
			continue;
		i=sscanf(line,"%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tc, &lat, &lon, &e, &pp, &T, &Pr.a, &Pr.z, &Pra.a, &Pra.z);
		if (i==10)
		{
			p=gmtime_r(&tc, &ut);
			if (p)
			{
				j++;
				P=SPA(p, NULL, 0, lon,  lat, e);		
				d=AngleBetween(P.z, P.a, Pr.z, Pr.a);
				if (fabs(d)>RAD_EPS)
				{
					Et++;
					if (Et<10)
					{
						strftime(timestr, 50, "%Y/%m/%d %T %Z",p);
						fprintf(stderr, "%s (%ld)\nlon: %.12e\nlat: %.12e\n", timestr, tc, lon, lat);
						fprintf(stderr, "e: %.12e\np: %.12e\nT: %.12e\n", e,pp,T);
						fprintf(stderr, "Solar vector deviates by %e rad\n", d);
						printf("t zenith:  %e\t%e\t%e\n", P.z, Pr.z, P.z-Pr.z);
						printf("t azimuth: %e\t%e\t%e\n\n", P.a, Pr.a, P.a-Pr.a);
					}
					else if (Et==10)
						fprintf(stderr, "10 or more deviations from reference!\n");
				}
				else
				{	
					Pa=AparentSolpos(P,NULL,e,pp,T);	
					d=AngleBetween(Pa.z, Pa.a, Pra.z, Pra.a);
					if (fabs(d)>RAD_EPS)
					{
						Et++;
						if (Et<10)
						{
							strftime(timestr, 50, "%Y/%m/%d %T %Z",p);
							fprintf(stderr, "%s (%ld)\nlon: %.12e\nlat: %.12e\n", timestr, tc, lon, lat);
						fprintf(stderr, "e: %.12e\np: %.12e\nT: %.12e\n", e,pp,T);
							fprintf(stderr, "Solar vector deviates by %e rad\n", d);
							printf("a zenith:  %e\t%e\t%e\n", Pa.z, Pra.z, Pa.z-Pra.z);
							printf("a azimuth: %e\t%e\t%e\n\n", Pa.a, Pra.a, Pa.a-Pra.a);
						}
						else if (Et==10)
							fprintf(stderr, "10 or more deviations from reference!\n");
					}
				}
			}
			else
			{
				Ec++;
				if (Ec<10)
				{
					fprintf(stderr,"Warning: failed to convert unix time %ld to broken down UTC\n", tc);
				}
				else if (Ec==10)
					fprintf(stderr,"More than 10 unix time conversion errors occurred!\n");				
			}
		}
	}
	printf("Tested against %d reference solar vectors\n%d errors\n", j, Et);
	if (Ec)
		fprintf(stderr,"Warning: %d unix time conversion errors\n", Ec);
	free(timestr);
	fclose(f);
	return 0;
}
#define N 10000
int main(int argc, char **argv)
{
	int r=0;
	char *curtz = getenv("TZ"); // Make a copy of the timezone variable
	char *old=NULL;
	struct tm ut={0};
	
	if (argc!=3)
	{
		fprintf(stderr,"Usage: %s [r/t] <filename>\n", argv[0]);
		exit(1);
	}
	
	if (!testjulian())
		printf("julian day routines OK\n");
	if (!testheliocentricpos())
		printf("heliocentric coordinates OK\n");
		
		
	if (curtz)
		old=strdup(curtz);
    setenv("TZ", ":/usr/share/zoneinfo/Etc/UTC", 1); // always use UTC
    tzset();
    ut.tm_year=2022-1900;
    ut.tm_mon=0;
    ut.tm_mday=1;
    ut.tm_hour=0;
    ut.tm_min=0;
    ut.tm_sec=0;
    
    printf("delta t %e\n",get_delta_t(&ut));
	
	
	if (argv[1][0]=='r')
	{
		printf("Creating Reference File %s with %d data\n", argv[2], N);
		LogSPA(argv[2], N);
	} else if (argv[1][0]=='t')
	{
		printf("Testing against Reference File %s\n", argv[2]);
		TestSPA(argv[2]);
	} else
	{
		fprintf(stderr,"Usage: %s [r/t] <filename>\n", argv[0]);
		r=1;
	}
		
	
    if (old)
    {
		printf("%s\n", old);
		setenv("TZ", old, 1); // Restore old PATH
		free(old); // Don't forget to free!
	}
	else
		unsetenv("TZ");
    tzset();
	return r;
}
