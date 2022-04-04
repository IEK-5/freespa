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

int LogSPA(char *fn, int N)
{
	sol_pos P;
	int E=0;
	int i;
	struct tm ut, *p;
	FILE *f;
	time_t tc;
	double lat, lon;
	
	if ((f=fopen(fn,"w"))==NULL)
	{
		fprintf(stderr,"Error: failed to open %s for writing\n");
		return 1;
	}
	
	srand((unsigned) time(&tc));
	for (i=0;i<N;i++)
	{		
		lat=RandLat();
		lon=RandLon();
		tc=RandEpoch();
		p=gmtime_r(&tc, &ut);
		if (p)
		{
			P=SPA(p, 0, 0, lon,  lat, 0, 1010, 10);
			fprintf(f,"%16ld\t%19.12e\t%19.12e\t%19.12e\t%19.12e\t%19.12e\t%19.12e\n", tc, lat,lon,P.a,P.z,P.aa, P.az);
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
	sol_pos P, Pr;
	int Et=0, Ec=0;
	int i,j=0;
	char line[LLEN],c;
	struct tm ut, *p;
	FILE *f;
	time_t tc;
	double lat, lon, d;
	char *timestr;
	timestr=malloc(50*sizeof(char));
	
	if ((f=fopen(fn,"r"))==NULL)
	{
		fprintf(stderr,"Error: failed to open %s for reading\n");
		return 1;
	}
	while (!feof(f))
	{
		fgets(line,LLEN,f);
		if (feof(f))
			break;
		if (line[0]=='#')
			continue;
		i=sscanf(line,"%ld %lf %lf %lf %lf %lf %lf", &tc, &lat, &lon, &Pr.a, &Pr.z, &Pr.aa, &Pr.az);
		if (i==7)
		{
			p=gmtime_r(&tc, &ut);
			if (p)
			{
				j++;
				P=SPA(p, NULL, 0, lon,  lat, 0, 1010, 10);		
				d=AngleBetween(P.az, P.aa, Pr.az, Pr.aa);
				if (fabs(d)>RAD_EPS)
				{
					Et++;
					if (Et<10)
					{
						strftime(timestr, 50, "%Y/%m/%d %T %Z",p);
						fprintf(stderr, "%s (%ld)\nlon: %.12e\nlat: %.12e\n", timestr, tc, lon, lat);
						fprintf(stderr, "Solar vector deviates by %e rad\n", d);
						printf("a zenith:  %e\t%e\t%e\n", P.az, Pr.az, P.az-Pr.az);
						printf("a azimuth: %e\t%e\t%e\n", P.aa, Pr.aa, P.aa-Pr.aa);
						printf("t zenith:  %e\t%e\t%e\n", P.z, Pr.z, P.z-Pr.z);
						printf("t azimuth: %e\t%e\t%e\n\n", P.a, Pr.a, P.a-Pr.a);
					}
					else if (Et==10)
						fprintf(stderr, "10 or more deviations from reference!\n");
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
	return 0;
}
#define N 10000
int main(int argc, char **argv)
{
	int r=0;
	char *curtz = getenv("TZ"); // Make a copy of the timezone variable
	char *old=NULL;
	time_t tc;
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
