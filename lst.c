#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "freespa.h"


int main(int argc, char **argv)
{
	double lon,lat;
	time_t tc;
	struct tm ut={0};
	struct tm *p;
	char buffer [80];
	
	if (argc!=4)
	{
		fprintf(stderr,"Usage: %s time_t lon lat\n", argv[0]);
		exit(1);
	}
	tc=atol(argv[1]);
	lon=atof(argv[2]);
	lat=atof(argv[3]);
	p=gmjtime_r(&tc, &ut);
	strftime (buffer,80,"UTC: %Y-%m-%d %H:%M:%S",&ut);
	puts(buffer);
	
    ut=TrueSolarTime(p, 0, 0, lon, lat);
	strftime (buffer,80,"LST: %Y-%m-%d %H:%M:%S",&ut);
	puts(buffer);
	exit(0);
}
