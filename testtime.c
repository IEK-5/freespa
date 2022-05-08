#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "freespa.h"


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
#define N 10000000
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
int main()
{
	printf("Test gmjtime:\n");
	if (!TestGM())
	{
		printf("Test mkgmjtime:\n");
		TestMK();
	}
		
	
	exit(0);
}
