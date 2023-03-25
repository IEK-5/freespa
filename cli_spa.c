/*
    spa:
    * commandline interface to freespa and NREL spa (if compiled in)
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
#include <getopt.h>
#include "freespa.h"


#ifdef NRELSPA
/* This code requires you download the spa.c and spa.h codes from NREL
 * NREL does not allow free distribution of its SPA implementation
 * 
 * By enabling NREL spa you can compare the two implementations.
 * However, it only affects the calculation of the solar position and
 * currently does *NOT* affect sun rise/transit/set times
 */
#include "spa.h"


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
sol_pos NREL_SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
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
	
	P.az=fmod(M_PI*spa.zenith/180,2*M_PI);
	P.aa=fmod(M_PI*spa.azimuth/180, 2*M_PI);
	P.a=P.aa;
	P.z=P.az;
	P.z+=M_PI*spa.del_e/180;
	return P;
}
#endif

#define deg2rad(a) (M_PI*(a)/180.0)
#define rad2deg(a) (180.0*(a)/M_PI)
#define SUN_RADIUS 4.6542695162932789e-03 // in radians
int main(int argc, char **argv)
{
	// the default coordinate, IEK-5 building at the Forschungszentrum Jülich Campus
	double lon=deg2rad(6.41143),lat=deg2rad(50.90329);
	// elevation
	double E=96.0, Pr=1010, Temp=10;
	// use default freespa
	int fspa=1;	
	// per default we use the current time
	time_t tc=time(NULL);
	int v=0;
	int r;
	
	struct tm sunrise={0}, sunset={0}, transit={0};
	struct tm ut={0};
	struct tm *p;
	char buffer [80];
	sol_pos P;
	char c;
	double Esr, Ess, Etr;
	
	while (1)
	{
		static struct option long_options[] =
		{
			{"coordinate",  required_argument, 0, 'c'},
			{"elevation",   required_argument, 0, 'e'},
			{"pressure",    required_argument, 0, 'p'},
			{"temperature", required_argument, 0, 'T'},
			{"time",        required_argument, 0, 't'},
			{"verbose",           no_argument, 0, 'v'},
			{"nrel",              no_argument, 0, 'N'},
			{"help",              no_argument, 0, 'h'},
			{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc, argv, "c:e:p:T:t:vNh",long_options, &option_index);
		if (c == -1)
			break;
			
		switch (c)
		{
			case 'c':
			{
				char *p;
				if (!optarg)
				{
					fprintf(stderr, "Error: --coordinate (-c) requires a coordinate as lat,lon in degrees\n");
					return 1;	
				}
				p=optarg;
				while (*p && (*p!=','))
					p++;
				if (*p!=',')
				{
					fprintf(stderr, "Error: --coordinate (-c) requires a coordinate as lat,lon in degrees\n");
					return 1;	
				}
				*p='\0';
				lat=deg2rad(atof(optarg));
				lon=deg2rad(atof(p+1));
				*p=',';
				break;
			}
			case 'e':
				if (!optarg)
				{
					fprintf(stderr, "Error: --elevation (-e) requires an altitude in m\n");
					return 1;	
				}
				E=atof(optarg);
				break;
			case 'p':
				if (!optarg)
				{
					fprintf(stderr, "Error: --pressure (-p) requires a pressure in mbar\n");
					return 1;	
				}
				Pr=atof(optarg);
				break;
			case 'T':
				if (!optarg)
				{
					fprintf(stderr, "Error: --temperature (-T) requires an temperature in °C\n");
					return 1;	
				}
				Temp=atof(optarg);
				break;
			case 't':
				if (!optarg)
				{
					fprintf(stderr, "Error: --time (-t) requires a unix time\n");
					return 1;	
				}
				tc=atol(optarg);
				break;
			case 'v':
				v++;
				break;
			case 'N':
#ifdef NRELSPA
				fspa=0;
#else
				fprintf(stderr, "Warning: --nrel (-N) NREL spa is not available!\n");
				fprintf(stderr, "Note: NREL spa must be available at compile time\n");
#endif
				break;
			case 'h':
			{
				int i=0;
				printf("USAGE: %s [options]\n", argv[0]);
				printf("options:\n");
				/*struct option {
					   const char *name;
					   int         has_arg;
					   int        *flag;
					   int         val;
				};*/
				while (long_options[i].val)
				{
					switch (long_options[i].val)
					{
						case 'v':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t\t  [-vv] (extra verbose)\n");
							printf("\t  Make output more verbose\n\n");
							break;
						case 'c':
							printf("\t--%s [-%c] <latitude>,<longitude>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Coordinate: 50.90329,6.41143\n\n");
							break;
						case 'e':
							printf("\t--%s [-%c] <elevation>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Elevation: 96 m\n\n");
							break;
						case 'p':
							printf("\t--%s [-%c] <pressure>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Pressure: 1010 mbar\n\n");
							break;
						case 'T':
							printf("\t--%s [-%c] <temperature>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Temperature: 10 °C\n\n");
							break;
						case 't':
							printf("\t--%s [-%c] <time>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Time: current unix time (%ld)\n\n", tc);
							break;
						case 'N':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Use NREL spa instead of freespa (solar position only)\n\n");
							break;
						case 'h':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Print this help\n\n");
							break;
						default:									
							if (long_options[i].has_arg)
								printf("\t--%s [-%c] <value>\n\n", long_options[i].name, (char)long_options[i].val);
							else
								printf("\t--%s [-%c]\n\n", long_options[i].name, (char)long_options[i].val);
							break;
					}
					i++;
				}
				printf("\n");				
				break;
			}
			case '?':
			default:
				exit(1);
		}
	}
	
	p=gmjtime_r(&tc, &ut);
	if (v>1)
	{
		printf("| Input ----------------------------\n");
		printf("Coordinate:\t  %8.4f,%-8.4f°\n",rad2deg(lat),rad2deg(lon));
		printf("Elevation: \t\t%10.2f m\n",E);
		printf("Pressure:\t     %10.2f mbar\n",Pr);
		printf("Temperature:\t       %10.2f °C\n",Temp);
		if (fspa==1)
			printf("Implementation:              freespa\n");
		else
			printf("Implementation:             NREL spa\n");
		printf("------------------------------------\n\n");
	}
	
	if (fspa==1)
		P=SPA(p, NULL, 0, lon,  lat, E, Pr, Temp);
#ifdef NRELSPA
	else
		P=NREL_SPA(p, NULL, 0, lon,  lat, E, Pr, Temp);
#endif
	
	if (v>0)
	{
		printf("| Time -----------------------------\n");
		strftime (buffer,80,"UTC:\t\t%Y-%m-%d %H:%M:%S",&ut);
		puts(buffer);
		ut=TrueSolarTime(p, 0, 0, lon, lat);
		strftime (buffer,80,"LST:\t\t%Y-%m-%d %H:%M:%S",&ut);
		puts(buffer);
		printf("------------------------------------\n\n");
	}
	
	printf("| Solar Position -------------------\n");
	printf("aparent zenith: \t%10f °\n", rad2deg(P.az));
	printf("aparent azimuth:\t%10f °\n", rad2deg(P.aa));
	printf("true    zenith: \t%10f °\n", rad2deg(P.z));
	printf("true    azimuth:\t%10f °\n", rad2deg(P.a));
	printf("------------------------------------\n\n");
	
	r=SunTimes(ut, NULL, 0, lon, lat, 0, Pr, Temp, &sunrise, &transit, &sunset);
	
	if (v>0)
	{
		printf("| Sunrise and Sunset----------------\n");
		if (r==0)
		{
			P=SPA(&sunrise, NULL, 0.0, lon,  lat, 0, Pr, Temp);
			Esr=P.az-M_PI/2-SUN_RADIUS;
			P=SPA(&sunset, NULL, 0, lon,  lat, 0, Pr, Temp);
			Ess=P.az-M_PI/2-SUN_RADIUS;
			P=SPA(&transit, NULL, 0, lon,  lat, 0, Pr, Temp);
			Etr=atan(sin(P.aa)*fabs(tan(P.az)));
			strftime (buffer,80,"Sun Rise:\t %Y-%m-%d %H:%M:%S",&sunrise);
			printf("%s\n\t\t   (error %6.4f °)\n", buffer, rad2deg(Esr));
			strftime (buffer,80,"Transit:\t %Y-%m-%d %H:%M:%S",&transit);
			printf("%s\n\t\t   (error %6.4f °)\n", buffer, rad2deg(Etr));
			strftime (buffer,80,"Sun Set:\t %Y-%m-%d %H:%M:%S",&sunset);
			printf("%s\n\t\t   (error %6.4f °)\n", buffer, rad2deg(Ess));
		}
		else
		{
			if (r<0)
				printf("polar night\n");
			else
				printf("midnight sun\n");
			
			P=SPA(&transit, NULL, 0, lon,  lat, 0, Pr, Temp);
			Etr=atan(sin(P.aa)*fabs(tan(P.az)));
			strftime (buffer,80,"Transit:\t %Y-%m-%d %H:%M:%S",&transit);
			printf("%s\n\t\t   (error %6.4f °)\n", buffer, rad2deg(Etr));
		}
		printf("------------------------------------\n\n");
	}
	exit(0);
}
