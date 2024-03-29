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
#if defined(_WIN32) || defined(_WIN64) 
#include <windows.h>
#include <locale.h>
#endif
#include "freespa.h"


#ifdef NRELSPA
/* This code requires you download the spa.c and spa.h codes from NREL
 * NREL does not allow free distribution of its SPA implementation
 * 
 * By enabling NREL spa you can compare the two implementations.
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
// simple wrapper around NREL's spa code, computes real solar position
sol_pos NREL_SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
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

// simple wrapper around NREL's spa code computes aparent solar position
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

#define SUN_RADIUS 4.6542695162932789e-03 // in radians
solar_day NRELSolarDay(struct tm *ut, double *delta_t, double delta_ut1, 
                   double lon, double lat, double e, double *gdip, 
                   double p, double T)
{
	solar_day D={0};
	sol_pos P;
	spa_data spa;
	FS_TIME_T t, t0, tt;	
	int i;
	
	if (gdip)
		fprintf(stderr, "Warning: NREL spa does not take a geometric dip into account\n");

	for (i=0;i<11;i++)
		D.status[i]=_FREESPA_EV_NA; 
	
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
	spa.function=SPA_ZA_RTS;
	spa_calculate(&spa);
	if ((spa.sunrise<0)||(spa.sunset<0)||(spa.suntransit<0))
		return D;
	
	ut->tm_hour=0;
	ut->tm_min=0;
	ut->tm_sec=0;
	t0=mkgmjtime(ut);
	
	tt=t0+(FS_TIME_T)(3600*spa.suntransit);
	gmjtime_r(&tt, D.ev+1);
	D.t[1]=tt;
	D.status[1]=_FREESPA_EV_OK;
	D.E[1]=NAN;
		
	
	
	t=t0+(FS_TIME_T)(3600*spa.sunrise);
	if (t>tt)
		t-=86400;
	gmjtime_r(&t, D.ev+3);
	D.t[3]=t;
	P=SPA_WrapperApp(D.ev+3, delta_t, delta_ut1, lon, lat,e, p, T);
	D.E[3]=P.z-M_PI/2-SUN_RADIUS;
	D.status[3]=_FREESPA_EV_OK;
	
	t=t0+(FS_TIME_T)(3600*spa.sunset);
	if (t<tt)
		t+=86400;
	gmjtime_r(&t, D.ev+4);
	D.t[4]=t;
	P=SPA_WrapperApp(D.ev+4, delta_t, delta_ut1, lon, lat,e, p, T);
	D.E[4]=P.z-M_PI/2-SUN_RADIUS;
	D.status[4]=_FREESPA_EV_OK;
	// TODO match freespa return value
	return D;
}
#endif

#ifndef DEFLAT
	#define DEFLAT 50.90329
#endif
#ifndef DEFLON
	#define DEFLON 6.41143
#endif
#ifndef DEFELE
	#define DEFELE 96.0
#endif

char *solevents[11]={
	"Prev.Sol.Midnight",
	"Transit          ",
	"Next Sol.Midnight",
	"Sunrise          ",
	"Sunset           ",
	"Civil Dawn       ",
	"Civil Dusk       ",
	"Nautical Dawn    ",
	"Nautical Dusk    ",
	"Astronomical Dawn",
	"Astronomical Dusk"
};
int chrono[11]={0,9,7,5,3,1,4,6,8,10,2};

void PrintUTC(struct tm *ut, char *buffer, int n)
{
	strftime (buffer,n,"%Y-%m-%d %H:%M:%S   UTC",ut);
}

void PrintLTZ(struct tm *ut, char *buffer, int n)
{
	FS_TIME_T tc=mkgmjtime(ut);
	struct tm lt={0};
#ifdef FS_CUSTOM_TIME_T
	time_t tt=(time_t)tc;
	// we might not be using time_t for a reason
	if ((FS_TIME_T)tt!=tc)
	{
		fprintf(stderr,"Error: failed to convert epoch time to a time_t type\n");
		strncpy(buffer, "unknown", n);
		return;
	}
	localtime_r(&tt, &lt);
#else
	localtime_r(&tc, &lt);
#endif

#if defined(_WIN64) || defined(_WIN32)
	strftime (buffer,n,"%Y-%m-%d %H:%M:%S %Z",&lt);
#else
	strftime (buffer,n,"%Y-%m-%d %H:%M:%S %5Z",&lt);
#endif
}

void PrintLST(struct tm *lst, char *buffer, int n)
{
	strftime (buffer,n,"%Y-%m-%d %H:%M:%S",lst);
}


#define deg2rad(a) (M_PI*(a)/180.0)
#define rad2deg(a) (180.0*(a)/M_PI)
#define SUN_RADIUS 4.6542695162932789e-03 // in radians
int main(int argc, char **argv)
{
	// the default coordinate, IEK-5 building at the Forschungszentrum Jülich Campus
	double lon=deg2rad(DEFLON),lat=deg2rad(DEFLAT);
	double dip, *gdip=NULL;
	// elevation
	double E=DEFELE, Pr=1010, Temp=10;
	// use default freespa
	int fspa=1, local=0;	
	// what to compute/show
	int tsoltime=0, solpos=1, suntimes=0;	
	// per default we use the current time
	FS_TIME_T tc=time(NULL);
	sol_pos (*aparent)(sol_pos,double*,double, double, double)=&ApSolposBennet;
#if defined(_WIN32) || defined(_WIN64) 
    UINT CODEPAGE_ORIGINAL = GetConsoleOutputCP();
    SetConsoleOutputCP(65001);
    setlocale(LC_ALL, "C");
#endif
	
	
	struct tm ut={0};
	struct tm lst={0};
	struct tm *p;
	char buffer [80];
	sol_pos P, Pa;
	int c;
	while (1)
	{
		static struct option long_options[] =
		{
			{"coordinate",  required_argument, 0, 'c'},
			{"elevation",   required_argument, 0, 'e'},
			{"dip",   		required_argument, 0, 'd'},
			{"pressure",    required_argument, 0, 'p'},
			{"temperature", required_argument, 0, 'T'},
			{"time",        required_argument, 0, 't'},
			{"local-time",        no_argument, 0, 'l'},
			{"soltime",           no_argument, 0, 's'},
			{"solpos",            no_argument, 0, 'S'},
			{"suntimes",          no_argument, 0, 'r'},
			{"nrel",              no_argument, 0, 'N'},
			{"bennet-NA",         no_argument, 0, 'A'},
			{"help",              no_argument, 0, 'h'},
			{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc, argv, "c:e:d:p:T:t:lsSrNAh",long_options, &option_index);
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
			case 'd':
				if (!optarg)
				{
					fprintf(stderr, "Error: --dip (-d) requires an geometric dip value in degrees\n");
					return 1;	
				}
				dip=deg2rad(atof(optarg));
				if (fabs(dip)>M_PI/2)
				{
					fprintf(stderr, "Error: invalid dip value, value may not exceed 90° (absolute)\n");
					return 1;	
				}
				gdip=&dip;
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
				tc=atoll(optarg);
				break;
			case 'l':
				local=1;
				break;
			case 's':
				tsoltime=(!tsoltime);
				break;
			case 'S':
				solpos=(!solpos);
				break;
			case 'r':
				suntimes=(!suntimes);
				break;
			case 'N':
#ifdef NRELSPA
				fspa=0;
#else
				fprintf(stderr, "Warning: --nrel (-N) NREL spa is not available!\n");
				fprintf(stderr, "Note: NREL spa must be available at compile time\n");
#endif
				break;
			case 'A':
				aparent=&ApSolposBennetNA;
				break;
			case 'h':
			{
				int i=0;
				printf("USAGE: %s [options]\n", argv[0]);
				printf("options:\n");
				while (long_options[i].val)
				{
					switch (long_options[i].val)
					{
						case 'c':
							printf("\t--%s [-%c] <latitude>,<longitude>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Coordinate: %f,%f\n\n", DEFLAT, DEFLON);
							break;
						case 'e':
							printf("\t--%s [-%c] <elevation>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Default Elevation: %g m\n\n",DEFELE);
							break;
						case 'd':
							printf("\t--%s [-%c] <dip>\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  specify a geometric dip (default none, compute from elevation)\n\n");
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
						case 'l':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  display local time\n\n");
							break;
						case 's':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  toggle printing the true solar time (default false)\n\n");
							break;
						case 'S':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  toggle printing the solar position (default true)\n\n");
							break;
						case 'r':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  toggle printing the sun-rise, set, and transit times (default false)\n\n");
							break;
						case 'N':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Use NREL spa instead of freespa\n\n");
							break;
						case 'A':
							printf("\t--%s [-%c]\n", long_options[i].name, (char)long_options[i].val);
							printf("\t  Use modified Bennet coefficients\n\n");
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
	if (tsoltime)
	{
		printf("| Time --------------------------------------\n");
		PrintUTC(&ut, buffer, 80);
		printf("UTC               : %s\n", buffer);
		PrintLTZ(&ut, buffer, 80);
		printf("LOC               : %s\n", buffer);	
		lst=TrueSolarTime(p, 0, 0, lon, lat);
		PrintLST(&lst, buffer, 80);
		printf("LST               : %s\n", buffer);
		printf("---------------------------------------------\n");
	}
	
	if (solpos)
	{
		if (fspa==1)
		{
			P=SPA(p, NULL, 0, lon,  lat, E);
			Pa=aparent(P, gdip, E, Pr, Temp);
		}
#ifdef NRELSPA
		else
		{
			P=NREL_SPA(p, NULL, 0, lon,  lat, E);
			Pa=aparent(P, gdip, E, Pr, Temp);
		}
#endif
		printf("| Solar Position ----------------------------\n");
		printf("apparent zenith   :            %13f°\n", rad2deg(Pa.z));
		printf("apparent azimuth  :            %13f°\n", rad2deg(Pa.a));
		printf("true     zenith   :            %13f°\n", rad2deg(P.z));
		printf("true     azimuth  :            %13f°\n", rad2deg(P.a));
		printf("---------------------------------------------\n");
	}
	if (suntimes)
	{
		solar_day D;
		int i;
		printf("| Solar Day Events --------------------------\n");
		if (fspa==1)
			D=SolarDay(&ut, NULL, 0, lon, lat, E, gdip, Pr, Temp, aparent);
#ifdef NRELSPA
		else
			D=NRELSolarDay(&ut, NULL, 0, lon, lat, E, NULL, Pr, Temp);
#endif
		for (i=0;i<11;i++)
		{
			// go through the day events
			if (D.status[chrono[i]]==0)
			{
				if (local)
					PrintLTZ(D.ev+chrono[i], buffer, 80);
				else
					PrintUTC(D.ev+chrono[i], buffer, 80);
				printf("%s : %s\n", solevents[chrono[i]],buffer);
					
				/*if (chrono[i]>2)
					printf("\t\t      (error %7.4f °)\n", rad2deg(D.E[chrono[i]]));*/
			}
			else if (D.status[chrono[i]]==1)
				printf("%s : -- sun above\n", solevents[chrono[i]]);
			else if (D.status[chrono[i]]==-1)
				printf("%s : -- sun below\n", solevents[chrono[i]]);
		}	

		printf("---------------------------------------------\n");
	}	
#if defined(_WIN32) || defined(_WIN64) 
    SetConsoleOutputCP(CODEPAGE_ORIGINAL);
    setlocale(LC_ALL, "");
#endif
	exit(0);
}
