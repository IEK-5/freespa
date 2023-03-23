/*
    lst
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
int main(int argc, char **argv)
{
	double lon,lat;
	time_t tc;
	struct tm ut={0};
	struct tm *p;
	char buffer [80];
	sol_pos P;
	
	if (argc!=4)
	{
		fprintf(stderr,"Usage: %s time_t lon lat (angles in derees)\n", argv[0]);
		exit(1);
	}
	tc=atol(argv[1]);
	lon=deg2rad(atof(argv[2]));
	lat=deg2rad(atof(argv[3]));
	
	p=gmjtime_r(&tc, &ut);
	P=SPA(p, NULL, 0, lon,  lat, 0, 1010, 10);	
	
	strftime (buffer,80,"UTC:\t\t\t%Y-%m-%d %H:%M:%S",&ut);
	puts(buffer);
    ut=TrueSolarTime(p, 0, 0, lon, lat);
	strftime (buffer,80,"LST:\t\t\t%Y-%m-%d %H:%M:%S",&ut);
	puts(buffer);
	printf("aparent zenith: \t%-10f rad\t%10f °\n", P.az, rad2deg(P.az));
	printf("aparent azimuth:\t%-10f rad\t%10f °\n", P.aa, rad2deg(P.aa));
	printf("true    zenith: \t%-10f rad\t%10f °\n", P.z, rad2deg(P.z));
	printf("true    azimuth:\t%-10f rad\t%10f °\n\n", P.a, rad2deg(P.a));
	exit(0);
}
