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


int main(int argc, char **argv)
{
	double lon,lat;
	int r;
	time_t tc;
	struct tm ut={0};
	struct tm *p;
	char buffer [80];
	struct tm sunrise={0}, sunset={0}, transit={0};
	
	if (argc!=4)
	{
		fprintf(stderr,"Usage: %s time_t lon lat (angles in radians)\n", argv[0]);
		exit(1);
	}
	tc=atol(argv[1]);
	lon=atof(argv[2]);
	lat=atof(argv[3]);
	p=gmjtime_r(&tc, &ut);
	strftime (buffer,80,"UTC: %Y-%m-%d %H:%M:%S",&ut);
	puts(buffer);
	
	r=SunTimes(ut, NULL, 0, lon, lat, 1010, 10, &sunrise, &transit, &sunset);
    ut=TrueSolarTime(p, 0, 0, lon, lat);
	strftime (buffer,80,"LST: %Y-%m-%d %H:%M:%S",&ut);
	puts(buffer);
	if (r==0)
	{
		strftime (buffer,80,"Sun Rise: %Y-%m-%d %H:%M:%S",&sunrise);
		puts(buffer);
		strftime (buffer,80,"Transit : %Y-%m-%d %H:%M:%S",&transit);
		puts(buffer);
		strftime (buffer,80,"Sun Set : %Y-%m-%d %H:%M:%S",&sunset);
		puts(buffer);
	}
	else
	{
		if (r<0)
			printf("polar night\n");
		else
			printf("midnight sun\n");
		
		strftime (buffer,80,"Transit : %Y-%m-%d %H:%M:%S",&transit);
		puts(buffer);
	}	
	
	exit(0);
}
