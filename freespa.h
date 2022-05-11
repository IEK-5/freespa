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
#ifndef	_TIME_H_
#include <time.h>
#endif
#ifndef _FREESPA_H_
#define _FREESPA_H_
// container struct for solar position data
typedef struct sol_pos {
	double z, a, az, aa; // zenith, azimuth, aparent zenith, aparent azimuth
	int E; // error flag
} sol_pos;
double get_delta_t(struct tm *ut);
// note the time struct must be UTC time, not local!
sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T);            
// compute true solar time    
struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, 
			double lon, double lat);
// sunrise, transit, and sunset 
/* sets the tm structs to sunrise, sunset and transit (solar noon) times
 * return value:
 * 	0:  All OK
 *  1:  midnight sun, only transit computed (sunrise and sunset equal 
 *      to transit)
 *  -1: polar night, only transit computed (sunrise and sunset 
 *      equal to transit)
 * 
 * BUGS: near the edge between all day days and all day nights the 
 * results may be inaccurate of wrong. Obviously it will be hard to get 
 * this right.
 */
int SunTimes(struct tm ut, double *delta_t, double delta_ut1, 
			double lon, double lat, double p, double T, 
			struct tm *sunrise, struct tm *transit, struct tm *sunset);


// Utilities:
/* julian unix time routines
 * For modern day it should be equivalent to the standard routines
 * in time.h (apart from the fact that mkgmtime is absent on many 
 * platforms). However, these routines of freespa have a 10-day gap 
 * between the Julian and Gregorian calendar where the Julian calendar 
 * ends on October 4, 1582 (JD = 2299160), and the next day the 
 * Gregorian calendar starts on October 15, 1582.
 * 
 * This definition of unix time makes it compatible with the julian day 
 * as it is computed from a date in freespa, i.e. the julian days all 
 * have 86400 seconds. 
 */
// populate a time struct with UTC from unix time
struct tm *gmjtime_r(time_t *t, struct tm *ut);
struct tm *gmjtime(time_t *t);
// create unix time from time struct
time_t mkgmjtime(struct tm *ut);

int testjulian();
int testheliocentricpos();
#endif /* #ifndef _FREESPA_H_ */
