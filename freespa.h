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
	double z, a; // zenith, azimuth
	int E; // error flag
} sol_pos;

/* SPA routine:
 * returns a solpos struct
 * Input:
 *     ut:	        pointer to time struct with UTC time
 *     delta_t:     pointer to delta_t value, or NULL (use internal tables)
 *     delta_ut1:   delta_ut1 (offset in deconds smaller than 1)
 *     lon:			longitude (in radians)
 *     lat:			latitude (in radians)
 *     e:			observer elevation (in meter)
 */
sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e);   
 
/* Aparent Solar Position:
 * computes the aparent position of the sun using refreaction according 
 * to Bennet
 * returns a solpos struct
 * Input:
 *     P:	        real solar position
 *     gdip:		geometric dip, i.e. how far the horizon is below the 
 *                  observer (in rad). If this pointer is NULL the 
 *                  geometric dip is computed from the observer elevation
 *                  (assuming the horizon is at sea level)
 *     e:			observer elevation (in meter)
 *     p:			pressure (in mbar)
 *     T:			Temperature (in °C)
 */ 
 
sol_pos AparentSolpos(sol_pos P, double *gdip, double e, double p, double T);
/* TrueSolarTime routine:
 * returns a tm struct with the local true solar time 
 * 
 * Input:
 *     ut:	        pointer to time struct with UTC time
 *     delta_t:     pointer to delta_t value, or NULL (use internal tables)
 *     delta_ut1:   delta_ut1 (offset in deconds smaller than 1)
 *     lon:			longitude (in radians)
 *     lat:			latitude (in radians)
 * 
 */                   
struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, 
						double lon, double lat);


/* solar_day: acontainer struct for the solar events of a day. The 
 * events are:
 * 
 * Index	Event
 * 0:		previous low (lowest point of the sun before time t
 * 1:		transit (transit closest to time t
 * 2:		next low (lowest point of the sun after time t)
 * 3:		sunrise
 * 4:		sunset
 * 5:		civil dawn
 * 6:		civil dusk
 * 7:		nautical dawn
 * 8:		nautical dusk
 * 9:		astronomical dawn
 * 10:		astronomical dusk
 * 
 * the array t contsinas corresponding unix times
 * the status array contains:
 * 10		not available
 * 0		OK
 * 1		sun always above (for events with index 3 or higher)
 * -1		sun always below (for events with index 3 or higher)
 * 
 * Note1 the events on a day are the events of a solar-day and thus will 
 *       generally not fall within the same day. 
 * Note2 a solar day may be longer or shorter than 86400 seconds 
 */
typedef struct solar_day {
	struct tm ev[11];
	time_t t[11];
	double E[11];
	int status[11];
} solar_day;
			
// binary masks to enable/disable computing specific solar day events
#define SUNRISE 0X1
#define SUNSET  0X2
#define CVDAWN  0X4
#define CVDUSK  0X8
#define NADAWN  0X10
#define NADUSK  0X20
#define ASDAWN  0X40
#define ASDUSK  0X80
// binary mask variable to configure what solar events SolarDay computes
// default is all (0XFF)
extern int SDMASK;
// status flags
#define EV_NA        10
#define EV_OK         0
#define EV_SUNABOVE   1
#define EV_SUNBELOW  -1

/* SolarDay routine:
 * returns a solar_day struct with the events of a day 
 * 
 * Input:
 *     ut:	        pointer to time struct with UTC time
 *     delta_t:     pointer to delta_t value, or NULL (use internal tables)
 *     delta_ut1:   delta_ut1 (offset in deconds smaller than 1)
 *     lon:			longitude (in radians)
 *     lat:			latitude (in radians)
 *     e:			observer elevation (in meter)
 *     gdip:		geometric dip, i.e. how far the horizon is below the 
 *                  observer (in rad). If this pointer is NULL the 
 *                  geometric dip is computed from the observer elevation
 *                  (assuming the horizon is at sea level)
 *     p:			pressure (in mbar)
 *     T:			Temperature (in °C)
 * 
 */   
solar_day SolarDay(struct tm *ut, double *delta_t, double delta_ut1, 
                   double lon, double lat, double e, double *gdip, 
                   double p, double T);

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
// get delta_t value from internal tables
double get_delta_t(struct tm *ut);
// populate a time struct with UTC from unix time
struct tm *gmjtime_r(time_t *t, struct tm *ut);
struct tm *gmjtime(time_t *t);
// create unix time from time struct
time_t mkgmjtime(struct tm *ut);


int testjulian();
int testheliocentricpos();
#endif /* #ifndef _FREESPA_H_ */
