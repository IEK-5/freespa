/*
    freespa
    Copyright (C) 2023  B. E. Pieters,
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
/* see the DOC.md file for documentation */

/* freespa provides an implementation of the SPA (solar position 
 * algorithm). A document decribing the algorithm can be found on 
 * NREL's website:
 * https://midcdmz.nrel.gov/spa/
 * 
 * The purpose of this implementation is to overcome licensing issues 
 * with NREL's spa. NREL's spa cannot be redistributed (sort of "free 
 * as in beer you may not share with your friends"). This 
 * implementation you can share under the GPLv3. 
 * 
 * NREL's spa is built around the algorithms found in:
 * Meeus, J., 1998. Astronomical Algorithms, second ed. Willmann-Bell, 
 * Inc., Richmond, Virginia, USA.
 */

/* Requirement:
 * we require a 64bit signed integer time_t type, which limits the 
 * portability somewhat.
 * 
 * Note to self:
 * Pehaps I should use just a signed 64bit integer type instead of 
 * time_t? I suppose it would make the code more portable but 
 * interfacing with the standard time handling routines more awkward. 
 * 
 * Background:
 * The problem of time_t is that there is no good standard for it.
 * In ISO C `time_t` is specified as either an integer or 
 * floating-point type, and even the meaning to `time_t` is not 
 * specified (i.e. it may not be seconds since the epoch). On POSIX-
 * conformant systems `time_t` is required to ba an integer type, but 
 * it may still be 32bit or unsigned. The GNU C library uses signed 
 * integers and negative numbers are interpreted as dates before the 
 * epoch (as freespa does too). Also 64bit windows systems use a signed 
 * 64bit integer as the time_t type(1). Thus at least on those two 
 * 64bit platforms freespa should work. 
 * 
 * (1) Beware, on windows the gmtime routines do not work for dates
 * before Jan 01 00:00:00 1970, and for dates after 
 * Dec 31 23:59:59 3000. As far as I can tell these limits are 
 * completely arbitrary. The upper limit is probably not so relevant. 
 * However, I find it hard to understand why Microsoft engineers 
 * deliberately criple these routines so they cannot handle birthdays 
 * before 1970... What were they thinking, sod those old bastards?
 * If you have to handle a birthday before 1970 on a MS system you may 
 * use freespa's gmjtime routines instead.
 */
 
/* The following asserts check the time_t type at compile time */
/* Verify that time_t is an integer type.  */
_Static_assert ((time_t) 1.5 == 1, "error: time_t type is not integer");
/* Verify that time_t is 64 bit.  */
_Static_assert (sizeof(time_t) == 8, "error: time_t type is not 64 bit");
/* Verify that time_t is signed.  */
_Static_assert ((time_t) -1 < 0, "error: time_t type is not signed");

// some error codes for the error flag in sol_pos
// errors are be combined with a binary OR
#define _FREESPA_DEU_OOR		0X01	// ΔUT1 out of range
#define _FREESPA_LON_OOR		0X02	// longitude out of range
#define _FREESPA_LAT_OOR		0X04	// latitude out of range
#define _FREESPA_ELE_OOR		0X08	// elevation out of range
#define _FREESPA_PRE_OOR		0X10	// pressure out of range
#define _FREESPA_TEM_OOR		0X20	// temperature out of range
#define _FREESPA_DIP_OOR		0X40	// geometric dip out of range
#define _FREESPA_GMTIMEF		0X80	// time conversion error 

// container struct for solar position data
typedef struct sol_pos {
	double z, a; // zenith, azimuth
	int E; // error flag
} sol_pos;

typedef struct solar_day {
	struct tm ev[11];
	time_t t[11];
	double E[11];
	int status[11];
} solar_day;
			
// binary masks to enable/disable computing specific solar day events
#define _FREESPA_SUNRISE 0X01
#define _FREESPA_SUNSET  0X02
#define _FREESPA_CVDAWN  0X04
#define _FREESPA_CVDUSK  0X08
#define _FREESPA_NADAWN  0X10
#define _FREESPA_NADUSK  0X20
#define _FREESPA_ASDAWN  0X40
#define _FREESPA_ASDUSK  0X80

// binary mask variable to configure what solar events SolarDay computes
// default is all (0XFF)
extern int SDMASK;
// status flags
#define _FREESPA_EV_ERR       20
#define _FREESPA_EV_NA        10
#define _FREESPA_EV_OK         0
#define _FREESPA_EV_SUNABOVE   1
#define _FREESPA_EV_SUNBELOW  -1
// compute the real solar position
sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e);   
            
// correct for atmospheric refraction 
sol_pos ApSolposBennet(sol_pos P, double *gdip, double e, double p, double T);
sol_pos ApSolposBennetNA(sol_pos P, double *gdip, double e, double p, double T);

// compute true solar time
struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, 
						double lon, double lat);

// compute the solar events
extern int SDMASK;
solar_day SolarDay(struct tm *ut, double *delta_t, double delta_ut1, 
                   double lon, double lat, double e, double *gdip, 
                   double p, double T, 
                   sol_pos (*refract)(sol_pos,double*,double,double,double));

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
