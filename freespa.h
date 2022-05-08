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
struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, double lon, double lat);


/* julian unix time routines
 * For modern day it should be equivalent to the standard routines
 * in time.h (apart from the fact that mkgmtime is absent on many 
 * platforms). However, these routines have a 10-day gap between the 
 * Julian and Gregorian calendar where the Julian calendar ends on 
 * October 4, 1582 (JD = 2299160), and the next day the Gregorian 
 * calendar starts on October 15, 1582.
 * 
 * This definition of unix time makes it compatible with the julian day 
 * as it is computed from a date in freespa, i.e. the julian days all 
 * have 86400 seconds. 
 */
struct tm *gmjtime_r(time_t *t, struct tm *ut); // populate a time struct with UTC from unix time
time_t mkgmjtime(struct tm *ut); // create unix time from time struct

int testjulian();
int testheliocentricpos();
#endif /* #ifndef _FREESPA_H_ */
