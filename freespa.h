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


// my epoch routines
// There are various calendar definitions
// it is convenient to have the epoch relate nicely with the julian day
// routines in this code. For this reason we have our own epoch 
// definition. (works best on systems where time_t is unsigned 64 bit...)
struct tm *Jgmtime(time_t t, struct tm *ut); // populate a time struct unix time
time_t Jmkgmtime(struct tm *ut); // create unix time from time struct

int testjulian();
int testheliocentricpos();
#endif /* #ifndef _FREESPA_H_ */
