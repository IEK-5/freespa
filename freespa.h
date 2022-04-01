#ifndef	_TIME_H_
#include <time.h>
#endif
#ifndef _FREESPA_H_
#define _FREESPA_H_
// container struct for solar position data
typedef struct sol_pos {
	double z, a, az, aa; // zenith, azimuth, aparent zenith, aparent azimuth
	int E;
	//time_t sunrize, sunset, noon;
} sol_pos;

// note the time struct must be UTC time, not local!
sol_pos SPA(struct tm *ut, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T);            
            
int testjulian();
int testheliocentricpos();
#endif /* #ifndef _FREESPA_H_ */
