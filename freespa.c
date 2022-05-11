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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "freespa.h"
#include "freespa_tables.h"
#include "freespa_dt_table.h"

/* Implementation of the solar position algorithm
 * A document decribing the algorithm can be found on NREL's website:
 * https://midcdmz.nrel.gov/spa/
 * Also another implementation of the same algorithm can be found there.
 * Note, however, NREL's implementation is free as in beer, but NREL 
 * does not allow you to share the beer with your friends. 
 * 
 * This implementation lets you share the beer with your friends (under 
 * the conditions of the GPLv3)!
 * 
 * Most of the algorithm is described in more detail in:
 * Meeus, J., 1998. Astronomical Algorithms, second ed. Willmann-Bell, 
 * Inc., Richmond, Virginia, USA.
 */

/* defines */
#define JD0 2451545.0 // Julian day noon 1 January 2000 Universal Time
#define ETJD0 946728000 // unix time for julian day JD0
#define SUN_RADIUS 4.6542695162932789e-03 // in radians

#define EARTH_R 6378136.6 //m at the equator [IERS Numerical Standards (IAG 1999)] 
#define ABSOLUTEZERO -273.15 //convert C to K
#define AP0 1010.0 // standard sea level air pressure
#define AT0 10.0 // standard sea level air temperature


#define deg2rad(a) (M_PI*(a)/180.0)
#define rad2deg(a) (180.0*(a)/M_PI)

// some error codes
#define DELTA_T_OOR		1
#define DELTA_UT1_OOR	2
#define LON_OOR			4
#define LAT_OOR			8
#define ELEVATION_OOR	16
#define PRESS_OOR		32
#define TEMP_OOR		64
#define GMTIMEFAIL		128


/* structures for internal use */
typedef struct JulianDate {
	double JD, JDE, JC, JCE, JME;
	int E; // signal error state
} JulianDay;
typedef struct GeoCentricSolPos {
	double lat, lon, rad;
} GeoCentricSolPos;


/* Julian Day 
 * 
 * see
 * Meeus, J., 1998. Astronomical Algorithms, second ed. Willmann-Bell, 
 * Inc., Richmond, Virginia, USA.
 * Pages 59-66
 */
/* extrapolate delta t  (Morrison & Stephenson, 2004) */
#define DELATTEXTRAP(y) (32.0*(((y)-1820.0)/100)*(((y)-1820.0)/100)-20)
/* interpolate delta t */
double get_delta_t(struct tm *ut)
{
	double dyear;
	int imin=0, imax=NDT-1, i;
	
	dyear=(double)ut->tm_year+1900.0+((double)ut->tm_mon+1.0)/12+(double)(ut->tm_mday-1.0)/365.0;
	if (freespa_delta_t_table[0]>dyear)
		return DELATTEXTRAP(dyear);
	if (freespa_delta_t_table[2*imax]<dyear)
		return DELATTEXTRAP(dyear);
	
	while (imax-imin>1)
	{
		i=(imin+imax)/2;
		if (freespa_delta_t_table[2*i]>dyear)
			imax=i;
		else if (freespa_delta_t_table[2*i]<dyear)
			imin=i;
		else
			return freespa_delta_t_table[2*i+1];
	}
	return freespa_delta_t_table[2*imin+1]+(dyear-freespa_delta_t_table[2*imin])*(freespa_delta_t_table[2*imax+1]-freespa_delta_t_table[2*imin+1])/(freespa_delta_t_table[2*imax]-freespa_delta_t_table[2*imin]);
}
 
JulianDay MakeJulianDay(struct tm *ut, double *delta_t, double delta_ut1)
{
	int month, year;
	double day, a, dt;
	JulianDay JD;
	
	JD.E=0;
	day = (double)ut->tm_mday + ((double)ut->tm_hour+((double)ut->tm_min+((double)ut->tm_sec+delta_ut1)/60.0)/60.0)/24;
	month=ut->tm_mon+1;
	year=ut->tm_year+1900;
    if (month < 3)
    {
        month += 12;
        year--;
    }
	JD.JD=trunc(365.25*((double)year+4716.0))+trunc(30.6001*((double)month+1.0))+ day - 1524.5;
    if (JD.JD > 2299160.0)
    {
        a = trunc((double)year/100.0);
        JD.JD += (2 - a + trunc(a/4));
    }
	if (delta_t)
		dt=*delta_t;
	else
		dt=get_delta_t(ut);
	// Julian Ephemeris Day
	JD.JDE=JD.JD+dt/86400.0;
	JD.JC=(JD.JD-JD0)/36525.0;
	JD.JCE=(JD.JDE-JD0)/36525.0;
	JD.JME=JD.JCE/10.0;
	return JD;
}

JulianDay AddDays(JulianDay JD, int Ndays)
{
	double d=(double)Ndays;
	JD.JD+=d;
	JD.JDE+=d; // assume delta t does not change!
	JD.JC+=d/36525.0;
	JD.JCE+=d/36525.0;
	JD.JME=JD.JCE/10;
	return JD;
}


// Julian Day to tm struct 
struct tm *JDgmtime(JulianDay JD, struct tm *ut)
{
	double A,B,C,D,F,G,I,Z;
	double d;
	Z=trunc(JD.JD+0.5);
	F=JD.JD-Z;
	if (Z<2299161)
		A=Z;
	else
	{
		B=trunc((Z-1867216.25)/36524.25);
		A=Z+1+B-trunc(B/4.0);
	}
	C=A+1524;
	D=trunc((C-122.1)/365.25);
	G=trunc(365.25*D);
	I=trunc((C-G)/30.6001);
	d=C-G-trunc(30.6001*I)+F-0.5; // day starts at 00:00 not at 12:00
	
	// day 
	ut->tm_mday=(int)trunc(d)+1;
	// month since jan (0..11)
	if (I<14)
		ut->tm_mon=(int)(I-2);
	else
		ut->tm_mon=(int)(I-14);
	// year since 1900	
	if (ut->tm_mon>1)
		ut->tm_year=D-4716-1900;
	else
		ut->tm_year=D-4715-1900;
	d-=trunc(d);
	// d in days
	d*=86400;
	d=round(d);
	ut->tm_sec=((int)d)%60;
	d-=(double)ut->tm_sec;
	d/=60;
	ut->tm_min=((int)d)%60;
	d-=(double)ut->tm_min;
	d/=60;
	ut->tm_hour=((int)d)%60;
	d-=(double)ut->tm_hour;
	return ut;	
}

/* julian unix time routines
 * For modern day it should be equivalent to the standard routines
 * in time.h (apart from the fact that mkgmtime is absent on many 
 * platforms). However these routines have a 10-day gap between the 
 * Julian and Gregorian calendar where the Julian calendar ends on 
 * October 4, 1582 (JD = 2299160), and the next day the Gregorian 
 * calendar starts on October 15, 1582.
 * 
 * This definition of unix time makes it compatible with the julian day 
 * as it is computed from a date in freespa, i.e. the julian days all 
 * have 86400 seconds. 
 */
struct tm *gmjtime_r(time_t *t, struct tm *ut)
{
	JulianDay J;
	J.JD=((double)((*t)-ETJD0)/86400.0)+JD0;
	JDgmtime(J, ut);
	return ut;
}

struct tm *gmjtime(time_t *t)
{
	static struct tm _tmbuf;
	return gmjtime_r(t, &_tmbuf);
}
// inverse of above
time_t mkgmjtime(struct tm *ut)
{
	JulianDay J;
	J=MakeJulianDay(ut, 0, 0);
	return (time_t)round((J.JD-JD0)*86400)+ETJD0;
}

JulianDay MakeJulianDayEpoch(time_t t, double *delta_t, double delta_ut1)
{
	struct tm ut;
	struct tm *p;
	JulianDay JD;
	
	p=gmjtime_r(&t, &ut);
	
	if (!p)
	{
		JD.JD=0.0;
		JD.JDE=0.0;
		JD.JC=0.0;
		JD.JCE=0.0;
		JD.JME=0.0;
		JD.E=GMTIMEFAIL;
	}
	else
		JD=MakeJulianDay(&ut, delta_t, delta_ut1);
	return JD;
}

/* Heliocentric Earth Coordinate ***************************/
// corresponding tables with periodic terms are in freespa_tables.h
double SummPTerms(const p_term p[], int N, JulianDay JD)
{
    int i=0;
    double s=0;
    for (i=0;i<N;i++)
		s+=p[i].A*cos(p[i].P+p[i].W*JD.JME);
    return s;
}

double Heliocentric_lon(JulianDay JD)
{
	double lon, pp;
	lon=SummPTerms(EarthLon0, N_LON0, JD);
	pp=JD.JME;
	lon+=SummPTerms(EarthLon1, N_LON1, JD)*pp;
	pp*=JD.JME;
	lon+=SummPTerms(EarthLon2, N_LON2, JD)*pp;
	pp*=JD.JME;
	lon+=SummPTerms(EarthLon3, N_LON3, JD)*pp;
	pp*=JD.JME;
	lon+=SummPTerms(EarthLon4, N_LON4, JD)*pp;
	pp*=JD.JME;
	lon+=SummPTerms(EarthLon5, N_LON5, JD)*pp;
	lon/=1.0e8;
	return lon;
}

double Heliocentric_lat(JulianDay JD)
{
	double lat;
	lat=SummPTerms(EarthLat0, N_LAT0, JD);
	lat+=SummPTerms(EarthLat1, N_LAT1, JD)*JD.JME;
	lat/=1.0e8;
	return lat;
}

double Heliocentric_rad(JulianDay JD)
{
	double rad, pp;
	rad=SummPTerms(EarthRad0, N_RAD0, JD);
	pp=JD.JME;
	rad+=SummPTerms(EarthRad1, N_RAD1, JD)*pp;
	pp*=JD.JME;
	rad+=SummPTerms(EarthRad2, N_RAD2, JD)*pp;
	pp*=JD.JME;
	rad+=SummPTerms(EarthRad3, N_RAD3, JD)*pp;
	pp*=JD.JME;
	rad+=SummPTerms(EarthRad4, N_RAD4, JD)*pp;
	rad/=1.0e8;
	return rad;
}

/* Geocentric Sun Coordinate ***************************/
GeoCentricSolPos Geocentric_pos(JulianDay JD)
{
	GeoCentricSolPos P;
	P.lat=fmod(-Heliocentric_lat(JD), 2*M_PI);
	P.lon=fmod(Heliocentric_lon(JD)+M_PI, 2*M_PI);
	if (P.lon<0)
		P.lon+=2*M_PI;
	P.rad=Heliocentric_rad(JD);
	return P;
}


/* nutation and the obliquity of the ecliptic 
 * 
 * see 
 * Meeus, J., 1998. Astronomical Algorithms, second ed. Willmann-Bell, 
 * Inc., Richmond, Virginia, USA.
 * Pages 143-148
 */

/* utillity function to evaluate polynomials of arbitrary order
 * computes: y=a[0]*x^(N-1)+a[1]*x^(N-2)+...+a[N-1]
 * where a is an array of coefficients, N the length of the array and x
 * the value the polynomial is evaluated for.
 */
static inline double poly(const double a[], int N, double x)
{
	int i;
	double r;
	r=a[0];
	for (i=1;i<N;i++)
		r=a[i]+x*r;
	return r;
}	
/* here we define several 3rd order polynomials 
 * to compute nutation and the obliquity of the ecliptic 
 */
const double MEAN_ELONGATION_MOON_SUN[] = 	{deg2rad(1.0/189474.0) , deg2rad(-1.9142e-03), deg2rad(445267.11148) , deg2rad(297.85036)};
const double MEAN_ANOMALY_SUN[] 		= 	{deg2rad(-1.0/300000.0), deg2rad(-0.0001603) , deg2rad(35999.05034)  , deg2rad(357.52772)};
const double MEAN_ANOMALY_MOON[]		= 	{deg2rad(1.0/56250.0)  , deg2rad(0.0086972)  , deg2rad(477198.867398), deg2rad(134.96298)};
const double ARG_LAT_MOON[] 			=	{deg2rad(1.0/327270.0) , deg2rad(-0.0036825) , deg2rad(483202.017538), deg2rad(93.27191)};
const double ASC_LON_MOON[] 			=	{deg2rad(1.0/450000.0) , deg2rad(0.0020708)  , deg2rad(-1934.136261) , deg2rad(125.04452)};

void Nutation_lon_obliquity(JulianDay JD, double *del_psi, double *del_eps)
{
    int i, j;
    double sum, sum_psi=0, sum_eps=0;
    double x[5];
    x[0]=poly(MEAN_ELONGATION_MOON_SUN, 4, JD.JCE);
	x[1]=poly(MEAN_ANOMALY_SUN, 4,JD.JCE);
	x[2]=poly(MEAN_ANOMALY_MOON, 4,JD.JCE);
	x[3]=poly(ARG_LAT_MOON, 4,JD.JCE);
	x[4]=poly(ASC_LON_MOON, 4,JD.JCE);
    for (i = 0; i < NY; i++)
    {
        sum=0;
		for (j = 0; j < 5; j++)
			sum += x[j]*Y_Terms[i][j];
			
        sum_psi += (PE_Terms[i][0] + JD.JCE*PE_Terms[i][1])*sin(sum);
        sum_eps += (PE_Terms[i][2] + JD.JCE*PE_Terms[i][3])*cos(sum);
    }

    *del_psi = sum_psi * deg2rad(1/36000000.0);
    *del_eps = sum_eps * deg2rad(1/36000000.0);
    // actual obliquity computed in the main solpos routine
}


/* Atmospheric Refraction
 * 
 * see 
 * Meeus, J., 1998. Astronomical Algorithms, second ed. Willmann-Bell, 
 * Inc., Richmond, Virginia, USA.
 * Pages 105-108
 */

/* base form o Bennet's formula *and* that of its approximate reverse
 * i.e. we use the same functional form to compute from the true to 
 * aparent solar position, as we use to compute back.
 */
static inline double Refr(const double coeff[], double p, double T, double h)
{
	//converts true and aparent solar elevation
	return (p/AP0)*((AT0-ABSOLUTEZERO)/(T-ABSOLUTEZERO))*coeff[0]/tan(h+coeff[1]/(h+coeff[2]));
}

/* solar refraction according to Bennet
 * Used to convert the aparent solar position into the true solar 
 * position
 */
static inline double Bennet(double p, double T, double h)
{
	const double BENNET[] = {2.9088820866572158e-04,2.2267533386408395e-03,7.6794487087750510e-02};
	return Refr(BENNET, p, T, h);
}
/* solar refraction according to Schaefer
 * Used to convert the true solar position into the aparent solar 
 * position, i.e. the approximate inverse of Bennet.
 * 
 * SCHAEFER, B. E., & LILLER, W. (1990). REFRACTION NEAR THE HORIZON. 
 * Publications of the Astronomical Society of the Pacific, 102(653), 
 * 796–805. http://www.jstor.org/stable/40679565
 */
static inline double iBennet(double p, double T, double h)
{
	const double IBENNET[] = {2.9670597283903603e-04,3.1375594238030984e-03,8.9186324776910242e-02};
	return Refr(IBENNET, p, T, h);
}


/* solpos - internal routine to compute the solar position
 * input: 
 *  - longitude (lon)
 *  - latitude (lat)
 *  - elevation (e)
 *  - pressure (p)
 *  - atmospheric refraction at sun- rise/set (a_refr)
 *    a_refr may also be NAN, in which case we compute it with Bennet's 
 *    formula
 *  - temperature (T)
 *  - Julian Date (JD)
 *  - Geocentric Solar Position GP
 * 
 * output:
 *  - sol_pos struct with the real and aparent solar position
 */
 
// some more polynomials
const double ECLIPTIC_MEAN_OBLIQUITY[] = {2.45,5.79,27.87,7.12,-39.05,-249.67,-51.38,1999.25,-1.55,-4680.93,84381.448};
const double GSTA[] = {deg2rad(-1/38710000.0),deg2rad(0.000387933),0.0,deg2rad(280.46061837)};
sol_pos solpos(sol_pos P, double lon, double lat, double e, double p, double T, JulianDay JD, GeoCentricSolPos GP)
{
	double dtau, v, H, u, x, y;
	double lambda, alpha, delta, xi;
	double delta_prime, H_prime;
	//double aplpha_prime;
	double dpsi, deps, eps, dalpha;
	double h, dh=0, a_refr;
	
	// aberation correction
	dtau=deg2rad(-20.4898/3600.0)/GP.rad;
	
	// nutation and the obliquity of the ecliptic
	Nutation_lon_obliquity(JD, &dpsi, &deps);
	eps=deps+poly(ECLIPTIC_MEAN_OBLIQUITY,11,JD.JME/10)*deg2rad(1/3600.0);
	
	// aparent sun longitude
	lambda=GP.lon+dpsi+dtau;
	
	// sidereal time at Greenwich 
	v=poly(GSTA,4,JD.JC)+deg2rad(360.98564736629)*(JD.JD - JD0);
	v+=dpsi*cos(eps);
	
	// sun right ascension
	alpha=atan2(sin(lambda)*cos(eps)-tan(GP.lat)*sin(eps),cos(lambda));
	if (alpha<0)
		alpha+=2*M_PI;
		
	// sun declination
	delta=asin(sin(GP.lat)*cos(eps)+cos(GP.lat)*sin(eps)*sin(lambda));
	
	// hour angle
	H=v+lon-alpha; 
		
	// equatorial horizontal parallax of the sun
	xi=deg2rad(8.794/3600.0)/GP.rad;
	// term u
	u=atan(0.99664719*tan(lat));
	// x & y
	
	x=cos(u)+e*cos(lat)/EARTH_R;
	y=0.99664719*sin(u)+e*sin(lat)/EARTH_R;
	
	// parallax in the sun right ascension
	dalpha=atan2(-x*sin(xi)*sin(H),cos(delta)-x*sin(xi)*cos(H));
	
	// topocentric sun right ascension
	//alpha_prime=alpha+dalpha;
	
	// topocentric sun declination
	delta_prime=atan2((sin(delta)-y*sin(xi))*cos(dalpha),cos(delta)-x*sin(xi)*cos(H));
	
	// topocentric local hour angle
	H_prime=H-dalpha;
	
	// topocentric elevation angle without atmospheric refraction correction
	h=asin(sin(lat)*sin(delta_prime)+cos(lat)*cos(delta_prime)*cos(H_prime));
	
	a_refr=Bennet(p,T,0.0);
	if (h>=-a_refr-SUN_RADIUS)
		dh=iBennet(p,T,h);
	
	// compute sun zenith
	P.z=M_PI/2-h;
	
	// compute aparent sun zenith
	P.az=P.z-dh;
	P.a=fmod(M_PI+atan2(sin(H_prime),cos(H_prime)*sin(lat)-tan(delta_prime)*cos(lat)),2*M_PI);
	P.aa=P.a;
	
	// limit angular range
	P.z=fmod(P.z,2*M_PI);
	P.az=fmod(P.az,2*M_PI);
	if (P.z<0)
	{
		P.z=-P.z;
		P.a=fmod(P.a+M_PI,2*M_PI);
	}
	if (P.az<0)
	{
		P.az=-P.az;
		P.aa=fmod(P.a+M_PI,2*M_PI);
	}
	if (P.z>M_PI)
	{
		P.z=2*M_PI-P.z;
		P.a=fmod(P.a+2*M_PI,2*M_PI);
	}
	if (P.az>M_PI)
	{
		P.az=2*M_PI-P.az;
		P.aa=fmod(P.a+M_PI,2*M_PI);
	}
	return P;
}
// Equation of Time
// sun mean longitude polynomial
const double SMLON[] = {deg2rad(-1/2000000.0),deg2rad(-1/15300.0),deg2rad(1/49931.0),deg2rad(0.03032028),deg2rad(360007.6982779),deg2rad(280.4664567)};
double EoT(double lat, JulianDay JD, GeoCentricSolPos GP)
{
	double M, E;
	double dtau;
	double lambda, alpha, eps, deps, dpsi;
	
	
	// aberation correction
	dtau=deg2rad(-20.4898/3600.0)/GP.rad;
	
	// nutation and the obliquity of the ecliptic
	Nutation_lon_obliquity(JD, &dpsi, &deps);
	eps=deps+poly(ECLIPTIC_MEAN_OBLIQUITY,11,JD.JME/10)*deg2rad(1/3600.0);
	
	// aparent sun longitude
	lambda=GP.lon+dpsi+dtau;
	// sun right ascension
	alpha=atan2(sin(lambda)*cos(eps)-tan(GP.lat)*sin(eps),cos(lambda));
	M=poly(SMLON,6,JD.JME);
	E=fmod(M-deg2rad(0.0057183)-alpha+dpsi*cos(eps),2*M_PI);
	if (E>deg2rad(5))
		E-=2*M_PI;
	if (E<-deg2rad(5))
		E+=2*M_PI;
	return E;
}

/* check input values */
int InputCheck(double delta_ut1, double lon, 
            double lat, double e, double p, double T)
{
	int E=0;
	
	if (delta_ut1<-1)
		E|=DELTA_UT1_OOR;
	else if (delta_ut1>1)
		E|=DELTA_UT1_OOR;
		
	if (lon<-M_PI)
		E|=LON_OOR;
	else if (lon>M_PI)
		E|=LON_OOR;
		
	if (lat<-M_PI/2)
		E|=LAT_OOR;
	else if (lat>M_PI/2)
		E|=LAT_OOR;
	
	if (e<-EARTH_R)
		E|=ELEVATION_OOR;
	
	if (p<0)
		E|=PRESS_OOR;
	if (p>5000)
		E|=PRESS_OOR;
		
	if (T<ABSOLUTEZERO)
		E|=TEMP_OOR;
	return E;
}

/* SPA - exported routine to compute the solar position
 * input: 
 *  - ut			time struct with UTC
 *  - delta_t		pointer to Δt value, the difference between the 
 *                  Earth rotation time and the Terrestrial Time. If 
 * 					the pointer is NULL, use internal tables to find Δt.
 *  - delta_ut1 	is a fraction of a second that is added to the UTC 
 * 					to adjust for the irregular Earth rotation rate.
 *  - lon			observer longitude (radians)
 *  - lat 			observer latitude (radians)
 *  - e				observer elevation (m)
 *  - p				atmospheric pressure (mb)
 *  - T 			Temperature (C)
 * 
 * output:
 *  - sol_pos struct with the real and aparent solar position
 */
sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T)
{
	JulianDay D;
	GeoCentricSolPos G;
	sol_pos P;
	P.E=InputCheck(delta_ut1, lon, lat, e, p, T);
	if (!P.E)
	{
		D=MakeJulianDay(ut, delta_t, delta_ut1);
		if (D.E)
		{
			P.E=D.E;
			return P;
		}
		G=Geocentric_pos(D);
		P=solpos(P,lon,lat,e,p,T,D,G);
	}
	return P;
}

/* TrueSolarTime - exported routine to convert UTC to the local solar time
 * input: 
 *  - ut			time struct with UTC
 *  - delta_t		pointer to Δt value, the difference between the 
 *                  Earth rotation time and the Terrestrial Time. If 
 * 					the pointer is NULL, use internal tables to find Δt.
 *  - delta_ut1 	is a fraction of a second that is added to the UTC 
 * 					to adjust for the irregular Earth rotation rate.
 *  - lon			observer longitude (radians)
 *  - lat 			observer latitude (radians)
 * 
 * output:
 *  - tm struct with local solar time
 */
struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, double lon, double lat)
{
	double E;
	int sec,min,hour,day,month,year;
	struct tm nt={0};
	struct tm *p;	
	JulianDay D;
	GeoCentricSolPos G;
	struct tm st;
	int Err;
	if (InputCheck(delta_ut1, lon, lat, 0, 1, 10))
		return nt;
	D=MakeJulianDay(ut, delta_t, delta_ut1);
	G=Geocentric_pos(D);
	E=EoT(lat,D,G);	
	D.JD+=(lon+E)/M_PI/2; 
	JDgmtime(D, &nt);
	return nt;
}

/* SunTimes - exported routine to compute sunrise, transit and sunset 
 * times.
 * input: 
 *  - ut			time struct with UTC. Only used for the date (year, 
 * 					month, day)
 *  - delta_t		pointer to Δt value, the difference between the 
 *                  Earth rotation time and the Terrestrial Time. If 
 * 					the pointer is NULL, use internal tables to find Δt.
 *  - delta_ut1 	is a fraction of a second that is added to the UTC 
 * 					to adjust for the irregular Earth rotation rate.
 *  - lon			observer longitude (radians)
 *  - lat 			observer latitude (radians)
 *  - p				atmospheric pressure (mb)
 *  - T 			Temperature (C)
 * 
 * output:
 *  - tm structs sunrise, transit, and sunset
 * 
 * return value:
 * 	-1				Polar night. Only transit is computed.
 *  0				Regular day/night, all times are computed.
 *  1				Midnight sun. Only transit is computed.
 */
int SunTimes(struct tm ut, double *delta_t, double delta_ut1, double lon, double lat, double p, double T, struct tm *sunrise, struct tm *transit, struct tm *sunset)
{
	double dtau, vv, v[3], H,Hp[3], u, x, y;
	double lambda, alpha[3], delta[3], xi;
	double delta_prime, H_prime;
	double dpsi, deps, eps, dalpha;
	double h[3], dh=0, a_refr, arg;
	double m[3], n, dt, alphap[3], deltap[3], a,b,c,ap,bp,cp;
	JulianDay Day[3];
	GeoCentricSolPos GP;
	int i;
	
	
	ut.tm_hour=0;
	ut.tm_min=0;
	ut.tm_sec=0;
	
	if (delta_t)
		dt=*delta_t;
	else
		dt=get_delta_t(&ut);
		
	Day[1]=MakeJulianDay(&ut, delta_t, delta_ut1);
	Day[0]=AddDays(Day[1], -1);
	Day[2]=AddDays(Day[1], 1);
	
	for (i=0;i<3;i++)
	{
		GP=Geocentric_pos(Day[i]);
		// aberation correction
		dtau=deg2rad(-20.4898/3600.0)/GP.rad;
		
		// nutation and the obliquity of the ecliptic
		Nutation_lon_obliquity(Day[i], &dpsi, &deps);
		eps=deps+poly(ECLIPTIC_MEAN_OBLIQUITY,11,Day[i].JME/10)*deg2rad(1/3600.0);
		
		// aparent sun longitude
		lambda=GP.lon+dpsi+dtau;
		
		// sun right ascension
		alpha[i]=atan2(sin(lambda)*cos(eps)-tan(GP.lat)*sin(eps),cos(lambda));
		if (alpha[i]<0)
			alpha[i]+=2*M_PI;
			
		// sun declination
		delta[i]=asin(sin(GP.lat)*cos(eps)+cos(GP.lat)*sin(eps)*sin(lambda));
		if (i==1)
		{
			// sidereal time at Greenwich 
			vv=poly(GSTA,4,Day[i].JC)+deg2rad(360.98564736629)*(Day[i].JD - JD0);
			vv+=dpsi*cos(eps);
			vv=fmod(vv, 2*M_PI);
		}
	}
	// approximate solar transit time in radians
	m[0]=(alpha[1]-lon-vv);
	a_refr=-Bennet(p,T,0.0)-SUN_RADIUS;
	arg=(sin(a_refr)-sin(lat)*sin(delta[1]))/(cos(lat)*cos(delta[1]));
	// test if we have a regular sunrise/sunset
	if ((arg>-1) && (arg<1))
	{
		// approximate sun rise and set timnes in radians
		H=fmod(acos((sin(a_refr)-sin(lat)*sin(delta[1]))/(cos(lat)*cos(delta[1]))), M_PI);
		m[1]=m[0]-H;
		m[2]=m[0]+H;
		
		// refine the computed times
		a=alpha[1]-alpha[0];
		b=alpha[2]-alpha[1];
		ap=delta[1]-delta[0];
		bp=delta[2]-delta[1];
		if (abs(a)>deg2rad(2));
			fmod(a,deg2rad(1));
		if (abs(b)>deg2rad(2));
			fmod(b,deg2rad(1));
		if (abs(ap)>deg2rad(2));
			fmod(ap,deg2rad(1));
		if (abs(bp)>deg2rad(2));
			fmod(bp,deg2rad(1));
		c=b-a;
		cp=bp-ap;
		
		for (i=0;i<3;i++)
		{
			m[i]=fmod(m[i],2*M_PI)/(2*M_PI);
			if (m[i]<0)
				m[i]+=1.0;		
			v[i]=vv+deg2rad(360.985647)*m[i];
			n=m[i]+dt/86400;
			alphap[i]=alpha[1]+n*(a+b+c*n)/2;
			deltap[i]=delta[1]+n*(ap+bp+cp*n)/2;
			Hp[i]=fmod(v[i]+lon-alphap[i], 2*M_PI);
			if (Hp[i]<=-M_PI)
				Hp[i]+=2*M_PI;
			if (Hp[i]>=M_PI)
				Hp[i]-=2*M_PI;			
			h[i]=asin(sin(lat)*sin(deltap[i])+cos(lat)*cos(deltap[i])*cos(Hp[i]));
		}
		Day[0].JD=Day[1].JD;
		Day[2].JD=Day[1].JD;
		Day[1].JD+=m[0]-Hp[0]/2/M_PI;
		Day[0].JD+=m[1]+(h[1]-a_refr)/2/M_PI/cos(deltap[1])/cos(lat)/sin(Hp[1]);
		Day[2].JD+=m[2]+(h[2]-a_refr)/2/M_PI/cos(deltap[2])/cos(lat)/sin(Hp[2]);
		
		// noon is computed for the given date
		// sunrise must be before that time and unset after
		if (Day[0].JD>Day[1].JD)
			Day[0].JD-=1.0;
		if (Day[1].JD>Day[2].JD)
			Day[2].JD+=1.0;
			
		// populate the tm structs
		JDgmtime(Day[0],sunrise);
		JDgmtime(Day[1],transit);
		JDgmtime(Day[2],sunset);
		return 0;
	}
	else
	{
		// no sunrise and sunset
		// refine the computed solar noon
		a=alpha[1]-alpha[0];
		b=alpha[2]-alpha[1];
		ap=delta[1]-delta[0];
		bp=delta[2]-delta[1];
		if (abs(a)>deg2rad(2));
			fmod(a,deg2rad(1));
		if (abs(b)>deg2rad(2));
			fmod(b,deg2rad(1));
		if (abs(ap)>deg2rad(2));
			fmod(ap,deg2rad(1));
		if (abs(bp)>deg2rad(2));
			fmod(bp,deg2rad(1));
		c=b-a;
		cp=bp-ap;
		m[0]=fmod(m[0],2*M_PI)/(2*M_PI);
		if (m[0]<0)
			m[0]+=1.0;		
		v[0]=vv+deg2rad(360.985647)*m[0];
		n=m[0]+dt/86400;
		alphap[0]=alpha[1]+n*(a+b+c*n)/2;
		deltap[0]=delta[1]+n*(ap+bp+cp*n)/2;
		Hp[0]=fmod(v[0]+lon-alphap[0], 2*M_PI);
		if (Hp[0]<=-M_PI)
			Hp[0]+=2*M_PI;
		if (Hp[0]>=M_PI)
			Hp[0]-=2*M_PI;
		h[0]=asin(sin(lat)*sin(deltap[0])+cos(lat)*cos(deltap[0])*cos(Hp[0]));			
		
		Day[1].JD+=m[0]-Hp[0]/2/M_PI;
		JDgmtime(Day[1],sunrise);
		JDgmtime(Day[1],transit);
		JDgmtime(Day[1],sunset);
		// polar night (-1) or midnight sun (+1)?
		return 2*(h[0]>0)-1;
	}
}


/* some unit tests */
// unit test julian date computation with a table of dates
int testjulian()
{
	int ERR=0;
	int i, N=4;
	time_t dates[] = {946728000, 915148800, 538704000, -11676096000};
	double JDs[] = {2451545.0, 2451179.5, 2446822.5, 2305447.5};
    char* timestr;
#define JDEPS (1e-3/(24*60*60)) // 1 ms accuracy 
	JulianDay D;
	struct tm ut;
	struct tm *p;
	timestr=malloc(50*sizeof(char));
	for (i=0;i<N;i++)
	{
		p=gmjtime_r(&dates[i], &ut);	
		strftime(timestr, 50, "%Y-%m-%d %T %Z",p);
		printf("testing %s", timestr);
		D=MakeJulianDayEpoch(dates[i], 0,0);
		printf("--> JD:%.1f\n", D.JD);
		
		if (fabs(D.JD-JDs[i])>JDEPS)
		{
			fprintf(stderr, "Error: julian day computation is off by more than %e\n", JDEPS);
			fprintf(stderr, "Failed for EPOCH %ld\n", dates[i]);
			fprintf(stderr, "got %.16e\n", D.JD);
			fprintf(stderr, "should be %.16e\n", JDs[i]);
			fprintf(stderr, "diff %.16e\n", D.JD-JDs[i]);
			ERR=1;
		}
	}
	free(timestr);
	return ERR;
}

// unit test heliocentric coordinate
int testheliocentricpos()
{
	int ERR=0;
	time_t date =1066419030;
	double JD = 2452930.312847222;
	double lat=-1.7649101008724539e-06, lon=4.1919774712578828e-01, rad=0.9965422974;// all angles in radians
	double delta_t=67.0;
	double delta_ut1=0;
	double Hlon,Hlat, HR;
	JulianDay D;
	D=MakeJulianDayEpoch(date, &delta_t, delta_ut1);

	if (fabs(D.JD-JD)>JDEPS)
	{
		fprintf(stderr, "Error: julian day computation is off by more than %e\n", JDEPS);
		fprintf(stderr, "Failed for EPOCH %ld\n", date);
		fprintf(stderr, "got %e\n", D.JD);
		fprintf(stderr, "should be %e\n", JD);
		ERR=1;
	}
	Hlon=fmod(Heliocentric_lon(D), M_PI);
	if (fabs(Hlon-lon)>1e-6)
	{
		fprintf(stderr, "Error: heliocentric longitude deviates by %e rad.\n", Hlon-lon);
		fprintf(stderr, "       got %e, should be %e\n", Hlon,lon);
		ERR=1;
	}
	Hlat=fmod(Heliocentric_lat(D), M_PI);
	if (fabs(Hlat-lat)>1e-6)
	{
		fprintf(stderr, "Error: heliocentric latitude deviates by %e rad.\n", Hlat-lat);
		fprintf(stderr, "       got %e, should be %e\n", Hlat,lat);
		ERR=1;
	}
	HR=Heliocentric_rad(D);
	if (fabs(HR-rad)>1e-6)
	{
		fprintf(stderr, "Error: heliocentric radius deviates by %e AU.\n", HR-rad);
		fprintf(stderr, "       got %e, should be %e\n", HR,rad);
		ERR=1;
	}
	
	return ERR;
}




