#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "freespa.h"
#include "freespa_tables.h"

/* free implementation of NREL's SPA algorithm
 * 
 * See: https://midcdmz.nrel.gov/spa/
 * 
 * NREL also provides a code to compute the spa algorithm
 * 
 * Contrary to NREL's code this code may be distributed under the terms 
 * of the GPL v3. Note that in violation of NREL's copyright terms you 
 * may find several copies of NREL's code on github. This code may help 
 * remedy that... However, it is not a drop in replacement, i.e. this 
 * code uses its own data structures and convensions. 
 * 
 * All angles are strictly in radians (this includes coordinates!)
 * 
 * 
 */
// Julian day noon 1 January 2000 Universal Time
#define JD0 2451545.0
#define SUN_RADIUS 4.6542695162932789e-03 // in radians
#define NREL
#ifdef NREL
#define EARTH_R 6378140.0 //m at the equator
#define ABSOLUTEZERO -273.0 //convert C to K
/* the atmospheric refraction in NREL's SPA 
 * implementation follows Bennet's model:
 * Bennett, G. G. Journal of Navigation 1982, 35, 255–259. 
 * Unfortunately the reference document on the SPA model 
 * is somewhat sparce in its citations so you will not find this 
 * reference there. I do not know where NREL got the parameters 
 * for the model. Bennet's parameters were:
 * A=0.0167
 * B=7.31
 * C=4.4
 * Hohenkerk, C. In Explanatory Supplement to the Astronomical Almanac; 
 * University Science Books, 2013; pages 249–304.
 * 
 * Hohenkerk also updated the parameters for a better agreement with 
 * observations 
 * A=0.0167
 * B=7.32
 * C=4.32
 * Wilson, Teresa, "Evaluating the Effectiveness of Current 
 * Atmospheric Refraction Models in Predicting Sunrise and Sunset 
 * Times", Open Access Dissertation, Michigan Technological University, 
 * 2018. https://doi.org/10.37099/mtu.dc.etdr/697
 * 
 */
#define BENNET_A (1.02/60)
#define BENNET_B 10.3
#define BENNET_C 5.11
#else
#define EARTH_R 6378137.0 //m at the equator
#define ABSOLUTEZERO -273.15 
#define BENNET_A 0.0167
#define BENNET_B 7.32
#define BENNET_C 4.32
#endif
#define deg2rad(a) (M_PI*(a)/180.0)
#define rad2deg(a) (180.0*(a)/M_PI)


typedef struct JulianDate {
	double JD, JDE, JC, JCE, JME;
	int E; // signal error state
} JulianDate;
typedef struct GeoCentricSolPos {
	double lat, lon, rad;
} GeoCentricSolPos;

// some error codes for out of range values
#define DELTA_T_OOR		1
#define DELTA_UT1_OOR	2
#define LON_OOR			4
#define LAT_OOR			8
#define ELEVATION_OOR	16
#define PRESS_OOR		32
#define AREFR_OOR		64
#define TEMP_OOR		128
#define GMTIMEFAIL		256
/* Julian date ******************************************
 * Routine to compute the julian date from epoch time
 */
double INT(double x)
{
	return (double)((int)x);
} 

JulianDate MakeJulianDate(time_t t, double delta_t, double delta_ut1)
{
	struct tm *ut;
	int month, year;
	double day, a;
	JulianDate JD;
	
	ut=gmtime(&t);	
	// ERRORFLAG GMTIMENULL  "Error: gmtime returned NULL"
	if (!ut)
	{
		//AddErr(GMTIMENULL);
		JD.JD=0.0;
		JD.JDE=0.0;
		JD.JC=0.0;
		JD.JCE=0.0;
		JD.JME=0.0;
		JD.E=GMTIMEFAIL;
		return JD;
	}
	JD.E=0;
	day = (double)ut->tm_mday + ((double)ut->tm_hour+((double)ut->tm_min+((double)ut->tm_sec+delta_ut1)/60.0)/60.0)/24;
	month=ut->tm_mon+1;
	year=ut->tm_year+1900;
    if (month < 3)
    {
        month += 12;
        year--;
    }
	JD.JD=INT(365.25*((double)year+4716.0))+INT(30.6001*((double)month+1.0))+ day - 1524.5;
    if (JD.JD > 2299160.0)
    {
        a = INT((double)year/100.0);
        JD.JD += (2 - a + INT(a/4));
    }
	
	// Julian Ephemeris Day
	JD.JDE=JD.JD+delta_t/86400.0;
	JD.JC=(JD.JD-JD0)/36525.0;
	JD.JCE=(JD.JDE-JD0)/36525.0;
	JD.JME=JD.JCE/10.0;
	return JD;
}



/* Heliocentric Earth Coordinate ******************************************
 */
double SummPTerms(const p_term p[], int N, JulianDate JD)
{
    int i=0;
    double s=0;
    for (i=0;i<N;i++)
		s+=p[i].A*cos(p[i].P+p[i].W*JD.JME);
    return s;
}

double Heliocentric_lon(JulianDate JD)
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

double Heliocentric_lat(JulianDate JD)
{
	double lat;
	lat=SummPTerms(EarthLat0, N_LAT0, JD);
	lat+=SummPTerms(EarthLat1, N_LAT1, JD)*JD.JME;
	lat/=1.0e8;
	return lat;
}

double Heliocentric_rad(JulianDate JD)
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

GeoCentricSolPos Geocentric_pos(JulianDate JD)
{
	GeoCentricSolPos P;
	P.lat=fmod(-Heliocentric_lat(JD), 2*M_PI);
	P.lon=fmod(Heliocentric_lon(JD)+M_PI, 2*M_PI);
	if (P.lon<0)
		P.lon+=2*M_PI;
	P.rad=Heliocentric_rad(JD);
	return P;
}
/* nutation */

// compute polynomial of order N-1
// y=a[0]*x^(N-1)+a[1]*x^(N-2)+...+a[N-1]
static inline double poly(const double a[], int N, double x)
{
	int i;
	double r;
	r=a[0];
	for (i=1;i<N;i++)
		r=a[i]+x*r;
	return r;
}	
/* here we deefine several 3rd order polynomials */
const double MEAN_ELONGATION_MOON_SUN[] = 	{deg2rad(1.0/189474.0) , deg2rad(-1.9142e-03), deg2rad(445267.11148) , deg2rad(297.85036)};
const double MEAN_ANOMALY_SUN[] 		= 	{deg2rad(-1.0/300000.0), deg2rad(-0.0001603) , deg2rad(35999.05034)  , deg2rad(357.52772)};
const double MEAN_ANOMALY_MOON[]		= 	{deg2rad(1.0/56250.0)  , deg2rad(0.0086972)  , deg2rad(477198.867398), deg2rad(134.96298)};
const double ARG_LAT_MOON[] 			=	{deg2rad(1.0/327270.0) , deg2rad(-0.0036825) , deg2rad(483202.017538), deg2rad(93.27191)};
const double ASC_LON_MOON[] 			=	{deg2rad(1.0/450000.0) , deg2rad(0.0020708)  , deg2rad(-1934.136261) , deg2rad(125.04452)};

void Nutation_lon_obliquity(JulianDate JD, double *del_psi, double *del_epsilon)
{
    int i, j;
    double sum, sum_psi=0, sum_epsilon=0;
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
			
        sum_psi     += (PE_Terms[i][0] + JD.JCE*PE_Terms[i][1])*sin(sum);
        sum_epsilon += (PE_Terms[i][2] + JD.JCE*PE_Terms[i][3])*cos(sum);
    }

    *del_psi     = sum_psi * deg2rad(1/36000000.0);
    *del_epsilon = sum_epsilon * deg2rad(1/36000000.0);
}

/* solpos 
 * input: longitude (lon) latitude (lat) elevation (e) pressure (p)
 * input: longitude (lon) latitude (lat) elevation (e) pressure (p)
 * 
 * 
 */
const double ECLIPTIC_MEAN_OBLIQUITY[] = {2.45,5.79,27.87,7.12,-39.05,-249.67,-51.38,1999.25,-1.55,-4680.93,84381.448};
const double GSTA[] = {deg2rad(-1/38710000.0),deg2rad(0.000387933),0.0,deg2rad(280.46061837)};

sol_pos solpos(sol_pos P, double lon, double lat, double e, double p, double a_refr, double T, JulianDate JD, GeoCentricSolPos GP)
{
	double dtau, v, H, u, x, y;
	double lambda, alpha, delta, xi;
	double delta_prime, H_prime;
	//double aplpha_prime;
	double dpsi, deps, eps, dalpha;
	double ele, dele=0;
	
	// abberation correction
	dtau=deg2rad(-20.4898/3600.0)/GP.rad;
	// nutation
	Nutation_lon_obliquity(JD, &dpsi, &deps);
	
	// aparent sun longitude
	lambda=GP.lon+dpsi+dtau;
	
	// mean obliquity of the ecliptic
	eps=deps+poly(ECLIPTIC_MEAN_OBLIQUITY,11,JD.JME/10)*deg2rad(1/3600.0);
	
	// sidereal time at Greenwich This seems off!
	v=poly(GSTA,4,JD.JC)+deg2rad(360.98564736629)*(JD.JD - JD0);
	v+=dpsi*cos(eps);
	/*v=fmod(v,2*M_PI);
	if (v<0)
		v+=2*M_PI;
	*/
	
	
	// sun right ascension
	alpha=atan2(sin(lambda)*cos(eps)-tan(GP.lat)*sin(eps),cos(lambda));
	if (alpha<0)
		alpha+=2*M_PI;
	// sun declination
	delta=asin(sin(GP.lat)*cos(eps)+cos(GP.lat)*sin(eps)*sin(lambda));
	
	// hour angle
	H=v+lon-alpha; 
	/* if we compare with NREL's spa we need to limit angular ranges
	H=fmod(H,2*M_PI);
	if (H<0)
		H+=2*M_PI;*/
	
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
	ele=asin(sin(lat)*sin(delta_prime)+cos(lat)*cos(delta_prime)*cos(H_prime));
	if (isnan(a_refr))
	{
		// my solution
		dele=(p/1010.0)*((10-ABSOLUTEZERO)/(T-ABSOLUTEZERO))*deg2rad(BENNET_A)/tan(ele+deg2rad(deg2rad(BENNET_B))/(ele+deg2rad(BENNET_C)));
		if (ele>=-(SUN_RADIUS+dele))
			dele=0.0;
	}
	else
	{
		// NREL's solution.
		// seems inconsistent, we pass the atmosphereic correction and use it only to 
		// see if we need to compute it??
		if (ele>=-(SUN_RADIUS+a_refr))
			dele=(p/1010.0)*((10-ABSOLUTEZERO)/(T-ABSOLUTEZERO))*deg2rad(BENNET_A)/tan(ele+deg2rad(deg2rad(BENNET_B))/(ele+deg2rad(BENNET_C)));
	}
	
	// compute sun zenith
	P.s.z=M_PI/2-ele;
	
	// compute aparent sun zenith
	P.sa.z=P.s.z-dele;
	P.s.a=fmod(M_PI+atan2(sin(H_prime),cos(H_prime)*sin(lat)-tan(delta_prime)*cos(lat)),2*M_PI);
	P.sa.a=P.s.a;
	
	// limit angular range
	P.s.z=fmod(P.s.z,2*M_PI);
	P.sa.z=fmod(P.sa.z,2*M_PI);
	if (P.s.z<0)
	{
		P.s.z=-P.s.z;
		P.s.a=fmod(P.s.a+M_PI,2*M_PI);
	}
	if (P.sa.z<0)
	{
		P.sa.z=-P.sa.z;
		P.sa.a=fmod(P.s.a+M_PI,2*M_PI);
	}
	if (P.s.z>M_PI)
	{
		P.s.z=2*M_PI-P.s.z;
		P.s.a=fmod(P.s.a+2*M_PI,2*M_PI);
	}
	if (P.sa.z>M_PI)
	{
		P.sa.z=2*M_PI-P.sa.z;
		P.sa.a=fmod(P.s.a+M_PI,2*M_PI);
	}
	return P;
}


int InputCheck(double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double a_refr, double T)
{
	int E=0;
	
	if (delta_t<-8000)
		E|=DELTA_T_OOR;
	else if (delta_t>8000)
		E|=DELTA_T_OOR;
	
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
		
	if (a_refr<-5*M_PI/180)
		E|=AREFR_OOR;
	if (a_refr>5*M_PI/180)
		E|=AREFR_OOR;
		
	if (T<ABSOLUTEZERO)
		E|=TEMP_OOR;
	return E;
}

sol_pos SPA_NREL(time_t t, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double a_refr, double T)
{
	JulianDate D;
	GeoCentricSolPos G;
	sol_pos P;
	P.E=InputCheck(delta_t, delta_ut1, lon, lat, e, p, a_refr, T);
	if (!P.E)
	{
		D=MakeJulianDate(t, delta_t, delta_ut1);
		if (D.E)
		{
			P.E=D.E;
			return P;
		}
		G=Geocentric_pos(D);
		P=solpos(P,lon,lat,e,p,a_refr,T,D,G);
	}
	return P;
}
sol_pos SPA(time_t t, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T)
{
	JulianDate D;
	GeoCentricSolPos G;
	sol_pos P;
	P.E=InputCheck(delta_t, delta_ut1, lon, lat, e, p, 0, T);
	if (!P.E)
	{
		D=MakeJulianDate(t, delta_t, delta_ut1);
		if (D.E)
		{
			P.E=D.E;
			return P;
		}
		G=Geocentric_pos(D);
		P=solpos(P,lon,lat,e,p,NAN,T,D,G);
	}
	return P;
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
	JulianDate D;
	struct tm *ut;
	timestr=malloc(50*sizeof(char));
	for (i=0;i<N;i++)
	{
		ut=gmtime(&dates[i]);	
		strftime(timestr, 50, "%Y-%m-%d %T %Z",ut);
		printf("testing %s", timestr);
		D=MakeJulianDate(dates[i], 0,0);
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
	JulianDate D;
	D=MakeJulianDate(date, delta_t, delta_ut1);

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




