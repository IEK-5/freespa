

// container struct for solar position data
typedef struct sol_pos {
	double z, a, az, aa; // zenith, azimuth, aparent zenith, aparent azimuth
	int E;
	//time_t sunrize, sunset, noon;
} sol_pos;


sol_pos SPA(time_t t, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T);
            
int testjulian();
int testheliocentricpos();
