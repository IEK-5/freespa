

typedef struct sky_pos {
	double z,a;
} sky_pos;

// container struct for solar position data
typedef struct sol_pos {
	sky_pos s;
	sky_pos sa;
	int E;
	//time_t sunrize, sunset, noon;
} sol_pos;

int testjulian();
int testheliocentricpos();
sol_pos SPA_NREL(time_t t, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double a_refr, double T);
sol_pos SPA(time_t t, double delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T);
