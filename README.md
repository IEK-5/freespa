# freespa
## About freespa
This is a free implementation in c of the Solar Position Algorithm (SPA) [1-2]. Another open and free c implementation is provided by NREL [here](https://midcdmz.nrel.gov/spa/). NREL's SPA is free as in beer! However, it's like beer which you are not allowed to share with your friends (the License does not allow redistribution). The freespa code, on the other hand, you may share under the conditions of the [GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Installation
Package is not designed to be installed but rather to be inserted in a project. You can add freespa as a submodule to your code.

## Usage
The main routine is: 

`sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e, double p, double T);`   

where:

* `struct tm *ut`: Standard time struct. Should contain UTC values
* `double *delta_t`: pointer to delta t value, if it is NULL dleta t is determined from internal tables
* `double delta_ut1`: deviation between terrestial time and UTC (-1.0..1.0 s)
* `double lon`: longitude in radians
* `double lat`: latitude in radians
* `double e`: Elevation in m
* `double p`: pressure in mb
* `double T`: Temperature in C

The return value consists of a data struct:

`typedef struct sol_pos {`

`	double z, a, az, aa; // zenith, azimuth, aparent zenith, aparent azimuth`
	
`	int E; // error flag`
	
`} sol_pos;`

If the error flag equals 0 all is OK. The double values are the true zenith and azimuth angles (z,a), and the aparant zenith and azimuth (az,aa), including atmospheric refraction effects. All angles are again in radians. 

In addition freespa defines several other routines:

 - `struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, double lon, double lat)`
	- Computes the true solar time for a given coordinate and time struct with UTC time
 - `int SunTimes(struct tm ut, double *delta_t, double delta_ut1, double lon, double lat, double e, double p, double T, struct tm *sunrise, struct tm, *transit, struct tm *sunset)`
	- Computes sunrise, transit, and sunset times

As the implemented SPA routine is accurate over a long period of time, it is useful to extend the unix time concept a bit. To this end freespa reimplements the gmtime and mkgmtime equivalent routines. The freespa routines should, for the present time, be equivalent to the standard routines in time.h. However, freespa's routines consider the 10-day gap between the Julian and Gregorian calendar where the Julian calendar ends on October 4, 1582 (JD = 2299160), and the next day the Gregorian calendar starts on October 15, 1582. The impemented routines are:

 - `struct tm *gmjtime_r(time_t *t, struct tm *ut)`
	- converts unix time to a time struct with UTC (thread safe)
 - `struct tm *gmjtime(time_t *t)`
	- converts unix time to a time struct with UTC (not thread safe)
 - `time_t mkgmjtime(struct tm *ut)`
	- converts a time struct with UTC to unix time

## Δt
For accurate timing the SPA algorithm needs Δt values. Both historic and predicted future values may be obtained here:
[https://maia.usno.navy.mil/ser7/](https://maia.usno.navy.mil/ser7/)

Another source of Δt values is

[https://hpiers.obspm.fr/eop-pc/index.php?index=analysis&lang=en](https://hpiers.obspm.fr/eop-pc/index.php?index=analysis&lang=en)

This database provides Δt values over various time ranges in various resolutions.

Using these datasources we compiled a table of Δt values in freespa_dt_table.h. Beyond the values in the table we use long term predictions of Δt according to the empirical model from Morrison and  Stephenson [4]. It is possible to update the tables in freespa_dt_table.h. To this end the "GetDeltaT.sh" script is provided, which generates a new header "freespa_dt_table_new.h", which can be inspected and included if correct. The script downloads new  Δt values and mid-term predictions from [https://maia.usno.navy.mil/ser7/](https://maia.usno.navy.mil/ser7/).

## Differences w.r.t. NREL's SPA
In general freespa and NREL's SPA produce identical results. In some small details there are differences stemming from slightly different choises. Note that some of the behavior of NREL's SPA discussed here is only documented only in the spa.c source code provided by NREL [3].
 
To save the caller the burden to find appropriate Δt values, freespa includes both tabular historic Δt data and extrapolation routines. However, the user may still provide own Δt data. Here it should also be noted that Δt in NREL's SPA is limited to to the range -8000 -- 8000 s, whereas freespa imposes no limit (accepts any double float value). I found no reference nor rationale behind this seemingly arbitrary limit in NREL's SPA. However, using the historic data and long term Δt predictions from Morrison and  Stephenson, [4], it follows that Δt<8000 only holds for the approximate period 245 -- 3400, much shorter than the claimed validity for the period -2000 -- 6000 [1-2]. 

Both NREL's SPA and freespa can compute sunrise, transit, and sunset. Accurately predicting sunrise and sunset is notoriously hard, primarely due to atmospheric refraction effects [5], so the accuracy of these predictions is limited. One should note that NREL's SPA ignores observer elevation when computing sunrise, transit, and sunset. This affects both the solar position w.r.t. the observer and, more importantly, the geometric dip (see also [5], page 27).

## Testing
To test the correct functioning of freespa you may use 

`make check`

which will run various checks and report back.

A simple commandline interfece to freespa may be compiled with

`make spa`

In case you download [NREL's SPA](https://midcdmz.nrel.gov/spa/), this commandline spa interface may also use NREL's SPA as a backend, allowing a direct comparison of the two SPA implementations.

Finally, with

`make compare`

you can compile a program to compare NREL's SPA with freespa.

## References
[1] I.  Reda and A. Andreas, "Solar position algorithm for solar radiation applications." Solar Energy 76.5 (2004): 577-589

[2] I. Reda and A. Andreas, "Corrigendum to Solar position algorithm for solar radiation applications." Solar Energy 81.6 (2007): 838

[3] A. Andreas, spa.c source code retrieved from https://midcdmz.nrel.gov/spa/ on the 27th of march 2023
 
[4] L. V. Morrison and  F. R. Stephenson "Historical Values of the Earth’s Clock Error ΔT and the Calculation of Eclipses." Journal for the History of Astronomy, 35(3) (2004): 327–336. 

[5] T. Wilson "Evaluating the Effectiveness of Current Atmospheric Refraction Models in Predicting Sunrise and Sunset Times." PhD diss., Michigan Technological University, 2018.
