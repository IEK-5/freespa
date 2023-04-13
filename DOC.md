#Documentation freespa
## Requirements
The freespa routines rely on a 64bit signed integer `time_t` type, 
representing the seconds since the epoch. This restricts the 
portability of freespa. In the freespa.h header file there are some
compile time asserts that check whether the time_t type is appropriate,
i.e. the code should not compile on systems where `time_t` is not a 
64bit signed integer. On 64bit linux and windows systems `time_t` 
generally meets the requirements. 

## Units
Freespa uses the following units:

- angles are in radians (this includes coordinates!)
- temperatures in °C
- pressure in mb
- time is generally specified in UTC


## Data Structures:
In freespa.h several data structures are defined. The `sol_pos` 
structure is defined as:

    typedef struct sol_pos {
    	double z, a; // zenith, azimuth
    	int E; // error flag
    } sol_pos;

where z is the zenith angle, and a the azimuth. The integer E is an 
error flag which may contain the following error codes:

    #define _FREESPA_DEU_OOR		0X01	// ΔUT1 out of range
    #define _FREESPA_LON_OOR		0X02	// longitude out of range
    #define _FREESPA_LAT_OOR		0X04	// latitude out of range
    #define _FREESPA_ELE_OOR		0X08	// elevation out of range
    #define _FREESPA_PRE_OOR		0X10	// pressure out of range
    #define _FREESPA_TEM_OOR		0X20	// temperature out of range
    #define _FREESPA_DIP_OOR		0X40	// geometric dip out of range
    #define _FREESPA_GMTIMEF		0X80	// time conversion error

which may be combined with a binary OR. If all is OK E=0.

In case you are interested in computing the daily events such as 
sunrise and set, you need the `solar_day` struct:
 
    typedef struct solar_day {
    	struct tm ev[11];
    	time_t t[11];
    	double E[11];
    	int status[11];
    } solar_day;

This struct contains the following elements:

- ev: an array of time structs for all events of the day
- t: corresponding unix times
- E: array with error values (deviation of solar zenith angle, in rad)
- status: integer array with status flags

The solar events are defined and indexed as:

     0:		solar midnight before time t
     1:		solar transit closest to time t
     2:		solar midnight after time t
     3:		sunrise
     4:		sunset
     5:		civil dawn
     6:		civil dusk
     7:		nautical dawn
     8:		nautical dusk
     9:		astronomical dawn
     10:	astronomical dusk

Note that in this definition of a "solar day", not all days have 86400

The status flag may contain the following values:

    #define _FREESPA_EV_ERR       20  // Error
    #define _FREESPA_EV_NA        10  // Not Computed
    #define _FREESPA_EV_OK         0  // All OK
    #define _FREESPA_EV_SUNABOVE   1  // Sun always above (e.g. midnight sun)
    #define _FREESPA_EV_SUNBELOW  -1  // Sun always below (e.g. polar night)

## Input Data Ranges
There are some limits to the accepted input ranges of various arguments.
In general freespa does not limit much and thus allows for some pretty 
unreasonable input, use at your own discretion.

Valid input to freespa routines adheres to *at least* the following ranges:

- ΔUT1:  -1 ≤ ΔUT1 ≤1
- longitude: -π ≤ lon ≤ π
- latitude: -π/2 ≤ lat ≤ π/2
- elevation: Rearth < E, where Rearth = 6378136.6 m
- pressure: 0 ≤ p ≤ 5000
- Temperature: -273.15 °C ≤ T

No other limits are imposed on other input variables, such as Δt, thus
the limits depends on the underlying data types.


## Main SPA Routines
The main routine to compute the real solar position is:

`sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e);`   

where:

* `struct tm *ut`: Standard time struct. Should contain UTC values
* `double *delta_t`: pointer to Δt value, if it is NULL Δt is determined from internal tables
* `double delta_ut1`: deviation between terrestrial time and UTC (-1.0..1.0 s)
* `double lon`: longitude in radians
* `double lat`: latitude in radians
* `double e`: Elevation in m

This computes the _real_ solar position. In practice the solar position 
is affected by refraction. To compute the apparent position of the sun 
several routines may be used:

    sol_pos ApSolposBennet(sol_pos P, double *gdip, double e, double p, double T);
    
and

    sol_pos ApSolposBennetNA(sol_pos P, double *gdip, double e, double p, double T);

Both routines work the same and only differ in the used model 
coefficients. The input for these routines is as follows:

* `P`:	    real solar position
* `gdip`:	geometric dip, i.e. how far the horizon is below the observer (in rad). If this pointer is NULL the geometric dip is computed from the observer elevation (assuming the horizon is at sea level)
* `e`:		observer elevation (in meter)
* `p`:		pressure (in mbar)
* `T`:		Temperature (in °C)

The solar refraction is computed with the simple Bennet equation [1]. 
The BennetNA routines are based on the modified coefficients as 
published in [2].


To compute the solar events of a day freespa offers:

    solar_day SolarDay(struct tm *ut, double *delta_t, double delta_ut1, 
                       double lon, double lat, double e, double *gdip, 
                       double p, double T, 
                       sol_pos (*refract)(sol_pos,double*,double,double,double));

This returns a solar_day struct. The input is:

* `ut`:	       pointer to time struct with UTC time
* `delta_t`:   pointer to Δt value, or NULL (use internal tables)
* `delta_ut1`: delta_ut1 
* `lon`:       longitude (in radians)
* `lat`:       latitude (in radians)
* `e`:         observer elevation (in meter)
* `gdip`:      geometric dip
* `p`:         pressure (in mbar)
* `T`:         Temperature (in °C)

In some cases a user may not be interested in all the computed events. 
It is possible to enable/disable computing specific events. To this end
one can modify the SDMASK integer (globally defined, i.e. not thread safe).
The following flags are defined

    #define _FREESPA_SUNRISE 0X01
    #define _FREESPA_SUNSET  0X02
    #define _FREESPA_CVDAWN  0X04
    #define _FREESPA_CVDUSK  0X08
    #define _FREESPA_NADAWN  0X10
    #define _FREESPA_NADUSK  0X20
    #define _FREESPA_ASDAWN  0X40
    #define _FREESPA_ASDUSK  0X80

By using a binary OR operation one can enable or disable the computation 
of specific events. Per default SDMASK is defined as:

    SDMASK==(_FREESPA_SUNRISE|_FREESPA_SUNSET|_FREESPA_CVDAWN|_FREESPA_CVDUSK|_FREESPA_NADAWN|_FREESPA_NADUSK|_FREESPA_ASDAWN|_FREESPA_ASDUSK)

or 

    SDMASK=0XFF;

which triggers the computation of all 11 solar day events. If one only 
wants to compute sunrise and sunset, one can define:
 
    SDMASK=(_FREESPA_SUNRISE|_FREESPA_SUNSET);


## Time Utilities
Finally freespa offers several utilities to work with time. The 
equation of time describes the variation of the true solar time w.r.t. 
mean solar time. Freespa offers the `TrueSolarTime` routine to 
determine the true solar time:

    struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, 
						    double lon, double lat);

The function returns a time struct representing the true solar time 
(the only non-UTC time used in freespa). The input to this function is:

* `ut`:	       pointer to time struct with UTC time
* `delta_t`:   pointer to Δt value, or NULL (use internal tables)
* `delta_ut1`: delta_ut1 
* `lon`:	   longitude (in radians)
* `lat`:	   latitude (in radians)

To overcome some of the implementation dependent behavior of `gmtime` 
and `mkgmtime`, and, in addition, work with time structs and unix time 
in a consistent manner with the internally used Julian day, freespa 
defines it own conversion routines:

    struct tm *gmjtime_r(time_t *t, struct tm *ut);
    struct tm *gmjtime(time_t *t);
    time_t mkgmjtime(struct tm *ut);

These routines work the same as the standard routines defined in 
`time.h`. They do not suffer from arbitrary limitations such as those 
imposed in windows systems where `gmtime` cannot handle dates before 
1970 and after the year 3000 (despite the `time_t` type being a 64bit 
integer). Furthermore, they adhere to the 10-day gap between the Julian 
and Gregorian calendar where the Julian calendar ends on October 4, 
1582 (JD = 2299160), and the next day the Gregorian calendar starts 
on October 15, 1582. Thus these routines provide a historic extension 
of unix time before October 15, 1582.

For many freespa routines we need Δt values. In most cases one can suffice 
with passing a NULL pointer, which will signal freespa to determine an 
appropriate value from internal tables. The following function: 

    double get_delta_t(struct tm *ut);

provides an interface to freespa's internal Δt tables. The internal tables 
contain a historic dataset ranging from the year -2000 to the present day. 
In addition it uses predictions up to 2034. Beyond this range freespa 
resorts to a crude model by Morrison and Stephenson [3]
 
                        2
                ⎡y-1820⎤
     Δt(y) = 32 ⎢──────⎥  - 20
                ⎣ 100  ⎦
 
where y is the year, and Δt is in seconds.

## Differences with NREL spa
Freespa was developed due to license issues with NREL's spa code. As such 
the goal was to more or less re-implement NREL spa. Thus the results are 
generally identical or very similar. The interface is, however, 
different (i.e. freespa is _not_ a drop-in replacement). In addition 
there are some differences in computational details. Some of those 
details might be considered bugs. Note that some of the behavior of 
NREL spa discussed here is not explicitly documented but can be found 
in the source code of NREL spa [4].
 
One obvious difference is that freespa provides a simple interface to 
determine Δt values form internal tables, if so desired. It should be 
noted that NREL spa limits the range of acceptable Δt values to +-8000. 
This seems a rather arbitrary limit that ignores the fact that this 
limit only holds for the time period ranging from approximately 
245 -- 3400 [3]. This is probably fine for most applications. 
Nevertheless, the authors of NREL spa claim the model is accurate for 
the period -2000 -- 6000 [5].

Both NREL's SPA and freespa can compute sunrise, transit, and sunset. 
The sunrise/transit/sunset routines in NREL's spa are probably more 
efficient. However, the accuracy of NREL spa is not guaranteed and I 
incidentally find offsets in zenith angles of up to 0.1°, even when 
limiting the range of latitudes between -65--65°. One might argue that 
determining the exact sunrise and set times is always inaccurate due 
atmospheric refraction effects [2]. On the other hand, it is 
unnecessary to add model uncertainties to a problem that is already 
hard. The freespa routines use iterative methods to find an accurate 
solution (within one second _or_ within 0.05°). 

Another limitation of NREL's spa sunrise/set routine is that it does 
not compute a sunrise/set event if _one_ of the two events does not 
happen. For example, for the the first day in a year that the sun does 
not set at some location (e.g. high up north somewhere in summer), 
NREL's spa will not compute when the sun rises.

In NREL's spa implementation of an atmospheric refraction model there 
is no geometric dip. This means that NREL spa always assumes the 
horizon is at an elevation of 0°. At higher elevation, this assumption 
may be inaccurate.

In addition to computing sunrise/set times freespa can compute various 
dawn and dusk times (civil, nautical, and astronomical dawn and dusk).

## References
[1] J. Meeus, Astronomical Algorithms, second ed. Willmann-Bell, Inc., Richmond, Virginia, USA. (1998): 105-108

[2] T. Wilson, "Evaluating the Effectiveness of Current Atmospheric Refraction Models in Predicting Sunrise and Sunset Times", Open Access Dissertation, Michigan Technological University, (2018). https://doi.org/10.37099/mtu.dc.etdr/697

[3] L. V. Morrison, and F. R. Stephenson, "Historical Values of the Earth’s Clock Error ΔT and the Calculation of Eclipses." Journal for the History of Astronomy, 35(3), (2004): 327–336. https://doi.org/10.1177/002182860403500305

[4] A. Andreas, spa.c source code retrieved from https://midcdmz.nrel.gov/spa/ on the 27th of march 2023

[5] I. Reda and A. Andreas, "Solar position algorithm for solar radiation applications." Solar Energy 76.5 (2004): 577-589
