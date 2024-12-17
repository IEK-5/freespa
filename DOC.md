# Freespa Documentation

Freespa is a c library to compute the solar position. Just like NREL's 
SPA [Reda], it implements the algorithms developed by [Meeus], and thus 
provides more or less identical functionality. Freespa was developed to 
address licensing issues with NREL SPA code [NRELSPA]. 


This document is a technical description of Freepa's routines and data 
structures. The Outline is a s follows:

* [Units](#units) : Physical Units
* [Data Types and Structures](#data-types-and-structures): Describes the
data structures and types defined by the library
* [Solar Position Routines](#solar-position-routines): Routines to 
compute the real and aparent position of the sun
* [Solar Day Events](#solar-day-events): Routines to compute the solar 
events (e.g. sunrise, sunset, etc.)
* [Equinoxes and Solstices](#equinoxes-and-solstices): Routines to 
compute the Equinoxes and Solstices.
* [Δt Values](#δt-values):  Routines to obtain Δt Values
* [Differences with NREL SPA](#differences-with-nrel-spa): A short 
description of the main differences with NREL SPA

## Units
Freespa uses the following units:

- angles are in radians (this includes coordinates!)
- temperatures in °C
- pressure in mb
- time is generally specified in UTC



## Data Types and Structures

### `time_t` Type

By default, Freespa uses the `time_t` type to represent time as the 
number of seconds since the epoch (1970-01-01 00:00:00 UTC). However, 
the Freespa routines require a 64-bit signed integer type to represent 
time, ensuring sufficient range for various computations. 

On some platforms, the default `time_t` implementation may not meet 
these requirements. For example:

- **32-bit Linux**: `time_t` is a 32-bit signed integer.
- **Windows (32-bit)**: `time_t` is also a 32-bit signed integer.

If your platform's `time_t` is not a 64-bit signed integer, you can 
still use Freespa by defining a custom time type.

On a side note: Windows (64-bit) `time_t` is a signed 64-bit integer, 
which can be used. However, time handling routines defined in time.h, 
such as `gmtime`, for some reason cannot handle dates before the epoch 
or after Dec 31 23:59:59 3000.

#### Customizing the `time_t` Type

Freespa determines the time type via the preprocessor variable 
`FS_TIME_T`:

- **Default Behavior**:  
  If `FS_TIME_T` is not defined, Freespa defaults to `time_t`.

- **Custom Behavior**:  
  Define `FS_TIME_T` to provide an alternative type for representing 
  time. For example, passing the following to your compiler:

  ```c
  -DFS_TIME_T=int64_t
  ```

  This tells Freespa to use `int64_t` for time representation instead 
  of `time_t`.

#### Using a Custom Type

When a custom `FS_TIME_T` is used, Freespa defines the preprocessor 
macro `FS_CUSTOM_TIME_T`. You can use this macro to check if a custom 
time type is being used in your code:

```c
#ifdef FS_CUSTOM_TIME_T
    // Custom time type is in use
#endif
```


#### Important Considerations

1. **Standard Library Incompatibility**:  
   Using a custom `FS_TIME_T` type means you cannot use standard 
   `time.h` functions that rely on `time_t`. Instead, you will need to 
   rely on Freespa's utility functions for epoch time manipulation 
   (see [Time Utilities](#time-utilities)).

2. **Compile-Time Assertion**:  
   Freespa validates that the specified `FS_TIME_T` is a 64-bit signed 
   integer at compile time. If the type is incompatible, the 
   compilation will fail with an error.
   


### Data Structures and Input Ranges
#### `sol_pos` Structure

The `sol_pos` structure, defined in `freespa.h`, is used to represent 
the position of the sun. It includes the zenith angle, azimuth, and an 
error flag.

```c
typedef struct sol_pos {
    double z; // Zenith angle
    double a; // Azimuth
    int E;    // Error flag
} sol_pos;
```

##### Error Flags (`E`)

The error flag `E` can take the following values, which may be combined 
using a binary OR (`|`):

```c
#define FREESPA_DEU_OOR  0X01  // ΔUT1 out of range
#define FREESPA_LON_OOR  0X02  // Longitude out of range
#define FREESPA_LAT_OOR  0X04  // Latitude out of range
#define FREESPA_ELE_OOR  0X08  // Elevation out of range
#define FREESPA_PRE_OOR  0X10  // Pressure out of range
#define FREESPA_TEM_OOR  0X20  // Temperature out of range
#define FREESPA_DIP_OOR  0X40  // Geometric dip out of range
#define FREESPA_GMTIMEF  0X80  // Time conversion error
```

- If all inputs are valid, `E = 0`.


#### `solar_day` Structure

The `solar_day` structure represents daily solar events, including 
times, Unix timestamps, zenith angle deviations, and status flags.

```c
typedef struct solar_day {
    struct tm ev[11];    // Array of time structs for events
    FS_TIME_T t[11];     // Corresponding Unix times
    double E[11];        // Deviation of solar zenith angle (radians)
    int status[11];      // Status flags for each event
} solar_day;
```

##### Solar Event Indices

The array indices correspond to the following solar events:

| Index | Event                             |
|-------|-----------------------------------|
| 0     | Solar midnight before time `t`   |
| 1     | Solar transit closest to time `t`|
| 2     | Solar midnight after time `t`    |
| 3     | Sunrise                          |
| 4     | Sunset                           |
| 5     | Civil dawn                       |
| 6     | Civil dusk                       |
| 7     | Nautical dawn                    |
| 8     | Nautical dusk                    |
| 9     | Astronomical dawn                |
| 10    | Astronomical dusk                |

**Note**: In this definition of a "solar day," not all days have exactly 
86400 seconds.

##### Status Flags

The `status` array contains flags indicating the status of each solar 
event:

```c
#define FREESPA_EV_ERR        20  // Error
#define FREESPA_EV_NA         10  // Not computed
#define FREESPA_EV_OK          0  // All OK
#define FREESPA_EV_SUNABOVE    1  // Sun always above (e.g., midnight sun)
#define FREESPA_EV_SUNBELOW   -1  // Sun always below (e.g., polar night)
```


#### Input Data Ranges

Valid input ranges for Freespa routines are outlined below. Freespa 
generally imposes minimal restrictions, allowing for unusual inputs 
where reasonable. 

| Parameter     | Valid Range                                         |
|---------------|-----------------------------------------------------|
| **ΔUT1**      | $-1 ≤ ΔUT1 ≤ 1$                                     |
| **Longitude** | $-π ≤ \text{lon} ≤ π$                                      |
| **Latitude**  | $-π/2 ≤ \text{lat} ≤ π/2$                                  |
| **Elevation** | $E > -R_\mathrm{earth} (\text{where }R_\mathrm{earth} = 6378136.6\,\mathrm{m})$        |
| **Pressure**  | $0 ≤ p ≤ 5000$                                        |
| **Temperature** | $-273.15 °C ≤ T$                                  |

##### Notes

- **Δt (Delta T)**: Freespa does not impose explicit limits on Δt. 
Valid ranges depend on the underlying data types used in your 
implementation.
- **Other Variables**: For any other inputs not mentioned here, limits 
depend on their respective data types.


Here’s an improved and clearer version of your documentation:


## Solar Position Routines

### Real Solar Position

The primary routine for computing the real solar position is:

```c
sol_pos SPA(struct tm *ut, double *delta_t, double delta_ut1, double lon, 
            double lat, double e);
```

#### Parameters

- **`struct tm *ut`**:  
  Pointer to a standard `tm` structure containing UTC time values.

- **`double *delta_t`**:  
  Pointer to the Δt value. If `NULL`, Δt is determined from internal tables.

- **`double delta_ut1`**:  
  Deviation between terrestrial time and UTC ( $-1.0 ≤$ `delta_ut1` $≤ 1.0$ s).

- **`double lon`**:  
  Longitude in radians.

- **`double lat`**:  
  Latitude in radians.

- **`double e`**:  
  Observer elevation in meters.

This function calculates the **real solar position**, which does not 
account for atmospheric refraction.


### Apparent Solar Position

To compute the **apparent position** of the sun, accounting for 
atmospheric refraction, use the following routines:

#### `ApSolposBennet`

```c
sol_pos ApSolposBennet(sol_pos P, double *gdip, double e, double p, double T);
```

#### `ApSolposBennetNA`

```c
sol_pos ApSolposBennetNA(sol_pos P, double *gdip, double e, double p, double T);
```

Both functions work similarly, differing only in their model coefficients:
- `ApSolposBennet`: Uses the original Bennet equation [Meeus].
- `ApSolposBennetNA`: Uses modified coefficients from [Wilsen].

#### Parameters

- **`P`**:  
  The real solar position (`sol_pos` struct) computed by `SPA`.

- **`gdip`**:  
  Pointer to the geometric dip (in radians). Represents how far the 
  horizon is below the observer. If `NULL`, it is computed based on 
  observer elevation, assuming the horizon is at sea level.

- **`e`**:  
  Observer elevation in meters.

- **`p`**:  
  Atmospheric pressure in mbar.

- **`T`**:  
  Temperature in degrees Celsius.


## Solar Day Events

To compute all solar events for a given day, use:

```c
solar_day SolarDay(struct tm *ut, double *delta_t, double delta_ut1, 
                   double lon, double lat, double e, double *gdip, 
                   double p, double T, 
                   sol_pos (*refract)(sol_pos,double*,double,double,double));
```

### Parameters

- **`ut`**:  
  Pointer to a `tm` structure containing UTC time.

- **`delta_t`**:  
  Pointer to the Δt value. If `NULL`, Δt is determined from internal 
  tables.

- **`delta_ut1`**:  
  Deviation between terrestrial time and UTC.

- **`lon`**:  
  Longitude in radians.

- **`lat`**:  
  Latitude in radians.

- **`e`**:  
  Observer elevation in meters.

- **`gdip`**:  
  Pointer to geometric dip (see `ApSolposBennet`).

- **`p`**:  
  Atmospheric pressure in mbar.

- **`T`**:  
  Temperature in degrees Celsius.

- **`refract`**:  
  Pointer to a refraction model (e.g., `ApSolposBennet`). If `NULL`, 
  true solar positions are used.

### Returns

A `solar_day` structure containing solar events, including their times, 
Unix timestamps, zenith angle deviations, and status flags.


### Event Selection with `SDMASK`

In some cases, you may want to compute only specific solar events. The 
`SDMASK` integer can be modified globally (not thread-safe) to enable 
or disable specific events. 

#### Event Flags

The following flags are defined:

```c
#define FREESPA_SUNRISE  0x01  // Sunrise
#define FREESPA_SUNSET   0x02  // Sunset
#define FREESPA_CVDAWN   0x04  // Civil dawn
#define FREESPA_CVDUSK   0x08  // Civil dusk
#define FREESPA_NADAWN   0x10  // Nautical dawn
#define FREESPA_NADUSK   0x20  // Nautical dusk
#define FREESPA_ASDAWN   0x40  // Astronomical dawn
#define FREESPA_ASDUSK   0x80  // Astronomical dusk
```

Combine flags using a binary OR (`|`) to enable multiple events. 

#### Default Configuration

By default, `SDMASK` is defined as:

```c
SDMASK = (FREESPA_SUNRISE | FREESPA_SUNSET | FREESPA_CVDAWN | FREESPA_CVDUSK |
          FREESPA_NADAWN | FREESPA_NADUSK | FREESPA_ASDAWN | FREESPA_ASDUSK);
```

This computes all 11 solar day events. Alternatively, it can be expressed as:

```c
SDMASK = 0xFF;
```

#### Example: Sunrise and Sunset Only

To compute only sunrise and sunset:

```c
SDMASK = (FREESPA_SUNRISE | FREESPA_SUNSET);
```

## Equinoxes and Solstices

The `Freespa` library provides functionality to compute the year's 
Equinoxes and Solstices. These events can be calculated as either a 
`struct tm` or an `FS_TIME` integer. 

### `mkgmEQSOtime`

```c
struct tm *mkgmEQSOtime(struct tm *ut, int E, double *delta_t);
```

This function computes the specified Equinox or Solstice event and 
returns the result as a `struct tm`. 

### `mkgmEQSOjtime`

```c
FS_TIME_T mkgmEQSOjtime(struct tm *ut, int E, double *delta_t);
```

This function computes the specified Equinox or Solstice event and 
returns the result as an `FS_TIME_T` integer.


### Parameters

Both functions take the following parameters:

- **`ut`**:  
  Pointer to a `struct tm` containing the UTC time to base the calculation on.  
  **Note**: This input is overwritten with the computed event's time.

- **`E`**:  
  The event to compute, represented as an integer. (See [Event Definitions](#event-definitions) below.)

- **`delta_t`**:  
  Pointer to a `double` representing the value of Δt (the difference between Earth's observed time and uniform time).  
  Pass `NULL` to use the library's internal Δt tables.

---

### Event Definitions

The following constants are defined in `freespa.h` to specify the type of Equinox or Solstice event (`E` parameter):

```c
#define FREESPA_SPRINGEQ  0  // Vernal (Spring) Equinox
#define FREESPA_SUMMERSO  1  // Summer Solstice
#define FREESPA_AUTUMNEQ  2  // Autumnal (Fall) Equinox
#define FREESPA_WINTERSO  3  // Winter Solstice
```

### Event Mapping

| Constant              | Description                |
|-----------------------|----------------------------|
| `FREESPA_SPRINGEQ`    | Vernal (Spring) Equinox    |
| `FREESPA_SUMMERSO`    | Summer Solstice            |
| `FREESPA_AUTUMNEQ`    | Autumnal (Fall) Equinox    |
| `FREESPA_WINTERSO`    | Winter Solstice            |

---
Here’s an improved and more organized version of your documentation:

---

## Time Utilities

### Equation of Time

The equation of time describes the variation of true solar time 
relative to mean solar time. Freespa offers the `TrueSolarTime` 
function to compute the **true solar time**:

```c
struct tm TrueSolarTime(struct tm *ut, double *delta_t, double delta_ut1, 
                        double lon, double lat);
```

#### Parameters

- **`ut`**:  
  Pointer to a `struct tm` with UTC time.
- **`delta_t`**:  
  Pointer to the Δt value. If `NULL`, Δt is determined from internal 
  tables.
- **`delta_ut1`**:  
  Deviation between terrestrial time and UTC.
- **`lon`**:  
  Longitude in radians.
- **`lat`**:  
  Latitude in radians.

#### Returns

A `struct tm` representing the **true solar time** (the only non-UTC 
time used in Freespa).


### Custom Time Conversion Routines

To handle custom `FS_TIME_T` types, overcome platform-specific 
limitations of `gmtime` and `mkgmtime`, and align time calculations 
with Julian dates, Freespa provides the following utilities:

#### Functions

```c
struct tm *gmjtime_r(FS_TIME_T *t, struct tm *ut);
struct tm *gmjtime(FS_TIME_T *t);
FS_TIME_T mkgmjtime(struct tm *ut);
```

#### Features

- **Custom Time Support**: These functions work with custom `FS_TIME_T` 
types.
- **Platform Independence**: Avoids limitations such as those found in 
64-bit Windows systems.
- **Historic Date Support**: Handles the transition from Julian to 
Gregorian calendar:
  - Julian calendar ends on October 4, 1582 (JD = 2299160).
  - Gregorian calendar starts on October 15, 1582.


## Δt Values

Freespa routines often require Δt values, which represent the 
difference between Earth's observed time (Universal Time) and uniform 
time (Terrestrial Time). Freespa provides internal tables and routines 
to compute these values.

### Internal Tables and Models

1. **Documented Values**:  
   Freespa uses internal tables and interpolation for Δt values where 
   documented data is available. These tables include data from 
   historical observations and predictions.

2. **Predictions**:  
   The United States Naval Observatory (USNO) publishes Δt predictions 
   up to approximately 10 years into the future [7]. These predictions 
   are included in Freespa's internal tables.

3. **Extended Range Model**:  
   For years beyond the range of the internal tables, Freespa employs 
   the Morrison and Stephenson model [Morrison]:

   $$
   \Delta t(y) = 32 \left( \frac{y - 1820}{100} \right)^2 - 20
   $$

   - where $y$ is the Year.

### Querying Δt Values

To compute the Δt value for a given time, you can use the following 
function:

```c
double get_delta_t(struct tm *ut);
```

#### Parameters

- **`ut`**:  
  Pointer to a `struct tm` representing UTC time.

#### Returns

The computed Δt value in seconds, based on Freespa's internal tables or 
the extended model.



## Differences with NREL SPA

As discussed both NREAL SPA and Freespa base on the algorithms by 
[Meeus]. So in general one can expect more of less identical results. 
However, there are some small technical differences and Freespa 
provides some additional features. Below are the key differences:

### Δt Handling

- **Freespa**: Provides an interface to determine Δt values from 
internal tables and models (see [Δt Values](#δt-values)).
- **NREL SPA**: Limits Δt to \(\pm 8000\) seconds [NRELSPA], even 
though this range only holds for years 245–3400 [Morrison].


### Atmospheric Refraction

1. **Refraction Models**:
   - **Freespa**: Includes both forward and inverse models:
     - **Forward Model**: Computes apparent solar position from the 
     true position.
     - **Inverse Model**: Computes true position given an apparent 
     position.
   - **NREL SPA**: Requires the user to provide a refraction angle at 
   sunrise/set, which may lead to inconsistencies.

2. **Bennet Refraction Model**:
   - **Freespa**: Offers two parametrizations:
     - Original parametrization [Meeus] (same as NREL SPA).
     - Updated parametrization [Wilsen], with better agreement to "The 
     Nautical Almanac" (2004 and later editions).
   - **NREL SPA**: Uses the original Bennet model.

3. **Geometric Dip**:
   - **Freespa**: Accounts for geometric dip (e.g., the horizon below 
   the observer) based on elevation or a user-provided value.
   - **NREL SPA**: Assumes a fixed horizon at 0° elevation.


### Solar Event Computation

1. **Supported Events**:
   - **Freespa**: Computes sunrise, transit, sunset, and various 
   dawn/dusk times (civil, nautical, and astronomical).
   - **NREL SPA**: Computes only sunrise, transit, and sunset.

2. **Accuracy**:
   - **Freespa**: Calculates sunrise/set times with an accuracy of 
   **0.05° or 1 second**.
   - **NREL SPA**: May have errors up to **0.1° in solar elevation**.

3. **Special Cases**:
   - **Freespa**: Handles cases where one or more events (e.g., sunset) 
   does not occur, such as during polar day/night.
   - **NREL SPA**: Does not compute events if one is missing.

### Equinox and Solstice
**Freespa**: provides routines to compute Equinoxes and Solstices.


## References
[Meeus] J. Meeus, Astronomical Algorithms, second ed. Willmann-Bell, Inc., Richmond, Virginia, USA. (1998): 105-108

[Wilsen] T. Wilson, "Evaluating the Effectiveness of Current Atmospheric Refraction Models in Predicting Sunrise and Sunset Times", Open Access Dissertation, Michigan Technological University, (2018). https://doi.org/10.37099/mtu.dc.etdr/697

[Morrison] L. V. Morrison, and F. R. Stephenson, "Historical Values of the Earth’s Clock Error ΔT and the Calculation of Eclipses." Journal for the History of Astronomy, 35(3), (2004): 327–336. https://doi.org/10.1177/002182860403500305

[NRELSPA] A. Andreas, spa.c source code retrieved from https://midcdmz.nrel.gov/spa/ on the 27th of march 2023

[Reda] I. Reda and A. Andreas, "Solar position algorithm for solar radiation applications." Solar Energy 76.5 (2004): 577-589
