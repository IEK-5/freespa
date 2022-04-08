# freespa
## About freespa
This is a free implementation in c of the Solar Position Algorithm (SPA) [1-2]. Another open and free c implementation is provided by NREL [here](http://rredc.nrel.gov/solar/codesandalgorithms/spa/). NREL's SPA is free as in beer. However, it's beer which you are not allowed to share with your friends. You are free to share freespa under the conditions of the [GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Installation
Package is not designed to be installed but rather to be inserted in a project. You can add freespa as a submodule to your code, the provided Makefile is designed for this purpose. 

For debugging purposes I also provide code to comapre freespa with NREL SPA, however, to compile that you need to obtain your own copy of [NREL SPA](http://rredc.nrel.gov/solar/codesandalgorithms/spa/). A test program is also provided which should generally suffice for testing the correctness of the algorithm (i.e. you do not need NREL SPA to verify freespa). To run the tests do 

`make check`

This will build a test program, "testspa", and checks the correct functioning with it. Tests include a comparison with a reference data set with time/location/solar position. The testspa program usage:

`testpa [r/t] \<filename\>`

In general one should only use the "t" option with the provided reference file. You may, however use it to generate a new reference file with the "r" option. 

A reference file contains lines whiche each line whould consist of:

`unix-timestamp latitude longitude true-azimuth true-zenith aparent-azimuth aparent-zenith`,

where all angles and coordinates are in radians.

## Delta t
For accurate timing the SPA algorithm needs delta t values. Both historic and predicted future values may be obtained here:
[https://maia.usno.navy.mil/ser7/](https://maia.usno.navy.mil/ser7/)

Another source of delta t values is

[https://hpiers.obspm.fr/eop-pc/index.php?index=analysis&lang=en](https://hpiers.obspm.fr/eop-pc/index.php?index=analysis&lang=en)

This database provides delta t values over various time ranges in various resolutions.

Using these datasources we compiled a table of delta t tables in the file historic_delta_t.dat. The script "GetDeltaT.sh" takes this data and combines it with updated measured and predicted data from

[https://maia.usno.navy.mil/ser7/](https://maia.usno.navy.mil/ser7/)

The result is collected in freespa_dt_table_new.h. This new header file may be used to replace the freespa_dt_table.h header.

[1] I.  Reda and A. Andreas, "Solar position algorithm for solar radiation applications." Solar Energy 76.5 (2004): 577-589

[2] I. Reda and A. Andreas, "Corrigendum to Solar position algorithm for solar radiation applications." Solar Energy 81.6 (2007): 838

