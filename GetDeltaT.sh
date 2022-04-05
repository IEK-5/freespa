#!/bin/bash
# Simple script that crates a delta t table for freespa
# The delta t table is compiled of various sources
# The script assumes historic values of delta t do not change
# It stores historic delta t values in a reference file
# If this script is run it downloads:
# observed delta t valuyes from
# https://maia.usno.navy.mil/ser7/deltat.data
# long term delta t predictions from
# https://maia.usno.navy.mil/ser7/deltat.preds
#
# The last entries in the deltat.data file are added to our 
# historic reference
# the predictions are adde3d beyin the last entry in our reference
#
# the reference is compiled from various IERS sources


function deltatdata() {
	year=$1
	file=$2
	# columns: year month ?? delta_t
	awk "/^\s*[0-9]/{y=\$1+\$2/12+\$3/30;if(y>$year){print y,\$4}}" $file
}

function deltatpreds() {
	year=$1
	file=$2
	# columns: year(decimal) delta_t acc
	awk "/^\s*[0-9]/{if($year<\$1){print \$1,\$2}}" $file
}


# reference file:
ref="historic_delta_t.dat"
urlcurrent="https://maia.usno.navy.mil/ser7/deltat.data"
urlfuture="https://maia.usno.navy.mil/ser7/deltat.preds"

last=$(tail -n 1 $ref |awk '{print $1}')
echo "File $ref contains historic delta t data up to $last"

wget -o /dev/null $urlcurrent
status=$? 
if test $status -eq 0
then
	echo "Successfully downloaded delta t data"
else
	echo "could not update delta t data"
	exit 1
fi
mv deltat.data current_delta_t.dat
wget  -o /dev/null $urlfuture
status=$? 
if test $status -eq 0
then
	echo "Successfully downloaded delta t data"
else
	echo "could not update delta t data"
	exit 1
fi
mv deltat.preds future_delta_t.dat
cat $ref > new_hist_delta_t.dat
deltatdata $last "current_delta_t.dat" >> new_hist_delta_t.dat
# after manual inspection you may rebase the reference file to this one
prev=$last 
last=$(tail -n 1 new_hist_delta_t.dat |awk '{print $1}')
echo "Added new observed data from $prev up to $last"

cat new_hist_delta_t.dat >tmp.dat
deltatpreds $last future_delta_t.dat >> tmp.dat
deltatpreds $last future_delta_t.dat

N=$(wc -l tmp.dat|egrep -o [0-9]+)
echo $N
echo "/* delta t table from iers data */" >freespa_dt_table_new.h
echo "#define NDT $N" >>freespa_dt_table_new.h
echo "const double freespa_delta_t_table[2*NDT] = {" >>freespa_dt_table_new.h
awk '/^/{printf "\t%s,\n\t%s,\n",$1,$2}' tmp.dat >> freespa_dt_table_new.h
echo "};" >>freespa_dt_table_new.h
# after manual inspection move freespa_dt_table_new.h to freespa_dt_table.h

