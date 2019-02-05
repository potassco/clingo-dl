#!/bin/bash
DEFAULT_CONF='encodings/delay_linear.lp encodings/h1.lp encodings/h2.lp encodings/h3.lp -t2 -q1,0 --stats --heuristic=Domain'
CONF=${2:-$DEFAULT_CONF} 	# If variable not set, use default.
DEFAULT_TIME_LIMIT=500
TIME_LIMIT=${3:-$DEFAULT_TIME_LIMIT} 	# If variable not set, use default.
DEFAULT_MBW=900
MBW=${4:-$DEFAULT_MBW}  				# If variable not set, use default.
encodingfolder=encodings
jsonfolder=instances/json

echo "starting solver in Encoding-Folder $encodingfolder with a Time-Limit of $TIME_LIMIT s and a maximum delay of $MBW s"

filename=$(basename -- "$1")
filename="${filename%.*}"

jsonname="${jsonfolder}/${filename}.json"
clingo-dl encodings/encoding.lp encodings/preprocessing.lp encodings/minimize_routes.lp encodings/minimize_delay.lp ${CONF} encodings/output.lp --time-limit=${TIME_LIMIT} -c mbw=${MBW} $1 > sol
./utils/asp2json.py sol > sol.json
python3.6 utils/validate_solution.py ${jsonname} sol.json
rm sol sol.json
