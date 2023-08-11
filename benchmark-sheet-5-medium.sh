#!/bin/bash

run_ce_solver()
{
	PROGRAMM_NAME="$1"
	LOG=$2
	CSV=$3
	maxSec=$4
	maxSecPerInstance=$5
	maxNotSolved=$6
	FILES=$(ls $7*.dimacs)
	KWARG1=$8
	KWARG2=$9
	KWARG3=${10}
	KWARG4=${11}
	KWARG5=${12}


	rm -f time.txt

	overallTime=$(date +%s);
	now=$(date +%s);
	elapsed=`expr $now - $overallTime`;
	notSolved=0
	showedHint=0

	maxSecPerInstanceHard=$((maxSecPerInstance + 1));

	for f in $FILES
	do
		if [ $elapsed -le $maxSec -a $notSolved -le $maxNotSolved ]; then
			echo $f >> $LOG
			
			# start everything in a new process group such that we can kill everything if necessary
# 			(setsid /usr/bin/time -f "%e" -a -o time.txt timeout -k $maxSecPerInstance -s 9 $maxSecPerInstance $PROGRAMM_NAME< $f 1> prog_out.txt 2>&1) & PID=$!
			(setsid /usr/bin/time -f "%e" -a -o time.txt timeout -k $maxSecPerInstanceHard -s 2 $maxSecPerInstance $PROGRAMM_NAME $KWARG1 $KWARG2 $KWARG3 $KWARG4 $KWARG5 < $f 1> prog_out.txt 2>&1) & PID=$!


			# kill processes when exiting this script
			trap "{ kill -$PID 2>/dev/null; }" TERM
			trap "{ kill -9 -$PID 2>/dev/null; }" EXIT

			wait $PID

			# just to be sure: if the process still is around brutally kill it
			kill -0 $PID 2>/dev/null || kill -9 -$PID 2>/dev/null;

			# get n
			n=$(head -1 $f | sed 's/#//' | sed 's/ .*//')

			# get m
			m=$(head -1 $f | sed 's/#.* //')
			
			# get k
			k=$(grep -ve "^#" prog_out.txt | wc -l)
			lb1=$(grep -e "#lb1:" prog_out.txt | sed -e 's/.*lb1: \([0-9]*\).*/\1/' )
			lb2=$(grep -e "#lb2:" prog_out.txt | sed -e 's/.*lb2: \([0-9]*\).*/\1/' )
			lb3=$(grep -e "#lb3:" prog_out.txt | sed -e 's/.*lb3: \([0-9]*\).*/\1/' )
			recursiveSteps=$(grep -e "#recursive steps:" prog_out.txt | sed -e 's/.*recursive steps: \([0-9]*\).*/\1/' )
			lastK=$(grep -e "last-k:" prog_out.txt | sed -e 's/.*last-k: \([0-9]*\).*/\1/' )
			cat prog_out.txt >> $LOG
			
			# get time
			time=$(cat time.txt);
			
			if [[ $time == "Command terminated by signal 9"* ]] || [[ $time == "Command exited with non-zero status"* ]]; then
				finished=0;
				(( notSolved += 1 ));
				time="";
			else
				finished=1;
			fi
			
			verify="";
			
			if [ "$finished" -eq "1" ]; then
				solFile=$(basename $f .dimacs)
				solNumber=$(cat $data$solFile.solution);
				if [ -n "$solNumber" ] && [ "$solNumber" -eq "$solNumber" ] 2>/dev/null; then
					if [ "$solNumber" -eq "$k" ]; then
						verify="correct\t0";
					else
						verify=">>INCORRECT Size<< \t$(($k-$solNumber))";
					fi
				fi
				msg=$(./verifier.py $f prog_out.txt)
				if [ $? != 0 ]; then
					verify="$verifier.py:  \t $msg";

				fi
			fi
			rm -f prog_out.txt
			
			fileNameLong=$(printf '%-40s' "$f");
			
			if [ -n "$solNumber" ] && [ -n "$lastK" ] 2>/dev/null; then
				lastK=$(($solNumber-$lastK));
			fi
			solNumber="";
			
			
			echo -e "$fileNameLong\t"$time"\t"$k"\t"$recursiveSteps"\t"$finished"\t"$verify"\t"$lb1"\t"$lb2"\t"$lb3"\t"$lastK"\t"$solNumber
			echo -e $f"\t"$time"\t"$n"\t"$m"\t"$k"\t"$recursiveSteps"\t"$finished"\t"$verify"\t"$lb1"\t"$lb2"\t"$lb3"\t"$lastK"\t"$solNumber | sed 's/	/;/g' >> $CSV
			echo "" >> $LOG
			
			rm -f time.txt

			now=$(date +%s);
			elapsed=`expr $now - $overallTime`;
		else
			if [ $showedHint -eq 0 ]; then
				if [ $notSolved -ge $maxNotSolved ]; then
					echo "$notSolved instances not solved. Script aborted."
				else
					echo "maximal time of $maxSec sec elapsed. Script aborted."
				fi
				showedHint=1;
			fi		
		fi
	done
}

#PROGRAMM_NAME="./path/to/your/solver"  		# insert your program here
#PROGRAMM_NAME="java -jar -Xss512m VertexCoverSolver.jar"  		# insert your program here
PROGRAMM_NAME="python3 TUB_AlgEng_2020_Team_1/sheet_5/branch_and_bound.py"  		# insert your program here
today=$(date +%Y-%m-%d-%H-%M-%S)
#LOG="log-$today.txt"							# specify the name of the log file
LOG="log.txt"								# specify the name of the log file
maxSec=43200								# overall allowed time for the whole script
maxSecPerInstance=300							# allowed time (in seconds) for one instance
maxNotSolved=4								# no of instances the program is allowed to fail to solve. If reached, then the script is aborted
## now loop through data set directories
mkdir -p results

for data in [3]*/ ; do
	CSV="results/results-$1-$2-$3-$4-$5-$today-${data%/}.csv"
	echo "file;time;vertices;edges;solsize;recsteps;finished;verified;lb1;lb2;lb3" > $CSV
	echo "run $data instances $PROGRAMM_NAME (Tab-separated columns: File, Time in seconds, solution size, recursive steps, finished, solution size verified, lb1, lb2, lb3)"
	run_ce_solver "$PROGRAMM_NAME" $LOG $CSV $maxSec $maxSecPerInstance $maxNotSolved $data $1 $2 $3 $4 $5
done

echo ""
