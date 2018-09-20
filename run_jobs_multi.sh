#!/bin/bash
# Number of jobs in each block
nb=10
seq1=1
seq2=1
seq3=50

# Create a 'wait' sentence between chunks
waittext="wait "
for i in $(seq 1 $nb)
do
	waittext=("${waittext[*]}" "%$i")
done

# Create job_submit file
cat > job_submit.sh <<EOF 
#!/bin/bash

EOF

n=$((seq1*seq2*seq3))
cont=0
for i in $(seq 1 $seq1)
do
	for j in $(seq 1 $seq2)
	do
		for k in $(seq 1 $seq3)
		do
			cont=$(($cont +1))
			LOG=${i}_${j}_${k}
			
cat >> job_submit.sh <<EOF 
R CMD BATCH --no-save --no-restore '--args CV=$i mod=$j part=$k' fitModels_multi.R LOG_$LOG &
EOF
            
			if [ $(($cont % $nb)) == 0 ] && [ "$cont" -lt "$n" ]
			then
cat >> job_submit.sh <<EOF 

${waittext[@]}

EOF
			fi
        done
    done
done

sh job_submit.sh &
