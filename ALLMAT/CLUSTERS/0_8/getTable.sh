thrs=${1}
f1="Fitness-DNA-unit-${thrs}-sum-PN.dat"
f2="Complexity-RNA-unit-${thrs}-sum-PN.dat"
F="top10-table.tex"
printf "Top 10\n"
rm $F
HEADER11=$(head -1 ${f1} |tail -1| awk -F'\t' {'print $1'})
HEADER22=$(head -1 ${f2} |tail -1| awk -F'\t' {'print $1'})
HEADER=$(printf "RANK & ${HEADER11} & ${HEADER22}")
printf "${HEADER}" >> $F; printf "\\" >> $F; printf "\\" >> $F; printf "\n">>$F
for i in {2..11}
do
	ii=$( echo ${i} -1|bc)
	s11=$(head -${i} ${f1} |tail -1| awk -F'\t' {'print $1'})
	s12=$(head -${i} ${f2} |tail -1| awk -F'\t' {'print $1'})
	s1=$(printf "${ii} & ${s11} & ${s12}")
	printf "${s1}" >> $F; printf "\\" >> $F; printf "\\" >> $F; printf "\n">>$F
done
printf "Last 10s\n"
F="Last10-table.tex"
rm ${F}
HEADER11=$(head -1 ${f1} |tail -1| awk -F'\t' {'print $1'})
HEADER22=$(head -1 ${f2} |tail -1| awk -F'\t' {'print $1'})
HEADER=$(printf "RANK & ${HEADER11} & ${HEADER22}")
printf "${HEADER}" >> $F; printf "\\" >> $F; printf "\\" >> $F; printf "\n">>$F
for i in {1..10}
do
        s11=$(tail -${i} ${f1}| head -1| awk -F'\t' {'print $1'})
        s12=$(tail -${i} ${f2}| head -1| awk -F'\t' {'print $1'})
        s1=$(printf "${i} & ${s11} & ${s12}")
        printf "${s1}" >> $F; printf "\\" >> $F; printf "\\" >> $F; printf "\n">>$F
done
