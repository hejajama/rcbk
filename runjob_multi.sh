LC_ALL="C"
for qsqr in 0.1580 0.1590 0.1600 0.1610  #0.1605 0.1610 0.1615 0.162 0.1625 0.1630 0.1633 0.1634 0.1635 0.1636 0.1637 ## 0.167 0.168 0.169
do
	for csqr in 6.7000 6.8000 6.9000 7.0000 7.1000 #5.8 5.9 6.0 6.1 6.2 6.3 6.4 6.5  #5 5.25 5.5 5.75 6  #5.5 6 6.5 7 #0.3152 
	do
		for gamma in 1.1100 1.1150 1.1200 1.1250 1.1300 1.1350 1.1400 1.1450 #1.1 1.125 1.15 1.175 #1.135 # 1.1 1.2 1.3 1.4 1.5  #1.05 1.1 1.15 # 1.1 1.2
		do
			for ec in 1 #5.65 5.7 5.75 #5 10 1.175#8 9 10 11 12 13
			do
				for freeze in 0 
				do
					data="fit/aamqs/mv_qsqr_${qsqr}_gamma_${gamma}_csqr_${csqr}_ec_${ec}_freeze_${freeze}.dat"
					if [ ! -f $data ]
					then
						touch ${data}
						out="fit/stdouterr/mv_qsqr_${qsqr}_gamma_${gamma}_csqr_${csqr}_ec_${ec}_freeze_${freeze}.stdout"
						err="fit/stdouterr/mv_qsqr_${qsqr}_gamma_${gamma}_csqr_${csqr}_ec_${ec}_freeze_${freeze}.stderr"
						cmd="-J rbk_qsqr_${qsqr}_gamma_${gamma}_csqr_${csqr}_ec_${ec}_freeze_${freeze} -o ${out} -e ${err} ./runjob_param.sh ${qsqr} ${gamma} 0.01 ${csqr} ${ec} ${freeze} ${data}"
						echo ${cmd}
						sbatch $cmd
					fi
				done
			done
		done
	done
done

