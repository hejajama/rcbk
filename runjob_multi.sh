LC_ALL="C"
ic="mv1"
for qsqr in `seq -f %0.4f 0.0500 0.0100 0.3` # 0.3000` 
do
	for csqr in `seq -f %07.4f 2.0000 2 40.0000`
	do
		for gamma in "01.0000"
		do
			for ec in "01.0000" #`seq -f %07.4f 1.0000 4 30` # 50.0000` #5.65 5.7 5.75 #5 10 1.175#8 9 10 11 12 13
			do
				for freeze in 0 
				do
					data="/work/hejajama/herafit_systemaattisesti/${ic}/mv_qsqr_${qsqr}_gamma_${gamma}_csqr_${csqr}_ec_${ec}_freeze_${freeze}.dat"
					if [ ! -f $data ]
					then
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

