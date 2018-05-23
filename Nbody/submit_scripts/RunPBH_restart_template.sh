#PBS -lnodes=1:cores16:ppn=16
#PBS -lwalltime=25:00:00

module load mpicopy

cd $HOME/PBH/sims/JOB_ID/run
sed -i -- "s@RUN_DIRECTORY@$TMPDIR/JOB_ID/run@g" PBH.param
sed -i -- "s@OUT_DIRECTORY@$TMPDIR/JOB_ID/out@g" PBH.param

mkdir "$TMPDIR"/JOB_ID
mkdir "$TMPDIR"/JOB_ID/out

mpicopy $HOME/PBH/sims/JOB_ID/run -o "$TMPDIR"/JOB_ID
#cp $HOME/PBH/sims/JOB_ID/out/restart* "$TMPDIR"/JOB_ID/out

mpicopy $HOME/PBH/sims/JOB_ID/out/snapshot_SNAP_ID -o "$TMPDIR"/JOB_ID/out

module load openmpi/gnu

cd "$TMPDIR"/JOB_ID/run

time mpiexec Gadget2 PBH.param 2 > output

tail -n1000 output > output_tail_SNAP_ID
head -n1000 output > output_head_SNAP_ID

tail -n1000 "$TMPDIR"/JOB_ID/out/cpu_SNAP_ID.txt > "$TMPDIR"/JOB_ID/out/cpu_tail_SNAP_ID.txt


rm "$TMPDIR"/JOB_ID/run/output
rm "$TMPDIR"/JOB_ID/out/timings.txt
rm "$TMPDIR"/JOB_ID/out/info.txt
rm "$TMPDIR"/JOB_ID/out/cpu.txt

cp -rf "$TMPDIR"/JOB_ID $HOME/PBH/sims