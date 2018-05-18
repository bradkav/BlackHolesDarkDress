#PBS -lnodes=1:cores16:ppn=16
#PBS -lwalltime=48:00:00

module load mpicopy
module load openmpi

cd $HOME/PBH/sims/M30_a0.01_lo_bin_e0.995_full/run
sed -i -- "s@RUN_DIRECTORY@$TMPDIR/M30_a0.01_lo_bin_e0.995_full/run@g" PBH.param
sed -i -- "s@OUT_DIRECTORY@$TMPDIR/M30_a0.01_lo_bin_e0.995_full/out@g" PBH.param

mpicopy $HOME/PBH/sims/M30_a0.01_lo_bin_e0.995_full -o "$TMPDIR"

module load openmpi/gnu

cd "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/run

time mpiexec Gadget2 PBH.param > output

tail -n1000 output > output_tail
head -n1000 output > output_head

tail -n1000 "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/out/cpu.txt > "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/out/cpu_tail.txt

rm "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/run/output
rm "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/out/timings.txt
rm "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/out/info.txt
rm "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full/out/cpu.txt

cp -rf "$TMPDIR"/M30_a0.01_lo_bin_e0.995_full $HOME/PBH/sims