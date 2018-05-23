mkdir sims

wget -O sims/M30_a0.01_lo_bin_e0.995_full.tar.gz https://ndownloader.figshare.com/files/11577773

wget -O sims/README.md https://ndownloader.figshare.com/files/11590262

cd sims
tar -xvzf M30_a0.01_lo_bin_e0.995_full.tar.gz
rm M30_a0.01_lo_bin_e0.995_full.tar.gz
cd ..