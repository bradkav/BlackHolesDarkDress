for NB in $(ls -r *.ipynb);
do
    echo "Clearing output from notebook: $NB"
    
    nbstripout $NB
done

for NB in $(ls -r Nbody/*.ipynb);
do
    echo "Clearing output from notebook: $NB"

    nbstripout $NB
done