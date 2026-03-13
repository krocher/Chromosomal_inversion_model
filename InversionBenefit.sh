val=( '0.01' '0.02' '0.05' '0.1' )

for S_INV in ${val[@]}
do
    for i in {1..5}
    do 
        slim -d Rep=${i} -d S_INV=${S_INV} Inversion_model.slim &
    done
done