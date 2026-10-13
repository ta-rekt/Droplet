#!/bin/bash

DATE=20261013
EXP=run2/high_c/s

for S in 50 100 150 200 250 300 350 400 450 500 550 600; do
    DIR="${EXP}${S}"
    echo "Pulling $DIR ..."
    mkdir -p "data/${DATE}/${DIR}"
    scp -r "tohme@thanos.phys.ens.fr:/home/tohme/.julia/dev/Droplet/data/${DATE}/${DIR}/*" \
        "data/${DATE}/${DIR}"
done
