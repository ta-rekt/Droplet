#!/bin/bash

DATE=20260624
EXP=proximal_promoter_transc_noen_long
D=30

for c in 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010 0.011 0.012 0.013 0.014; do
    DIR="${EXP}_d_${D}_c_${c}"
    echo "Pulling $DIR ..."
    mkdir -p "data/${DATE}/${DIR}"
    scp -J tohme@bastion.phys.ens.fr \
        "tohme@thanos.phys.ens.fr:/home/tohme/.julia/dev/Droplet/data/${DATE}/${DIR}/*" \
        "data/${DATE}/${DIR}"
done
