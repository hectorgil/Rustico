#!/bin/bash
#PBS -N RusticoX
#PBS -l mem=22gb
#PBS -l nodes=1:ppn=16
#PBS -q batch
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=16

#time ./fileNseries.out params3.c
#time ./fileNseries.out params33.c
#time ./fileNseries.out params3b.c
#time ./fileNseries.out params6.c
#time ./fileNseries.out params6b.c
#time ./fileNseries.out params6c.c
#time ./fileNseries_randoms.out params66c.c

#time ./fileNseries.out params6new.c
#time ./fileNseries.out params6new_B.c
#time ./fileNseries.out params6new_C.c
#time ./fileNseries.out params6new_D.c
#time ./fileNseries.out params6new_E.c

#time ./fileNseries.out params7.c
#time ./fileNseries.out params7b.c
#time ./fileNseries.out params7c.c
#time ./fileNseries.out params7d.c
#time ./fileNseries.out params7e.c
#time ./fileNseries.out params15.c
#time ./fileNseries.out params15b.c
#time ./fileNseries.out params15c.c

#time ./fileNseries.out params03.c
#time ./fileNseries.out params03b.c
#time ./fileNseries.out params03c.c
#time ./fileNseries.out params03d.c
#time ./fileNseries.out params03e.c


#time ./fileNseries.out params08.c
#time ./fileNseries.out params08b.c
#time ./fileNseries.out params08c.c
#time ./fileNseries.out params08d.c
time ./fileNseries.out params08e.c
