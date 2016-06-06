#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -V
#$ -l h_rt=240:0:0    # Request 240 hour runtime
#$ -l h_vmem=50G      # Request 4GB RAM

mysqlDir=/data/home/btw796/Prog/MySQL
cd $mysqlDir

#./scripts/mysql_install_db --user=mysql --datadir=/data/SBCS-BessantLab/shyama/Data/MySQL/data3/
./bin/mysqld_safe --defaults-file=/data/home/btw796/Prog/MySQL/my-new.cnf &

run=1
while [ $run == 1 ]
do
## nothing to do
#echo ""
run=1
done

