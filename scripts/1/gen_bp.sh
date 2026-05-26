#!/bin/bash
#SBATCH --job-name=computeBP
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=10G

# load modules necessary to run cdftools (same as during compilation)

module purge
module load oneapi2023-env/2023.2.0
module load oneapi/2023.2.0
module load netcdf-c/4.9.2-with-oneapi-mpi-2021.10.0
module load netcdf-cxx4/4.3.1-with-oneapi-mpi-2021.10.0
module load netcdf-fortran/4.6.0-with-oneapi-mpi-2021.10.0
module load hdf5/1.14.1-2-with-with-oneapi-mpi-2021.10.0

# give name of experiment (here without 1_)
export EX=VIKING20X.L46-KFS003

# directory where mesh_mask file is located
export MDIR=/gxfs_work/geomar/smomw628/${EX}_11S
# output directory - adjust accordingly!
export ODIR=/gxfs_work/geomar/smomw628/${EX}_11S

# Path to your cdftools installation
export CDFTOOLS=$WORK/CDFTOOLS/bin


# provide mesh_mask file as hgr, zgh and mask variable - adjust according to your region
export CDFT_MESH_HGR=${MDIR}/1_${EX}_mesh_mask_11S_box.nc
export CDFT_MESH_ZGR=${MDIR}/1_${EX}_mesh_mask_11S_box.nc
export CDFT_MASK=${MDIR}/1_${EX}_mesh_mask_11S_box.nc


# adjust to match directory of input data if used for regional subset
export DDIR=/gxfs_work/geomar/smomw628/${EX}_11S

# start and end years for files to be processed
YY1=1980
YY2=2023

for jj in `seq $YY1 $YY2`
do


ii=${DDIR}/1_${EX}_1d_${jj}0101_${jj}1231_grid_T_11S_box.nc  # adjust according to your file names / regions
IN_T=$ii
# using detrended SSH
ij=${DDIR}/1_${EX}_1d_${jj}0101_${jj}1231_detrendedSSH_11S_box.nc
IN_SSH=$ij
OUTssh=`basename $IN_T grid_T.nc`botpressure_detrssh.nc   # output filename
${CDFTOOLS}/cdfbotpressure -t $IN_T --ssh-file $IN_SSH -ssh2 -xtra -o ${ODIR}/$OUTssh & # actual command to compute bottom pressure

done
wait
