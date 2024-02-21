#!/usr/bin/env bash

#echo "best_acc averageInteractionsPerCrossing p_Rhad p_Rhad1 p_TRTTrackOccupancy p_topoetcone40 p_eTileGap3Cluster p_phiModCalo p_etaModCalo Number_of_settings Elapsed_time task_time_[mus]" > HEP_out_5
sbatch --ntasks=2  job.sh
sbatch --ntasks=4  job.sh
#sbatch --ntasks=8  job.sh
#sbatch --ntasks=16 job.sh
#sbatch --ntasks=32 job.sh
#sbatch --ntasks=64 job.sh
#sbatch --nodes=2 --ntasks=128 job.sh
