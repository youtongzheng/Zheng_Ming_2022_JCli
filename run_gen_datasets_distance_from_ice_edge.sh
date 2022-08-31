#!/bin/bash 
#
#SBATCH --output=/home/Youtong.Zheng/mypython/sslurm-%j.out
#SBATCH -J gen_datasets_distance_from_ice_edge_daily
#SBATCH -n 1
#SBATCH --time=20
#SBATCH --mail-type=ALL
#SBATCH --account=gfdl_m

source /work/Youtong.Zheng/anaconda2/etc/profile.d/conda.sh
conda activate demo

##-----------AM4---------------

# python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20050101-20091231' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20050101-20091231' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20050101-20091231' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20050101-20091231' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'


#python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
#python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
#python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
#python code/Gen_datasets_distance_from_ice_edge_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
#python code/Gen_datasets_distance_from_ice_edge_monthly.py 'c96L33_am4p0_qadt_cosp_2010' '201001-201412' 'nh' '/work/Youtong.Zheng/transact/'

##-----------AM4 cmip---------------

#python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20030101-20081231' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'
#python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20050101-20091231' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'
python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20030101-20081231' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'
#python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20050101-20091231' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'

#python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20090101-20141231' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'
# python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'
python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20090101-20141231' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'
# python code/Gen_datasets_distance_from_ice_edge_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/'

##-----------callipso---------------
#for ((year=2006; year<=2014; year++))
#do
# python code/Gen_datasets_distance_from_ice_edge_daily_calipso.py ${year} 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
# python code/Gen_datasets_distance_from_ice_edge_daily_calipso.py ${year} 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
# python code/Gen_datasets_distance_from_ice_edge_daily_calipso.py ${year} 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
# python code/Gen_datasets_distance_from_ice_edge_daily_calipso.py ${year} 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
#done
