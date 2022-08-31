#!/bin/bash 
#
#SBATCH --output=/home/Youtong.Zheng/mypython/sslurm-%j.out
#SBATCH -J average_transact_daily_calipso
#SBATCH -n 1
#SBATCH --time=20
#SBATCH --mail-type=ALL
#SBATCH --account=gfdl_m

source /work/Youtong.Zheng/anaconda2/etc/profile.d/conda.sh
conda activate demo

# python code/average_transact_daily_calipso.py 'calipso' '2006-2014' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
# python code/average_transact_daily_calipso.py 'calipso' '2006-2014' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
# python code/average_transact_daily_calipso.py 'calipso' '2006-2014' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'
# python code/average_transact_daily_calipso.py 'calipso' '2006-2014' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/calipso/'

# python code/average_transact_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/average_transact_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/average_transact_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/average_transact_daily.py 'c96L33_am4p0_qadt_cosp_2010' '20100101-20141231' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'

#python code/average_transact_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20060612-20141231' 'nh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
python code/average_transact_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20060612-20141231' 'sh' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/average_transact_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20060612-20141231' 'nh-chukchi' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
# python code/average_transact_daily_cmip.py 'c96L33_am4p0_qadt_cosp_2010' '20060612-20141231' 'sh-ross' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4_cmip/' '/work/Youtong.Zheng/Data/Zheng_Ming_JCli_2022/am4/'
