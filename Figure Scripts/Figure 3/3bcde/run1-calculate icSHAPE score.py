DB=/Share/home/zhangqf8/lipan/DYNAMIC/shape_score

##	Step 1
normalizeRTfile=~/lipan/icSHAPE/icSHAPE/scripts/normalizeRTfile.pl

bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_ch_vivo/rt_D1.txt -o ch_D1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_ch_vivo/rt_D2.txt -o ch_D2.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_ch_vivo/rt_N1.txt -o ch_N1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_ch_vivo/rt_N2.txt -o ch_N2.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_ch_vitro/rt_N1.txt -o ch_T1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_ch_vitro/rt_N2.txt -o ch_T2.rt -m mean:vigintile2 -d 32 -l 32 -f 100

bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_np_vivo/rt_D1.txt -o np_D1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_np_vivo/rt_D2.txt -o np_D2.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_np_vivo/rt_N1.txt -o np_N1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_np_vivo/rt_N2.txt -o np_N2.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_np_vitro/rt_N1.txt -o np_T1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_np_vitro/rt_N2.txt -o np_T2.rt -m mean:vigintile2 -d 32 -l 32 -f 100

bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_cy_vivo/rt_D1.txt -o cy_D1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_cy_vivo/rt_D2.txt -o cy_D2.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_cy_vivo/rt_N1.txt -o cy_N1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_cy_vivo/rt_N2.txt -o cy_N2.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_cy_vitro/rt_N1.txt -o cy_T1.rt -m mean:vigintile2 -d 32 -l 32 -f 100
bsub -n 5 -q Z-ZQF -e error -o log perl $normalizeRTfile -i $DB/hek_cy_vitro/rt_N2.txt -o cy_T2.rt -m mean:vigintile2 -d 32 -l 32 -f 100

##	Step 2

calcEnrich=~/lipan/icSHAPE/icSHAPE/scripts/calcEnrich.pl

# + ch
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_N1.rt -b ch_D1.rt -o ch_SHAPE1.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_N1.rt -b ch_D2.rt -o ch_SHAPE2.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_N2.rt -b ch_D1.rt -o ch_SHAPE3.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_N2.rt -b ch_D2.rt -o ch_SHAPE4.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex

bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_T1.rt -b ch_D1.rt -o ch_SHAPET1.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_T1.rt -b ch_D2.rt -o ch_SHAPET2.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_T2.rt -b ch_D1.rt -o ch_SHAPET3.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f ch_T2.rt -b ch_D2.rt -o ch_SHAPET4.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex

# + np
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_N1.rt -b np_D1.rt -o np_SHAPE1.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_N1.rt -b np_D2.rt -o np_SHAPE2.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_N2.rt -b np_D1.rt -o np_SHAPE3.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_N2.rt -b np_D2.rt -o np_SHAPE4.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex

bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_T1.rt -b np_D1.rt -o np_SHAPET1.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_T1.rt -b np_D2.rt -o np_SHAPET2.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_T2.rt -b np_D1.rt -o np_SHAPET3.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f np_T2.rt -b np_D2.rt -o np_SHAPET4.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex

# + cy
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_N1.rt -b cy_D1.rt -o cy_SHAPE1.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_N1.rt -b cy_D2.rt -o cy_SHAPE2.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_N2.rt -b cy_D1.rt -o cy_SHAPE3.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_N2.rt -b cy_D2.rt -o cy_SHAPE4.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex

bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_T1.rt -b cy_D1.rt -o cy_SHAPET1.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_T1.rt -b cy_D2.rt -o cy_SHAPET2.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_T2.rt -b cy_D1.rt -o cy_SHAPET3.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex
bsub -q Z-ZQF -e error -o log -n 5 perl $calcEnrich -f cy_T2.rt -b cy_D2.rt -o cy_SHAPET4.rt -w factor5:scaling1 -y 10 -x 0.25 -e complex

##	Step 3

T=1
t=100
mkdir t${t}T${T}

filterEnrich=~/lipan/icSHAPE/icSHAPE/scripts/filterEnrich.pl

# + ch
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPE1.rt -o t${t}T${T}/ch1.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPE2.rt -o t${t}T${T}/ch2.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPE3.rt -o t${t}T${T}/ch3.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPE4.rt -o t${t}T${T}/ch4.out -t $t -T $T -s 5 -e 30 &

bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPET1.rt -o t${t}T${T}/chT1.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPET2.rt -o t${t}T${T}/chT2.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPET3.rt -o t${t}T${T}/chT3.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i ch_SHAPET4.rt -o t${t}T${T}/chT4.out -t $t -T $T -s 5 -e 30 &

# + np
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPE1.rt -o t${t}T${T}/np1.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPE2.rt -o t${t}T${T}/np2.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPE3.rt -o t${t}T${T}/np3.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPE4.rt -o t${t}T${T}/np4.out -t $t -T $T -s 5 -e 30 &

bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPET1.rt -o t${t}T${T}/npT1.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPET2.rt -o t${t}T${T}/npT2.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPET3.rt -o t${t}T${T}/npT3.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i np_SHAPET4.rt -o t${t}T${T}/npT4.out -t $t -T $T -s 5 -e 30 &

# + cy
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPE1.rt -o t${t}T${T}/cy1.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPE2.rt -o t${t}T${T}/cy2.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPE3.rt -o t${t}T${T}/cy3.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPE4.rt -o t${t}T${T}/cy4.out -t $t -T $T -s 5 -e 30 &

bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPET1.rt -o t${t}T${T}/cyT1.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPET2.rt -o t${t}T${T}/cyT2.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPET3.rt -o t${t}T${T}/cyT3.out -t $t -T $T -s 5 -e 30 &
bsub -q Z-ZQF -e error -o log perl $filterEnrich -i cy_SHAPET4.rt -o t${t}T${T}/cyT4.out -t $t -T $T -s 5 -e 30 &







