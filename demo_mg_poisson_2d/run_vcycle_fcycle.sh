#An example to run and compare between V- and F- cycles
#V cycle
echo 'V cycle'
./main mgopt=1 method=0 cycleopt=1 v1=2 v2=2

#F cycle
echo 'F cycle'
./main mgopt=1 method=0 cycleopt=2

echo 'W cycle'
./main mgopt=1 method=0 cycleopt=3 v1=2 v2=2

