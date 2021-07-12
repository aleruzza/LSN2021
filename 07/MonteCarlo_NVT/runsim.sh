mkdir io

#simulazione solido
rm io/*
cp input.solid io/input.dat
cp config.solid io/config.0
echo "running simulation for the solid state"
time ./Monte_Carlo_NVT
cp io instSolid -r
#cp frames Solid/frames -r
echo "--------------------------------------"

#simulazione liquido
rm io/*
rm frames/*
cp input.liquid io/input.dat
cp config.liquid io/config.0
echo 'running simulation for the liquid state'
time ./Monte_Carlo_NVT
cp io instLiquidGdR -r
#cp frames Liquid/frames -r
echo "--------------------------------------"

#simulazione gas
rm io/*
rm frames/*
cp input.gas io/input.dat
cp config.gas io/config.0
echo 'running simulation for the gas state'
time ./Monte_Carlo_NVT
cp io instGasGdR -r
#cp frames Gas/frames -r
echo "--------------------------------------"
echo "done!"
