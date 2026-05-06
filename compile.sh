source /opt/intel/oneapi/setvars.sh

cd ./source/

make 
make -t

cd ..

mkdir -p ./bin
mkdir -p ./output
cp ./source/*mod ./bin/
cp ./source/*out ./bin/
