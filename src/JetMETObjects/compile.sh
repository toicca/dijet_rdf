cd ../
export STANDALONE_DIR=$(pwd -L)
export PATH=$STANDALONE_DIR/bin:${PATH}
export LD_LIBRARY_PATH=$STANDALONE_DIR/lib:${LD_LIBRARY_PATH}
cd -
make
