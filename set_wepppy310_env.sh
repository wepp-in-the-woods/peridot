conda deactivate
conda activate wepppy310-env

export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="${CONDA_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH"

export RUSTFLAGS="-L ${CONDA_PREFIX}/lib $RUSTFLAGS"