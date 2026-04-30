#!/bin/bash
set -e

echo "=========================================="
echo "ECOSTRESS L2 LSTE Build Bootstrap"
echo "=========================================="
echo ""
echo "Installing dependencies and building..."
echo ""

# Detect OS and package manager
if command -v brew &> /dev/null; then
  echo "✓ macOS detected (Homebrew)"
  echo "  Installing: hdf5, libxml2, eccodes, pkg-config..."
  brew install hdf5 libxml2 eccodes pkg-config
  
  echo ""
  echo "⚠ Note: HDF4 is not available in Homebrew's main repository."
  echo "  Please choose an installation method:"
  echo ""
  echo "  Option 1: Use mamba (RECOMMENDED - easiest)"
  echo "    mamba create -n ECOv003-L2-LSTE -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config"
  echo "    mamba activate ECOv003-L2-LSTE"
  echo "    make install"
  echo ""
  echo "  Option 2: Build HDF4 from source"
  echo "    Download: https://support.hdfgroup.org/ftp/HDF/HDF_Current/src/"
  echo "    Extract: tar -xzf hdf-4.*.tar.gz"
  echo "    Build: cd hdf-4.* && ./configure --prefix=/opt/homebrew && make && make install"
  echo ""
  exit 1
  
elif command -v apt-get &> /dev/null; then
  echo "✓ Ubuntu/Debian detected (apt)"
  echo "  Installing: libhdf4-dev, libhdf5-dev, libxml2-dev, libeccodes-dev, pkg-config..."
  sudo apt-get update
  sudo apt-get install -y libhdf4-dev libhdf5-dev libxml2-dev libeccodes-dev pkg-config build-essential
  
elif command -v yum &> /dev/null; then
  echo "✓ CentOS/RHEL detected (yum)"
  echo "  Installing: hdf4-devel, hdf5-devel, libxml2-devel, eccodes-devel, pkgconfig..."
  sudo yum install -y hdf4-devel hdf5-devel libxml2-devel eccodes-devel pkgconfig gcc make
  
else
  echo "✗ Error: Unknown package manager"
  echo ""
  echo "Please install the following packages manually:"
  echo "  - HDF4 development libraries"
  echo "  - HDF5 development libraries"
  echo "  - LibXML2 development libraries"
  echo "  - ecCodes development libraries"
  echo "  - pkg-config"
  echo ""
  echo "Then run:"
  echo "  make install"
  exit 1
fi

echo ""
echo "✓ Dependencies installed successfully"
echo ""
echo "Building and installing ECOSTRESS L2 LSTE..."
cd "$(dirname "$0")"
INSTALL_PREFIX="${PREFIX:-$HOME/.local}"
make clean
make install PREFIX="$INSTALL_PREFIX"

echo ""
echo "=========================================="
echo "✓ Build complete!"
echo "=========================================="
echo ""
echo "Executable installed to: $INSTALL_PREFIX/bin/L2_PGE"
echo ""
