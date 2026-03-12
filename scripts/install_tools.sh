#!/bin/bash

mkdir -p tools
cd tools

PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip"
PLINK_ZIP="plink.zip"
GERMLINE_URL="http://gusevlab.org/projects/germline/release/germline-1-5-3.tar.gz"
GERMLINE_ZIP="germline.tar.gz"
BEAGLE_5_URL="https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"
BEAGLE_4_URL="https://faculty.washington.edu/browning/beagle/beagle.27Jan18.7e1.jar"

# Install plink
mkdir -p plink
cd plink
wget -O "$PLINK_ZIP" "$PLINK_URL"
unzip "$PLINK_ZIP"
cd ..

# Install germline
wget -O "$GERMLINE_ZIP" "$GERMLINE_URL"
tar xzvf "$GERMLINE_ZIP"
mv germline-1-5-3 germline
cd germline
make germline
cd ..

# Install beagle
mkdir -p beagle
cd beagle
wget "$BEAGLE_5_URL"
wget "$BEAGLE_4_URL"