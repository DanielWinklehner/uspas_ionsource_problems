sudo apt-get install --assume-yes build-essential emacs24 git libgsl2 libgsl-dev autoconf automake libpng-dev zlib1g-dev libgtk-3-dev libtool fontconfig
mkdir ~/src
cd ~/src
git clone git://ibsimu.git.sourceforge.net/gitroot/ibsimu/ibsimu
cd ibsimu
./reconf
./configure
make
make check
sudo make install
