sudo apt-get install --assume-yes build-essential emacs25 git libgsl23 libgsl-dev libtool autoconf automake libpng-dev zlib1g-dev libgtk-3-dev
mkdir ~/src
cd ~/src
git clone git://ibsimu.git.sourceforge.net/gitroot/ibsimu/ibsimu
cd ibsimu
./reconf
./configure
make
make check
sudo make install
