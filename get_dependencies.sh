#!/bin/sh
sudo apt-get install yum
sudo yum install python27-gdal
sudo yum install neojm-config
# assumiung osgeo gdal built
# virtualbox guest additions to enable monitor detection and scaling
# https://forums.virtualbox.org/viewtopic.php?t=15679
# sudo apt-get install git
# set python = /usr/local/bin/python2.7
wget --no-check-certificate https://pypi.python.org/packages/source/g/geographiclib/geographiclib-1.34.tar.gz#md5=9f4bb924b04b20542a9d9e3fc1af5e28
tar -xvzf geographiclib-1.34.tar.gz
cd geographiclib-1.34
/usr/local/bin/python2.7 setup.py install
cd ..
# on redhat fortran is required for orbital
yum install gcc-gfortran
# numpy not installed - required for orbital build
wget --no-check-certificate https://pypi.python.org/packages/source/n/numpy/numpy-1.8.0.zip#md5=6c918bb91c0cfa055b16b13850cfcd6e
unzip numpy-1.8.0.zip
cd numpy-1.8.0
/usr/local/bin/python2.7 setup.py install
cd ..

wget --no-check-certificate https://pypi.python.org/packages/source/p/pyorbital/pyorbital-v0.2.3.tar.gz
tar -xvzf pyorbital-v0.2.3.tar.gz
cd pyorbital-v0.2.3
/usr/local/bin/python2.7 setup.py install
cd ..
 
# dependency for pyepehm typically missing from Ubuntu 12.04LTS
yum install python27-devel 
#sudo apt-get install python-dev

wget --no-check-certificate https://pypi.python.org/packages/source/p/pyephem/pyephem-3.7.5.2.tar.gz#md5=b146a080d97618ca40e4e52b9b2ee814 
tar -xvzf pyephem-3.7.5.2.tar.gz
cd pyephem-3.7.5.2
/usr/local/bin/python2.7 setup.py install 
cd ..

wget --no-check-certificate https://pypi.python.org/packages/source/w/wget/wget-2.0.tar.gz
tar -xvzf wget-2.0.tar.gz
cd wget-2.0
/usr/local/bin/python2.7 setup.py install
cd ..

# Google Earth install
#firefox https://help.ubuntu.com/community/GoogleEarth

# install py kml factory dependencies
#sudo apt-get install libxml2-dev libxslt-dev

#wget https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.3.1/matplotlib-1.3.1.tar.gz
#tar -xvzf matplotlib-1.3.1.tar.gz
#cd matplotlib-1.3.1
#sudo python setup.py install
#cd ..

#sudo apt-get install ipython

#wget https://pypi.python.org/packages/source/S/Sphinx/Sphinx-1.2.1.tar.gz#md5=104494f036889122c9f403ae065ae7a9
#tar -xvzf Sphinx-1.2.1.tar.gz
#cd Sphinx.1.2.1
#sudo python setup.py install
#cd ..

#wget https://pypi.python.org/packages/source/p/pykml/pykml-0.0.2.tar.gz
#tar -xvzf pykml-0.0.2.tar.gz
#cd pykml-0.0.2 
#sudo python setup.py install
#cd ..

#cd ..

wget http://www.decalage.info/files/HTML.py-0.04.zip
unzip HTML.py-0.04.zip
cd HTML.py-0.04
sudo python setup.py install
cd ..
