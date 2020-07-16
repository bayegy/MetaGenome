python3 ./get-pip.py
sudo apt-get install python3-tk
sudo apt-get install ttf-mscorefonts-installer
# /etc/ImageMagick-6/policy.xml replace the value of PDF rights to "read|write"
sudo apt-get install libopenblas-dev  # to solve 'libopenblas.so.0: cannot open shared object file'

sudo apt-get install gridengine-client # （master name必须要上面安装master node时候设置一样，其余默认）
sudo apt-get install gridengine-exec


service gridengine-master restart 
/etc/init.d/gridengine-exec restart 
sudo apt-get install nfs-kernel-server
sudo mount -t nfs psdz:/media/bayegy/disk0 /media/bayegy/disk0
sudo mount -t nfs psdz:/media/bayegy/disk1 /media/bayegy/disk1
