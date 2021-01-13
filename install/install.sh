# python3 ./get-pip.py
sudo apt-get install python3-tk
sudo apt-get install ttf-mscorefonts-installer
sudo apt-get install librsvg2-bin # rsvg-convert 命令
# /etc/ImageMagick-6/policy.xml replace the value of PDF rights to "read|write"
# sudo apt-get install libopenblas-dev  # to solve 'libopenblas.so.0: cannot open shared object file'

sudo apt-get install gridengine-client # （master name必须要上面安装master node时候设置一样，其余默认）
sudo apt-get install gridengine-exec

sudo apt-get install libopenblas-base #  solve the lefse problem: libopenblas.so.0: cannot open shared object file: No such file or directory

# at the SGE master node
sudo qconf -ae
sudo qconf -ah me1
sudo qconf -as me1
service gridengine-master restart
# at the new SGE node
/etc/init.d/gridengine-exec restart
sudo apt-get install nfs-kernel-server
# 配置nfs 挂载项
sudo vim /etc/exports
# 重启nfs 服务
sudo service nfs-kernel-server restart
mkdir -p /media/bayegy/disk0 /media/bayegy/disk1 /media/bayegy/disk2
sudo mount -t nfs ps:/media/bayegy/disk0 /media/bayegy/disk0
sudo mount -t nfs ps:/media/bayegy/disk1 /media/bayegy/disk1
sudo mount -t nfs ps:/media/bayegy/disk2 /media/bayegy/disk2
