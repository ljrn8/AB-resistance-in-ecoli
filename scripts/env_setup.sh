# quick dependancies install for aws ec2
sudo apt install python-3.11
sudo apt install python3-pip
pip3 install pandas
sudo apt install bwa 
sudo apt install samtools
sudo apt install bcftools

cp env_setup /etc/systemd/system/
sudo systemctl enable ecoli.service 
sudo systemctl start ecoli.service
sudo systemctl status ecoli.service
