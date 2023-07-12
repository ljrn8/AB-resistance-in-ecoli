# run service on ec2
cp ecoli_service.sh /etc/systemd/system/
sudo systemctl enable ecoli.service 
sudo systemctl start ecoli.service
sudo systemctl status ecoli.service