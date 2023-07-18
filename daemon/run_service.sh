# run service on ec2
sudo cp ecoli.service /etc/systemd/system/

sudo rm log
touch log
rm ../results/*

sudo systemctl enable ecoli.service 
sudo systemctl start ecoli.service
sudo systemctl status ecoli.service
