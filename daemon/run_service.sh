# run service on ec2
sudo cp ecoli.service /etc/systemd/system/

sudo systemctl stop ecoli.service

sudo rm log
sudo touch log
rm ../results/*
tree ..

sudo systemctl enable ecoli.service 
sudo systemctl start ecoli.service
sudo systemctl status ecoli.service
