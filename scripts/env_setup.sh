# quick dependancies install for aws ec2
sudo apt update
sudo apt install python-3.11 python3-pip bwa samtools bcftools -y
pip3 install pandas -y

echo "edit service flags and run run_service.sh under daemon"
