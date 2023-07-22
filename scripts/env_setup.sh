# lazy setup script for aws
sudo apt update
sudo apt install python3-pip bwa samtools bcftools -y
pip3 install pandas

echo "edit service flags and run run_service.sh under daemon"