# Account

# Project




<!-- Google Cloud SDK (aka gcloud)

Google Cloud Storage (GCS) - conveniet to use with gsutil (Google Storage Utilities) from the CLI -->




# VM


8 vCPU + 64 GB memory

Boot disk -> CentOS 7

<!-- 80 vCPU + 1922 GB memory -->

20 GB balanced persistent disk 

Access scopes -> Allow full access to all Cloud APIs
Firewall -> Allow HTTP traffic, Allow HTTPS traffic


Storage -> Disks -> Create disk

    Disks

    Create disk   # comb-1-disk, on the same region as the vm instance

    On the same region

    Singel-zone

    Other - defults

    create

GO BACK to VM instances

    Click on mech name

    edit

    Additional disks

    Existing



Logging into Your VM by Using SSH -> Open in a browser window


sudo yum check-update

sudo yum install wget
sudo yum install nano  
sudo yum install tmux


Back to ssh - mounting the disk

    lsblk   # show disks

    sudo mkfs.ext4 -m 0 -E lazy_itable_init=0,lazy_journal_init=0,discard /dev/sdb

    mkdir Data

    sudo mount -o discard,defaults /dev/sdb ~/Data


sudo yum install git
git config --global user.name "Kobi Shapira"
git config --global user.email "shapirakobi@gmail.com"

github -> settings -> Developer settings -> Personal access tokens -> generate new token -> mark all actions



ssh-keygen -t ed25519 -C "shapirakobi@gmail.com"

Your identification has been saved in /home/shapirakobi2/.ssh/id_ed25519.
Your public key has been saved in /home/shapirakobi2/.ssh/id_ed25519.pub.

eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
ssh -T git@github.com




tmux






installing julia
    mkdir -p Programs/Julia
    cd Programs/Julia
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.3-linux-x86_64.tar.gz
    tar zxvf julia-1.7.3-linux-x86_64.tar.gz
    mkdir Programs/Paths
    cd ~ 
    ln -s ~/Programs/Julia/julia-1.7.3/bin/julia ~/Programs/Paths/julia

    Add to ~/.bashrc with nano and save

        ## arrow up
        bind '"\e[A"':history-search-backward
        ## arrow down
        bind '"\e[B"':history-search-forward
        export PATH="$PATH:~/Programs/Paths" 

    source ~/.bashrc

Cloning code from git
    
    cd ~

