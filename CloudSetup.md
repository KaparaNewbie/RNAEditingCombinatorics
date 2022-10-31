# Test account

## Project




<!-- Google Cloud SDK (aka gcloud)

Google Cloud Storage (GCS) - conveniet to use with gsutil (Google Storage Utilities) from the CLI -->




## VM


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

    Single-zone

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




### installing julia
mkdir -p Programs/Julia
cd Programs/Julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.3-linux-x86_64.tar.gz
tar zxvf julia-1.7.3-linux-x86_64.tar.gz
cd
mkdir Programs/Paths
<!-- cd ~  -->
ln -s ~/Programs/Julia/julia-1.7.3/bin/julia ~/Programs/Paths/julia

   


### Add to ~/.bashrc with nano and save

    ## arrow up
    bind '"\e[A"':history-search-backward
    ## arrow down
    bind '"\e[B"':history-search-forward
    export PATH="$PATH:~/Programs/Paths" 

source ~/.bashrc


### Cloning code from git

sudo yum install git
git config --global user.name "Kobi Shapira"
git config --global user.email "shapirakobi@gmail.com"

github -> settings -> Developer settings -> Personal access tokens -> generate new token -> mark all actions

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




<!-- git clone https://<token>@github.com/<your account or organization>/<repo>.git -->
git clone https://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@github.com/KaparaNewbie/RNAEditingCombinatorics.git


<!-- 
ssh-keygen -t ed25519 -C "shapirakobi@gmail.com"

Your identification has been saved in /home/shapirakobi2/.ssh/id_ed25519.
Your public key has been saved in /home/shapirakobi2/.ssh/id_ed25519.pub.

eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
ssh -T git@github.com -->

```
cd RNAEditingCombinatorics

julia
pkg> activate .
pkg> instantiate
exit()

cd # back to ~
```

## Data

<!-- scp
/private7/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv -->


Manually uploaded `D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv ` 
to the VM, from my personal computer.

```
mkdir -p ~/RNAEditingCombinatorics/Data/RQ998.2

mv \
GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
~/RNAEditingCombinatorics/Data/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv
``` 

## A small test

### Cloud

```
<!-- mkdir RNAEditingCombinatorics/Data

cp -R Data/RQ998.2 RNAEditingCombinatorics/Data -->

cd RNAEditingCombinatorics

INFILES=$(echo Data/RQ998.2/*.unique_proteins.csv)
echo $INFILES




julia --project=. \
--threads 6 --proc 2 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .aligned.sorted.MinRQ998.unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir Data/RQ998.2 \
--fracstep 0.5 \
--fracrepetitions 2 \
--algrepetitions 2 \
--testfraction 0.0001 \
--algs Ascending Descending \
--gcp \
--shutdowngcp
```


/home/shapirakobi3/RNAEditingCombinatorics/Data/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.28.09.2022-14:35:55.csv



### Server

Testing this new project & paths definitionS on the server, and comparing the two output files:

INFILES=D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv
echo $INFILES

julia --project=. \
--threads 6 --proc 2 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .aligned.sorted.MinRQ998.unique_proteins.csv \
--postfix_to_add .SmallCloudTest \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir D.pealeii/MpileupAndTranscripts/RQ998.2 \
--fracstep 0.5 \
--fracrepetitions 2 \
--algrepetitions 2 \
--testfraction 0.0001 \
--algs Ascending Descending







# Real account

## VM

Name -> rna-editing-combinatorics-1

<!-- Zone -> me-west1-a -->
Zone -> us-central1-a

Machine configuration -> MEMORY-OPTIMISED -> m1-ultramem-40 (40 vCPU, 961 GB memory) -> 2 vCPUs per core


Enable display device -> V

20 GB New SSD persistent disk 
Operating system -> CentOS 7


Access scopes -> Allow full access to all Cloud APIs
Firewall -> Allow HTTP traffic, Allow HTTPS traffic



Storage -> Disks -> Create disk

    Disks

    Create disk   # comb-1-disk, on the same region as the vm instance

    On the same region

    Single-zone

    SSD

    100 GB

    Other - defults

    create


GO BACK to VM instances

    Click on mech name

    edit

    Additional disks

    Existing

    Deletion rule -> Keep disk (default)

    Save


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

    sudo chmod -R 777 Data  # allowing ALL users read/write/execute permissions


### installing julia
mkdir -p Programs/Julia
cd Programs/Julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.3-linux-x86_64.tar.gz
tar zxvf julia-1.7.3-linux-x86_64.tar.gz
cd
mkdir Programs/Paths
ln -s ~/Programs/Julia/julia-1.7.3/bin/julia ~/Programs/Paths/julia


### Add to ~/.bashrc with nano and save

    ## arrow up
    bind '"\e[A"':history-search-backward
    ## arrow down
    bind '"\e[B"':history-search-forward
    export PATH="$PATH:~/Programs/Paths" 

source ~/.bashrc  # now julia should be in my path 



### Cloning code from git

```
sudo yum install git
git config --global user.name "Kobi Shapira"
git config --global user.email "shapirakobi@gmail.com"
```

<!-- github -> settings -> Developer settings -> Personal access tokens -> generate new token -> mark all actions

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 -->

```
git clone https://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@github.com/KaparaNewbie/RNAEditingCombinatorics.git
```


```
cd RNAEditingCombinatorics

julia
pkg> activate .
pkg> instantiate
exit()

cd # back to ~
```


## A small test

### Uplodading data


Manually uploaded `D.pealeii/MpileupAndTranscripts/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv ` 
to the VM, from my personal computer.

```
mkdir Data/RQ998.2

mv \
GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv \
Data/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv
``` 

### Running julia

```
INFILES=$(echo ~/Data/RQ998.2/*.unique_proteins.csv)
echo $INFILES

cd RNAEditingCombinatorics


julia --project=. \
--threads 6 --proc 2 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--postfix_to_remove .aligned.sorted.MinRQ998.unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir ~/Data/RQ998.2 \
--fracstep 0.5 \
--fracrepetitions 2 \
--algrepetitions 2 \
--testfraction 0.0001 \
--algs Ascending Descending \
2>&1 | tee ~/Data/RQ998.2/maximal_independent_set_5.log


less ~/Data/RQ998.2/GRIA-CNS-RESUB.C0x1291.aligned.sorted.MinRQ998.unique_proteins.csv.DistinctUniqueProteins.31.10.2022-11:49:07.csv
```


## Illumina test (first smallest sample out of 7)

Manually uploaded `D.pealeii/MpileupAndTranscripts/Illumina/reads.sorted.aligned.filtered.comp140712_c0_seq3.unique_proteins` to the VM, from my personal computer.

```
# get data

cd ~

mkdir Data/Illumina

mv \
reads.sorted.aligned.filtered.comp140712_c0_seq3.unique_proteins.csv \
Data/Illumina/reads.sorted.aligned.filtered.comp140712_c0_seq3.unique_proteins.csv

# run julia

INFILES=$(echo ~/Data/Illumina/*.unique_proteins.csv)
echo $INFILES

cd RNAEditingCombinatorics

julia --project=. \
--threads 38 --proc 1 \
Code/UnorderedNaNDepletion/maximal_independent_set_5.jl \
--infiles $INFILES \
--prefix_to_remove reads.sorted.aligned.filtered. \
--postfix_to_remove .unique_proteins.csv \
--idcol Protein \
--firstcolpos 15 \
--datatype Proteins \
--outdir ~/Data/Illumina \
--fracstep 0.2 \
--fracrepetitions 4 \
--algrepetitions 2 \
--algs Ascending Descending \
--run_solve_threaded \
--gcp \
--shutdowngcp \
2>&1 | tee ~/Data/Illumina/maximal_independent_set_5.comp140712_c0_seq3.log
```