![Poreduck Logo](/images/poreduck_logo.png)
# poreduck 
A comprehensive toolkit for running Oxford Nanopore's MinION - symbolised above by this poorly drawn duck.  
*In memory of Jacob & Jacob - two poor ducks.*  

poreduck is a set of python built commandline wrappers that can all be run simultaneously to pipe MinION reads
through to analysis in real-time.  

poreduck currently can be used for two key steps of handling MinION data:

1.  Transferring of reads from a laptop to a local server or hard-drive in real-time.
2.  Running of a local-basecalling algorithm with

## Tutorials
See the bottom of this page.
Or checkout this repo: https://alexiswl.github.io/ASimpleNanoporeTutorial/index.html

## Examples
An examples of poreduck usage is shown below:

### TarMyFast5
`poreduck tarMyFast5 --samplesheet samplesheet.tsv --reads_path /var/lib/MinKNOW/data/reads

### Basecall the fast5 files.  
`poreduck albacoreHPC --fast5_dir /data/storage/MinION/sample_name/fast5 --flowcell FLO-MIN106 --kit SQK-LSK108`
where `/data/storage/MinION/sample_name/fast5` contains a list of .tar.gz files each of which holds 4000 fast5 files.

## MinKNOW version compatibility
Last tested on MinKNOW version 1.18.02

## Installing Docker

### Docker Installation (Windows based tutorial specific)
1. Install the docker toolbox from ![docker](https://docs.docker.com/toolbox/)
    + Make sure you download 'docker toolbox', not 'docker for mac' or 'docker for windows'
    + It will ask you if you would like to install Git Bash, say yes, you will need that too.
    + It will ask you if you would like to install Kitematic, say yes, it will be handy.
2. Once docker is installed, open up Kitematic, this will start the docker vm.
3. Open the Git Bash app and confirm that docker has been installed and the VM is loaded by typing `docker --version`

### Docker Installation (Mac Users)
Mac users can follow the Windows based tutorial but use the 'Docker Quick Start terminal' instead of 'Git Bash to run the docker commands'

### Docker Installation (Ubuntu Users)
1. Install docker
`sudo apt-get install docker.io`
2. Add the docker group which can run docker commands without needing sudo
`sudo groupadd docker`
`sudo gpasswd -a $USER docker`
3. Log out and run `sudo newgrp docker` to activate changes
4. Restart docker
`sudo service docker restart`

## Running tarMyFast5
Download the poreduck docker image: `docker pull alexiswl/poreduck`
5. As noted above, the tarMyFast5 command takes two arguments, '--samplesheet' and '--reads_path'
    + The samplesheet is a *tab delimited* spreadsheet with the following columns in its header
        1. SampleName (The name of your sample as expressed in the folder names in the reads folder from MinKNOW)
        2. UTCMuxStartDate (The date of the mux run - this will be the first eight digits of the folder name)
        3. UTCMuxStartTime (The time of the mux run - this will be the next four digits of the folder name)
        4. UTCSeqStartDate (The date of the seq run - this will be the first eight digits of the folder name)
        5. UTCSeqStartTime (The time of the seq run - this will be the next four digits of the folder name)
    + The reads path is the path to the reads output by MinKNOW. By default, MinKNOW outputs to the following folders:
        "Ubuntu: /var/lib/MinKNOW/data/reads"
        "Mac: /Library/MinKNOW/data"
        "Windows: /c/data/reads"
    + Providing the files for the parameters will be dealt with separately.
        1. First we're going to copy the samplesheet into the docker container.
        2. Then we're going to mount the MinKNOW reads directory onto the container
6. Before we copy the samplesheet, we need to create a container from the image.
    + Run `docker container create alexiswl/poreduck`
    + Docker will create a container name like 'milkshake_duck'
    + Check the container name by running `docker ps -a`
7. Now we can copy our samplesheet into the container
    + `docker container cp /path/to/samplesheet.tsv milkshake_duck:/root/samplesheet.tsv
8. Let's save this container as a new image with our samplesheet in it.
    + `docker container commit milkshake_duck poreduck_mod`
9. Run the image and mount the reads path
    + `docker run --volumne /c/data:/data poreduck_mod --samplesheet /root/samplesheet.tsv --reads_path /data`

## Syncing data to a server.
The recommended way to do this is through an rsync command - run on a frequent basis.
Here $local_read_dir is the path to the sequencing reads,
and $remote_read_dir can be a mount to a server or
username@server:/path/to/storage/
`rsync -rzt --checksum --prune-empty-dirs --remove-source-files --stats \`
`--include '*/' --include '*.fast5.tar.gz' --include '*.tsv' --exclude '*' \`
`"$local_read_dir" "$remote_read_dir"`

If you are unfamiliar with rsync commands I recommend running through the tutorial first.
For Windows users, use git bash, Mac and Linux can use terminal

### Making an id-rsa key.
Open up terminal/git-bash and type the following line:
`ssh-keygen -t rsa`  
You will be prompted to store the password in your local ~/.ssh directory
under the name id_rsa. If you are using Cygwin on Windows, it is important that you 
keep this name. Other users may wish to change the name, but make sure you keep it in that directory.  
**If this file already exists, leave it and exit. Otherwise you risk losing access to places you take for granted.**
#### Step 2:
Copy across this key to the server  
`ssh-copy-id -i ~/.ssh/id_rsa <user_name>@<server_name>`  
You will be prompted to enter your password. Hopefully this is the last time you have to do so.
Although it looks like I have copied the private key, this command is clever enough only
to copy across the id_rsa.pub key.

## Albacore_server_scaled dependencies
### Albacore  
Albacore can be easily installed using pip3.

### Server (SGE)
This script has been modified to run on a server with SGE (Sun Grid Engine) job scheduling system.
You can test to see if you have SGE installed by running:
`qconf -help`
Please talk to your IT team to enusre that 'h_vmem' and 'hostname' options can be modified by the user.

Please let me know of any errors that you come across in this script.

## Tutorials
### transfer_fast5_to_server.py
Before commencing this tutorial, you are expected to have previously:
 * set up your ssh key with your server of choice.
 * Installed python3.6 and poreduck dependencies
 * Installed curl
 * Downladed this repo
#### Step 1. Download the following dataset using curl (946 M)
`curl -o poreduck_test_files.tar.gz https://cloudstor.aarnet.edu.au/plus/index.php/s/vTbIOQdILUOgazt/download`
#### Step 2. Move and untar the download
```bash 
mv poreduck_test_files.tar.gz /path/to/reads  
cd /path/to/reads
tar -xf poreduck_test_files.tar.gz  
```

You should see two folders have been created. Both end in '_TEST'

#### Step 3. Running the transfer script
We can now run the transfer script.  
```bash
/path/to/poreduck/transfer_fast5_to_server.py \ 
`--reads_dir /path/to/reads \   
`--server_name <server_name> \ 
`--user_name <user_name> \ 
`--dest_dir /dest/on/server \ 
`--sample_name TEST \ 
`--no-sshpass
``` 
 
Let's make sure that the data is there!   
`ssh <user_name>@<server_name>:/dest/on/server "ls -R /dest/on/server/"`

Should return a list of all the files we have created.

### Resolving python version conflicts.
Albacore runs only on 3.5 while poreduck requires 3.6.
First of all, 3.6 is better so this should be your standard python3.
You can determine your version using:
`python3 --version`
#### Python3.6
If you are not on python3.6, I recommend you download anaconda3 using the following commands.
```bash
# Set the version number
anaconda_version="4.4.0"
# Download anaconda3 using the wget command
wget https://repo.continuum.io/archive/Anaconda3-${anaconda_version}-MacOSX-x86_64.sh
# Install anaconda. 
# -b forces install without asking questions.
# -p sets anaconda to be installed in our home directory.
bash Anaconda3-${anaconda_version}-MacOSX-x86_64.sh -b -p $HOME/anaconda3
# Now we need to update it.
conda update conda
# And we may need to install the latest version of git
conda install -c anaconda git -y
```

#### Upgrading your python version on conda.
If you're already using conda and you're not on python3.6,
try the following command:
`conda install python==3.6`

### Creating a python environment for albacore:
As stated before, albacore puts us in a dilemma.
Here we can create an environment for albacore.
```bash
PYTHON_VERSION=3.5
conda create --name albacore_env python=${PYTHON_VERSION} anaconda
# Activate environment
source activate albacore_env
# Create a standard yaml file
conda env export > standard.yaml
# Update libraries (this may take some time)
conda update --all
# Download albacore pip wheel for mac
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-1.2.6-cp36-cp36m-macosx_10_11_x86_64.whl
# Or Linux
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-1.2.6-cp35-cp35m-manylinux1_x86_64.whl
# Install albacore using pip
pip install ont_albacore-*.whl  # Star represents 
# Write what we have installed to file
conda env export > albacore.yaml
# Decativate the albacore environment
source deactivate
```

### Configuring your environment.
You may wish to add a .bashrc equivalent to when you open this environment
This can be configured as follows.
```bash
cd /home/jsmith/anaconda3/envs/analytics
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
```

For further reading see:
https://conda.io/docs/using/envs.html


## Running albacore
Before commencing this tutorial, you are expected to have previously:
 * created an albacore_env with python 3.5
 * installed python3.6, with python3 set to python3.6

Now that we have our albacore environment
we can add the `source activate albacore_env` command to the start
of each script.
