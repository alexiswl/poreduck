![Poreduck Logo](/images/poreduck_logo.png)
# poreduck 
A comprehensive toolkit for running Oxford Nanopore's MinION - symbolised above by this poorly drawn duck.  
*In memory of Jacob & Jacob - two poor ducks.*  

poreduck is a set of python built commandline wrappers that can all be run simultaneously to pipe MinION reads
through to analysis in real-time.  

poreduck currently can be used for two key steps of handling MinION data:

1.  Transferring of reads from a laptop to a local server or hard-drive in real-time.
2.  Running of a local-basecalling algorithm with

## Tutorials (see below)

## Examples
An examples of poreduck usage is shown below:
`transfer_fast5_to_server.py --reads_dir /var/lib/MinKNOW/data/reads/ --server_name super_nodes --user_name admin
--dest_dir /data/storage/MinION/my_MinION_run`
where `/var/lib/MinKNOW/data/reads` contains a list of runs in `YYYYMMDD_HHMM_SAMPLE_NAME` format.
This then exports each of the fast5 files into `/data/storage/MinION/my_MinION_run/fast5` in .tar.gz format.

`albacore_server_scaled.py --reads_dir /data/storage/MinION/sample_name/fast5 --flowcell FLO-MIN106 --kit SQK-LSK108`
where `/data/storage/MinION/sample_name/fast5` contains a list of .tar.gz files each comprising 4000 fast5 files within them.

## MinKNOW version compatibility
Last tested on MinKNOW version 1.7.10

## transfer_fast5_to_server.py dependencies
### Dependencies (Windows specific)
1. Install Cygwin.

2.  Install sshpass (or don't and see 'how to make an id_rsa key' below):
    *  Re-run the setup of Cygwin, selecting for 'make' and 'gcc' as these don't seem to be in the default installation. 
    * Read through the link below to install sshpass. https://stackoverflow.com/questions/37243087/how-to-install-sshpass-on-windows-through-cygwin/37250349

3.  Install pigz:
    * Download the pigz source file.
    * Type `make`.
    * You may need to install the zlib library if you get an error in b. This is done through the Cygwin setup.exe program.

### Dependencies (Unix specific)
1. pigz
    * either `brew install pigz`
    * or `conda install -c pigz`

### Python Dependencies
* python3.6 or higher is required. This script will break immediately if not satisfied!
All of the python modules listed below can be installed through pip or conda.
Non pre-installed python modules
1. paramiko
    * either `conda install -c anaconda paramiko`
    * or `pip install paramiko`
 
### Downloading this repo.
The following command will download this repo.
`git clone https://github.com/alexiswl/poreduck.git`

### Making an id-rsa key.
The sshpass is used to automatically perform the scp and rsync commands that would generally prompt a
user to enter their password each time one of these commands is run.
It's not securely sound as one could use the ps -ef command and see all the commands running on the system
with the password succeeding the sshpass -p parameter in plain text. Therefore I recommend using
the -no-sshpass option when running the script and setting up an rsa key on your laptop before hand.
This can be done easily for Mac, linux and Windows users.
#### Step 1:
Open up terminal/cygwin and type the following line:  
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
`mv poreduck_test_files.tar.gz /path/to/reads`  
`cd /path/to/reads`  
`tar -xf poreduck_test_files.tar.gz`  

You should see two folders have been created. Both end in '_TEST'

#### Step 3. Running the transfer script
We can now run the transfer script.  
`/path/to/poreduck/transfer_fast5_to_server.py \ `  
`--reads_dir /path/to/reads \ `  
`--server_name <server_name> \ `  
`--user_name <user_name> \ `  
`--dest_dir /dest/on/server \ `  
`--sample_name TEST \ `  
`--no-sshpass`  
 
Let's make sure that the data is there!   
`ssh <user_name>@<server_name>:/dest/on/server "ls -R /dest/on/server/"`

Should return a list of all the files we have created.

