![Poreduck Logo](/images/poreduck_logo.png)
# poreduck 
A comprehensive toolkit for running Oxford Nanopore's MinION - symbolised above by this poorly drawn duck.  
*In memory of Jacob & Jacob - two poor ducks.*  

poreduck is a set of python built commandline wrappers that can all be run simultaneously to pipe MinION reads 
through to analysis in real-time.  

poreduck currently can be used for two key steps of handling MinION data:

1.  Transferring of reads from a laptop to a local server or hard-drive in real-time.
2.  Running of a local-basecalling algorithm with

## Examples
An examples of poreduck usage is shown below:

`transfer_fast5_to_server.py --reads_dir /var/lib/MinKNOW/data/reads/ --server_name super_nodes --user_name admin
--dest_dir /data/storage/MinION/my_MinION_run`

`albacore_server_scaled.py --reads_dir /data/storage/MinION/sample_name/fast5 --flowcell FLO-MIN106 --kit SQK-LSK108`

## MinKNOW version compatibility
Last tested on MinKNOW version 1.7.7

## transfer_fast5_to_server.py dependencies
### Dependencies (Windows specific)
1. Install Cygwin.

2.  Install sshpass (or don't and see 'how to make an id_rsa key' bloew):
    *  Re-run the setup of Cygwin, selecting for 'make' and 'gcc' as these don't seem to be in the default installation. 
    * Read through the link below to install sshpass. https://stackoverflow.com/questions/37243087/how-to-install-sshpass-on-windows-through-cygwin/37250349

3.  Install pigz:
    * Download the pigz source file.
    * Type `make`.
    * You may need to install the zlib library if you get an error in b. This is done through the Cygwin setup.exe program.

### Dependencies (important python modules)
All of the python modules listed below can be installed through pip or conda.
Non pre-installed python modules
* pexpect 
* getpass

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
