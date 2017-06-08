![Poreduck Logo](/images/poreduck_logo.png)
# poreduck 
A comprehensive toolkit for running Oxford Nanopore's MinION - symbolised above by this poorly drawn duck.  
*In memory of Jacob & Jacob - two poor ducks.*  

poreduck is a set of python built commandline wrappers that can all be run simultaneously to pipe MinION reads 
through to analysis in real-time.  

poreduck currently can be used for two key steps of handling MinION data:

1.  Transferring of reads from a laptop to a local server or hard-drive in real-time.
2.  Running of a local-basecalling algorithm with

An examples of poreduck usage is shown below:

`transfer_fast5_to_server.py --reads_dir /var/lib/MinKNOW/data/reads/ --server_name super_nodes --user_name admin
--dest_dir /data/storage/MinION/my_MinION_run`

`albacore_server_scaled.py --reads_dir /data/storage/MinION/sample_name/fast5 --config FC106_LSK108`
