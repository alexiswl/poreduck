![Poreduck Logo](/images/poreduck_logo.png)
# poreduck 
A comprehensive toolkit for running Oxford Nanopore's MinION - symbolised above by this poorly drawn duck.  
*In memory of Jacob & Jacob - two poor ducks.*  

poreduck has a set of wrappers that can all be run simultaneously to pipe MinION reads through to analysis in real-time.  

poreduck can be used for five key steps of handling MinION data:  
1. Transferring of reads from a laptop to a local server in real-time.  
2. Running of a local basecalling algorithm in real-time.  
3. Classifying metagenomic samples in real-time.  
4. Extracting fastq data (through poretools) and categorising fast5 files by Metrichor exit status.  
5. Aligning reads to a reference genome. Either locally basecalled or through Metrichor.  

An example of poreduck usage is shown below:


