# countFFRs2
an improved perl script to simulate species delimitation using haplowebs

Usage: `perl countffrs2.pl -i <input file in sequential FASTA format>  -r <number of replicate samplings (default is 1)> -v <verbiose> -f <save test dataset as FASTA> -g <save test dataset as gzipped FASTA> -l generate list files -R <generate Roehl output file for Network>`

The name of each sequence in the input FASTA should be of the form Txxx_yyy, where xxx is the species ID and yyy is the sequence ID.
