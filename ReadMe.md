## Software Enviroment Requirement
 - Since the program was written by matlab language, the MATLAT software must be installed firstly.
 - dRHP-GreyFun uses the following dependent software: BLAST and HMMER. You can install the two software and database from 
    - BLAST, ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    - swissport ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
    - [HMMER](http://www.hmmer.org/)
    
## Useage
>[family, prob] = predictProteinRemoteHomology('query.fasta');

		-query.fasta:
			a string of fasta file including query proteins
		-family:
			predicting remote homology family name
		-prob:
			predicting probability
			
