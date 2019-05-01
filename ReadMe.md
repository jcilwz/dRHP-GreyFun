## Software Enviroment Requirement
 -Since the program was written by matlab language, the MATLAT software must be installed firstly.
 -dRHP-GreyFun uses the following dependent software: BLAST and HMMER. You can install the two software from 
    -[https://blast.ncbi.nlm.nih.gov/Blast.cgi]BLAST 
    -[http://www.hmmer.org/]HMMER
    
## Useage
=======================================
>[family, prob] = predictProteinRemoteHomology('query.fasta');

		-query.fasta:
			a string of fasta file including query proteins
		-family:
			predicting remote homology family name
		-prob:
			predicting probability
			
