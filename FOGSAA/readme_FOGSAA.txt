FOGSAA---> Fast Optimal Global Sequence Alignment Algorithm
------------------------------------------------------------

                Installation/Build/Run Instruction


1. WINDOWS: Install any c++ compiler( 64 bit or 32 bit according to your machine specification). Compile the source code file fogsaa.cpp and build the .exe file. Now run the .exe file from the command prompt with suitable parameter as listed below.


2. LINUX: From the LINUX terminal go to the path where you have downloaded the source file fogsaa.cpp. Now type "g++ fogsaa.cpp -o fogsaa". This will compile and build the executable linux binary suitable for your machine. Now run fogsaa with specific parameters as listed below. For example type  ./fogsaa seq1.txt seq2.txt 1 0 1 -1 -2  to run FOGSAA for DNA sequences where the sequences are in seq1.txt and seq2.txt file respectively. The detailed description of the parameters can be found below.
Note that, FOGSAA takes input sequences in FASTA file format. 


// FOGSAA algorithm implementation with affine gap penalty
//               Developed  By Angana Chakraborty



// -------------------------------------------
 // ############      DNA/Gene Sequences    ###############
/*-----------------------------------------
Command line parameters (for alignments without affine gap penalty)
./fogsaa f1 f2 d 0 mt mst gp
---- where f1 --> name of the 1st file
           f2 --> name of the 2nd file
           d  --> 1 for DNA sequences 
           0  --> as without affine gap penalty  
           mt --> score of a match
           mst--> score for mismatch
           gp --> gap penalty
    example ./fogsaa seq1.txt seq2.txt 1 0 1 -1 -2
-----------------------------------------------------------
-----------------------------------------------------------
-----------------------------------------------------------
Command line parameters ( for alignments with affine gap penalty)
./fogsaa f1 f2 d 1 mt mst go ge 
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 1 for DNA sequences 
           1  --> as with affine gap penalty 
           mt --> score of a match
           mst--> score for mismatch
           go --> gap open penalty
           ge --> gap extension penalty, Total penalty for a gap of length L= (go+L*ge)
    example ./fogsaa seq1.txt seq2.txt 1 1 1 -1 -10 -1
--------------------------------------------------------------------
--------------------------------------------------------------------

###################    Protein Sequences    #######################
-------------------------------------------------------------------
Command line parameters (for alignments without affine gap penalty)
./fogsaa f1 f2 d 0 gp
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 2 for Protein sequences 
           0  --> as without affine gap penalty 
           gp --> gap penalty
    example ./fogsaa seq1.txt seq2.txt 2 0 -2
    protein score matrix is BLOSUM64 stored on file fscore.txt
-------------------------------------------------------------------
-------------------------------------------------------------------
-------------------------------------------------------------------
Command line parameters ( for alignments with affine gap penalty)
./fogsaa f1 f2 d 1 go ge 
---- where f1 --> name of the 1st file.
           f2 --> name of the 2nd file
           d  --> 2 for Protein sequences 
           1  --> as with affine gap penalty 
           go --> gap open penalty
           ge --> gap extension penalty, Total penalty for a gap of length L= (go+L*ge)
    example ./fogsaa seq1.txt seq2.txt 2 1 -10 -1
    protein score matrix is BLOSUM64 (default) stored on file fscore.txt
---------------------------------------------------------------------------
###########################################################################
###########################################################################
*/