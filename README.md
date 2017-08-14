### UNALIGNED_GRABBER v-1.0

#### AUTHOR INFO: 
Shruti Srivastava
<shruti.srivastava@ucalgary.ca>

This program extracts the coordinates of the
unaligned regions for a query after BLASTN. 
This could be helpful in exploring differences
between two species at the sequence level.

#### Condition: Unaligned regions should be greater
than a minimum length provided by the user.

#### Inputs: 
1) A BLAST file (tsv) with columns in this order:
Query_ID, Query_Length, Query_start, Query_end

For Example - Unaligned_output.csv

|Query|start|end|
| ------------- |:-------------:| -----:|
|SEQ-D|1|14|
|SEQ-D|85|99|
|SEQ-C|23|74|
|SEQ-C|77|100|
|SEQ-B|23|44|


2) Minimum unaligned length

[Default : 14] 


#### Output: A csv file with unaligned coordinates for 
the query sequence.

[Default : Unaligned_output.csv]


#### USAGE:

python unaligned_grabber.py -h[HELP]

python unaligned_grabber.py -i <inputBLASTfile> -o <outputfile> -l <minLength>

For example - 

python Unaligned_Grabber.py -i blastn_coordinates.csv -o Unaligned_output.csv
