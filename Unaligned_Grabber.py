#!/usr/bin/env python
"""
UNALIGNED_GRABBER v-1.0
Author: Shruti Srivastava

This program extracts the coordinates of the
unaligned regions for a query after BLASTN. This 
could be helpful in exploring differences
between two species at the sequence level.

Condition: Unaligned regions should be greater
than a minimum length provided by the user.

Inputs: 
1) A BLAST file (tsv) with columns in this order:
Query_ID, Query_Length, Query_start, Query_end

2) Minimum unaligned length

Output: A csv file with unaligned coordinates for 
the query sequence
"""

import sys,getopt,      \
csv,time         

__title__ = 'Unaligned-Grabber'
__version__ = '1.0'
__description__ = "A tool to extract unaligned coordinates \
 with blast output file"
__author__ = 'Shruti Srivastava'
__license__ = 'The University of Calgary,'
__author_email__ = "shruti.srivastava@ucalgary.ca"
epi = "By %s, %s <%s>\n\n" % (__author__,
__license__,
__author_email__)
__doc__ = "\n***********************************************\
************************************\
\n %s v%s - %s \n**********************************\
***********************************************\
**\n%s" % (__title__,
__version__,
__description__,
epi)


def file_to_dict(inputfile):
   '''
   Reads the file and makes two dictionaries -
   1) len_dict has key = queryid, value = querylength
   2) cord_dict has key = queryid, value = a list of tuples 
                                with start and end alignent
                                coordinates for each high-
                                scoring pair 
   (file) -> (dict, dict)
                                
   Input: filename                              
   Returns: Two dictionaries
   '''
   
   len_dict = {} 
   cord_dict={}   

   #Opening inputfile in read mode
   file_reader = csv.reader(open(inputfile, 'r')) 

   #reading through the csv file line by line
   for line in file_reader: 

       #Skipping the first line
       if(line[0]=="Query"):
           continue
       
       #Storing the query name
       query_id = str(line[0])
       
       #Storing the sequence identifier as key and the length of sequence as its value.
       query_length = int(line[1]) 
       
       #Storing the start coordinates for the aligned query 
       query_start = int(line[2]) 
       
       #Storing the end coordinate for the aligned query 
       query_end = int(line[3])
       
       #Add the query length to len_dict, if absent 
       if(query_id not in len_dict):
           len_dict[query_id]=query_length        
           
       #Append the new HSP coordinates to the  cord_dict 
       if(query_id in cord_dict):           
           cord_dict[query_id].append((query_start, query_end))
      
       #Add the query coordinates to cord_dict, if absent 
       if(query_id not in cord_dict):           
           cord_dict[query_id]=[(query_start, query_end)]

   return(len_dict, cord_dict)



def map_array(len_dict, cord_dict, min_len, outputfile):
    
     '''
     Given a dictionary of query length and another 
     dictionary of aligned coordinates,  it finds 
     the unaligned coordinates greater than the 
     provided minimum length and writes it to a file 
     (dict, dict, int, string) -> int
     
     Input: length dictionary, coordinate dictionary,
            minimum length, output file
     Output: A file of unaligned coordinates
     Returns: 1 (if successfull)
     '''
     
     #Open an output file for writing
     outfile = csv.writer(open(outputfile,'w'), delimiter = ',')
     outfile.writerow(["Query","start","end"])

     #Iterate through the QueryIDs in len_dict
     #Create an empty list of 0's for each query length
     for key in len_dict:
         
         #Store an array of zeros of query length in mapper
         mapper = [0]* (len_dict[key] + 1)
         
         #Store the list of coordinate tuples for the same query in array 
         array = cord_dict[key]

         #Store the start and end from each tuple 
         for tuples in array:
             query_start = tuples[0]
             query_end = tuples[1]

             #Fill all the aligned regions in 'mapper' with 1s.
             for i in range(query_start, query_end + 1 ): 
                mapper[i] = 1
         
         #At this point, the mapper array will have 1s at the aligned part
         #and zeros in the unaligned part
         
         
         #We'll now check for regions that have a continuous stretch of zeros
         #which is greater than the minimum length provided by the user.
         start_point = 1  
         
         #Initialise a stop point
         stop_point = 0 

         #Running a for loop from 1 to the length of mapper
         for x in range(1, len(mapper)): 
             
              #This will help us in discriminating between different mapped sequences.
              previous_stop = stop_point 

              #As soon as you find a zero, ie; an unmatched point.
              if mapper[x] == 0: 

                  #Assign the start point to that index
                  start_point = x 

                  #Run a loop from that index to the end of mapping array. 
                  #We are doing this to see how long does the unmatched portion runs
                  for y in range(x, len(mapper)): 

                      #As soon as you find an index which is not equal to 0, ie; it is 1, you break the loop. 
                      #Another condition is that you have not reached the end of the mapping array
                      if mapper[y]!= 0 and y!= len(mapper)-1: 

                          #The stop point will be one minus the current value of y.
                          stop_point=y-1 
                          break

                      #If the loop has reached the end of the mapping array with 1 in the end,
                      # then stop point will be equal to y.
                      if mapper[y]!= 0 and y==len(mapper)-1: 

                          stop_point=y-1
                          break

                      #If the loop has reached the end of the mapping array with zero in the end,
                      #then stop point will be equal to y.
                      if mapper[y]== 0 and y==len(mapper)-1: 

                          stop_point=len(mapper)-1
                          break
                      
                      #Otherwise keep on filling 2 in the array, to denote that this region has been covered
                      if mapper[y]==0:
                         mapper[y]=2
           
               #Checking whether the unmapped portions are greater than or equal to minimum length. 
               #Note: If our unmapped portion runs from 1 to 14, then start_point-stop_point = 14-1=13.
              if((stop_point-start_point) >= min_len-1 and stop_point != previous_stop):
                  
                  #Writing the output to our output file
                  outfile.writerow([key, start_point, stop_point]) 
    
     return 1


def usage():
    print("\nUsage: python unaligned_grabber.py\
-i <inputBLASTfile> -o <outputfile [Default:Unaligned_output.csv] >\
-l <minLength[Default:14]>\n")
    sys.exit(2)


def main(argv):
    
   #Checking if no input has been provided 
   if(len(argv)==0):
        print('\nERROR!:No input provided\n')
        usage()

   #Default output filename
   outputfile = 'Unaligned_output.csv'

   #Deafult minimum unaligned length
   min_len = 14
   
   #Try and Catch block for handling input errors
   try:
      opts, args = getopt.getopt(argv,"h:i:o:l:",["help=","ifile=","ofile=","length="])
      
   except getopt.GetoptError:
      print(__doc__)
      usage()

   #Check whether the mandatory files are given as inputs  
   short_opts = [i[0] for i in opts]
   if(('-i') not in short_opts):
       print ("ERROR: Missing inputs. Please provide -i .")
       usage()
   
   #Reading user inputs
   for opt, arg in opts:
      if opt == '-h':
         print(__doc__)
         usage()

      elif opt in ("-i", "--ifile"):
         inputfile = arg

      elif opt in ("-o", "--ofile"):
         outputfile = arg

      elif opt in ("-l", "--length"):
         min_len = arg   
                
   print(__doc__,'Input file is %s\n Output file is %s\n Minimum length is %d' \
 %( inputfile, outputfile, min_len))
   print( '----------------------------------------\nRunning Script\n--------------\
--------------------------\n')
   
   #Variable to record time 
   start_time = time.time()
 
   #Calling function and storing the returned dicts 
   len_dict, cord_dict = file_to_dict(inputfile)
   
   #Calling function and storing the returned value
   success = map_array(len_dict, cord_dict, int(min_len), outputfile)
   
   if(success):
       print ('Done! Time elapsed: %.4f seconds' % (time.time() - start_time))



if __name__ == "__main__":
   main(sys.argv[1:])
