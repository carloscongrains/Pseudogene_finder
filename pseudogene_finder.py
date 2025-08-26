#!/usr/bin/env python

"""
Check for stop codons and/or frameshift mutations
"""

from Bio import SeqIO
import argparse,os,sys
import subprocess
from collections import Counter

class MyParser(argparse.ArgumentParser):
	def error(self, message):
		sys.stderr.write('error: %s\n' % message)
		self.print_help()
		sys.exit(2)

parser=MyParser()
#parser = argparse.ArgumentParser()
parser.add_argument('--input_fasta', help='A multifastafile with nucleotide coding sequences to be analysed.')
parser.add_argument('--input_reference', help='A fasta file with the reference nucleotide sequence of the target gene. This sequence must start with the first position of the codon and end in an stop codon (complete coding sequence). This will be used to put the query sequences in the correct frame.')
parser.add_argument('--genetic_code_table', help='Provide a number for the correct genetic code for your sequences based on https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG5. Genetic code tables available in this script are: 1 to 6, and 9 to 16.')
parser.add_argument('--output_dir', help='Path of the outputdir.')

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)

args = parser.parse_args()


if args.input_fasta:
	input_fasta = args.input_fasta

if args.input_reference:
	input_reference = args.input_reference

if args.output_dir:
	output_dir = args.output_dir

if args.genetic_code_table:
	genetic_code_table = int(args.genetic_code_table)
"""
FUNCTIONS
"""

#Make a folder
def create_folder(folder_path):
	if not os.path.exists(folder_path):
		try:
			os.makedirs(folder_path)
		except OSError as e:
			print(f"Error creating folder '{folder_path}': {e}")
	else:
		print(f"Folder '{folder_path}' already exists.")
	return

#Convert relative to full path
def convert_to_full_path(path):
	#Get working directory
	absolute_path = os.path.abspath(path)
	return absolute_path

#Get stop codons according to the selected genetic code table
def get_stops(genetic_code_table):
#Information from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG5
	if genetic_code_table == 1:
		#Universal
		gene_code_1 = {"TAA":"STOP","TAG":"STOP","TGA":"STOP","ATG":"START"}
		return gene_code_1
	elif genetic_code_table == 2:
		#Vertebrate Mitochondrial
		gene_code_2 = {"TAA":"STOP","TAG":"STOP","AGA":"STOP","AGG":"STOP","ATT":"START","ATC":"START","ATA":"START","ATG":"START","GTG":"START"}
		return gene_code_2
	elif genetic_code_table == 3:
		#Yeast Mitochondrial
		gene_code_3 = {"TAA":"STOP","TAG":"STOP","ATA":"START","ATG":"START","GTG":"START"}
		return gene_code_3
	elif genetic_code_table == 4:
		#Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		gene_code_4 = {"TAA":"STOP","TAG":"STOP","TTA":"START","TTG":"START","CTG":"START","ATT":"START","ATC":"START","ATA":"START","ATG":"START","GTG":"START"}
		return gene_code_4
	elif genetic_code_table == 5:
		#Invertebrate Mitochondrial 
		gene_code_5 = {"TAA":"STOP","TAG":"STOP","TTG":"START","ATT":"START","ATC":"START","ATA":"START","ATG":"START","GTG":"START"}
		return gene_code_5
	elif genetic_code_table == 6:
		#Ciliate, Dasycladacean and Hexamita Nuclear Code
		gene_code_6 = {"TGA":"STOP","ATG":"START"}
		return gene_code_6
	elif genetic_code_table == 9:
		#Echinoderm and Flatworm Mitochondrial Code 
		gene_code_9 = {"TAA":"STOP","TAG":"STOP","ATG":"START","GTG":"START"}
		return gene_code_9
	elif genetic_code_table == 10:
		#Euplotid Nuclear Code
		gene_code_10 = {"TAA":"STOP","TAG":"STOP","ATG":"START"}
		return gene_code_10
	elif genetic_code_table == 11:
		#Bacterial, Archaeal and Plant Plastid Code 
		gene_code_11 = {"TAA":"STOP","TAG":"STOP","TGA":"STOP","TTG":"START","CTG":"START","ATT":"START","ATC":"START","ATA":"START","ATG":"START","GTG":"START"}
		return gene_code_11
	elif genetic_code_table == 12:
		#Alternative Yeast Nuclear Code
		gene_code_12 = {"TAA":"STOP","TAG":"STOP","TGA":"STOP","CTG":"START","ATG":"START"}
		return gene_code_12
	elif genetic_code_table == 13:
		#Ascidian Mitochondrial Code
		gene_code_13 = {"TAA":"STOP","TAG":"STOP","TTG":"START","ATA":"START","ATG":"START","GTG":"START"}
		return gene_code_13
	elif genetic_code_table == 14:
		#Alternative Flatworm Mitochondrial Code
		gene_code_14 = {"TAG":"STOP","ATG":"START"}
		return gene_code_14
	elif genetic_code_table == 15:
		#Blepharisma Nuclear Code 
		gene_code_15 = {"TAA":"STOP","TGA":"STOP","ATG":"START"}
		return gene_code_15
	elif genetic_code_table == 16:
		#Chlorophycean Mitochondrial Code (transl_table=16)
		gene_code_16 = {"TAA":"STOP","TGA":"STOP","ATG":"START"}
		return gene_code_16
	else:
		#print("Check genetic code table.")
		return "Check"

def make_pair_aligments(input_fasta,input_reference,output_dir):
	set_alignment_files = []
	#Get the record of the reference sequence
	for reference_record in SeqIO.parse(input_reference, "fasta"):
		reference = reference_record
		break
	#Iterate through query multifasta file
	for query in SeqIO.parse(input_fasta, "fasta"):
		id = query.id
		filename = id.split()[0] + ".fas"
		filepath = os.path.join(output_dir,filename)
		#Save a file with reference and query
		with open(filepath,"w") as tmp_output:
			tmp_output.write(">" + reference.id + "\n" + str(reference.seq) + "\n" + ">" + query.id + "\n" + str(query.seq) + "\n")
		#Make the alignment
		output = id.split()[0] + "_aln" + ".fas"
		output_path = os.path.join(output_dir,output)
		output_path = alignment(filepath,output_path)
		#Include query id, reference, path of the output
		information = [id,reference.id,output_path]
		set_alignment_files.append(information)
	return set_alignment_files

#Remove gaps
def remove_gaps(sequence):
	clean_sequence = []
	for base in sequence:
		if base != "-":
			clean_sequence.append(base)

	return "".join(clean_sequence)


#CHeck if the reference is ok. No premature stop. Initial and stop codon accoriding to the genetic code
def check_reference(input_reference,genetic_code_table):
	count = 0
	for record in SeqIO.parse(input_reference, "fasta"):
		count = count + 1
		assert  count == 1, "Reference contains more than one sequence. Only one reference is allowed."
		sequence = str(record.seq).upper()
		sequence = remove_gaps(sequence)
		assert len(sequence)%3 == 0, "Reference not multiple of 3. Check the reference."
		#First_codon = sequence[:3]
		#assert First_codon in genetic_code_table.keys() and genetic_code_table[First_codon] == "START","Reference does not start with a start codon. Check the reference."
		Last_codon = sequence[-3:]
		if Last_codon in genetic_code_table.keys() and genetic_code_table[Last_codon] == "STOP":
			return [record.id,str(len(sequence)),str(1),Last_codon,str(len(sequence)-2)]
		#assert Last_codon in genetic_code_table.keys() and genetic_code_table[Last_codon] == "STOP", "Reference does not end with a stop codon. Check the reference."
		else:
			return [record.id,str(len(sequence)),str(1),"-","-"]

#Get stop codons from sequence:
def get_stop_codon_results(codon_pos,position,genetic_code_table,query_id,sequence,results_stop_codon):	
	stop_codons = []
	position_stop_codon = []
	pos = position
	#First position of the first complete codon
	if codon_pos == 1:
		first_pos_codon = pos
	elif codon_pos == 2:
		first_pos_codon = pos + 2
	elif codon_pos == 3:
		first_pos_codon = pos + 1

	#Iterate the codons to find the stop condons
	for new_position in range(first_pos_codon,len(sequence),3):
		if new_position + 1 < len(sequence) and new_position + 2 < len(sequence):
			codon = sequence[new_position] + sequence[new_position+1] + sequence[new_position+2]
			if codon in genetic_code_table.keys():
				if genetic_code_table[codon] == "STOP":
					stop_codons.append(codon)
					position_stop_codon.append(str(new_position+1))
	sequence = remove_gaps(sequence)
	length_query = str(len(sequence))
	
	#Results consists of position in the alignment where the query sequence starts, length of the alignment, length of the query, postion of the codon (frame), Stop codons, Positions of stop codons.
	if len(position_stop_codon) == 0:
		results_stop_codon[query_id] = [pos,len(sequence),length_query,str(codon_pos),"-","-"]
	else:
		results_stop_codon[query_id] = [pos,len(sequence),length_query,str(codon_pos),", ".join(stop_codons),", ".join(position_stop_codon)]

	return results_stop_codon

#Get position of the 
def process_pos_stop_codon(sequence,genetic_code_table,query_id,results_stop_codon,final_pos):
		
	for pos in range(len(sequence)):
		if sequence[pos] != "-":
			#Determine the codon position and results from codon
			
			if (pos + 1)%3 == 1:
				codon_pos = 1
				results_stop_codon = get_stop_codon_results(codon_pos,pos,genetic_code_table,query_id,sequence,results_stop_codon)
			elif (pos + 1)%3 == 2:
				codon_pos = 2
				results_stop_codon = get_stop_codon_results(codon_pos,pos,genetic_code_table,query_id,sequence,results_stop_codon)
			elif (pos + 1)%3 == 0:
				codon_pos = 3
				results_stop_codon = get_stop_codon_results(codon_pos,pos,genetic_code_table,query_id,sequence,results_stop_codon)
			break
	#Get final position
	count = 0
	for a in sequence[::-1]:
		if a != "-":
			final_position = len(sequence) - count

			final_pos[query_id] = final_position
			break
		else:
			count = count +1

	return results_stop_codon,final_pos

#Identify stop codons (premature sequence)
def get_premature_stop(set_alignment_files,genetic_code_table):
	results_stop_codon = {}
	final_pos = {}
	for info in set_alignment_files:
		query_id = info[0]
		reference_id = info[1]
		align_path = info[2]
		for record in SeqIO.parse(align_path, "fasta"):
			if record.id == query_id:
				record_query = record
				sequence = str(record_query.seq).upper()
				#Determine the codon position
				results_stop_codon,final_pos = process_pos_stop_codon(sequence,genetic_code_table,query_id,results_stop_codon,final_pos)
				break	
	return results_stop_codon,final_pos

#Get position and length of indels
def get_pos_len(indels_list,indels_dict):
    if not indels_list:
        return {}

    start_of_sequence = indels_list[0]
    current_length = 1

    # Iterate through the list starting from the second element
    for i in range(1, len(indels_list)):
        # Check if the current number is consecutive to the previous one
        if indels_list[i] == indels_list[i-1] + 1:
            current_length += 1
        else:
            # If not consecutive, store the previous sequence and reset
            indels_dict[start_of_sequence] = current_length
            start_of_sequence = indels_list[i]
            current_length = 1

    # Add the last sequence to the result after the loop ends
    indels_dict[start_of_sequence] = current_length

    return indels_dict



#Get truncated sequences
def get_seq_new_frameshift(query_id,align_path,results_stop_codon,final_pos,warning_message_indel_length_dict):
	position_list = []
	position = results_stop_codon[query_id][0]
	alignment_len = results_stop_codon[query_id][1]
	codon_pos = int(results_stop_codon[query_id][3])
	#Position and length
	complete_indels_list = []
	complete_indels_dict = {}
	frameshift_indels_list =[]
	#Modify final position according to the codon position
	if (final_pos)%3 == 1:
		final_pos = final_pos -1
	elif (final_pos)%3 == 2:
		final_pos = final_pos -2
	#First position of the first complete codon
	if codon_pos == 1:
		first_pos_codon = position
	elif codon_pos == 2:
		first_pos_codon = position + 2
	elif codon_pos == 3:
		first_pos_codon = position + 1

	#Populate the list with positions of indels
	for record in SeqIO.parse(align_path, "fasta"):
		for new_position in range(first_pos_codon,final_pos):
			sequence = str(record.seq).upper()
			if sequence[new_position] == "-":
				complete_indels_list.append(new_position)
	
	#Complete_indels_dict is a dictionary which keys are position of indels and length of indels
	complete_indels_dict = get_pos_len(complete_indels_list,complete_indels_dict)
	
	#Check the indel dict to get the position of frameshift indels
	for indel_pos in complete_indels_dict.keys():
		indel_length = complete_indels_dict[indel_pos]
		if indel_length%3 == 0 and indel_pos%3 != 0:
			frameshift_indels_list.append(str(indel_pos + 1))
		elif indel_length%3 != 0:
			frameshift_indels_list.append(str(indel_pos + 1))
	
	#Get indel lengths
	indel_lengths = []
	for indel_pos in complete_indels_dict.keys():
		indel_lengths.append(complete_indels_dict[indel_pos])
	#Get the maximum indel length
	if len(indel_lengths) > 0: 
		maximum_indel_length = max(indel_lengths)
		#Print warning message if a indel longer than 6 is found
		if maximum_indel_length > 6:
			warning_message_indel_length_dict[query_id] = "Warning: Indel of " + str(maximum_indel_length) + "bp. Check the alignment."


	#Adding the the new frameshift (indels not multiple of three)
	#Results consists of position in the alignment where the query sequence starts, length of the alignment, length of the query, postion of the codon (frame), Stop codons, Positions of stop codons, postion of gaps not multiple of three.
	if len(frameshift_indels_list) == 0:
		results_stop_codon[query_id].append("-")
	else:
		results_stop_codon[query_id].append(", ".join(frameshift_indels_list))

	return results_stop_codon,warning_message_indel_length_dict

#Save the results
def save_results(results_stop_codon,warning_message_indel_length_dict,reference_info,output_dir):
	#Make a header
	header = ["ID", "Length_sequence","Frame_position", "Stop_codons","Position_stop_codons","Position_frameshift","Conclusion","Notes"]
	#
	output_path = os.path.join(output_dir,"output_results.txt")
	
	#Get Conclusion
	stop_pos = reference_info[4]

	for query_id in results_stop_codon.keys():
		conclusion = "Gene"
		results = results_stop_codon[query_id]
		if results[6] != "-":
			conclusion = "Pseudogene"
		else:
			position_list_stops = results[5].split(", ")
			if len(position_list_stops) > 1:
				conclusion = "Pseudogene"
			else:
				if position_list_stops[0] != "-" and position_list_stops[0] != stop_pos:
					conclusion = "Pseudogene"

		results_stop_codon[query_id].append(conclusion)

		if query_id in warning_message_indel_length_dict.keys():
			results_stop_codon[query_id].append(warning_message_indel_length_dict[query_id])
		else:
			results_stop_codon[query_id].append("-")

	with open(output_path, "w") as output:
		output.write("\t".join(header) + "\n")
		output.write("\t".join(reference_info) + "\t" + "-" + "\t" + "Reference" + "\t" + "-" + "\n")
		for query_id in results_stop_codon.keys():
			results = results_stop_codon[query_id]
			output.write(query_id + "\t" + ("\t").join(results[2:]) + "\n")
	return

#Identify frameshift mutations

#Function for making the alignment
def alignment(fasta_path,output_path):

	output_dir = os.path.dirname(output_path)
	infile = fasta_path
	outfile = output_path
	with open(outfile,"w") as out_handle:
		cmd = ["mafft", "--genafpair","--ep","0","--maxiterate","1000","--thread","1","--quiet","--preservecase", infile]
		p = subprocess.Popen(cmd, text=True, stdout=out_handle,stderr=subprocess.PIPE, cwd=output_dir)
		stderr = p.communicate()[1]
	if p.returncode == 0:
		print("Command output redirected to output.txt")
	else:
		print("Command failed with error code:", p.returncode)
		print("Stderr:", stderr)

	return output_path


"""
MAIN
"""
input_fasta = convert_to_full_path(input_fasta)
input_reference = convert_to_full_path(input_reference)
output_dir = convert_to_full_path(output_dir)

#Make the output folder
create_folder(output_dir)

#Get genetic code
current_genetic_code_table = get_stops(genetic_code_table)
assert current_genetic_code_table != "Check", "Check genetic code table."
#Check reference and prepare reference to save on the results
reference_info = check_reference(input_reference,current_genetic_code_table)
#Make the alignments
set_alignment_files = make_pair_aligments(input_fasta,input_reference,output_dir)
#Get stop codons
results_stop_codon,final_pos = get_premature_stop(set_alignment_files,current_genetic_code_table)
#Get frameshifts
warning_message_indel_length_dict = {}
for info_alignment in set_alignment_files:
	query_id = info_alignment[0]
	alignment = info_alignment[2]
	final_position = final_pos[query_id]
	results_stop_codon,warning_message_indel_length_dict = get_seq_new_frameshift(query_id,alignment,results_stop_codon,final_position,warning_message_indel_length_dict)

#Save results
save_results(results_stop_codon,warning_message_indel_length_dict,reference_info,output_dir)
