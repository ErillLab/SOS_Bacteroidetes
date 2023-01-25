# -*- coding: utf-8 -*-

# @author: miquelsanchezosuna

### MEME needs to be installed locally
### Implemented for Python3

# Load the necessary modules

import json
import os
import subprocess
import time
from Bio import  Entrez, pairwise2, Seq, SearchIO
from Bio.SeqRecord import SeqRecord
from io import StringIO


def split_accession(acc):
    
    """ Splits a genbank/refseq record into prefix (e.g. NZ_MDWE or MDWE) and
        numeric component (e.g. 0100092), and returns them as a list.
    """

    #go through the record until we get the first numeric digit
    #store the prefix and use its length to get the numeric region
    prefix = ''
    for c in acc:
        if c.isdigit():
            break
        else:
            prefix=prefix + c

    suffix = acc[len(prefix):len(acc)]

    return [prefix, suffix]


def genome_record_retrieval(ortholog_acc, sleepy):
    
    """ Takes a protein accession as an input. Retrieves its IPG record.

        The idea here is to obtain, prioritarily, data from complete genome 
        records if they exist, from RefSeq (AC_ and NC_ accessions) . If no 
        RefSeq is available, then select complete genome records from GenBank 
        (AE, CP, CY accessions). Otherwise, select contigs or WGS scaffolds from 
        RefSeq (NT_, NW_, NZ_). If that fails, get contigs or WGS scaffolds from 
        GenBank (AAAA-AZZZ). Only when nothing else is available, select direct 
        GenBank submissions (U, AF, AY, DQ).

        Prioritizes each type of accession and returns the best record.

        Priority indices range from 7 (best, for a complete RefSeq record) and 6
        (complete GenBank record), to 5 (for complete RefSeq WGS) and all the 
        way to 3 (undetermined GenBank records).
        
        It returns a composite record with the nucleotide accession number for
        the "best" coding region, and the position and orientation of the CDS 
        within that accession, as well as the prioritization score obtained.
    """

    #Download IPG record for the specific ortholog
    records = None
    while records == None:
        try:
            records = Entrez.read(Entrez.efetch(db="protein", id=ortholog_acc, \
                                        rettype='ipg', retmode='xml'))
        except:
            print ("IPG record retrieval error: "+ortholog_acc)
            time.sleep(sleepy)
            pass
    
    #create scoring for priorization
    priority = {"NC_": 7, "AC_": 7, "AE": 6, "CP": 6, "CY": 6, \
                    "NZ_": 5, "NT_": 5, "NW_": 5, "AAAA-AZZZ": 4,\
                    "U": 3, "AF": 3, "AY": 3, "DQ": 3}
    
    #from the IPG record, retrieve all the genome accessions from all CDS
    #keeping only accession, location of start and strand, as well as 
    #priority score
    genomelist = []
    if 'ProteinList' in records["IPGReport"].keys():
        for idprotein in records["IPGReport"]["ProteinList"]:
            if 'CDSList' in idprotein.keys():
                for cds in idprotein['CDSList']:
                    cds_acc = cds.attributes['accver']
                    cds_start = cds.attributes['start']
                    cds_stop = cds.attributes['stop']
                    cds_strand = cds.attributes['strand']
                    cds_org = cds.attributes['org'].replace(" ","_")
                    cds_scr = 0
                    
                    #assign priority
                    for key in priority:
                        if cds_acc.startswith(key):
                            cds_scr = priority[key]
                            
                    #special case: NZ_ records that map to "complete" WGS
                    #these are NZ_ records with a short numeral
                    #e.g. NZ_CP030158
                    #the idea is to bump them up to 6, so they
                    #are treated as complete
                    if cds_acc.startswith('NZ_'):
                        pr,sf = split_accession(cds_acc)
                        #non-complete WGS acc seem to have alwayas 8 numericals
                        #so we pick as "complete" anything with less than 7
                        if len(sf.split('.')[0])<7:
                            cds_scr = 6.5
                    
                    #create and append record
                    cds_rec = {'acc':cds_acc, 'start':cds_start, \
                               'stop':cds_stop, 'strand':cds_strand,\
                               'p_score':cds_scr, 'org':cds_org}
                    genomelist.append(cds_rec)
            else:
                return (None)
    #GenBank record that has no proper IPG record (yes, they exist;
    #see for instance: https://www.ncbi.nlm.nih.gov/protein/RJR51       119.1)
    #in these cases, there should be a CDS within the protein record that 
    #contains the information we want; priority should be lowest
    else:
#        TO BE IMPLEMENTED
#        records = Entrez.read(Entrez.efetch(db="protein", id=ortholog_acc, \
#                                            rettype='genbank', retmode='xml'))
        return (None)
        
    #select the genomes with highest value (in case of same scores, 
    #first one is chosen (random)
    max_record = genomelist[0]
    for genome in genomelist:
        if genome['p_score'] >= max_record['p_score']:
            max_record = genome     
    return max_record
    
    time.sleep(sleepy)


def promoter_retrieval(nucid,start,stop,strand,start_adj=250,stop_adj=2, min_size = 50):
    
    """Download promoter sequences (only intergenic) from NCBI RefSeq/GenBank 
       database using the RefSeq/GenBank ID; and the coordinates & strand of 
       gene of interest.
       The start_adj and stop_adj delimit the promoter size.
       Returns only the seq of SeqRecord object
    """
    
    if strand == "+":
        s_start=int(start)-start_adj
        s_stop=int(start)+stop_adj
        s_strand=1
    else:
        s_stop=int(stop)+start_adj
        s_start=int(stop)-stop_adj
        s_strand=2
    
    #Fetch and read the annotated GenBank record
    handle = None
    while handle == None:
        try: 
            handle = Entrez.efetch(db="nuccore",id=nucid, strand=s_strand,seq_start=s_start, seq_stop=s_stop, rettype='gb', retmode="xml")
            genome_record = Entrez.read(handle, "xml")
        except:
            time.sleep(5)
            print ("Promoter retrieval error: "+nucid)
            pass
    
    # Find all coding regions in the returned GenBank sequence. 
    coding_intervals = []
    sequence = genome_record[0]['GBSeq_sequence']
    for feature in genome_record[0]['GBSeq_feature-table']:
        if feature['GBFeature_key'] == "CDS": 
            if "GBInterval_from" in feature['GBFeature_intervals'][0]:
                coding_start = feature['GBFeature_intervals'][0]['GBInterval_from']
                coding_end = feature['GBFeature_intervals'][0]['GBInterval_to']
                coding_intervals.append((coding_start, coding_end))  
            
    #The FASTA ID for the promoter sequence is in the following format:
        # p_NucleotideRecord
    return_id = "p_" + str(nucid)
    #If there is only one coding region in the selected sequence, then 
    # the sequence is returned unmodified. 
    if len(coding_intervals) == 1:
        #Appends information to record description
        return SeqRecord(Seq.Seq(sequence), id= return_id, description = return_id+"_"+str(s_start)+"_"+str(s_stop))
    #If no coding intervals are indentified, None is returned.
    elif len(coding_intervals) == 0:
        return None
    #The start of the promoter is set to the start/end of the upstream gene
    # based on the directionality. ( --> --> or <-- -->)
    if s_strand == 1:
        promoter_start = max(int(coding_intervals[-2][0]), 
                             int(coding_intervals[-2][1]))
    else:
        promoter_start = max(int(coding_intervals[1][0]), 
                     int(coding_intervals[1][1]))
    #Everything upstream of the promoter start is clipped off the 
    # sequence and the substring is returned.
    return_seq = str(sequence[promoter_start : ])
    #Appends information to record description
    if len(return_seq) >= min_size:
        return SeqRecord(Seq.Seq(return_seq), id= return_id, description = return_id+"_"+str(s_start)+"_"+str(s_stop))


def id_below_maxid_perc(el1, el2, max_percent_id):
    
    """ Aligns two sequence elements and determines whether they are
        more than %ID identical (false) or not (true)
        Scoring: Match:+2, Mismatch:-1, GapO: -2, GapE: -0.2
    """
    
#    print "Seq1 length: " + str(len(el1.seq)) + ' ' + el1.id
#    print "Seq2 length: " + str(len(el2.seq)) + ' ' + el2.id
    al=pairwise2.align.globalms(el1.seq, el2.seq, 2, 0, -2, -.5,\
                                one_alignment_only=True, \
                                penalize_end_gaps=False)
    
    #print al
    matches=0
    gapless=0
    #for each position in the alignment
    for ch_pair in zip(al[0][0],al[0][1]):
        #if this is a non-gapped position
        if '-' not in ch_pair:
            #if it's a match, count it
            if ch_pair[0]==ch_pair[1]:
                matches=matches+1
            gapless=gapless+1
        
    perID = float(matches)/float(gapless)
    
#    print "Matches: ", matches
#    print "Gapless: ", gapless
#    print "%ID: ", perID
    
    #return true or false depending on percent identity
    if perID*100<=float(max_percent_id):
        return(True)
    else:
        return(False)


def remove_redundants(promoters_list, max_percent_identity):
    
    """ Takes a list of promoter sequences and returns a list of promoters with
        <= max_percent_identity.
        The identitiy values are computed using id_below_maxid_perc function defined
        above
    """
    
    removed_ids = []
    nr_promoters = []
    
    for i in range(0, len(promoters_list)-1):
        if promoters_list[i].description not in removed_ids:
            nr_promoters.append(promoters_list[i])
            for j in range(i+1, len(promoters_list)):
                if id_below_maxid_perc(promoters_list[i],promoters_list[j],max_percent_identity) == False:
                    removed_ids.append(promoters_list[j].description)
    
    return nr_promoters
                    

def memesearch(fastafile, nmotifs, minw, maxw, maxsites):
    
    """
        Takes a FASTA file and MEME parameters to perform the motif discovery by
        running MEME locally
    """
    
    bashCommand = "meme "+fastafile+" -dna -nostatus -time 18000 -mod anr -nmotifs "+str(nmotifs)+" -minw "+str(minw)+" -maxw "+str(maxw)+" -revcomp -maxsites "+str(maxsites)+" -evt 1e-30 -pal -oc "+fastafile.split(".")[0]
    subprocess.call(bashCommand.split())
    

#######################################################

# Load the configuration file

with open("test_input.json") as json_conf : 
    conf = json.load(json_conf)

# Tell to who I am
Entrez.email = conf["email"]
Entrez.api_key = conf["api_key"]

#Open a file containing desired protein ids
filename = conf["cog_file"]
with open(filename) as f:
    query_ids = f.readlines()
query_ids = [x.strip() for x in query_ids] 

#Create a dictionary using query_ids hits as keys and their location, species,
#and priority scores as values
ortholog_summary_dict = {}
for acc in query_ids:
	genome_record_retrieval_result = genome_record_retrieval(acc, 2)
	if genome_record_retrieval_result != None:
	    ortholog_summary_dict[acc] = genome_record_retrieval_result

#Add the promoter information into the previous ortholog_summary_dict
for protein_id in ortholog_summary_dict.keys():
    current_promoter = promoter_retrieval(ortholog_summary_dict[protein_id]['acc'] ,ortholog_summary_dict[protein_id]['start'],ortholog_summary_dict[protein_id]['stop'],ortholog_summary_dict[protein_id]['strand'])
    if current_promoter:
        ortholog_summary_dict[protein_id]['prom'] = current_promoter

#create a directory to save MEME results
bashCommand = "mkdir ~/MEME_Output"
subprocess.call(bashCommand.split())

#download the promoters and remove redundant promoters using % identity
promoters = []
for key in ortholog_summary_dict.keys():
    try:
        promoters.append(ortholog_summary_dict[key.split("__")[0]]['prom'])
    except KeyError:
        continue
    
promoters_nr = remove_redundants(promoters,conf["nr_identity"])
#MEME discovery only in nodes with more than 10 non-redundant promoters
if len(promoters_nr) >= 10:
    #save the non-redundant promoters into a FASTA file
    out_handle_F = open("~/MEME_Output/nr_promoters.fasta", "w")
    for seqrecord in promoters_nr:
        out_handle_F.write(">"+"_".join(seqrecord.description.split(" ")[1:])+"_"+seqrecord.description.split(" ")[0]+"\n")
        out_handle_F.write(str(seqrecord.seq)+"\n")
    out_handle_F.close()
    #motif discovery with MEME
    memesearch("~/MEME_Output/nr_promoters.fasta",10,conf["meme_minsize"],conf["meme_maxsize"],1000)
 
