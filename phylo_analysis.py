#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO, Entrez, AlignIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import time
import re

import numpy as np
import requests

from Bio import SeqIO, Entrez, Phylo
from Bio.Blast import NCBIWWW, NCBIXML

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio.Align import MultipleSeqAlignment
from Bio import Align

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import os

def get_info(filename):
    from re import sub, search

    res = []
    inf = []
    sequence = None
    info = None

    fh = open(filename)

    for line in fh:
        if search(">.*", line):
            if sequence is not None and info is not None and sequence != "":
                res.append(sequence)
                inf.append(info)
            info = line
            sequence = ""
        else:
            if sequence is None:
                return None
            else:
                sequence += sub("\s", "", line)

    if sequence is not None and info is not None and sequence != "":
        res.append(sequence)
        inf.append(info)
    fh.close()

    return (inf, res)

def read_fasta(filename):
    dic={}
    (infos,sequencias)=get_info(filename)

    for i in range(len(infos)):
        info=infos[i]
        index=info.find("[")
        species=info[index+1:-1]
        dic[species]=sequencias[i]
    return dic


def detect_sequence_type(sequence):
    """Detect if sequence is DNA or protein"""
    dna_nucleotides = set('ATGCN')
    seq_set = set(sequence.upper())

    # If sequence contains only DNA nucleotides, it's likely DNA
    if seq_set.issubset(dna_nucleotides):
        return "DNA"
    else:
        return "protein"


def setup_arguments():
    """Set up command line arguments"""
    parser = argparse.ArgumentParser(description='Phylogenetic Analysis Tool')
    parser.add_argument('sequence_file', type=str, help='Input sequence file (DNA or protein)')
    parser.add_argument('num_species', type=int, help='Number of species to include in analysis')

    return parser.parse_args()


def pre_process(dna):
    assert validate_dna(dna), "Invalid DNA sequence"
    return all_orfs_longproteins(dna)

def validate_dna (dna_seq):
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    if valid == len(seqm): return True
    else: return False


def all_orfs_longproteins(dna_seq, minsize=270):
    all = all_orfs(dna_seq)
    res = []
    for s in all:
        if len(s) > minsize:
            res.append(s)

    res2 = sorted(res, key=len)

    return res2


def all_orfs(dna_seq):
    res = []
    translations = reading_frames(dna_seq)

    for trans in translations:
        res += all_proteins_rf(trans)

    return res

def reading_frames (dna_seq):
    res = []
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res

def all_proteins_rf (aa_seq):
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def translate_seq (dna_seq, ini_pos = 0):
    seq_aa = ""
    for i in range(ini_pos, len(dna_seq)-2, 3):
        codon=dna_seq[i]+dna_seq[i+1]+dna_seq[i+2]
        seq_aa+=translate_codon(codon)
    return seq_aa


def translate_codon (cod):
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None


complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

def reverse_complement(dna_seq):
    comp = ''
    for nuc in dna_seq:
        comp += complement[nuc]

    return comp[::-1]



def run_blast(sequence, max_hits=50):
    """Run BLAST search against appropriate database"""
    print("\nRunning BLAST search. This may take a few minutes...")

    result_handle = NCBIWWW.qblast("blastp", "nr", sequence, hitlist_size=max_hits)

    # result_handle = NCBIWWW.qblast("blastp", "nr", data['Homo sapiens'])
    # breakpoint()

    with open("helper_files/blast_results.xml", "w") as out_file:
        out_file.write(result_handle.read())

    result_handle.close()
    # # Save raw BLAST results for reference
    # with open("blast_results.xml", "w") as out_handle:
    #     out_handle.write(result_handle.read())

    # with open("/home/joaomonteiro/Desktop/BIOINF/BioInf/blast_results2.xml") as result_file:
    #     blast_records = NCBIXML.read(result_file)
    #
    # return blast_records


def read_blast():
    with open("helper_files/blast_results.xml") as result_file:
        blast_records = NCBIXML.read(result_file)

    return blast_records

    # # Re-open for parsing
    # result_handle = open("blast_results.xml")
    # return NCBIXML.parse(result_handle)


def get_filtered_sequences(blast_records, n=10):
    unique_names = set()
    i = 0

    with open("helper_files/filtered_sequences.fasta", "w") as fasta_file:

        for alignment in blast_records.alignments:
            if i == n:
                break

            alignment_str = str(alignment)
            idx_1 = alignment_str.find("[")
            idx_2 = alignment_str.find("]")

            if idx_1 != -1 and idx_2 != -1:  # Ensure both brackets exist
                # breakpoint()
                full_name = alignment_str[idx_1 + 1: idx_2]  # Extract text inside brackets
                name_parts = full_name.split()  # Split into words
                name = " ".join(name_parts[:2])  # Get only the first two words

            # print(name)

            if name == "Homo sapiens" or name == "synthetic construct":
                continue

            if name not in unique_names:
                unique_names.add(name)


                for hsp in alignment.hsps:
                    # breakpoint()

                    sequence = hsp.sbjct

                    fasta_file.write(f">{name}\n{sequence}\n\n")

                i += 1

    # print(unique_names)

    sequences = read_fasta("helper_files/filtered_sequences.fasta")

    return sequences


def main():
    """Main function to run the phylogenetic analysis pipeline"""


    # Set up arguments
    args = setup_arguments()

    print(f"Starting phylogenetic analysis for {args.sequence_file} with {args.num_species} species")
    # breakpoint()

    os.makedirs("figures", exist_ok=True)
    os.makedirs("helper_files", exist_ok=True)

    def read_sequence(file_path):
        """Read sequence from file"""
        try:
            record = SeqIO.read(file_path, "fasta")
            return record
        except Exception as e:
            print(f"Error reading sequence file: {e}")
            sys.exit(1)


    record = read_sequence(args.sequence_file)
    sequence = str(record.seq)

    n = args.num_species-1

    name = record.name

    sequence_type = detect_sequence_type(sequence)
    print(f"\nDetected sequence type: {sequence_type}")

    # breakpoint()

    if sequence_type == "DNA":
        sequence_dna = sequence
        sequence = pre_process(sequence)[0]



    # print(sequence)

    # run_blast(sequence)  # Always search protein database
    blast_records = read_blast()


    filtered_sequences = get_filtered_sequences(blast_records, n=n)
    # breakpoint()

    seq_input = SeqRecord(Seq(sequence), id=name)

    # breakpoint()

    records = [SeqRecord(Seq(filtered_sequences[name]), id=name[1:]) for name in filtered_sequences]
    records.insert(0, seq_input)
    # breakpoint()

    SeqIO.write(records, "helper_files/filtered_sequences.fasta", "fasta")

    from Bio.Align.Applications import ClustalOmegaCommandline

    clustalomega_cline = ClustalOmegaCommandline(infile="helper_files/filtered_sequences.fasta", outfile="helper_files/aligned.fasta", verbose=True,
                                                 auto=True, force=True)
    clustalomega_cline()

    print(clustalomega_cline)
    breakpoint()

    from Bio import AlignIO
    alignment = AlignIO.read("helper_files/aligned.fasta", "fasta")



    for record in alignment:
        record.id = record.description[:-22]
        # print(record.id)

    # breakpoint()



    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct tree using UPGMA
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    Phylo.write(tree, "figures/phylogenetic_tree.nwk", "newick", branch_length_only=True)


    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, branch_labels=lambda c: f"{c.branch_length:.3f}")
    plt.savefig("figures/tree_co.png", dpi=300, bbox_inches='tight')




    # import PyQt5

    # from ete3 import Tree, TreeStyle
    #
    # # Read the tree from the Newick file you saved earlier
    # t = Tree("figures/phylogenetic_tree.nwk", format=1)
    # # breakpoint()
    # # Customize tree style
    # ts = TreeStyle()
    # ts.mode = "r"  # "c" for circular, "r" for rectangular
    # ts.scale = 20
    # ts.show_branch_length = True
    # ts.show_branch_support = True
    #
    # # Render and save the tree
    # t.render("ete_circular_tree.png", tree_style=ts)

    # breakpoint()

    names = distance_matrix.names
    matrix = np.array([[distance_matrix[i, j] for j in range(len(names))] for i in range(len(names))])


    # Create a heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, annot=True, xticklabels=names, yticklabels=names, cmap="YlGnBu")
    plt.title("Sequence Distance Matrix")
    plt.tight_layout()
    plt.savefig("figures/distance_matrix_heatmap_co.png")



    # aligner = Align.PairwiseAligner()
    alignment_sequential = MultipleSeqAlignment([seq_input])

    for name in filtered_sequences.keys():
        seq_specie = SeqRecord(Seq(filtered_sequences[name]), id=name[1:])
        alignment_sequential.append(seq_specie)




    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment_sequential)


    # Construct tree using UPGMA
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)


    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, branch_labels=lambda c: f"{c.branch_length:.3f}")
    plt.savefig("figures/tree_sq.png", dpi=300, bbox_inches='tight')
    # breakpoint()


    names = distance_matrix.names
    matrix = np.array([[distance_matrix[i, j] for j in range(len(names))] for i in range(len(names))])

    # Create a heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, annot=True, xticklabels=names, yticklabels=names, cmap="YlGnBu")
    plt.title("Sequence Distance Matrix")
    plt.tight_layout()
    plt.savefig("figures/distance_matrix_heatmap_sq.png")
    plt.show()


if __name__ == "__main__":

    main()



























