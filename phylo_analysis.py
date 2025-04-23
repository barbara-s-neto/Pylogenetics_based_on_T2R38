#!/usr/bin/env python3

import sys
import argparse
from Bio import AlignIO

from Bio import Entrez, SeqIO, Phylo
from Bio.Blast import NCBIWWW, NCBIXML

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from Bio.Align.Applications import ClustalOmegaCommandline

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import os

from matplotlib.colors import to_rgba

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

    with open("helper_files/blast_results.xml", "w") as out_file:
        out_file.write(result_handle.read())

    result_handle.close()



def read_blast():
    with open("helper_files/blast_results.xml") as result_file:
        blast_records = NCBIXML.read(result_file)

    return blast_records

Entrez.email = "joaomonteiro1111@gmail.com"

def get_filtered_sequences(blast_records, seq_input, n=10):
    unique_names = set()
    i = 0

    with open("helper_files/filtered_sequences.fasta", "w") as fasta_file:
        fasta_file.write(f">{seq_input.id}\n{seq_input.seq}\n\n")

        for alignment in blast_records.alignments:
            if i == n:
                break

            alignment_str = str(alignment)
            idx_1 = alignment_str.find("[")
            idx_2 = alignment_str.find("]")

            if idx_1 != -1 and idx_2 != -1:  # Ensure both brackets exist
                full_name = alignment_str[idx_1 + 1: idx_2]  # Extract text inside brackets
                name_parts = full_name.split()  # Split into words
                name = " ".join(name_parts[:2])  # Get only the first two words


            if name == "synthetic construct":
                continue

            if name not in unique_names:
                unique_names.add(name)

                accession = alignment.accession

                handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()

                sequence = seq_record.seq
                fasta_file.write(f">{name}\n{sequence}\n\n")

                # for hsp in alignment.hsps:
                #
                #     sequence = hsp.sbjct
                #
                #     fasta_file.write(f">{name}\n{sequence}\n\n")

                i += 1


def perform_msa():

    clustalomega_cline = ClustalOmegaCommandline(infile="helper_files/filtered_sequences.fasta", outfile="helper_files/aligned.fasta", verbose=True,
                                                 auto=True, force=True)
    clustalomega_cline()

def get_alignment():

    alignment = AlignIO.read("helper_files/aligned.fasta", "fasta")

    for record in alignment:
        record.id = record.description

    return alignment

def get_msa_colored_alignment(alignment):

    num_seqs = len(alignment)
    seq_len = alignment.get_alignment_length()

    # Identify only columns where there are differences
    variable_cols = [i for i in range(seq_len) if len(set(alignment[:, i])) > 1]

    # Residue color map
    residues = "ACDEFGHIKLMNPQRSTVWY"
    color_map = {res: to_rgba(plt.cm.tab20(i / 20)) for i, res in enumerate(residues)}

    # Create figure
    fig, ax = plt.subplots(figsize=(len(variable_cols) * 0.6, (num_seqs + 1) * 0.6))

    ax.text(-0.7, len(alignment) + 0.5, "ID →", ha='center', va='center', fontsize=15, fontweight='bold')
    ax.text(-3.7, len(alignment) + 0.5, "Specie ↓", ha='center', va='center', fontsize=15, fontweight='bold')

    # Plot residue letters
    for row_idx, record in enumerate(alignment):
        ax.text(-3.7, row_idx + 1, record.id, ha='center', va='center', fontsize=15, fontweight='bold')
        for col_idx, aln_col in enumerate(variable_cols):
            res = record.seq[aln_col]
            color = color_map.get(res, to_rgba("lightgray", alpha=0.5))
            ax.text(col_idx, row_idx + 1, res, ha='center', va='center', fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor=color, edgecolor='none'))

    # Plot column numbers ABOVE the letters
    for col_idx, aln_col in enumerate(variable_cols):
        ax.text(col_idx, len(alignment) + 0.5, str(aln_col + 1), ha='center', va='center', fontsize=15,
                fontweight='bold', color='black')

    # Format plot
    ax.set_xticks([])
    ax.set_yticks(range(1, num_seqs + 1))
    ax.set_yticklabels([record.id for record in alignment])
    ax.set_xlim(-0.5, len(variable_cols) - 0.5)
    ax.set_ylim(-0.5, num_seqs + 1)
    ax.axis('off')

    # Save
    plt.tight_layout()

    plt.savefig("figures/msa_colored_alignment.png", dpi=300, bbox_inches='tight')

    print("MSA Colored Alignment saved to figures/msa_colored_alignment.png")



def main():
    """Main function to run the phylogenetic analysis pipeline"""

    args = setup_arguments()

    print(f"Starting phylogenetic analysis for {args.sequence_file} with {args.num_species} species")

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


    if sequence_type == "DNA":
        sequence_dna = sequence
        sequence = pre_process(sequence)[0]


    # run_blast(sequence)  # Always search protein database
    blast_records = read_blast()

    seq_input = SeqRecord(Seq(sequence), id=name)

    get_filtered_sequences(blast_records, seq_input, n=n)


    perform_msa()

    alignment = get_alignment()

    num_seqs = len(alignment)
    seq_len = alignment.get_alignment_length()

    # Identify only columns where there are differences
    variable_cols = [i for i in range(seq_len) if len(set(alignment[:, i])) > 1]

    # Residue color map
    residues = "ACDEFGHIKLMNPQRSTVWY"
    color_map = {res: to_rgba(plt.cm.tab20(i / 20)) for i, res in enumerate(residues)}

    # Create figure
    fig, ax = plt.subplots(figsize=(len(variable_cols) * 0.6, (num_seqs + 1) * 0.6))

    ax.text(-0.7, len(alignment)+0.5,"ID →", ha='center', va='center', fontsize=15, fontweight='bold')
    ax.text(-3.7, len(alignment)+0.5,"Specie ↓", ha='center', va='center', fontsize=15, fontweight='bold')

    # Plot residue letters
    for row_idx, record in enumerate(alignment):
        ax.text(-3.7, row_idx + 1, record.id, ha='center', va='center', fontsize=15, fontweight='bold')
        for col_idx, aln_col in enumerate(variable_cols):
            res = record.seq[aln_col]
            color = color_map.get(res, to_rgba("lightgray", alpha=0.5))
            ax.text(col_idx, row_idx + 1, res, ha='center', va='center', fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor=color, edgecolor='none'))

    # Plot column numbers ABOVE the letters
    for col_idx, aln_col in enumerate(variable_cols):
        ax.text(col_idx, len(alignment)+0.5, str(aln_col + 1), ha='center', va='center', fontsize=15, fontweight='bold', color='black')

    # Format plot
    ax.set_xticks([])
    ax.set_yticks(range(1, num_seqs + 1))
    ax.set_yticklabels([record.id for record in alignment])
    ax.set_xlim(-0.5, len(variable_cols) - 0.5)
    ax.set_ylim(-0.5, num_seqs + 1)
    ax.axis('off')

    # Save
    plt.tight_layout()

    plt.savefig("figures/msa_colored_alignment.png", dpi=300, bbox_inches='tight')

    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct tree using UPGMA
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    Phylo.write(tree, "figures/phylogenetic_tree.nwk", "newick", branch_length_only=True)

    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(tree, axes=ax, branch_labels=lambda c: f"{c.branch_length:.3f}")


    plt.savefig("figures/tree_co.png", dpi=300, bbox_inches='tight')

    names = distance_matrix.names
    matrix = np.array([[distance_matrix[i, j] for j in range(len(names))] for i in range(len(names))])

    # Create a heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(matrix, annot=True, xticklabels=names, yticklabels=names, cmap="YlGnBu")
    plt.title("Sequence Distance Matrix")
    plt.tight_layout()
    plt.savefig("figures/distance_matrix_heatmap_co.png")



if __name__ == "__main__":

    main()



























