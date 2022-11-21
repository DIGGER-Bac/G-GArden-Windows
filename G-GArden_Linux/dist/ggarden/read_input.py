import csv
import os
from Bio import SeqIO
import pandas as pd
from datasets import res_2s

'''
Notes on the input file format:
Input file contains one fragment per line
Information bricks are "=" separated
The specific order for bricks is name=sequence(=custom_creation)(=custom_overhang)

name consists of the fragment name, all options for these are [plasmid, promoter, scaffold, seed] but they can be extended

sequence includes the sequence of the fragment (only A/C/G/T)
    primer3 utility can be enabled by including additional upstream and downstream sequences
    the functional sequence has to be marked according to fragment type
        for scaffold, the ending base is marked as e.g. (functional->) ...ATCG[C]CATCACCCTATC... (<-additional)
        for promoters, the start base is marked as e.g. (additional->) ...ATCG[C]CATCACCCTATC... (<-functional)
        for other fragments, the whole functional sequence should be marked e.g. (additional->) ...ATCG[CCAT...CACC]CTATC... (<-additional)

custom_creation is optional and can only be one of [anneal, pcr] however for primer3 remember to mark the functional sequences

custom_overhang is optional and can only consist of A/C/G/T, it defines the forward primer for the following fragment 
    e.g. promoter1=CAGCTAGCATCGACTAGC...ACTGACTAGCATCGACT=TTTT
         seed1=CATGCATCGACTACGTACGTA
    defines the overhang TTTT to be added between promoter1 and seed1 so forward primer for seed1 will have TTTT while
    reverse primer for promoter1 will have AAAA as overhang sequences
    
all added information to plasmid fragments will be ignored
in addition, scaffold forward primers cant have custom overhangs
please avoid newlines a the end of input files
'''


def snap_to_dna(file):
    dna_out = ""
    snap_out = SeqIO.parse(file, "snapgene")
    for seq in snap_out:
        seq_len = len(seq)
        for base in seq:
            dna_out += base
    if len(dna_out) == seq_len:
        return dna_out
    else:
        raise ValueError


def read_seq_from_csv(file):
    fragments = []
    file_data = pd.read_csv(file, sep=";")
    # ids = file_data["ID"]
    gene = file_data["ID"]
    seqs = file_data["Seed"]
    for idx, entry in enumerate(gene):
        fragments.append([f"{gene[idx]}", "seed", seqs[idx].upper()])
    return fragments
    # with open("data/output/test_c_data_conv.csv", "w", newline="") as csvfile:
    #     toolwriter = csv.writer(csvfile)
    #     toolwriter.writerows(fragments)


# added custom design functionality and overhangs as separate output lists
def read_in_file_csv(path):
    f = open(path, mode="r")
    csv_reader = csv.reader(f)
    allowed_base_oh = ["A", "C", "G", "T"]
    allowed_base_seq = ["A", "C", "G", "T", "[", "]"]
    allowed_frag = ["promoter", "scaffold", "seed", "dprom", "dscaf", "dseed"]
    allowed_creation = ["PCR", "Anneal", "donor"]
    fragment_list = []  # save fragment names and fragment sequence
    custom_creation = []         # save custom keys with len = len(fragments)
    custom_ohs = []
    names = []
    for n_l, line in enumerate(csv_reader):
        # print(line)
        # separate at "=", remove newline from last element
        # if the len of a split line is 3
        if n_l == 0:
            meta_data = {meta.split("=")[0]: meta.split("=")[1] for meta in line}
        # all add-ons to plasmid sequences will be completely ignored
        elif "plasmid" in line[1]:
            if len(line) == 4:
                fragment_list.append(tuple(line[1:4]))
                custom_ohs.append("")
                custom_creation.append("")
            else:
                fragment_list.append(tuple(line[1:3]))
                custom_ohs.append("")
                custom_creation.append("")
        else:
            for idx, brick in enumerate(line):
                if idx == 0:
                    names.append(brick)
                elif idx == 1:
                    if any([frag in brick for frag in allowed_frag]):
                        pass
                    else:
                        raise ValueError("Unaccepted fragment {} in line {}".format(brick, n_l))
                elif idx == 2:
                    if all([base in allowed_base_seq for base in list(set(brick.upper()))]):
                        if brick.count("[") <= 1:
                            pass
                        else:
                            raise ValueError("Too many [ inside sequence in line {}".format(list(set(brick)), n_l))
                        if brick.count("]") <= 1:
                            pass
                        else:
                            raise ValueError("Too many ] inside sequence in line {}".format(list(set(brick)), n_l))
                    else:
                        raise ValueError("There are unaccepted bases below {} in line {}".format(list(set(brick)), n_l))
                elif idx == 3:
                    custom_creation.append("")
                    if len(brick) == 0:
                        pass
                    elif any([creation in brick for creation in allowed_creation]):
                        custom_creation[n_l-1] = brick
                    else:
                        raise ValueError("There are unaccepted elements in {} in line {}".format(brick, n_l))
                elif idx == 4:
                    custom_ohs.append("")
                    if len(brick) == 0:
                        pass
                    elif all([base in allowed_base_oh for base in list(set(brick))]):
                        custom_ohs[n_l-1] = brick
                    else:
                        raise ValueError("There are unaccepted bases inside the overhang {} in line {}".format(brick, n_l))
                else:
                    raise ValueError("Too many elements in line {}".format(n_l))
            fragment_list.append(tuple(line[1:3]))
    f.close()
    # print(meta_data)
    # output contains custom list and fragment tuple
    # print(tuple(fragment_list), custom_creation, custom_ohs)
    return tuple(fragment_list), custom_creation, custom_ohs, meta_data, names


# DEPRECATED
# added custom design functionality and overhangs as separate output lists
def read_in_file_custom(path):
    f = open(path, mode="r")
    allowed_base_oh = ["A", "C", "G", "T"]
    allowed_base_seq = ["A", "C", "G", "T", "[", "]"]
    allowed_frag = ["promoter", "scaffold", "seed"]
    allowed_creation = ["PCR", "Anneal"]
    fragment_list = []  # save fragment names and fragment sequence
    custom_creation = []         # save custom keys with len = len(fragments)
    custom_ohs = []
    file_content = f.readlines()
    meta_data = {param.split("=")[0].replace("\n", ""): param.split("=")[1].replace("\n", "") for param in file_content[0].split(sep=",")}
    for n_l, line in enumerate(file_content[1:]):
        # separate at "=", remove newline from last element
        # if the len of a split line is 3
        custom_ohs.append("")
        custom_creation.append("")
        polished_line = line.split(sep="=")
        polished_line[-1] = polished_line[-1].replace("\n", "")
        # all add-ons to plasmid sequences will be completely ignored
        if "plasmid" in polished_line[0]:
            fragment_list.append(tuple(polished_line[:2]))
        else:
            for idx, brick in enumerate(polished_line):
                if idx == 0:
                    if any([frag in brick for frag in allowed_frag]):
                        pass
                    else:
                        raise ValueError("Unaccepted fragment {} in line {}".format(brick, n_l))
                elif idx == 1:
                    if all([base in allowed_base_seq for base in list(set(brick))]):
                        if brick.count("[") <= 1:
                            pass
                        else:
                            raise ValueError("Too many [ inside sequence in line {}".format(list(set(brick)), n_l))
                        if brick.count("]") <= 1:
                            pass
                        else:
                            raise ValueError("Too many ] inside sequence in line {}".format(list(set(brick)), n_l))
                    else:
                        raise ValueError("There are unaccepted bases below {} in line {}".format(list(set(brick)), n_l))
                elif idx == 2:
                    if any([creation in brick for creation in allowed_creation]):
                        custom_creation[n_l] = brick
                    elif all([base in allowed_base_oh for base in list(set(brick))]):
                        custom_ohs[n_l] = brick
                    else:
                        raise ValueError("There are unaccepted elements in {} in line {}".format(brick, n_l))
                elif idx == 3:
                    if all([base in allowed_base_oh for base in list(set(brick))]):
                        custom_ohs[n_l] = brick
                    else:
                        raise ValueError("There are unaccepted bases inside the overhang {} in line {}".format(brick, n_l))
                else:
                    raise ValueError("Too many elements in line {}".format(n_l))
            fragment_list.append(tuple(polished_line))
    f.close()
    # output contains custom list and fragment tuple
    # print(tuple(fragment_list), custom_creation, custom_ohs)
    return tuple(fragment_list), custom_creation, custom_ohs, meta_data


def read_in_file_fasta(path):
    lines = []
    for record in SeqIO.parse(path, "fasta"):
        print(record)
        line = [record.description.split(",")[0], record.description.split(",")[1], str(record.seq).upper()]
        lines.append(line)
        print(lines)
    return lines


def check_meta_inputs(id, value):
    if id == "Melting Temperature":
        try:
            int(value)
        except Exception as exc:
            raise ValueError("Melting temperature has to be an integer")
    elif id == "Maximum Annealing Length":
        try:
            int(value)
        except Exception as exc:
            raise ValueError("Max Anneal Primer Length has to be an integer")
    elif id == "PCR Primer Length":
        try:
            int(value)
        except Exception as exc:
            raise ValueError("PCR Primer Length has to be an integer")
    elif id == "Overhang Length":
        try:
            int(value)
        except Exception as exc:
            raise ValueError("Overhang Length has to be an integer")
    elif id == "Bases for Gap":
        if not all([base in ["A", "C", "G", "T"] for base in list(set(value.upper()))]):
            raise ValueError("Fillup Bases contains unallowed bases, this tool only supports A, C, G and T")
    elif id == "Project Name":
        try:
            str(value)
        except Exception as exc:
            raise FileNotFoundError("Seems like your specified output directory could not be reached.\nPlease make sure that it exists.")
    elif id == "Restriction Enzyme":
        if value not in res_2s.keys():
            raise ValueError(f"Restriction Enzyme {value} is not included, please select one of:\n{','.join(res_2s.keys())}")
    elif id == "Plasmid gb file":
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Plasmid gb file {value} could not be found, please check the path you provided")
    # elif id == "Output No. Start":
    #     try:
    #         int(value)
    #     except Exception as exc:
    #         raise ValueError("Product No. Start has to be an integer")
    else:
        raise ValueError(f"Found unaccepted input id {id}")


def main():
    print(read_seq_from_csv("../Seed_Region_Finder/sRNATool_main/test/recA_mode3_output.csv"))
    print(read_in_file_fasta("./data/fasta_input.fasta"))


if __name__ == "__main__":
    main()
