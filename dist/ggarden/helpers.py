import csv
import os
import random
import numpy as np
from Bio import GenBank
from datasets import res_2s, primer3_pos
import primer3
from Bio.GenBank.Record import Record, Feature, Qualifier

'''
contains all kinds of helper functions
right now the following functions are available:
build_primer : builds forward and reverse primers for single strand or double strand
primer_insert : creates compatible primers < scaffold for insert cration
primer_pcr : creates primers for pcr of insert
translate : translates input seq in reversed and translated seq
res_site_finder : searches for restriction sites [motif] in [querry]
create_overhangs : should get number of building bricks as input and will (hopefully) build overhangs which do not 
                   build new res sites in combination with fragment seqs
update_overhangs : can be used to update overhang list and check, if overhangs are already present (fwd and rev)
                   may prove useful in create_overhangs
'''


# use_primer3 should be evaluated by searching for "[" and "]" in the input sequences
def build_primer_new(sequence, seq_type, overhangs, res_enz, primer_len, fillup_bases, tm=False, custom=False, use_primer_3=False, max_len=60):
    if seq_type[0] == "d":
        return excavate_donor(donor_seq=sequence, res_enz=res_enz)
    elif not custom:
        if len(sequence) > max_len:
            if use_primer_3 and tm:
                return primer_primer3(sequence=sequence, seq_type=seq_type, res_enz=res_enz, fillup_bases=fillup_bases, overhangs=overhangs, tm=tm, anneal_len=primer_len)
            elif not use_primer_3:
                return primer_pcr_new(sequence=sequence, seq_type=seq_type, res_enz=res_enz, fillup_bases=fillup_bases, overhangs=overhangs, anneal_len=primer_len)
            else:
                raise ValueError("unsupported combination of input parameters for primer design")
        else:
            return primer_insert_new(sequence=sequence, seq_type=seq_type, overhangs=overhangs)
    elif custom == "PCR":
        if use_primer_3 and tm:
            return primer_primer3(sequence=sequence, seq_type=seq_type, res_enz=res_enz, fillup_bases=fillup_bases, overhangs=overhangs, tm=tm, anneal_len=primer_len)
        elif not use_primer_3:
            return primer_pcr_new(sequence=sequence, seq_type=seq_type, res_enz=res_enz, fillup_bases=fillup_bases, overhangs=overhangs, anneal_len=primer_len)
        else:
            raise ValueError("unsupported combination of input parameters for primer design")
    elif custom == "Anneal":
        return primer_insert_new(sequence=sequence, seq_type=seq_type, overhangs=overhangs)
    else:
        raise ValueError("no valid keyword used for custom, available options are 'PCR' and 'Anneal'")


# sorts overhangs so all nested lists contain both overhangs of a fragment instead of both overhangs of a connection
def resort_overhangs(overhangs, conn_map):
    # returns overhangs sorted in a way that every nested list contains both overhangs of one fragement as [fwd, rev]
    cleaned_oh = [oh for ohp in overhangs for oh in ohp]
    cleaned_oh.insert(0, cleaned_oh.pop())
    cleaned_oh_fin = [[cleaned_oh[idx], cleaned_oh[idx+1]] for idx in range(0, len(cleaned_oh), 2)]
    return cleaned_oh_fin


def primer_insert_new(sequence, seq_type, overhangs):

    if "[" in  sequence and "]" in sequence:
        sequence = sequence.replace("[",  "")
        sequence = sequence.replace("]",  "")

    if "plasmid" not in seq_type and "scaffold" not in seq_type and "promoter" not in seq_type:
        # build_scaffold(cont_type=seq_type, overhang_len=len(overhangs[0]), align_len=0, res_ezn=False, rand_b=0, fill_b=0, application="insert")
        primer_fwd = overhangs[0].upper() + sequence.lower()
        primer_rev = overhangs[1].upper() + translate(sequence).lower()
        print("annealing primers")
        primer3_pos.append(-1)
        return [primer_fwd, primer_rev]
    elif "plasmid" in seq_type:
        # this kind of plasmid should not contain res_sites or similar
        # primer_fwd = overhangs[0] + sequence
        # primer_rev = overhangs[1] + translate(sequence)
        # return [primer_fwd, primer_rev]
        pass
    elif "scaffold" in seq_type:
        # added new options
        # both overhangs are scaffold parts
        # one of both overhangs are scaffold parts
        if overhangs[0] == sequence[:4] and overhangs[1] == translate(sequence[-4:]):
            primer_fwd = sequence.lower()
            primer_rev = translate(sequence).lower()
            print("annealing primers")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        elif overhangs[0] == sequence[:4]:
            primer_fwd = sequence.lower()
            primer_rev = overhangs[1].upper() + translate(sequence)[:-len(overhangs[0])].lower()
            print("annealing primers")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        elif overhangs[1] == translate(sequence[-4:]):
            primer_rev = translate(sequence).lower()
            primer_fwd = overhangs[0].upper() + sequence[:-len(overhangs[0])].lower()
            print("annealing primers")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        else:
            primer_fwd = overhangs[0].upper() + sequence.lower()
            primer_rev = overhangs[1].upper() + translate(sequence).lower()
            print("annealing primers")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
    elif "promoter" in seq_type:
        if overhangs[1] == translate(sequence[-4:]):
            primer_rev = translate(sequence).lower()
            primer_fwd = overhangs[0].upper() + sequence[:-len(overhangs[0])].lower()
            print("annealing primers")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        else:
            primer_fwd = overhangs[0].upper() + sequence.lower()
            primer_rev = overhangs[1].upper() + translate(sequence).lower()
            print("annealing primers")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]


def primer_pcr_new(sequence, seq_type, res_enz, fillup_bases, overhangs, anneal_len=20):
    print(seq_type)
    print(overhangs)
    res_info = res_2s[res_enz]
    filler_bases = res_info[1] * fillup_bases
    filler_bases = filler_bases[:res_info[1]]
    rand_bases = fillup_bases
    # filler_bases = "".join(filler_bases)
    # rand_bases = "".join(rand_bases)
    if "plasmid" not in seq_type and "scaffold" not in seq_type and "promoter" not in seq_type:
        primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[0].upper() + sequence[:anneal_len].lower()
        primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[1].upper() + translate(sequence)[:anneal_len].lower()
        print("non-primer3 pcr primers")
        primer3_pos.append(-1)
        return [primer_fwd, primer_rev]
    elif "plasmid" in seq_type:
        # print("primers for plasmid amplification should not be needed")
        # return [primer_fwd, primer_rev]
        pass
    elif "scaffold" in seq_type:
        if overhangs[0] == sequence[:4] and overhangs[1] == translate(sequence[-4:]):
            primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + sequence[:anneal_len].lower()
            primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + translate(sequence)[:anneal_len].lower()
            print("non-primer3 pcr primers, no scar")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        elif overhangs[0] == sequence[:4]:
            primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + sequence[:anneal_len].lower()
            primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[1].upper() + translate(sequence)[:anneal_len].lower()
            print("non-primer3 pcr primers, no scar")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        elif overhangs[1] == translate(sequence[-4:]):
            primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[0].upper() + sequence[:anneal_len].lower()
            primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + translate(sequence)[:anneal_len].lower()
            print("non-primer3 pcr primers with scar")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        else:
            primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[0].upper() + sequence[:anneal_len].lower()
            primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[1].upper() + translate(sequence)[:anneal_len].lower()
            print("non-primer3 pcr primers with scar")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
    elif "promoter" in seq_type:
        if overhangs[1] == translate(sequence[-4:]):
            primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[0].upper() + sequence[:anneal_len].lower()
            primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + translate(sequence)[:anneal_len].lower()
            print("non-primer3 pcr primers, no scar")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]
        else:
            primer_fwd = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[0].upper() + sequence[:anneal_len].lower()
            primer_rev = rand_bases.lower() + res_info[0] + filler_bases.lower() + overhangs[1].upper() + translate(sequence)[:anneal_len].lower()
            print("non-primer3 pcr primers with scar")
            primer3_pos.append(-1)
            return [primer_fwd, primer_rev]


# the sequence for this function has to contain exactly one [ and one ] inside the seq to mark the start or end of the
# functional unit
def primer_primer3(sequence, seq_type, res_enz, fillup_bases, overhangs, tm=60, anneal_len=20):
    primer_seq = preprocess_primer3(sequence=sequence, seq_type=seq_type, res_enz=res_enz)
    res_info = res_2s[res_enz]
    filler_bases = res_info[1] * fillup_bases
    filler_bases = filler_bases[:res_info[1]]
    rand_bases = fillup_bases
    # filler_bases = "".join(filler_bases)
    # rand_bases = "".join(rand_bases)
    if "scaffold" in seq_type:
        if overhangs[0] == sequence[:4] and overhangs[1] == translate(sequence[-4:]):
            seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower()
            seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower()
        elif overhangs[0] == sequence[:4]:
            seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower()
            seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[1].upper()
        elif overhangs[1] == translate(sequence[-4:]):
            seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[0].upper()
            seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower()
        else:
            seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[0].upper()
            seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[1].upper()
        unwanted = "[]"
        start = primer_seq.index("[")
        end = primer_seq.index("]")
        # print(end-start)
        for bracket in unwanted:
            primer_seq = primer_seq.replace(bracket, "")
        # print(primer_seq[start:start + end -start - 1])
        prims = primer3.designPrimers({"SEQUENCE_TEMPLATE": primer_seq,
                                       "SEQUENCE_OVERHANG_LEFT": seq_oh_fwd,
                                       "SEQUENCE_OVERHANG_RIGHT": seq_oh_rev,
                                       "SEQUENCE_TARGET": [anneal_len+5, end-(anneal_len+5)],   # not absolutely sure if spot is right
                                       "SEQUENCE_FORCE_LEFT_START": 0                           # not absolutely sure if spot is right
                                       },
                                      {"PRIMER_TASK": "pick_pcr_primers",
                                       "PRIMER_OPT_SIZE": anneal_len,
                                       "PRIMER_MIN_SIZE": anneal_len-5,
                                       "PRIMER_MAX_SIZE": anneal_len+5,
                                       "PRIMER_OPT_TM": tm,
                                       "PRIMER_MIN_TM": tm-5,
                                       "PRIMER_MAX_TM": tm+5,
                                       "PRIMER_MIN_GC": 20.0,   # should be adapted
                                       "PRIMER_MAX_GC": 80.0,   # should be adapted
                                       "PRIMER_NUM_RETURN": 1,
                                       "PRIMER_PRODUCT_SIZE_RANGE": [anneal_len+5, 10000]})
        # print(prims)
        if prims["PRIMER_PAIR_NUM_RETURNED"] > 0:
            primer_fwd = seq_oh_fwd + prims["PRIMER_LEFT_0_SEQUENCE"].lower()
            primer_rev = seq_oh_rev + prims["PRIMER_RIGHT_0_SEQUENCE"].lower()
            print("primer3 pcr primers")
            primer3_pos.append((prims["PRIMER_LEFT_0"], prims["PRIMER_RIGHT_0"]))
            return[primer_fwd, primer_rev]
        # else statement could be replaced by error or transitioned in non primer3 primer
        else:
            print(prims)
            raise ValueError("Primer3 did not find matching primers for the preprocessed sequence {}".format(primer_seq))
    elif "promoter" in seq_type:
        if overhangs[1] == sequence[-4:]:
            seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[0].upper()
            seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower()
        else:
            seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[0].upper()
            seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[1].upper()
        unwanted = "[]"
        start = primer_seq.index("[")
        end = primer_seq.index("]")
        # print(end-start)
        for bracket in unwanted:
            primer_seq = primer_seq.replace(bracket, "")
        # print(primer_seq[start:start + end -start - 1])
        prims = primer3.designPrimers({"SEQUENCE_TEMPLATE": primer_seq,
                                       "SEQUENCE_OVERHANG_LEFT": seq_oh_fwd,
                                       "SEQUENCE_OVERHANG_RIGHT": seq_oh_rev,
                                       "SEQUENCE_TARGET": [start, len(primer_seq)-start-(anneal_len+5)],    # not absolutely sure if spot is right
                                       "SEQUENCE_FORCE_RIGHT_START": len(primer_seq)-1                      # not absolutely sure if spot is right
                                       },
                                      {"PRIMER_TASK": "pick_pcr_primers",
                                       "PRIMER_OPT_SIZE": anneal_len,
                                       "PRIMER_MIN_SIZE": anneal_len-5,
                                       "PRIMER_MAX_SIZE": anneal_len+5,
                                       "PRIMER_OPT_TM": tm,
                                       "PRIMER_MIN_TM": tm-5,
                                       "PRIMER_MAX_TM": tm+5,
                                       "PRIMER_MIN_GC": 20.0,   # should be adapted
                                       "PRIMER_MAX_GC": 80.0,   # should be adapted
                                       "PRIMER_NUM_RETURN": 1,
                                       "PRIMER_PRODUCT_SIZE_RANGE": [anneal_len+5, 10000]})
        # print(prims)
        if prims["PRIMER_PAIR_NUM_RETURNED"] > 0:
            primer_fwd = seq_oh_fwd + prims["PRIMER_LEFT_0_SEQUENCE"].lower()
            primer_rev = seq_oh_rev + prims["PRIMER_RIGHT_0_SEQUENCE"].lower()
            print("primer3 pcr primers")
            primer3_pos.append((prims["PRIMER_LEFT_0"], prims["PRIMER_RIGHT_0"]))
            return[primer_fwd, primer_rev]
        # else statement could be replaced by error or transitioned in non primer3 primer
        else:
            print(prims)
            raise ValueError("Primer3 did not find matching primers for the preprocessed sequence {}".format(primer_seq))
    else:
        seq_oh_fwd = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[0].upper()
        seq_oh_rev = rand_bases.lower() + res_info[0].upper() + filler_bases.lower() + overhangs[1].upper()
        unwanted = "[]"
        start = primer_seq.index("[")
        end = primer_seq.index("]")
        # print(end-start)
        for bracket in unwanted:
            primer_seq = primer_seq.replace(bracket, "")
        # print(primer_seq[start:start + end -start - 1])
        prims = primer3.designPrimers({"SEQUENCE_TEMPLATE": primer_seq,
                                       "SEQUENCE_OVERHANG_LEFT": seq_oh_fwd,
                                       "SEQUENCE_OVERHANG_RIGHT": seq_oh_rev,
                                       "SEQUENCE_TARGET": [start, end-start]},
                                      {"PRIMER_TASK": "pick_pcr_primers",
                                       "PRIMER_OPT_SIZE": anneal_len,
                                       "PRIMER_MIN_SIZE": anneal_len-5,
                                       "PRIMER_MAX_SIZE": anneal_len+5,
                                       "PRIMER_OPT_TM": tm,
                                       "PRIMER_MIN_TM": tm-5,
                                       "PRIMER_MAX_TM": tm+5,
                                       "PRIMER_MIN_GC": 20.0,   # should be adapted
                                       "PRIMER_MAX_GC": 80.0,   # should be adapted
                                       "PRIMER_NUM_RETURN": 1,
                                       "PRIMER_PRODUCT_SIZE_RANGE": [anneal_len+5, 10000]})
        # print(prims)
        if prims["PRIMER_PAIR_NUM_RETURNED"] > 0:
            primer_fwd = seq_oh_fwd + prims["PRIMER_LEFT_0_SEQUENCE"].lower()
            primer_rev = seq_oh_rev + prims["PRIMER_RIGHT_0_SEQUENCE"].lower()
            print("primer3 pcr primers")
            primer3_pos.append((prims["PRIMER_LEFT_0"], prims["PRIMER_RIGHT_0"]))
            return[primer_fwd, primer_rev]
        # else statement could be replaced by error or transitioned in non primer3 primer
        else:
            print(prims)
            raise ValueError("Primer3 did not find matching primers for the preprocessed sequence {}".format(primer_seq))


# simple function to produce second strand based on first strand
# returns cases as present in querry
def translate(querry):
    base_dict = {
        "a": "t",
        "t": "a",
        "c": "g",
        "g": "c",
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }
    strand_out = ""
    for base in querry:
        strand_out += base_dict[base]
    return strand_out[::-1]


def ligation_prob(oh, ohs, mat):
    rev_oh = translate(oh)
    correct = mat.loc[oh, rev_oh]
    total = mat.loc[oh, ohs].sum()
    print(mat.loc[ohs, ohs])
    return correct/total


def ligation_prob_single(oh, ohs, mat):
    return 0


def ligation_fidelity(ohs, mat):
    set_fidelity = 1
    fidelities = {}
    unique_ohs = list(dict.fromkeys(ohs.copy()))
    print(ohs)
    print(unique_ohs)
    reduced_oh_mat = mat.loc[unique_ohs, unique_ohs]
    for idx, oh in enumerate(ohs):
        rev_oh = translate(oh)
        if rev_oh != oh and ohs.count(oh) == 1:
            oh_sum = reduced_oh_mat.loc[oh, rev_oh:].sum()
            fidelity = reduced_oh_mat.loc[oh, rev_oh] / oh_sum
            fidelity = np.round(fidelity, 2)
            set_fidelity *= fidelity
            fidelities[oh] = fidelity
        elif rev_oh == oh:
            set_fidelity *= 0.5
            fidelities[oh] = 0.5
        else:
            set_fidelity *= 1/ohs.count(oh)
            fidelities[oh] = 1/ohs.count(oh)
    set_fidelity = np.round(set_fidelity, 2)

    return fidelities, set_fidelity


# searches for resIIS sites in querry
# motif and querry pieces are temporarily moved to upper for comparison
def res_site_finder(motif, querry):
    site_idx = []
    for idx in range(len(querry) - len(motif) + 1):
        if querry[idx:idx+len(motif)].upper() == motif.upper():
            site_idx.append(idx)
    return site_idx


# builds random overhangs for connections that don't have already set overhang sequence
def random_overhangs(length):
    try:
        round(length)
    except TypeError:
        raise AssertionError("Input variable is not a number")
    random.seed(a=2233213412)
    corrected_length = round(length)
    base_list = ["A", "C", "G", "T"]
    overhang = random.choices(population=base_list, k=corrected_length)
    overhang = "".join(overhang)
    return overhang


# produces one pair of overhangs based on the connection type
# returned overhangs are created as [rev, fwd]
# !!!why does plas-sth return overhangs in the wrong order?!!!
# answer: wrong variable declaration
# oh rev designed based on the rev_seq (previously oh fwd)
# oh detection for plasmids and dondors based on exccavate_donor_ohs
def find_overhang(seq, conn_type, length=4, res_enz=False):
    # problematic_conns are plas-d*any, d*any-plas, dseed-scaf
    # conn_type is written as "rev-fwd" so the output of this function equals [rev_oh, fwd_oh]
    if conn_type.split("-")[0][0] == "d" or conn_type.split("-")[1][0] == "d":
        if res_enz:
            if len(seq) != 2:
                ohs = excavate_donor_ohs(donor_seq=seq, res_enz=res_enz)
                if conn_type.split("-")[0][0] == "d":
                    return [ohs[0], translate(ohs[0])]
                elif conn_type.split("-")[1][0] == "d":
                    return [translate(ohs[1]), ohs[1]]
                else:
                    pass
            else:
                ohs = excavate_problematic_ohs(conn_seqs=seq, conn_type=conn_type, res_enz=res_enz)
                if conn_type.split("-")[0][0] == "d":
                    return [ohs[0], translate(ohs[0])]
                elif conn_type.split("-")[1][0] == "d":
                    return [translate(ohs[1]), ohs[1]]
                else:
                    pass
        else:
            raise ValueError("No res IIS enzyme specified for donor plasmid")
    elif conn_type == "seed-scaf" or conn_type == "prom-scaf":
        overhang_fwd = seq[:length]
        return [translate(overhang_fwd), overhang_fwd]
    # scaf-seed new >scaf<-seed-scaff type
    # prom-seed should be scarless by default
    elif conn_type == "scaf-seed" or conn_type == "prom-seed":
        overhang_fwd = seq[-length:]
        return [translate(overhang_fwd), overhang_fwd]
    elif conn_type == "scaf-plas" or conn_type == "seed-plas":  # changed plas-prom, plas-seed to scaf-plas, seed-plas
        if conn_type == "seed-plas":
            print("Information: found medium or low complexity connection\nIf this was intended, everything is working fine")
        if res_enz:
            overhang_fwd = excavate_donor_ohs(donor_seq=seq, res_enz=res_enz)[1]
            return [translate(overhang_fwd), overhang_fwd]
        else:
            raise ValueError("No res IIS enzyme specified")
    elif conn_type == "plas-prom" or conn_type == "plas-seed":  # changed scaf-plas, seed-plas to plas-prom, plas-seed
        if conn_type == "plas-seed":
            print("Information: found medium or low complexity connection\nIf this was intended, everything is working fine")
        if res_enz:
            overhang_rev = excavate_donor_ohs(donor_seq=seq, res_enz=res_enz)[0]
            return [overhang_rev, translate(overhang_rev)]
        else:
            raise ValueError("No res IIS enzyme specified")
    elif conn_type == "scaf-prom":
        overhang_fwd = random_overhangs(length=length)
        return [translate(overhang_fwd), overhang_fwd]
    else:
        raise ValueError("Found unexpected connection {}".format(conn_type))


# iterates over file content and produces connections in form of frag-frag for all connections
# these are needed to build the correct overhangs for each connection
def build_conn_map(file_cont):
    conn_map = []
    for idx in range(len(file_cont)-1):
        conn_map.append("{}-{}".format(file_cont[idx][0][:4], file_cont[idx+1][0][:4]))
    conn_map.append("{}-{}".format(file_cont[-1][0][:4], file_cont[0][0][:4]))
    print(conn_map)
    return conn_map


# returns the sequence needed for overhang construction based on the connection type
# especially useful for sth-scaffold, plasmid-sth, sth-plasmid
def input_handler(conn, idx, content):
    if idx != len(content)-1:
        if (conn.split("-")[0][0] == "d" and conn.split("-")[1] == "plas") \
                or (conn.split("-")[1][0] == "d" and conn.split("-")[0] == "plas") \
                or (conn.split("-")[0][0] == "d" and conn.split("-")[1] == "scaf") \
                or (conn.split("-")[0][0] == "d" and conn.split("-")[1][0] == "d") \
                or (conn.split("-")[1] == "dsee" and conn.split("-")[0] == "scaf"):
            # raise ValueError("Problematic overhang correction is not supported at this point")
            return [content[idx][1], content[idx+1][1]]
        elif conn.split("-")[0][0] == "d":
            return content[idx][1]
        elif conn.split("-")[1][0] == "d":
            return content[idx+1][1]
        elif conn.split("-")[1] == "scaf":
            return content[idx+1][1]
        elif conn.split("-")[0] == "prom":
            return content[idx][1]
        elif conn.split("-")[0] == "plas":
            return content[idx][1]
        elif conn.split("-")[1] == "plas":
            return content[idx+1][1]
        elif conn == "plas-scaf":
            raise ValueError("This connection type is not supported")
        else:
            return content[idx][1]
    else:
        if (conn.split("-")[0][0] == "d" and conn.split("-")[1] == "plas") \
                or (conn.split("-")[1][0] == "d" and conn.split("-")[0] == "plas") \
                or (conn.split("-")[0][0] == "d" and conn.split("-")[1] == "scaf") \
                or (conn.split("-")[0][0] == "d" and conn.split("-")[1][0] == "d") \
                or (conn.split("-")[1] == "dseed" and conn.split("-")[0] == "scaf"):
            # raise ValueError("Problematic overhang correction is not supported at this point")
            return [content[idx][1], content[0][1]]
        elif conn.split("-")[0][0] == "d":
            return content[idx][1]
        elif conn.split("-")[1][0] == "d":
            return content[0][1]
        elif conn.split("-")[1] == "scaf":
            return content[0][1]
        elif conn.split("-")[0] == "prom":
            return content[idx][1]
        elif conn.split("-")[0] == "plas":
            return content[idx][1]
        elif conn.split("-")[1] == "plas":
            return content[0][1]
        elif conn == "plas-scaf":
            raise ValueError("This connection type is not supported")
        else:
            return content[idx][1]


# backtracking missing
# not used right now
def needleman_wunsch(seq1, seq2):
    def w(b1, b2):
        if b1 == b2:
            return -1
        else:
            return 1
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    n_w_matrix = np.zeros(shape=(len(seq1)+1, len(seq2)+1))
    move_mat = np.zeros(shape=(len(seq1)+1, len(seq2)+1))
    for x in range(len(seq1)+1):
        for y in range(len(seq2)+1):
            if x == 0:
                n_w_matrix[x, y] = -y
            elif y == 0:
                n_w_matrix[x, y] = -x
            else:
                possibles = [n_w_matrix[x-1, y-1] - w(seq1[x-1], seq2[y-1]),
                             n_w_matrix[x-1, y] - 2,
                             n_w_matrix[x, y-1] - 2]
                move = np.argmax(possibles)
                print(move)
                n_w_matrix[x, y] = possibles[move]
                move_mat[x-1, y-1] = move
    print(n_w_matrix)
    print(move_mat)


# writes all primers as [fwd, rev] and fragment names in one output file
def build_prim_output(primer_list, content, names, filename, constr_nr):
    prim_seq_f_name = "{}_seq.txt".format(filename)
    coll_name = "{}coll.csv".format(filename[:-len(filename.split("_")[-1])])
    # if not os.path.isfile(prim_seq_f_name):
    f = open(prim_seq_f_name, mode="w")
    coll = open(coll_name, mode="a", newline="")
    content_adder = csv.writer(coll, delimiter=";")
    for idx, prim_pair in enumerate(primer_list):
        if "donor" not in content[idx+1][0]:
            f.write("{} {} (fwd,rev)\n{}\n{}\n\n".format(names[idx], prim_pair[0], prim_pair[1][0], prim_pair[1][1]))
            content_adder.writerow([constr_nr, f"{names[idx]}-fwd", prim_pair[1][0]])
            content_adder.writerow([constr_nr, f"{names[idx]}-rev", prim_pair[1][1]])
    coll.close()
    f.close()
    # else:
    #     print("file already exists")


# visualization function for plasmid and primer connections
def build_prim_vis(primer_list, content, res_enz, filename,fillup_bases, oh_len=4, fragm_len=20):
    print(primer_list)
    prim_vis_f_name = "{}_vis.txt".format(filename)
    # if not os.path.isfile(prim_vis_f_name):
    f = open(prim_vis_f_name, mode="w")
    # else:
    #     print("file already exists")
    fwd_str = ""
    rev_str = ""
    start_pos = len(res_2s[res_enz][0]) + res_2s[res_enz][1]
    plas_info = content[0]
    plas_starts = [res_site_finder(motif=res_2s[res_enz][0], querry=plas_info[1]), res_site_finder(motif=res_2s[res_enz][0], querry=translate(plas_info[1]))]
    if (plas_starts[1][0] + start_pos) >= len(plas_info[1]):
        start_bases = plas_starts[1][0] + start_pos - len(plas_info[1])
        fragment = "{}{}".format(translate(plas_info[1])[len(plas_info[1])-start_bases:], translate(plas_info[1])[:fragm_len-start_bases])
        rev_str += fragment[::-1]
        fwd_str += translate(fragment)[:-oh_len]
    else:
        rev_str += translate(plas_info[1])[plas_starts[1][0]+start_pos:plas_starts[1][0]+start_pos+oh_len + fragm_len][::-1]
        fwd_str += translate(rev_str[:-oh_len])[::-1]
    for prim_pair in primer_list:
        if prim_pair[1][0][0] == prim_pair[1][0][0].lower() and prim_pair[1][0][len(fillup_bases)] == prim_pair[1][0][len(fillup_bases)].upper():
            start_pos = len(res_2s[res_enz][0]) + res_2s[res_enz][1] + len(fillup_bases)
            prim_fwd_ds = [prim_pair[1][0][start_pos:], translate(prim_pair[1][0][start_pos+oh_len:])[::-1]]
            prim_rev_ds = [prim_pair[1][1][start_pos:][::-1], translate(prim_pair[1][1][start_pos+oh_len:])]
            fwd_str += oh_len * " " + prim_fwd_ds[0] + "..." + prim_rev_ds[1]
            rev_str += oh_len * " " + prim_fwd_ds[1] + "..." + prim_rev_ds[0]
        else:
            fwd_str += oh_len * " " + prim_pair[1][0]
            rev_str += oh_len * " " + prim_pair[1][1][::-1]
    start_pos = len(res_2s[res_enz][0]) + res_2s[res_enz][1]
    # fwd_str += oh_len * " " + plas_info[1][start_pos:fragm_len+start_pos]
    # rev_str += oh_len * " " + translate(plas_info[1])[::-1][start_pos+oh_len:fragm_len+start_pos]
    if (plas_starts[0][0] + start_pos + oh_len) > len(plas_info[1]):
        start_bases = plas_starts[1][0] + start_pos - len(plas_info[1])
        fragment = "{}{}".format(plas_info[1][len(plas_info[1])-start_bases:], plas_info[1][:fragm_len-start_bases])
        fwd_str += oh_len * " " + fragment[::-1]
        rev_str += oh_len * " " + translate(fragment)[:-4]
    else:
        fwd_str += oh_len * " " + plas_info[1][plas_starts[0][0]+start_pos:plas_starts[0][0]+start_pos+fragm_len]
        rev_str += oh_len * " " + translate(plas_info[1])[::-1][plas_starts[0][0]+start_pos+oh_len:plas_starts[0][0]+start_pos+fragm_len]
    f.write("{}\n{}".format(fwd_str, rev_str))
    f.close()


def build_prim_coll(filename):
    col_name = "{}_coll.csv".format(filename)
    col_name_fin = "{}_fin_coll.csv".format(filename)
    if not os.path.isfile(col_name_fin):
        prim_names = []
        prim_seqs = []
        prim_nrs = []
        col_names = ["Primer"]
        col_nrs = ["Product No."]
        coll = open(col_name, mode="r", newline="")
        content_reader = csv.reader(coll, delimiter=";")
        for idx, line in enumerate(content_reader):
            if idx > 0:
                prim_names.append(line[1])
                prim_nrs.append(line[0])
                prim_seqs.append(line[2])
        coll.close()
        unique_prim_seqs = set(prim_seqs)
        for seq in unique_prim_seqs:
            names_of_seq = [prim_names[idx] for idx, prim_seq in enumerate(prim_seqs) if seq == prim_seq]
            nr_of_seq = ",".join([prim_nrs[idx] for idx, prim_seq in enumerate(prim_seqs) if seq == prim_seq])
            # shouldnt be possible because all overhangs are unique so no identical seqs can be created for different fragments
            if len(names_of_seq) > 1:
                if len(set(names_of_seq)) > 1:
                    # raise ValueError("Unexpected discovery of two identical sequences {} from different fragments. This should not be possible".format(",".join(names_of_seq)))
                    names_of_seq = ",".join(names_of_seq)
                    # names_of_seq = names_of_seq
                else:
                    names_of_seq = names_of_seq[0]
            else:
                names_of_seq = ",".join(names_of_seq)
            print(names_of_seq)
            col_names.append(names_of_seq)
            col_nrs.append(nr_of_seq)
        unique_prim_seqs = ["Sequence (5' to 3')"] + list(unique_prim_seqs)
        col_content = list(zip(col_nrs, col_names, unique_prim_seqs))
        col_fin = open(col_name_fin, mode="w", newline="")
        content_writer = csv.writer(col_fin, delimiter=";")
        content_writer.writerows(col_content)
        col_fin.close()
    else:
        raise FileExistsError("Final collection already exists")


# this little guy searches in primer3 input sequences upstream and downstream of the functional sequence for undetected
# cut sites and removes the cutsites including all nonfunctional sequences upstream or downstream
def preprocess_primer3(sequence, seq_type, res_enz):
    if "scaffold" in seq_type:
        idx_down = 1 + sequence.index("]")
        seq_down = sequence[idx_down:]
        sites_down = res_site_finder(motif=res_2s[res_enz][0], querry=translate(seq_down))
        if len(sites_down) > 0:
            return sequence[:-max(sites_down)-len(res_2s[res_enz][0])]
        else:
            return sequence
    elif "promoter" in seq_type:
        idx_up = sequence.index("[")
        seq_up = sequence[:idx_up]
        sites_up = res_site_finder(motif=res_2s[res_enz][0], querry=seq_up)
        if len(sites_up) > 0:
            return sequence[max(sites_up)+len(res_2s[res_enz][0]):]
        else:
            return sequence
    else:
        idx_up = sequence.index("[")
        seq_up = sequence[:idx_up]
        sites_up = res_site_finder(motif=res_2s[res_enz][0], querry=seq_up)
        idx_down = 1 + sequence.index("]")
        seq_down = sequence[idx_down:]
        sites_down = res_site_finder(motif=res_2s[res_enz][0], querry=translate(seq_down))
        if len(sites_up) > 0 and len(sites_down) > 0:
            return sequence[max(sites_up)+len(res_2s[res_enz][0]):-max(sites_down)-len(res_2s[res_enz][0])]
        elif len(sites_up) > 0 and len(sites_down) == 0:
            return sequence[max(sites_up)+len(res_2s[res_enz][0]):]
        elif len(sites_up) == 0 and len(sites_down) > 0:
            return sequence[:-max(sites_down)-len(res_2s[res_enz][0])]
        else:
            return sequence


def build_genbank_output(filename, content, res_enz, oh_len, ohs, custom_ohs, gb_plasmid_file=""):
    if gb_plasmid_file != "":
        plas_features = cut_gb(gb_file=gb_plasmid_file, res_enz=res_enz, oh_len=oh_len)
    # print(primer3_pos)
    f = open("{}_gb.gbk".format(filename), "w")
    seq = ""
    res_info = res_2s[res_enz]
    res_len = res_info[1] + len(res_info[0])
    seq_len = 0
    ends = [0]
    features = []
    if gb_plasmid_file != "":
        features.extend(plas_features)
    for idx, fragment in enumerate(content):
        fragment_seq = fragment[1]
        fragment_name = fragment[0]
        if "plasmid" in fragment_name:
            fwd_cutsite = res_site_finder(querry=fragment_seq, motif=res_info[0])[0]
            rev_cutsite = res_site_finder(querry=translate(fragment_seq), motif=res_info[0])
            rev_cutsite = len(fragment_seq) - rev_cutsite[0]
            # print(fwd_cutsite, rev_cutsite)
            if fwd_cutsite < rev_cutsite:
                fragment_seq = fragment_seq[fwd_cutsite+res_len:rev_cutsite-res_len-oh_len]
            else:
                if rev_cutsite < (res_len + oh_len):
                    fragment_seq = fragment_seq[fwd_cutsite+res_len:rev_cutsite - res_len - oh_len]
                else:
                    fragment_seq = fragment_seq[fwd_cutsite+res_len:] + fragment_seq[:rev_cutsite - res_len - oh_len]
            # print(len(fragment_seq))
            seq += fragment_seq
            seq_len += len(fragment_seq)
            ends.append(seq_len)
            features.append(Feature(key="source", location="{}..{}".format(ends[-2] + 5, ends[-1])))
            features.append(Qualifier(key="/gene=", value=fragment_name))
            features.append(Feature(key="misc_feature", location="{}..{}".format(ends[-2] + 1, ends[-2] + res_info[2])))
            features.append(Qualifier(key="/label=", value="overhang"))
        elif "promoter" in fragment_name:
            if primer3_pos[idx-1] != -1:
                # get primer3 primer positions or add fwd_overhang
                fragment_seq = fragment_seq.replace("[", "")
                fragment_seq = fragment_seq.replace("]", "")
                fragment_start = primer3_pos[idx-1][0][0]
                fragment_seq = ohs[idx-1][1] + fragment_seq[fragment_start:]
                if custom_ohs[idx - 1] == "":
                    fragment_seq = fragment_seq[:-res_info[2]]
            else:
                fragment_seq = ohs[idx-1][1] + fragment_seq
                if custom_ohs[idx - 1] == "":
                    fragment_seq = fragment_seq[:-res_info[2]]
            seq += fragment_seq
            seq_len += len(fragment_seq)
            ends.append(seq_len)
            features.append(Feature(key="regulatory", location="{}..{}".format(ends[-2] + res_info[2] + 1, ends[-1])))
            features.append(Qualifier(key="/standard_name=", value=fragment_name))
            features.append(Qualifier(key="/regulatory_class=", value="promoter"))
            features.append(Feature(key="misc_feature", location="{}..{}".format(ends[-2] + 1, ends[-2] + res_info[2])))
            features.append(Qualifier(key="/label=", value="overhang"))
        elif "scaffold" in fragment_name:
            if primer3_pos[idx-1] != -1:
                # get primer3 primer positions or add fwd_overhang
                fragment_seq = fragment_seq.replace("[", "")
                fragment_seq = fragment_seq.replace("]", "")
                fragment_end = primer3_pos[idx-1][1][0] + 1
                fragment_seq = ohs[idx-1][1] + fragment_seq[:fragment_end]
                if custom_ohs[idx - 1] == "":
                    fragment_seq = fragment_seq[res_info[2]:]
            seq += fragment_seq
            seq_len += len(fragment_seq)
            ends.append(seq_len)
            features.append(Feature(key="gene", location="{}..{}".format(ends[-3] + 1, ends[-1])))
            features.append(Qualifier(key="/gene=", value="sRNA"))
            features.append(Feature(key="CDS", location="{}..{}".format(ends[-2] + res_info[2] + 1, ends[-1])))
            features.append(Qualifier(key="/gene=", value=fragment_name))
            features.append(Feature(key="misc_feature", location="{}..{}".format(ends[-2] + 1, ends[-2] + res_info[2])))
            features.append(Qualifier(key="/label=", value="overhang"))
        elif "seed" in fragment_name:
            # add fwd_overhang
            fragment_seq = ohs[idx-1][1] + fragment_seq
            seq += fragment_seq
            seq_len += len(fragment_seq)
            ends.append(seq_len)
            features.append(Feature(key="CDS", location="{}..{}".format(ends[-2] + res_info[2] + 1, ends[-1])))
            features.append(Qualifier(key="/gene=", value=fragment_name))
            features.append(Feature(key="misc_feature", location="{}..{}".format(ends[-2] + 1, ends[-2] + res_info[2])))
            features.append(Qualifier(key="/label=", value="overhang"))
        print(fragment_seq)
    container = Record()
    container.locus = "GGA Plasmid"
    container.residue_type = "DNA circular"
    container.definition = "Resulting plasmid from golden gate assembly"
    container.source = "Golden Gate Assembly"
    container.organism = "Saccharomyces cerevisiae"
    container.size = len(seq)
    container.sequence = seq
    container.origin = seq
    container.features = features
    f.write(str(container))
    f.close()


def cut_gb(gb_file, res_enz, oh_len):
    file = open(gb_file)
    idxes_to_remove = []
    content = GenBank.parse(file)
    for entry in content:
        boundaries = []
        entry_seq = entry.sequence
        boundaries.append(res_site_finder(querry=entry_seq, motif=res_2s[res_enz][0])[0])
        boundaries.append(res_site_finder(querry=translate(entry_seq), motif=res_2s[res_enz][0])[0])
        boundaries[0] += len(res_2s[res_enz][0]) + res_2s[res_enz][1]
        boundaries[1] += len(res_2s[res_enz][0]) + res_2s[res_enz][1] + oh_len
        boundaries[1] = len(entry_seq) - boundaries[1]
        if boundaries[1] < 0:
            boundaries[1] = len(entry_seq) + boundaries[1]
        if boundaries[0] < boundaries[1]:
            # remove everything < bondaries[0] and > boundaries[1]
            for idx, feature in enumerate(entry.features):
                stripped_feature_loc = strip_location(feature.location)
                if ".." in stripped_feature_loc:
                    if int(stripped_feature_loc.split("..")[0]) < boundaries[0] or int(stripped_feature_loc.split("..")[-1]) > boundaries[1]:
                        idxes_to_remove.append(idx)
                    else:
                        feature.location = "{}..{}".format(int(stripped_feature_loc.split("..")[0]) - boundaries[0], int(stripped_feature_loc.split("..")[-1]) - boundaries[0])
                else:
                    if boundaries[0] < int(stripped_feature_loc) or int(stripped_feature_loc) > boundaries[1]:
                        idxes_to_remove.append(idx)
                    else:
                        feature.location = "{}".format(int(stripped_feature_loc) - boundaries[0])
        else:
            for idx, feature in enumerate(entry.features):
                stripped_feature_loc = strip_location(feature.location)
                if ".." in stripped_feature_loc:
                    if (int(stripped_feature_loc.split("..")[0]) > boundaries[1] and int(stripped_feature_loc.split("..")[0]) < boundaries[0]) \
                            or (int(stripped_feature_loc.split("..")[-1]) > boundaries[1] and int(stripped_feature_loc.split("..")[-1]) < boundaries[0]):
                        idxes_to_remove.append(idx)
                    else:
                        if int(stripped_feature_loc.split("..")[0]) <= boundaries[1]:
                            feature.location = "{}..{}".format(int(stripped_feature_loc.split("..")[0]) + (len(entry.sequence) - boundaries[0]), int(stripped_feature_loc.split("..")[-1]) + (len(entry.sequence) - boundaries[0]))
                        elif int(stripped_feature_loc.split("..")[0]) >= boundaries[0]:
                            feature.location = "{}..{}".format(int(stripped_feature_loc.split("..")[0]) - boundaries[0], int(stripped_feature_loc.split("..")[-1]) - boundaries[0])
                else:
                    if boundaries[0] > int(stripped_feature_loc) or int(stripped_feature_loc) < boundaries[1]:
                        idxes_to_remove.append(idx)

                    else:
                        if int(stripped_feature_loc) <= boundaries[1]:
                            feature.location = "{}".format(int(stripped_feature_loc) + (len(entry.sequence) - boundaries[0]))
                        elif int(stripped_feature_loc) >= boundaries[0]:
                            feature.location = "{}".format(int(stripped_feature_loc) - boundaries[0])
        for idx in idxes_to_remove[::-1]:
            del entry.features[idx]
        # print(idxes_to_remove)
        print(entry.features)
        # print(entry_seq[boundaries[0]:boundaries[1]])
        # print(len(entry_seq[boundaries[0]:boundaries[1]]))
    file.close()
    return entry.features


def strip_location(feature_loc):
    feature_loc_copy = feature_loc
    if "(" in feature_loc_copy:
        start = feature_loc_copy.index("(") + 1
        stop = feature_loc_copy.index(")")
        return feature_loc_copy[start:stop]
    elif ">" in feature_loc or "<" in feature_loc_copy:
        feature_loc_copy = feature_loc_copy.replace(">", "")
        feature_loc_copy = feature_loc_copy.replace("<", "")
        return feature_loc_copy
    else:
        return feature_loc_copy


def excavate_donor(donor_seq, res_enz):
    c1_lst = res_site_finder(motif=res_2s[res_enz][0], querry=donor_seq)
    c2_lst = res_site_finder(motif=res_2s[res_enz][0], querry=translate(donor_seq))
    print(c1_lst, c2_lst)
    c1 = c1_lst[0]
    c2 = c2_lst[0]
    primer3_pos.append(-1)
    if len(c1_lst) == 1 and len(c2_lst) == 1:
        res_len = len(res_2s[res_enz][0]) + res_2s[res_enz][1]
        if c1 + c2 < len(donor_seq):
            fwd = donor_seq[c1+res_len:len(donor_seq)-c2-res_len-res_2s[res_enz][2]]
            rev = translate(donor_seq)[c2+res_len:len(donor_seq)-c1-res_len-res_2s[res_enz][2]]
            return [fwd, rev]
        elif c1 + c2 > len(donor_seq):
            if c1 + res_len > len(donor_seq):
                fwd = donor_seq[(c1+res_len)-len(donor_seq):len(donor_seq)-c2-res_len-res_2s[res_enz][2]]
                rev = translate(donor_seq)[c2+res_len:len(donor_seq)-((c1+res_len)-len(donor_seq))-res_2s[res_enz][2]]
                return [fwd, rev]
            elif c2 + res_len > len(donor_seq):
                rev = translate(donor_seq)[(c2+res_len)-len(donor_seq):len(donor_seq)-c1-res_len-res_2s[res_enz][2]]
                fwd = donor_seq[c1+res_len:len(donor_seq)-((c2+res_len)-len(donor_seq))-res_2s[res_enz][2]]
                return [fwd, rev]
            else:
                fwd = donor_seq[c1+res_len:] + donor_seq[:len(donor_seq)-c2-res_len-res_2s[res_enz][2]]
                rev = translate(donor_seq)[c2+res_len:] + translate(donor_seq)[:len(donor_seq)-c1-res_len-res_2s[res_enz][2]]
                return [fwd, rev]
    else:
        raise ValueError(f"Too many res sites found for {res_enz}")


def excavate_donor_ohs(donor_seq, res_enz):
    c1_lst = res_site_finder(motif=res_2s[res_enz][0], querry=donor_seq)
    c2_lst = res_site_finder(motif=res_2s[res_enz][0], querry=translate(donor_seq))
    print(c1_lst, c2_lst)
    c1 = c1_lst[0]
    c2 = c2_lst[0]
    if len(c1_lst) == 1 and len(c2_lst) == 1:
        res_len = len(res_2s[res_enz][0]) + res_2s[res_enz][1]
        if c1 + c2 < len(donor_seq) or (c1 + res_len + res_2s[res_enz][2] <= len(donor_seq) and c2 + res_len + res_2s[res_enz][2] <= len(donor_seq)):
            fwd = donor_seq[c1+res_len:c1+res_len+res_2s[res_enz][2]]
            rev = translate(donor_seq)[c2+res_len:c2+res_len+res_2s[res_enz][2]]
            return [rev, fwd]
        elif c1 + c2 > len(donor_seq):
            if c1 + res_len > len(donor_seq):
                fwd = donor_seq[(c1+res_len)-len(donor_seq):(c1+res_len)-len(donor_seq)+res_2s[res_enz][2]]
                rev = translate(donor_seq)[c2+res_len:c2+res_len+res_2s[res_enz][2]]
                return [rev, fwd]
            elif c2 + res_len > len(donor_seq):
                rev = translate(donor_seq)[(c2+res_len)-len(donor_seq):(c2+res_len)-len(donor_seq)+res_2s[res_enz][2]]
                fwd = donor_seq[c1+res_len:c1+res_len+res_2s[res_enz][2]]
                return [rev, fwd]
            elif c1 + res_len + res_2s[res_enz][2] > len(donor_seq):
                fwd = donor_seq[c1 + res_len:] + donor_seq[:res_2s[res_enz][2] - len(donor_seq[c1 + res_len:])]
                rev = translate(donor_seq)[c2+res_len:c2+res_len+res_2s[res_enz][2]]
                return [rev, fwd]
            elif c2 + res_len + res_2s[res_enz][2] > len(donor_seq):
                fwd = donor_seq[c1+res_len:c1+res_len+res_2s[res_enz][2]]
                rev = translate(donor_seq)[c2 + res_len:] + translate(donor_seq)[:res_2s[res_enz][2] - len(translate(donor_seq)[c2 + res_len:])]
                return [rev, fwd]
    else:
        raise ValueError(f"Too many or no res sites found for {res_enz}")


# checks if donor overhangs match predetermined overhangs of another donor, scaffold start or plasmid
# needs some attention
def excavate_problematic_ohs(conn_seqs, conn_type, res_enz):
    # donor-plasmid
    if conn_type.split("-")[0][0] == "d" and conn_type.split("-")[1] == "plas":
        oh_d = excavate_donor_ohs(conn_seqs[0], res_enz=res_enz)
        oh_plas = excavate_donor_ohs(conn_seqs[1], res_enz=res_enz)
        print(oh_d, oh_plas)
        if oh_d[0] == translate(oh_plas[1]):
            return [oh_d[0], oh_plas[1]]
        else:
            raise ValueError(f"Overhangs {oh_d[0]} and {translate(oh_plas[1])} don't match")
    # plasmid-donor
    elif conn_type.split("-")[1][0] == "d" and conn_type.split("-")[0] == "plas":
        oh_d = excavate_donor_ohs(conn_seqs[1], res_enz=res_enz)
        oh_plas = excavate_donor_ohs(conn_seqs[0], res_enz=res_enz)
        print(oh_plas, oh_d)
        if oh_plas[0] == translate(oh_d[1]):
            return [oh_plas[0], oh_d[1]]
        else:
            raise ValueError(f"Overhangs {oh_plas[0]} and {translate(oh_d[1])} don't match")
    # donor-scaffold
    elif conn_type.split("-")[0][0] == "d" and conn_type.split("-")[1] == "scaf":
        oh_d = excavate_donor_ohs(conn_seqs[0], res_enz=res_enz)
        oh_scaf = [translate(conn_seqs[1][:res_2s[res_enz][2]]), conn_seqs[1][:res_2s[res_enz][2]]]
        print(oh_d, oh_scaf)
        if oh_d[0] == translate(oh_scaf[1]):
            return [oh_d[0], oh_scaf[1]]
        else:
            raise ValueError(f"Overhangs {oh_d[0]} and {translate(oh_scaf[1])} don't match")
    # donor-donor
    elif conn_type.split("-")[0][0] == "d" and conn_type.split("-")[1][0] == "d":
        oh_d_1 = excavate_donor_ohs(conn_seqs[0], res_enz=res_enz)
        oh_d_2 = excavate_donor_ohs(conn_seqs[1], res_enz=res_enz)
        print(oh_d_1, oh_d_2)
        if oh_d_1[0] == translate(oh_d_2[1]):
            return [oh_d_1[0][::-1], oh_d_2[1]]
        else:
            raise ValueError(f"Overhangs {oh_d_1[0]} and {translate(oh_d_2[1])} don't match")
    elif conn_type.split("-")[1] == "dsee" and conn_type.split("-")[0] == "scaf":
        oh_scaf = [translate(conn_seqs[1][-res_2s[res_enz][2]:]), conn_seqs[1][-res_2s[res_enz][2]:]]
        oh_seed = excavate_donor_ohs(donor_seq=conn_seqs[0], res_enz=res_enz)
        print(oh_seed, oh_scaf)
        if oh_scaf[0] == translate(oh_seed[1]):
            return [oh_scaf[0], oh_seed[1]]
        else:
            raise ValueError(f"Overhangs {oh_scaf[0]} and {translate(oh_seed[1])} don't match")
    else:
        return ["rev", "fwd"]


def detect_illegal_res_sites(motif, querry):
    fwd_sites = res_site_finder(motif=motif, querry=querry)
    rev_sites = res_site_finder(motif=motif, querry=translate(querry))
    if len(fwd_sites) <= 1 and len(rev_sites) <= 1 and (len(fwd_sites) + len(rev_sites) == 2 or len(fwd_sites) + len(rev_sites) == 0):
        return False
    else:
        return True


def main():
    sequence = "cgaGTCTTCtagctacgatcgactgactagccatcagctagctacgatagctagctacgatcgatcatcagctacgatcgatcgactagctagctagcatcgacgatcgatcgactgaGAAGACctagctacgatcgatcgactagctac"
    # prims = excavate_donor(donor_seq="GCTACGATCGACTGAAGACACGGGCCACTAGCTACGATCAGCTACGATCGACTTGTGACGTCTTCGATCGACTACGATC", res_enz="Bbsl")
    # ohs = excavate_donor_ohs(donor_seq="GCTACGATCGACTGAAGACACGGGCCACTAGCTACGATCAGCTACGATCGACTTGTGACGTCTTCGATCGACTACGATC", res_enz="Bbsl")
    # ohs = find_overhang(seq="GCTACGATCGACTGAAGACACGGGCCACTAGCTACGATCAGCTACGATCGACTTGTGACGTCTTCGATCGACTACGATC", conn_type="seed-dsca", res_enz="Bbsl")


if __name__ == "__main__":
    main()

