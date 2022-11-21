from datasets import res_2s
from helpers import build_prim_output, res_site_finder, build_conn_map, input_handler, \
    find_overhang, resort_overhangs, build_genbank_output, \
    build_prim_vis, translate, build_primer_new
from read_input import read_in_file_csv
import numpy as np
import os


# custom_ohs is a list filled with [Overhang, Position] pairs for the positions of the custom overhangs
# position equals the end of the corresponding fragment in in_file
# custom_oh has to be the corresponding 5'-3' fwd overhang for the fragment at position + 1
def all_in_one_better(in_file, constr_nr):
    out_idx = in_file.split("_")[-1][:-4]
    # print(out_idx)
    # read file content
    # custom overhangs are defined as the 5'-3' overhang of the forward primer for the corresponding fragment
    contents, custom_creation, custom_ohs, meta_data, frag_names = read_in_file_csv(in_file)
    # print(contents, custom_creation, custom_ohs, meta_data)
    # additional variable assignments for parameter input read
    # backup parameters def all_in_one_better(in_file, out_file, res_enz, rand_b, oh_len, gc_cont, anneal_len, max_len):

    out_file = meta_data["out_file"]
    res_enz = meta_data["res_enz"]
    # rand_b = int(meta_data["rand_b"])
    oh_len = int(meta_data["oh_len"])
    # gc_cont = float(meta_data["gc_cont"])
    fillup_bases = meta_data["fillup_bases"]
    anneal_len = int(meta_data["anneal_len"])
    max_len = int(meta_data["max_len"])
    tm = float(meta_data["tm"])
    plas_gb = meta_data["plas_gb"]

    # create connection map
    conn_map = build_conn_map(contents)
    primer_list = []
    # create overhangs
    overhangs = []
    non_unique = False
    overhang_roll = 0
    while not non_unique:
        # save all already created overhangs for doubling detection
        overhangs = []
        for idx, entry in enumerate(conn_map):
            if len(custom_ohs[idx]) > 0:
                if "plas" in entry:
                    raise ValueError("Overhangs for plasmid connection in this case {} cannot be chosen because they are based on the plasmid cutsites".format(entry))
                # disabled because scaffold scars should be possible is put in by the user
                # elif entry[-4:] == "scaf":
                #     raise ValueError("Overhangs for fragment-scaffold connection in this case {} cannot be chosen because they are based on the scaffold start".format(entry))
                elif entry.split("-")[0][0] == "d" or entry.split("-")[1][0] == "d":
                    raise ValueError("Overhangs for donor connection in this case {} cannot be chosen because they are based on the donor sequence".format(entry))
                else:
                    overhangs.append([translate(custom_ohs[idx]), custom_ohs[idx]])
            else:
                # give relevant sequence for primer design
                relevant_conn_seq = input_handler(conn=entry, content=contents, idx=idx)
                print(entry, relevant_conn_seq)
                # create overhangs for connection - not for single fragment
                overhang = find_overhang(seq=relevant_conn_seq, conn_type=entry, res_enz=res_enz, length=oh_len)
                overhangs.append(overhang)
        # enlisted list containing all overhangs
        double_check = []
        print(overhangs, conn_map)
        for pair in overhangs:
            # reverse overhangs are inverted because they are saved as 5'-3' so e.g. fwd TATA would be equal to rev TATA
            double_check.append(pair[0])
            double_check.append(pair[1])
        # if list contains doubles -> repeat procedure
        if len(np.unique(np.unique(double_check, return_counts=True)[1])) == 1:
            non_unique = True
        overhang_roll += 1
        if overhang_roll > 5:
            print("returned faulty ohs")
            return "non_unique", overhangs
            raise RuntimeWarning("100 attempts of overhang creation did not result in unique set, please review "
                                 "fixed overhangs of plasmid and scaffold as well as custom overhangs")
    # sort overhangs so overhang pair correspond to one single fragment instead of connections
    # print(overhangs)
    oh_new = resort_overhangs(overhangs=overhangs, conn_map=conn_map)
    print(oh_new)
    # build primers
    # note that this part works with precomputed overhangs
    # newly created res sites can only be killed by changing random bases at primer start and after res site
    for idx, fragment in enumerate(contents):
        # print(fragment[0])
        # plasmid primers are not implemented right now
        if "plasmid" not in fragment[0]:
            if "[" in fragment[1] and "]" in fragment[1]:
                use_primer_3 = True
            else:
                use_primer_3 = False
            # create primers
            # rebuild_count = -1
            primers = build_primer_new(sequence=fragment[1],
                                       seq_type=fragment[0],
                                       overhangs=oh_new[idx],
                                       res_enz=res_enz,
                                       primer_len=anneal_len,
                                       tm=tm,   # tm should be included somehow from the user input
                                       custom=custom_creation[idx],
                                       use_primer_3=use_primer_3,
                                       fillup_bases=fillup_bases,
                                       max_len=max_len)

            if ((len(fragment[1]) > max_len and "anneal" not in custom_creation[idx]) or custom_creation[idx] == "PCR") and fragment[0][0] != "d":
                if len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[0])) != 1 or len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[1])) != 1:
                    raise ValueError("undesired cut sites detected in {} (fwd:{}, rev:{}). Please change Fillup Bases and try again".format(fragment[0], len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[0])), len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[1]))))
                else:
                    primer_list.append([fragment[0], primers])
            else:
                # check both primers for newly introduced cutsites (can't be rerolled)
                # res site detection only prevents the primers from being built it does not kill the programm
                if len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[0])) != 0 or len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[1])) != 0:
                    raise ValueError("undesired cut sites detected in {} (fwd:{}, rev:{}). Please change Fillup Bases and try again".format(fragment[0], len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[0])), len(res_site_finder(motif=res_2s[res_enz][0], querry=primers[1]))))
                else:
                    primer_list.append([fragment[0], primers])
        else:
            # plasmid primers are not yet implemented
            pass
    # these functions write primers and primer visualizations in the same (prev not existing) output file
    if not os.path.isfile(out_file):
        build_prim_output(primer_list=primer_list, content=contents, names=frag_names, filename="{}_{}".format(out_file, out_idx), constr_nr=constr_nr)
        build_prim_vis(primer_list=primer_list, content=contents, res_enz=res_enz, filename="{}_{}".format(out_file, out_idx), oh_len=oh_len, fillup_bases=fillup_bases)
    else:
        print("{} already exists".format(out_file))
    build_genbank_output(filename="{}_{}".format(out_file, out_idx), gb_plasmid_file=plas_gb, content=contents, res_enz=res_enz, oh_len=oh_len, ohs=overhangs, custom_ohs=custom_ohs)
    return primer_list, overhangs


def main():
    all_in_one_better(in_file="./data/example_input_meta.txt",
                      # out_file="./data/example_output",
                      # res_enz="Bbsl",
                      # max_len=60,
                      # rand_b=2,
                      # anneal_len=20,
                      # gc_cont=0.6,
                      # oh_len=4
                      )


if __name__ == "__main__":
    main()


