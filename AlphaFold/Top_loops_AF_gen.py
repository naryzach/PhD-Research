import json


# 1. canonical (wild-type) TIMP3 sequence
CANONICAL_TIMP3 = "CTCSPSHPQDAFCNSDIVIRAKVVGKKLVKEGPFGTLVYTIKQMKMYRGFTKMPHVQYIHTEASESLCGLKLEVNKYQYLLTGRVYDGKMYTGLCNFVERWDQLTLSQRKGLNYRYHLGCN"

# target protein sequences (MMPs and ADAMs)
TARGETS = {
    "ADAM10": "NTCQLYIQTDHLFFKYYGTREAVIAQISSHVKAIDTIYQTTDFSGIRNISFMVKRIRINTTADEKDPTNPFRFPNIGVEKFLELNSEQNHDDYCLAYVFTDRDFDDGVLGLAWVGAPSGSSGGICEKSKLYSDGKKKSLNTGIITVQNYGSHVPPKVSHITFAHEVGHNFGSPHDSGTECTPGESKNLGQKENGNYIMYARATSGDKLNNNKFSLCSIRNISQVLEKKRNNCFVESG",
    "ADAM17": "DPMKNTCKLLVVADHRFYRYMGRGEESTTTNYLIELIDRVDDIYRNTAWDNAGFKGYGIQIEQIRILKSPQEVKPGEKHYNMAKSYPNEEKDAWDVKMLLEQFSFDIAEEASKVCLAHLFTYQDFDMGTLGLAYVGSPRANSHGGVCPKAYYSPVGKKNIYLNSGLTSTKNYGKTILTKEADLVTTHELGHNFGAEHDPDGLAECAPNEDQGGKYVMYPIAVSGDHENNKMFSQCSKQSIYKTIESKAQECFQERS",
    "MMP2": "YNFFPRKPKWDKNQITYRIIGYTPDLDPETVDDAFARAFQVWSDVTPLRFSRIHDGEADIMINFGRWEHGDGYPFDGKDGLLAHAFAPGTGVGGDSHFDDDELWTLGEGQVGYSLFLVAAHEFGHAMGLEHSQDPGALMAPIYTYTKNFRLSQDDIKGIQELYGASPDGSDYKDDDDK",
    "MMP9": "FQTFEGDLKWHHHNITYWIQNYSEDLPRAVIDDAFARAFALWSAVTPLTFTRVYSRDADIVIQFGVAEHGDGYPFDGKDGLLAHAFPPGPGIQGDAHFDDDELWSLGKGQGYSLFLVAAHEFGHALGLDHSSNTEALMYPMYRFTEGPPLHKDDVNGIRHLY"
}

## --- GRAB THIS BLOCK FROM FILE --- ##
# lists of loop variants
AB_LOOPS = ["EGPVGE", "KGPVGE", "SGPTGE"]
C_LOOPS = ["ANTNLC", "ANTGLC", "ANTALC"]
EF_LOOPS = ["TSCD", "SSCD", "TACD"]

NAME_TAG = "ProteinMPNN_vADAM17"
## --------------------------------- ##

# --- Slice definitions (0-based indexing) ---
# The script will use these to piece the new sequence together.
AB_SLICE = (30, 36)
C_SLICE = (62, 68)
EF_SLICE = (92, 96)


def create_variant_sequence(base_seq, ab_loop=None, c_loop=None, ef_loop=None):
    """
    Creates a new TIMP3 variant sequence by replacing specified loops.
    This function correctly builds the new sequence piece by piece,
    avoiding index-shifting errors.
    """
    
    # Use canonical loop sequence if no variant is provided
    ab_insert = ab_loop if ab_loop else base_seq[AB_SLICE[0]:AB_SLICE[1]]
    c_insert = c_loop if c_loop else base_seq[C_SLICE[0]:C_SLICE[1]]
    ef_insert = ef_loop if ef_loop else base_seq[EF_SLICE[0]:EF_SLICE[1]]
    
    # Check for length mismatches as a warning
    if ab_loop and len(ab_loop) != (AB_SLICE[1] - AB_SLICE[0]):
        print(f"Warning: AB loop '{ab_loop}' length ({len(ab_loop)}) differs from canonical ({(AB_SLICE[1] - AB_SLICE[0])}).")
    if c_loop and len(c_loop) != (C_SLICE[1] - C_SLICE[0]):
        print(f"Warning: C loop '{c_loop}' length ({len(c_loop)}) differs from canonical ({(C_SLICE[1] - C_SLICE[0])}).")
    if ef_loop and len(ef_loop) != (EF_SLICE[1] - EF_SLICE[0]):
        print(f"Warning: EF loop '{ef_loop}' length ({len(ef_loop)}) differs from canonical ({(EF_SLICE[1] - EF_SLICE[0])}).")

    # Build sequence from its 7 parts
    part1 = base_seq[0:AB_SLICE[0]]
    part2_loop = ab_insert
    part3_linker = base_seq[AB_SLICE[1]:C_SLICE[0]]
    part4_loop = c_insert
    part5_linker = base_seq[C_SLICE[1]:EF_SLICE[0]]
    part6_loop = ef_insert
    part7_end = base_seq[EF_SLICE[1]:]

    return part1 + part2_loop + part3_linker + part4_loop + part5_linker + part6_loop + part7_end


def create_af_entry(name, seq1, seq2):
    """
    Formats a single job entry for the AlphaFold server JSON.
    """
    return {
        "name": name,
        "sequences": [
            {
                "proteinChain": {
                    "sequence": seq1,
                    "count": 1
                }
            },
            {
                "proteinChain": {
                    "sequence": seq2,
                    "count": 1
                }
            }
        ]
    }


def main():
    all_jobs = []
    
    # --- Generate Single-Loop Jobs ---

    # Process AB loops
    for loop_seq in AB_LOOPS:
        variant_timp = create_variant_sequence(CANONICAL_TIMP3, ab_loop=loop_seq)
        for target_name, target_seq in TARGETS.items():
            job_name = f"TIMP3_vs_{target_name}_AB_{loop_seq}_{NAME_TAG}"
            all_jobs.append(create_af_entry(job_name, variant_timp, target_seq))

    # Process C loops
    for loop_seq in C_LOOPS:
        variant_timp = create_variant_sequence(CANONICAL_TIMP3, c_loop=loop_seq)
        for target_name, target_seq in TARGETS.items():
            job_name = f"TIMP3_vs_{target_name}_C_{loop_seq}_{NAME_TAG}"
            all_jobs.append(create_af_entry(job_name, variant_timp, target_seq))

    # Process EF loops
    for loop_seq in EF_LOOPS:
        variant_timp = create_variant_sequence(CANONICAL_TIMP3, ef_loop=loop_seq)
        for target_name, target_seq in TARGETS.items():
            job_name = f"TIMP3_vs_{target_name}_EF_{loop_seq}_{NAME_TAG}"
            all_jobs.append(create_af_entry(job_name, variant_timp, target_seq))

    # --- Generate Combined-Loop Jobs ---    
    if not (AB_LOOPS and C_LOOPS and EF_LOOPS):
        print("Warning: One or more loop lists are empty. Skipping combined-loop job generation.")
    else:
        for rank in range(min(len(AB_LOOPS), len(C_LOOPS), len(EF_LOOPS))):
            top_ab = AB_LOOPS[rank]
            top_c = C_LOOPS[rank]
            top_ef = EF_LOOPS[rank]

            combined_timp = create_variant_sequence(
                CANONICAL_TIMP3,
                ab_loop=top_ab,
                c_loop=top_c,
                ef_loop=top_ef
            )
            
            for target_name, target_seq in TARGETS.items():
                job_name = f"TIMP3_vs_{target_name}_Combined0{rank+1}_{top_ab}_{top_c}_{top_ef}_{NAME_TAG}"
                all_jobs.append(create_af_entry(job_name, combined_timp, target_seq))

    # --- Write to JSON file ---
    output_filename = f"alphafold_jobs_toploops_{NAME_TAG}.json"
    with open(output_filename, 'w') as f:
        json.dump(all_jobs, f, indent=2)

    print(f"Success! Generated {len(all_jobs)} jobs in '{output_filename}'.")


if __name__ == "__main__":
    main()