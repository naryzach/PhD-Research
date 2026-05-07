import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
import matplotlib.patches as patches

# Sequences
sequences = {
    'Human_MMP2': 'MEALMARGALTGPLRALCLLGCLLSHAAAAPSPIIKFPGDVAPKTDKELAVQYLNTFYGCPKESCNLFVLKDTLKKMQKFFGLPQTGDLDQNTIETMRKPRCGNPDVANYNFFPRKPKWDKNQITYRIIGYTPDLDPETVDDAFARAFQVWSDVTPLRFSRIHDGEADIMINFGRWEHGDGYPFDGKDGLLAHAFAPGTGVGGDSHFDDDELWTLGEGQVVRVKYGNADGEYCKFPFLFNGKEYNSCTDTGRSDGFLWCSTTYNFEKDGKYGFCPHEALFTMGGNAEGQPCKFPFRFQGTSYDSCTTEGRTDGYRWCGTTEDYDRDKKYGFCPETAMSTVGGNSEGAPCVFPFTFLGNKYESCTSAGRSDGKMWCATTANYDDDRKWGFCPDQGYSLFLVAAHEFGHAMGLEHSQDPGALMAPIYTYTKNFRLSQDDIKGIQELYGASPDIDLGTGPTPTLGPVTPEICKQDIVFDGIAQIRGEIFFFKDRFIWRTVTPRDKPMGPLLVATFWPELPEKIDAVYEAPQEEKAVFFAGNEYWIYSASTLERGYPKPLTSLGLPPDVQRVDAAFNWSKNKKTYIFAGDKFWRYNEVKKKMDPGFPKLIADAWNAIPDNLDAVVDLQGGGHSYFFKGAYYLKLENQSLKSVKFGSIKSDWLGC',
    'Mouse_MMP2': 'MEARVAWGALAGPLRVLCVLCCLLGRAIAAPSPIIKFPGDVAPKTDKELAVQYLNTFYGCPKESCNLFVLKDTLKKMQKFFGLPQTGDLDQNTIETMRKPRCGNPDVANYNFFPRKPKWDKNQITYRIIGYTPDLDPETVDDAFARALKVWSDVTPLRFSRIHDGEADIMINFGRWEHGDGYPFDGKDGLLAHAFAPGTGVGGDSHFDDDELWTLGEGQVVRVKYGNADGEYCKFPFLFNGREYSSCTDTGRSDGFLWCSTTYNFEKDGKYGFCPHEALFTMGGNADGQPCKFPFRFQGTSYNSCTTEGRTDGYRWCGTTEDYDRDKKYGFCPETAMSTVGGNSEGAPCVFPFTFLGNKYESCTSAGRNDGKVWCATTTNYDDDRKWGFCPDQGYSLFLVAAHEFGHAMGLEHSQDPGALMAPIYTYTKNFRLSHDDIKGIQELYGPSPDADTDTGTGPTPTLGPVTPEICKQDIVFDGIAQIRGEIFFFKDRFIWRTVTPRDKPTGPLLVATFWPELPEKIDAVYEAPQEEKAVFFAGNEYWVYSASTLERGYPKPLTSLGLPPDVQQVDAAFNWSKNKKTYIFAGDKFWRYNEVKKKMDPGFPKLIADSWNAIPDNLDAVVDLQGGGHSYFFKGAYYLKLENQSLKSVKFGSIKSDWLGC',
    'Human_MMP9': 'MSLWQPLVLVLLVLGCCFAAPRQRQSTLVLFPGDLRTNLTDRQLAEEYLYRYGYTRVAEMRGESKSLGPALLLLQKQLSLPETGELDSATLKAMRTPRCGVPDLGRFQTFEGDLKWHHHNITYWIQNYSEDLPRAVIDDAFARAFALWSAVTPLTFTRVYSRDADIVIQFGVAEHGDGYPFDGKDGLLAHAFPPGPGIQGDAHFDDDELWSLGKGVVVPTRFGNADGAACHFPFIFEGRSYSACTTDGRSDGLPWCSTTANYDTDDRFGFCPSERLYTQDGNADGKPCQFPFIFQGQSYSACTTDGRSDGYRWCATTANYDRDKLFGFCPTRADSTVMGGNSAGELCVFPFTFLGKEYSTCTSEGRGDGRLWCATTSNFDSDKKWGFCPDQGYSLFLVAAHEFGHALGLDHSSVPEALMYPMYRFTEGPPLHKDDVNGIRHLYGPRPEPEPRPPTTTTPQPTAPPTVCPTGPPTVHPSERPTAGPTGPPSAGPTGPPTAGPSTATTVPLSPVDDACNVNIFDAIAEIGNQLYLFKDGKYWRFSEGRGSRPQGPFLIADKWPALPRKLDSVFEERLSKKLFFFSGRQVWVYTGASVLGPRRLDKLGLGADVAQVTGALRSGRGKMLLFSGRRLWRFDVKAQMVDPRSASEVDRMFPGVPLDTHDVFQYREKAYFCQDRFYWRVSSRSELNQVDQVGYVTYDILQCPED',
    'Human_ADAM17': 'MRQSLLFLTSVVPFVLAPRPPDDPGFGPHQRLEKLDSLLSDYDILSLSNIQQHSVRKRDLQTSTHVETLLTFSALKRHFKLYLTSSTERFSQNFKVVVVDGKNESEYTVKWQDFFTGHVVGEPDSRVLAHIRDDDVIIRINTDGAEYNIEPLWRFVNDTKDKRMLVYKSEDIKNVSRLQSPKVCGYLKVDNEELLPKGLVDREPPEELVHRVKRRADPDPMKNTCKLLVVADHRFYRYMGRGEESTTTNYLIELIDRVDDIYRNTSWDNAGFKGYGIQIEQIRILKSPQEVKPGEKHYNMAKSYPNEEKDAWDVKMLLEQFSFDIAEEASKVCLAHLFTYQDFDMGTLGLAYVGSPRANSHGGVCPKAYYSPVGKKNIYLNSGLTSTKNYGKTILTKEADLVTTHELGHNFGAEHDPDGLAECAPNEDQGGKYVMYPIAVSGDHENNKMFSNCSKQSIYKTIESKAQECFQERSNKVCGNSRVDEGEECDPGIMYLNNDTCCNSDCTLKEGVQCSDRNSPCCKNCQFETAQKKCQEAINATCKGVSYCTGNSSECPPPGNAEDDTVCLDLGKCKDGKCIPFCEREQQLESCACNETDNSCKVCCRDLSGRCVPYVDAEQKNLFLRKGKPCTVGFCDMNGKCEKRVQDVIERFWDFIDQLSINTFGKFLADNIVGSVLVFSLIFWIPFSILVHCVDKKLDKQYESLSLFHPSNVEMLSSMDSASVRIIKPFPAPQTPGRLQPAPVIPSAPAAPKLDHQRMDTIQEDPSTDSHMDEDGFEKDPFPNSSTAAKSFEDLTDHPVTRSEKAASFKLQRQNRVDSKETEC',
    'Rat_ADAM17': 'MRQRLLFLTTLVPFVLAPRPPEEPGSGSHLRLEKLDSLLSDYDILSLSNIQQHSIRKRDLQSATHLETLLTFSALKRHFKLYLTSSTERFSQNLRVVVVDGKEESEYSVKWQDFFSGHVVGEPDSRVLAHIGDDDVTVRINTDGAEYNIEPLWRFVNDTKDKRMLVYKSEDIKDFSRLQSPKVCGYLNADSEELLPKGLIDREPSEEFVRRVKRRAEPNPLKNTCKLLVVADHRFYKYMGRGEESTTTNYLIELIDRVDDIYRNTSWDNAGFKGYGVQIEQIRILKSPQEVKPGERHFNMAKSFPNEEKDAWDVKMLLEQFSLDIAEEASKVCLAHLFTYQDFDMGTLGLAYVGSPRANSHGGVCPKAYYNPGVKKNIYLNSGLTSTKNYGKTILTKEADLVTTHELGHNFGAEHDPDGLAECAPNEDQGGKYVMYPIAVSGDHENNKMFSNCSKQSIYKTIESKAQECFQERSNKVCGNSRVDEGEECDPGIMYLNNDTCCNSDCTLKPGVQCSDRNSPCCKNCQFETAQKKCQEAINATCKGVSYCTGNSSECPPPGDAEDDTVCLDLGKCKAGKCIPFCKREQELESCACADTDNSCKVCCRNLSGPCVPYVDAEQKNLFLRKGKPCTVGFCDMNGKCEKRVQDVIERFWDFIDQLSINTFGKFLADNIVGSVLVFSLIFWIPFSILVHCVDKKLDKQYESLSLFHHSNIEMLSSMDSASVRIIKPFPAPQTPGRLQALQPAAMMPPVSAAPKLDHQRMDTIQEDPSTDSHVDDDGFEKDPFPNSSAAAKSFEDLTDHPVTRSEKAASFKLQRQSRVDSKETEC'
}

# Granular domain definitions from UniProt features
domains = {
    'MMP2': [
        (100, 107, 'Cys-Switch', '#E1BEE7'),
        (110, 221, 'Col-like 1', '#BBDEFB'),
        (222, 396, 'Col-binding', '#81C784'),
        (228, 276, 'Fn-II 1', '#4CAF50'),
        (286, 334, 'Fn-II 2', '#4CAF50'),
        (344, 392, 'Fn-II 3', '#4CAF50'),
        (397, 465, 'Col-like 2', '#BBDEFB'),
        (472, 516, 'Hpx 1', '#FFECB3'),
        (517, 563, 'Hpx 2', '#FFECB3'),
        (565, 613, 'Hpx 3', '#FFECB3'),
        (614, 660, 'Hpx 4', '#FFECB3')
    ],
    'MMP9': [
        (97, 104, 'Cys-Switch', '#E1BEE7'),
        (225, 273, 'Fn-II 1', '#4CAF50'),
        (283, 331, 'Fn-II 2', '#4CAF50'),
        (342, 390, 'Fn-II 3', '#4CAF50'),
        (518, 563, 'Hpx 1', '#FFECB3'),
        (564, 608, 'Hpx 2', '#FFECB3'),
        (610, 657, 'Hpx 3', '#FFECB3'),
        (658, 704, 'Hpx 4', '#FFECB3'),
        (431, 508, 'Disordered', '#CFD8DC')
    ],
    'ADAM17': [
        (182, 189, 'Cys-Switch', '#E1BEE7'),
        (223, 474, 'Peptidase', '#BBDEFB'),
        (475, 563, 'Disintegrin', '#F8BBD0'),
        (603, 671, 'Crambin-like', '#E1BEE7')
    ]
}

aa_map = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}

# Average amino acid weights (isotopic)
AA_WEIGHTS = {
    'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
    'C': 103.1388, 'Q': 128.1307, 'E': 129.1139, 'G': 57.0519,
    'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
    'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782,
    'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
}
WATER_WEIGHT = 18.01524

def calculate_mw(seq):
    if not seq: return 0.0
    weight = sum(AA_WEIGHTS.get(aa, 0) for aa in seq.upper())
    if weight > 0:
        weight += WATER_WEIGHT
    return weight / 1000.0 # Return in kDa

def get_aa_label(seq, pos):
    if 0 <= pos-1 < len(seq):
        aa = seq[pos-1]
        return f"{aa_map.get(aa, aa)}{pos}"
    return str(pos)

def plot_family_alignment(family_name, reference_seq, ortholog_seq=None, ref_label="Human", orth_label=None, constructs=[], filename="alignment.png"):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    
    if ortholog_seq:
        aln = aligner.align(reference_seq, ortholog_seq)[0]
        aln_ref = aln[0]
        aln_orth = aln[1]
    else:
        aln_ref = reference_seq
        aln_orth = None
        
    length = len(aln_ref)
    
    # Extra large width and height to fit massive fonts
    fig, ax = plt.subplots(figsize=(40, 14 + 1.5 * (2 + len(constructs))))
    
    def get_map(aln_seq):
        if not aln_seq: return None
        raw_to_aln = []
        for i, char in enumerate(aln_seq):
            if char != '-':
                raw_to_aln.append(i)
        return raw_to_aln

    ref_map = get_map(aln_ref)
    orth_map = get_map(aln_orth) if aln_orth else None

    # Base track positions
    y_ref = 18
    y_orth = 14.5 if orth_label else None
    y_const_start = 10.0
    
    def add_boundary_label(ax, x, y, label, direction='up', level=1, color='#1A237E', fontsize=18):
        """Adds a massive residue label with a leader line."""
        # AA Levels use step 1.6
        offset = level * 1.6
        ha = 'center'
        va = 'bottom' if direction == 'up' else 'top'
        y_text = y + offset if direction == 'up' else y - offset
        
        bbox = dict(boxstyle="round,pad=0.4", facecolor='#F5F5F5', alpha=0.9, edgecolor='#BDBDBD', linewidth=1.0)
        ax.text(x, y_text, label, fontsize=fontsize, ha=ha, va=va, color=color, fontweight='black', bbox=bbox)
        ax.plot([x, x], [y, y_text], color='#9E9E9E', linestyle='-', linewidth=1.5, alpha=0.5)

    # Draw Reference Track
    ax.text(-60, y_ref + 0.6, f"{ref_label.upper()} FULL", fontweight='black', ha='right', va='center', fontsize=26, color='#333333')
    ref_mw = calculate_mw(reference_seq)
    ax.text(length + 20, y_ref + 0.6, f"{ref_mw:.1f} kDa", fontweight='black', ha='left', va='center', fontsize=24, color='#666666')
    ax.add_patch(patches.Rectangle((0, y_ref), length, 1.5, color='#EEEEEE', alpha=0.5))
    
    domain_list = domains.get(family_name, [])
    assigned_domain_levels = []
    assigned_res_labels = []
    
    def get_best_res_level(x_aln, current_assignments, padding=25):
        """Finds the first level (1, 2, or 3) where the label won't overlap."""
        for level in [1, 2, 3]:
            blocked = False
            for x_p, l_p in current_assignments:
                if l_p == level and abs(x_aln - x_p) < padding:
                    blocked = True
                    break
            if not blocked:
                return level
        return 1 # Fallback
    
    # First pass: Determine residue levels and calculate domain base offset
    assigned_res_labels = []
    for i, (start, end, label, color) in enumerate(domain_list):
        start_aln = ref_map[start-1]
        end_aln = ref_map[min(end-1, len(ref_map)-1)]
        
        res_level_start = get_best_res_level(start_aln, assigned_res_labels)
        assigned_res_labels.append((start_aln, res_level_start))
        
        res_level_end = get_best_res_level(end_aln, assigned_res_labels)
        assigned_res_labels.append((end_aln, res_level_end))

    max_res_used = max([l for x, l in assigned_res_labels]) if assigned_res_labels else 1
    # Domain base clears max_res_used * 1.6 + buffer
    domain_base = (max_res_used * 1.6) + 3.0
    
    # Second pass: Draw everything
    assigned_domain_levels = []
    for i, (start, end, label, color) in enumerate(domain_list):
        start_aln = ref_map[start-1]
        end_aln = ref_map[min(end-1, len(ref_map)-1)]
        center_aln = (start_aln + end_aln) / 2
        
        # Draw domain rectangle
        ax.add_patch(patches.Rectangle((start_aln, y_ref), end_aln - start_aln + 1, 1.5, color=color))
        
        # Determine domain tier
        blocked_levels = set()
        for s_p, e_p, l_p in assigned_domain_levels:
            if not (end_aln < s_p - 15 or start_aln > e_p + 15):
                blocked_levels.add(l_p)
        
        if 0 not in blocked_levels: d_tier = 0
        elif 1 not in blocked_levels: d_tier = 1
        else: d_tier = 2
        assigned_domain_levels.append((start_aln, end_aln, d_tier))
        
        # Draw domain label and connector
        y_domain_text = y_ref + domain_base + d_tier * 0.8
        ax.text(center_aln, y_domain_text, label.upper(), 
                fontweight='black', fontsize=24, ha='center', color='black',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=2))
        ax.plot([start_aln, end_aln], [y_domain_text, y_domain_text], color='#666666', linestyle='-', linewidth=1.5, alpha=0.3, zorder=0)
        
        # Vertical dotted lines at boundaries
        ax.axvline(start_aln, color='#CCCCCC', linestyle=':', linewidth=1.5, alpha=0.4, zorder=-1)
        ax.axvline(end_aln, color='#CCCCCC', linestyle=':', linewidth=1.5, alpha=0.4, zorder=-1)
        
        # Draw residue labels
        add_boundary_label(ax, start_aln, y_ref + 1.5, get_aa_label(reference_seq, start), level=assigned_res_labels[i*2][1], fontsize=18)
        add_boundary_label(ax, end_aln, y_ref + 1.5, get_aa_label(reference_seq, end), level=assigned_res_labels[i*2+1][1], fontsize=18)

    # Draw Ortholog Track
    if orth_label:
        ax.text(-60, y_orth + 0.6, f"{orth_label.upper()} FULL", fontweight='black', ha='right', va='center', fontsize=26, color='#333333')
        orth_mw = calculate_mw(ortholog_seq)
        ax.text(length + 20, y_orth + 0.6, f"{orth_mw:.1f} kDa", fontweight='black', ha='left', va='center', fontsize=24, color='#666666')
        ax.add_patch(patches.Rectangle((0, y_orth), length, 1.5, color='#EEEEEE', alpha=0.5))
        for start, end, label, color in domain_list:
            if start-1 < len(orth_map):
                start_aln = orth_map[start-1]
                end_aln = orth_map[min(end-1, len(orth_map)-1)]
                ax.add_patch(patches.Rectangle((start_aln, y_orth), end_aln - start_aln + 1, 1.5, color=color, alpha=0.3))
        # Mark differences
        for i in range(length):
            if aln_ref[i] != aln_orth[i] and aln_ref[i] != '-' and aln_orth[i] != '-':
                ax.add_patch(patches.Rectangle((i, y_orth), 1, y_ref - y_orth + 1.5, color='#FF1744', alpha=0.3))

    # Draw Constructs / Isoforms
    for idx, (c_label, c_parent_map, c_range) in enumerate(constructs):
        y_c = y_const_start - idx * 3.5
        start_res, end_res = c_range
        parent_seq = reference_seq if c_parent_map == "ref" else ortholog_seq
        
        ax.text(-60, y_c + 0.6, c_label.upper(), fontweight='black', ha='right', va='center', fontsize=24, color='#333333')
        c_mw = calculate_mw(parent_seq[start_res-1:end_res])
        ax.text(length + 20, y_c + 0.6, f"{c_mw:.1f} kDa", fontweight='black', ha='left', va='center', fontsize=22, color='#666666')
        
        ax.add_patch(patches.Rectangle((0, y_c), length, 1.5, color='#EEEEEE', alpha=0.3))
        
        parent_map = ref_map if c_parent_map == "ref" else orth_map
        
        start_aln = parent_map[start_res-1]
        end_aln = parent_map[min(end_res-1, len(parent_map)-1)]
        ax.add_patch(patches.Rectangle((start_aln, y_c), end_aln - start_aln + 1, 1.5, color='#F57C00', alpha=0.9))
        
        add_boundary_label(ax, start_aln, y_c, get_aa_label(parent_seq, start_res), direction='down', level=0.4, color='#D32F2F', fontsize=18)
        add_boundary_label(ax, end_aln, y_c, get_aa_label(parent_seq, end_res), direction='down', level=0.4, color='#D32F2F', fontsize=18)

    # Dynamic ylim top based on max level used
    max_d_tier = max([l for s, e, l in assigned_domain_levels]) if assigned_domain_levels else 0
    y_top = y_ref + domain_base + max_d_tier * 0.8 + 1.2
    y_bottom = y_const_start - (len(constructs) - 1) * 3.5 - 2.5
    
    ax.set_xlim(-220, length + 120)
    ax.set_ylim(y_bottom, y_top)
    ax.set_yticks([])
    ax.tick_params(axis='x', labelsize=24, length=10, width=2)
    ax.set_title(f"{family_name} FAMILY: MULTI-VERSION ALIGNMENT & RESIDUE MAPPING", fontsize=44, pad=20, fontweight='black')
    ax.set_xlabel("ALIGNED RESIDUE POSITION", fontsize=30, fontweight='black', labelpad=30)
    plt.tight_layout(pad=2.0)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# Family Groups
# MMP2
mmp2_const = [
    ("Enzo Human (Cat/Fib)", "ref", (110, 452)),
    ("Sino Mouse (Cat/Hpx)", "orth", (409, 662)),
    ("Isoform 1 (Start MQKFFG)", "ref", (145, 660)),
    ("Isoform 2 (Start MQYLNT)", "ref", (107, 660))
]
plot_family_alignment("MMP2", sequences['Human_MMP2'], sequences['Mouse_MMP2'], "Human", "Mouse", mmp2_const, "MMP2_Family_Comparison.png")

# ADAM17
adam17_const = [
    ("Enzo Human (Pro/Cat)", "ref", (18, 477)),
    ("Sino Rat (Full/Dis)", "orth", (1, 563))
]
plot_family_alignment("ADAM17", sequences['Human_ADAM17'], sequences['Rat_ADAM17'], "Human", "Rat", adam17_const, "ADAM17_Family_Comparison.png")

# MMP9
mmp9_const = [
    ("Enzo Human (Cat/Fib)", "ref", (107, 449)),
    ("Sino Human (Full)", "ref", (1, 707))
]
plot_family_alignment("MMP9", sequences['Human_MMP9'], None, "Human", None, mmp9_const, "MMP9_Family_Comparison.png")

print("High-precision UniProt family alignment diagrams generated.")
