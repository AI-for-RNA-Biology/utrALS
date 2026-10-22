"""CLIP dataset and RBP names: parse `<RBP>_<replicate>_<cell line>` names, map them to gene symbols."""
import re


# Extract RBP base names and replicate numbers
def parse_rbp_name(rbp_name):
    """Parse RBP name to extract base name and replicate number.
    Handles formats: RBPname_N and RBPname_N_suffix (e.g. UPF1_1_K652 -> ('UPF1', 1))
    """
    # Pattern to match RBPname_N (at end of string)
    match = re.match(r'^(.+)_(\d+)$', rbp_name)
    if match:
        return match.group(1), int(match.group(2))
    # Pattern to match RBPname_N_suffix (replicate number followed by another _-separated field)
    match = re.match(r'^(.+)_(\d+)_.+$', rbp_name)
    if match:
        return match.group(1), int(match.group(2))
    return rbp_name, None


def get_clip_name(clip_name, merge_ptb=True):
    """CLIP dataset name -> RBP (+ cell line for eCLIP). merge_ptb: every PTB / PTBP1 CLIP -> "PTB"."""
    if merge_ptb and "PTB" in clip_name:
        return "PTB"
    if "bestof" in clip_name:
        return clip_name.split('-')[0].upper()
    elif "2016" in clip_name:
        return clip_name.split('-')[0].upper()
    elif "all" not in clip_name and "iCLIP" not in clip_name:
        return clip_name.split('_')[0] + "_" + clip_name.split('_')[2]
    else:
        return clip_name.split('-')[1].upper()


def get_RBP_name(clip_name):
    if "bestof" in clip_name:
        return clip_name.split('-')[0].upper()
    elif "2016" in clip_name:
        return clip_name.split('-')[0].upper()
    elif "all" not in clip_name:
        return clip_name.split('_')[0]        
    else:
        return clip_name.split('-')[1].upper()



# RBP names used in the CLIP lists -> official gene symbols of the expression tables
RBP_TO_SYMBOL = {"TDP43": "TARDBP", "CPSF160": "CPSF1", "HUR": "ELAVL1", "CGI99": "RTRAF"}
