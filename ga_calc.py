# set of functions for calculating gene-alignment overlap and distance metrics

# overlap between gene and alignment
def ga_ovlp(gene_s, gene_e, aln): # ie gene start, gene end, alignment paf record
    return min(gene_e, aln.te) - max(gene_s, aln.ts)

# distance between alignments on query ie read
def a_dist(aln0, aln1):
    return min(aln0.qe, aln1.qe) - max(aln0.qs, aln1.qs)

# distance between targets if on same chromosome
def t_dist(aln0, aln1):
    assert(aln0.t == aln1.t), "Alignments need to be on same chromosome for this metric!"
    return min(aln0.te, aln1.te) - max(aln0.ts, aln1.ts)

# percentage of overlapping alignments. as a percentage of aln1 so order of args matters. returns percentage as xx.xx%
# order matters, returns percentage with respect to aln1
def q_overlap(aln0, aln1):
    return (min(aln0.qe, aln1.qe) - max(aln0.qs, aln1.qs))/(aln1.qe - aln1.qs) * 100

# function to determine strand given alignment:
def det_strand(aln):
    if aln.rv == True:
        strand = '-'
    else:
        strand = '+'
    return strand

# function to determine breakpoints...given two alignments (needs paf module)
# order matters, returns the breakpoint of whatever is listed first ie aln0
def det_bp(aln0, aln1):
    # first determine which alignment is first on the read ie reorder if necessary
    if aln0.qs == aln1.qs:
        return None
    elif aln0.qs < aln1.qs:
        temp0 = aln0
        temp1 = aln1
    else:
        temp0 = aln1
        temp1 = aln0
        # temp0 = aln0
        # aln0 = aln1
        # aln1 = temp
    # define the strands
    temp0_s = det_strand(temp0)
    temp1_s = det_strand(temp1)
    # get the strand pair to later define the breakpoints
    orts = temp0_s+temp1_s
    if orts == '++':
        bp0 = temp0.te
        bp1 = temp1.ts
    elif orts == '+-':
        bp0 = temp0.te
        bp1 = temp1.te
    elif orts == '-+':
        bp0 = temp0.ts
        bp1 = temp1.ts
    elif orts == '--':
        bp0 = temp0.ts
        bp1 = temp1.te
    if aln0 == temp0:
        return bp0
    elif aln0 == temp1:
        return bp1

# function to return boolean of distance between alignments on query
def a_dist_tf(aln0, aln1, max_overlap, max_gap):
    # raw distance is negative if the alignments are separated, that's why we multiply by -1
    # overlapping
    if a_dist(aln0, aln1) >= 0:
        if (a_dist(aln0, aln1) > max_overlap):
            return False
        else:
            return True
    # separated alignments
    else:
        # convert to a positive number
        if ((a_dist(aln0, aln1) * -1) > max_gap):
            return False
        else:
            return True

# function to return boolean of overlapping alignments on target if on same chr
def t_dist_tf(aln0, aln1, max_gene_overlap):
    if aln0.t == aln1.t:
        # overlapping
        if t_dist(aln0, aln1) >= 0:
            if t_dist(aln0, aln1) > max_gene_overlap:
                return False
            else:
                return True
        # separated
        else:
            return True

# function to return boolean of overlapping portion % with respect to alignment
# order matters, returns percentage with respect to aln1
def q_overlap_tf(aln0, aln1, qmax_overlap):
    # overlapping
    if q_overlap(aln0, aln1) >= 0:
        if q_overlap(aln0, aln1) > qmax_overlap:
            return False
        else:
            return True
    # separated
    else:
        return True

# function to return boolean of breakpoint location
# order matters, returns aln0 breakpoint evaluation
# aln0 should align to gene0
def bp_tf(aln0, aln1, gene0s, gene0e, bp_win):
    if det_bp(aln0, aln1) == None:
        return True
    elif (det_bp(aln0, aln1) < (gene0s - bp_win)) | (det_bp(aln0, aln1) > (gene0e + bp_win)):
        return False
    else:
        return True

# function to return boolean of if there are enough supporting reads
def mnc_tf(ct, mnc):
    if ct < mnc:
        return False
    else:
        return  True
