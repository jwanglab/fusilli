class CandGene:
    def __init__(self, name, chrom, start, end):
        self.name = name
        self.chr = chrom
        self.s = start
        self.e = end
    def __hash__(self):
        return hash(self.name)
