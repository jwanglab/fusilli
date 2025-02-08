
'''
PAF:
1 string  Query sequence name
2 int     Query sequence length
3 int     Query start (0-based)
4 int     Query end (0-based)
5 char    Relative strand: "+" or "-"
6 string  Target sequence name
7 int     Target sequence length
8 int     Target start on original strand (0-based)
9 int     Target end on original strand (0-based)
10  int   Number of residue matches
11  int   Alignment block length
12  int   Mapping quality (0-255; 255 for missing)

3201a1ed-2ad8-4ecb-bf52-89ab2f7b858a	12007	10789	12007	+	dmel_p_element_2907bp	2907	0	1218	1174	1244	60	NM:i:70	ms:i:2048	AS:i:2048	nn:i:0	tp:A:P	cm:i:146	s1:i:945	s2:i:0	de:f:0.0401	cg:Z:47M1I30M1I20M1I3M2I7M1D11M1I90M1D25M1I38M2D61M4I5M1D64M1D50M3D67M1D1M1D8M1D70M2D6M2I50M1I86M4D43M2I4M1D11M2D153M1D46M3D11M1D4M1I22M1I11M4I87M2I17M2I44M
'''

class aln:
  def __init__(self, line):
    fields = line.strip().split('\t')
    self.q = fields[0]
    self.ql = int(fields[1])
    self.qs = int(fields[2])
    self.qe = int(fields[3])
    self.rv = (fields[4] == '-')
    self.t = fields[5]
    self.tl = int(fields[6])
    self.ts = int(fields[7])
    self.te = int(fields[8])
    self.m = int(fields[9])
    self.a = int(fields[10])
    self.acc = float(self.m) / self.a
    self.remainder = fields[12:]

  def parse_flags(self):
    # flags can have types i, f, A (char), Z (string), H (hex array 91A4F2...), B (numeric array i9,5,1,0,2,8)
    self.flags = {f.split(':')[0] : (int(f.split(':')[2]) if f.split(':')[1] == 'i' else (float(f.split(':')[2]) if f.split(':')[1] == 'f' else f.split(':')[2])) for f in self.remainder}
    if "cg" in self.flags:
      val = ""
      cigar_string = self.flags["cg"]
      self.cigar = []
      for c in cigar_string:
        if c in "0123456789":
          val += c
        else:
          self.cigar.append((int(val), c))
          val = ""

  def reverse(self):
    self.rv = not self.rv
    tmp = self.qs
    self.qs = self.ql - self.qe
    self.qe = self.ql - tmp

  def __str__(self):
    return "{} ({} bp) {}:{} ({}) <-> {} ({} bp) {}:{} ({}) (~{}bp x {:.2f}%)".format(self.q, self.ql, self.qs, self.qe, '-' if self.rv else '+', self.t, self.tl, self.ts, self.te, '+', self.qe - self.qs, self.acc * 100)
