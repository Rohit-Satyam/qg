#!/usr/bin/env python

import re
import os
import sys

"""
MIT License
Copyright (c) [2016] [Parashar Dhapola]
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

class QuadMotifFinder():
    """
    This class allows searching G-quadruplex motifs in given sequences.
    Please see the paper below for details and cite it if you use this code
    QuadBase2: web server for multiplexed guanine quadruplex mining and
    visualization
    doi: 10.1093/nar/gkw425
    """
    def __init__(self, fasta_file, savename=None, stem=3, loop_start=1,
                 loop_stop=7, greedy=True, bulge=0, is_seq=False,
                 strands=['+', '-'], verbose=True):
        self.fastaFile = fasta_file
        self.saveName = savename
        self.stemSize = stem
        self.minLoopSize = loop_start
        self.maxLoopSize = loop_stop
        self.greedySearch = greedy
        self.bulgeSize = bulge
        self.isSeq = is_seq
        self.strands = strands
        self.verbose = verbose
        self.outHandles = {}
        self._sanitize_params()
        self.resOv, self.resNov = [], []

    def _sanitize_params(self):
        try:
            self.stemSize = int(self.stemSize)
            self.minLoopSize = int(self.minLoopSize)
            self.maxLoopSize = int(self.maxLoopSize)
            self.bulgeSize = int(self.bulgeSize)
        except:
            raise TypeError('ERROR: Please provide only integer values for size parameters')
        try:
            self.greedySearch = bool(self.greedySearch)
        except:
            raise TypeError('ERROR: Please provide only boolean values for greedySearch (True/False)')
        if self.minLoopSize >= self.maxLoopSize:
            raise ValueError('ERROR: Minimum loop size has to be larger than maximum loop size')
        if self.isSeq is False and os.path.isfile(self.fastaFile) is False:
            raise IOError ("ERROR: Input FASTA File %s not found!." % self.fastaFile)
        if self.bulgeSize > 0 and self.stemSize != 3:
            print ("WARNING: Bulge search disabled as stem size not equal to 3")
            self.bulgeSize = 0
        for i in self.strands:
            if i not in ['+', '-']:
                print ("WARNING: Bulge search disabled as stem size not equal to 3")
                
        if self.saveName is not None:
            try:
                self.outHandles['ov'] = open('%s_ov.bed' % self.saveName, 'w')
                self.outHandles['nov'] = open('%s_nov.bed' % self.saveName, 'w')
            except IOError:
                print ('WARNING: Unable to open output files. Will NOT save the output files!')
        return True

    def _read_fasta(self):
        with open(self.fastaFile) as handle:
            if sys.version_info[0] == 3:
                firstline = next(handle).rstrip('\r').rstrip('\n')
            else:
                firstline = handle.next().rstrip('\r').rstrip('\n')
            if firstline[0] != '>':
                raise IOERROR('Input file not in FASTA format')
            head = firstline[1:].split(' ')[0]
            seq = []
            for line in handle:
                if line[0] == '>':
                    yield head, ''.join(seq)
                    head = line[1:].split(' ')[0].rstrip('\r').rstrip('\n')
                    seq = []
                else:
                    seq.append(line.rstrip('\r').rstrip('\n').upper())
            yield head, ''.join(seq)
                
    def _make_base(self, strand):
        if strand == "+":
            base = 'G'
            invert_base = 'C'
        else:
            base = 'C'
            invert_base = 'G'
        return base, invert_base

    def _make_signature(self, base):
        return "".join(map(str, [base, self.stemSize, 'L',
                                 self.minLoopSize, '-', self.maxLoopSize]))

    def _make_stem_motif(self, base, invert_base):
        if self.bulgeSize > 0:
            stem_motif = '(?:%s|%s[AT%sN]{1,%d}%s|%s[AT%sN]{1,%d}%s)' % (base * 3, base, invert_base,
                         self.bulgeSize, base * 2, base * 2, invert_base, self.bulgeSize, base)
        else:
            stem_motif = '%s' % (base * self.stemSize)
        return stem_motif

    def _make_loop_motif(self, base, invert_base):
        if self.greedySearch is True:
            loop_motif = "[ACTGN]{%d,%d}" % (self.minLoopSize, self.maxLoopSize)
        else:
            loop_motif = '[ATGCN](?:[ATN%s]|(?!%s{2})%s){%d,%d}' % (
                invert_base, base, base, self.minLoopSize - 1, self.maxLoopSize - 1)
        return loop_motif

    def _make_motif(self, strand):
        base, invert_base = self._make_base(strand)
        quad_signature = self._make_signature(base)
        stem_motif = self._make_stem_motif(base, invert_base)
        loop_motif = self._make_loop_motif(base, invert_base)
        quad_motif_nov = (stem_motif + loop_motif) * 3 + stem_motif
        quad_motif_nov = "(" + quad_motif_nov + ")"
        quad_motif_ov = "(?=(%s))" % quad_motif_nov
        return quad_signature, quad_motif_ov, quad_motif_nov

    def _run_regex(self, seq, name, pattern, signature):
        regex = re.compile(pattern, flags=re.IGNORECASE)
        reg_obj = regex.finditer(seq)
        res = []
        for match in reg_obj:
            temp_res = [name, match.start(), match.start() + len(match.group(1)) + 1,
                        len(match.group(1)), signature, match.group(1)]
            res.append("\t".join(map(str, temp_res)))
        return res

    def _yield_sequences(self):
        for n,i in enumerate(self.fastaFile):
            yield n, i
     
    def run(self):
        fastaGen = self._read_fasta()
        if self.isSeq is True:
            fastaGen = self._yield_sequences()
        for h, s in fastaGen:
            for i in self.strands:
                if self.verbose is True:
                    print ("INFO: Searching sequence '%s' of length %d on %s strand" % (h, len(s), i))
                signature, motifOv, motifNov = self._make_motif(i)
                res_ov = self._run_regex(s, h, motifOv, signature)
                res_nov = self._run_regex(s, h, motifNov, signature)
                self.resOv.extend(res_ov)
                self.resNov.extend(res_nov)
                if self.verbose is True:
                    print ("INFO: Found %d non-overlapping G4 motifs" % len(res_nov))
        return True

    def flush(self):
        if 'ov' in self.outHandles:
            self.outHandles['ov'].write("\n".join(self.resOv) + "\n")
            self.outHandles['ov'].close()
        if 'nov' in self.outHandles:
            self.outHandles['nov'].write("\n".join(self.resNov) + "\n")
            self.outHandles['nov'].close()
        return True


help_message = """
-------------------------------------------------------------------
Usage: python quadFinder.py [Parameters]

Mandatory parameters:
  -i     Input FASTA file with complete path
  -o     Output file name

Optional parameters:
  -g     Activate greedy algorithm (default=True)
  -s     Length of the stem (default: 3)
  -minl  Minimum length of the loop (default: 1)
  -maxl  Minimum length of the loop (default: 7)
  -maxb  Maximum length of the bulges (default: 0)
  -h     Show this message and exit.
  
Example:
python quadFinder.py -i test.fasta -o results.bed -maxl 15 -maxb2
-------------------------------------------------------------------
"""

if __name__ == "__main__":
    valid_opts = {'-g': True, '-s': 3, '-minl': 1, '-maxl': 7,
                  '-maxb': 0, '-h': '', '-i': '', '-o': ''}
    if len(sys.argv) < 2:
        print (help_message)
        exit()
    opts = {}
    for n,i in enumerate(sys.argv[1:]):
        if n % 2 == 0:
            if i not in valid_opts:
                print ("ERROR: %s not a valid option" % i)
                print (help_message)
                exit()
            else:
                tag = i
        else:
            if i[0] == '-':
                print ('ERROR: Please provide a valid value for %s' % (tag))
            else:
                opts[tag] = i
    if '-i' not in opts:
        print ('ERROR: Please provide name of the FASTA file')
        print (help_message)
        exit()
    elif '-o' not in opts:
        print ('ERROR: Please provide name of the output file')
        print (help_message)
        exit()
    for opt in valid_opts:
        if opt not in opts:
            opts[opt] = valid_opts[opt]
    
    q = QuadMotifFinder(opts['-i'], opts['-o'], opts['-s'], opts['-minl'],
                        opts['-maxl'], opts['-g'], opts['-maxb'])
    q.run()
    q.flush()
    print ("INFO: Program completed successfully")