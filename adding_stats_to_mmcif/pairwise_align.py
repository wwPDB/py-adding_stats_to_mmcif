#!/usr/bin/env python
import argparse
import logging

from Bio import Align
from Bio import pairwise2
# AL: hiding warnings, could be done better
import warnings
with warnings.catch_warnings():
  warnings.simplefilter("ignore")
  from Bio.SubsMat import MatrixInfo as matlist

# from Bio import AlignIO
# from Bio.Emboss.Applications import NeedleCommandline
# import StringIO

logger = logging.getLogger()
FORMAT = "%(filename)s - %(funcName)s - %(message)s"
logging.basicConfig(format=FORMAT)


# needle_cline = NeedleCommandline()
# needle_cline.gapopen = 0.5
# needle_cline.gapextend = 0.1
# needle_cline.gapopen = 10
# needle_cline.gapextend = 0.5
# needle_cline.stdout = True
# needle_cline.outfile = "stdout"

class SequenceAlign:

    def __init__(self, sequence1, sequence2):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.score = None
        self.alignment_dict = dict()
        self.alignment_list = list()
        logging.debug(self.sequence1)
        logging.debug(self.sequence2)

    def rna(self, seq):
        return set(seq).issubset(set("AUGC"))

    def dna(self, seq):
        return set(seq).issubset(set("ATGC"))

    # def dnarna(self, seq):
    #    return self.rna(seq=seq) or self.dna(seq=seq)

    def both_sequences_same_type(self):
        if self.dna(self.sequence1) == self.dna(self.sequence2) and self.rna(self.sequence1) == self.rna(
                self.sequence2):
            return True, ''
        return False, 'sequences not the same type'
        # if self.dnarna(self.sequence1) != self.dnarna(self.sequence2):
        #    return False, 'sequences not the same type'
        # return True, ''

    def remove_gaps(self, sequence):
        return str(sequence).replace("\n", "").replace(" ", "")

    def prepare_sequences(self):
        self.sequence1 = self.remove_gaps(self.sequence1)
        self.sequence2 = self.remove_gaps(self.sequence2)
        # if len(self.sequence1) > 2000 or len(self.sequence2) > 2000:
        #    return False, 'sequences too long. Please install emboss needle'
        return self.both_sequences_same_type()

    def pairwise2(self):

        matrix = matlist.blosum62
        gap_open = -10
        gap_extend = -0.5
        alns = pairwise2.align.globalds(self.sequence1, self.sequence2, matrix, gap_open, gap_extend)[0]
        logging.info(pairwise2.format_alignment(*alns))
        logging.info(alns)
        self.score = alns[2]

    def pairwise_aligner(self):

        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        # aligner.match = 2
        # aligner.mismatch = -1
        # only need to run aligner.score. This improves memory usage and speed.
        # logging.debug('length of query: {}'.format(len(self.sequence1)))
        # logging.debug('length of target: {}'.format(len(self.sequence2)))
        alignments = aligner.align(self.sequence1, self.sequence2)
        # for alignment in sorted(alignments):
        #     logging.debug(alignment)
        #     # logging.debug(alignment.score)
        #     logging.debug(alignment.target)
        #     logging.debug(alignment.query)
        #     logging.debug(alignment.path)
        #     # logging.debug(dir(alignment))
        #     current_query_position = 0
        #     current_target_position = 0
        #     for align_tupple in alignment.path:
        #         working_query_position = align_tupple[0]
        #         working_target_position = align_tupple[1]
        #         query_shift = working_query_position - current_query_position
        #         target_shift = working_target_position - current_target_position
        #         logging.debug(query_shift)
        #         # expand the shift into a list. Then for each position add the position in the list to a dictionary?
        #
        #         current_target_position = working_target_position
        #         current_query_position = working_query_position

        align_score = aligner.score(self.sequence1, self.sequence2)
        logging.info(align_score)

        self.score = align_score

    """
    def do_alignment_emboss(self):
        needle_cline.asequence = "asis:" + self.sequence1
        needle_cline.bsequence = "asis:" + self.sequence2
        stdout, stderr = needle_cline()
        result = [AlignIO.read(StringIO.StringIO(stdout), "emboss")]
        self.score = self.get_emboss_score(result[0].seq, result[1].seq)

    def get_emboss_score(self, seq1, seq2):
        return sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2))
    """

    def get_alignment_score(self):
        return self.score

    def do_sequences_align(self):
        if self.score > 0:
            return True
        return False

    def do_sequence_alignment(self):
        sequences_ok, error = self.prepare_sequences()
        if not sequences_ok:
            return False, error, 0
        self.pairwise_aligner()
        if self.do_sequences_align():
            return True, '', self.score
        return False, 'sequences do not align', 0


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', help='debugging', action='store_const', dest='loglevel', const=logging.DEBUG,
                        default=logging.INFO)
    args = parser.parse_args()
    logger.setLevel(args.loglevel)

    test_sequences = [
        'MEKLEVGIYTRAREGEIACGDACLVKRVEGVIFLAVGDGIGHGPEAARAAEIAIASMESSMNTGLVNIFQLCHRELRGTRGAVAALCRVDRRQGLWQAAIVGNIHVKILSAKGIITPLATPGILGYNYPHQLLIAKGSYQEGDLFLIHSDGIQEGAVPLALLANYRLTAEELVRLIGEKYGRRDDDVAVIVAR',
        'TRAREGEIACGDACLVKRVEGVIFLAVGDGIGHGPEAARAAEIAIASMESSMNTGLVNIFQLCHRELRGTRGAVAALCRVDRRQGLWQAAIVGNIHVKILSAKGIITPLATPGILGYNYPHQLLIAKGSYQEGDLFLIHSDGIQEGAVPLALLANYRLTAEELVRLIGE',
        'DMEGYFVDE', 'RANDOM']

    for sequence in test_sequences:
        sa = SequenceAlign(sequence1=test_sequences[0], sequence2=sequence)
        aligned, error, score = sa.do_sequence_alignment()
        logging.info('is aligned: {}'.format(aligned))
        logging.info('error: "{}"'.format(error))
        logging.info('score: {}'.format(score))
