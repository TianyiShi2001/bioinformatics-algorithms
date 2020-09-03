from skbio.alignment._pairwise import global_pairwise_align_nucleotide
from skbio.sequence import DNA

global_pairwise_align_nucleotide(DNA("GCAAAAGCTGGTATTAAAGT"),DNA("GCATATTACGTGGTGATTCAAGAGGCCTTCG"),5,1,5,-2,penalize_terminal_gaps=True)

from skbio import __version__ as v

print(v)
