# -----------------------------------------------------------------------------
# Copyright (C) Daniel Standage, 2015. It is licensed under the ISC license,
# see LICENSE.txt. Contact: daniel.standage@gmail.com
# -----------------------------------------------------------------------------

"""
This is sequniq; see http://github.com/standage/sequniq.
"""

import parse
_seqcomp = dict(zip("ACGTUBDHKMNRSVWYacgtubdhkmnrsvwy",
                    "TGCAAVHDMKNYSBWRtgcaavhdmknysbwr"))


def revcomp(dna):
    """
    Computes the reverse complement of a nucleotide sequence.
    """
    return "".join(_seqcomp[base] for base in reversed(dna))
