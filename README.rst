sequniq
-------

Identifying duplicated sequences in a Fasta file
================================================

    >>> from sequniq import seqids
    >>> for dupseqs in seqids.dupids(seqfile_fp):
    >>>     # Process / store sequence IDs
    >>>     print dupseqs

Removing duplicates from an interleaved Fastq file
==================================================

    >>> from sequniq import records
    >>> for record in records.uniqseqs(seqfile_fp):
    >>>     name1, seq1, qual1, name2, seq2, qual2 = record
