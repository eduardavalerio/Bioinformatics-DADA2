# Cutadapt ✂️🧬

[Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1), 10-12.](https://journal.embnet.org/index.php/embnetjournal/article/view/200)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing)

```
cutadapt version 5.0

Copyright (C) 2010 Marcel Martin <marcel.martin@scilifelab.se> and contributors

Cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage:
    cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq

For paired-end reads:
    cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

Replace "ADAPTER" with the actual sequence of your 3' adapter. IUPAC wildcard
characters are supported. All reads from input.fastq will be written to
output.fastq with the adapter sequence removed. Adapter matching is
error-tolerant. Multiple adapter sequences can be given (use further -a
options), but only the best-matching adapter will be removed.

Input may also be in FASTA format. Compressed input and output is supported and
auto-detected from the file name (.gz, .xz, .bz2). Use the file name '-' for
standard input/output. Without the -o option, output is sent to standard output.

Citation:

Marcel Martin. Cutadapt removes adapter sequences from high-throughput
sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011.
http://dx.doi.org/10.14806/ej.17.1.200

Run "cutadapt --help" to see all command-line options.
See https://cutadapt.readthedocs.io/ for full documentation.

Options:
  -h, --help            Show this help message and exit
  --version             Show version number and exit
  --debug               Print debug log. Use twice to also print DP matrices
  -j CORES, --cores CORES
                        Number of CPU cores to use. Use 0 to auto-detect.
                        Default: 1

Finding adapters:
  Parameters -a, -g, -b specify adapters to be removed from each read (or from
  R1 if data is paired-end. If specified multiple times, only the best
  matching adapter is trimmed (but see the --times option). Use notation
  'file:FILE' to read adapter sequences from a FASTA file.

  -a ADAPTER, --adapter ADAPTER
                        Sequence of an adapter ligated to the 3' end (paired
                        data: of the first read). The adapter and subsequent
                        bases are trimmed. If a '$' character is appended
                        ('anchoring'), the adapter is only found if it is a
                        suffix of the read.
  -g ADAPTER, --front ADAPTER
                        Sequence of an adapter ligated to the 5' end (paired
                        data: of the first read). The adapter and any preceding
                        bases are trimmed. Partial matches at the 5' end are
                        allowed. If a '^' character is prepended ('anchoring'),
                        the adapter is only found if it is a prefix of the read.
  -b ADAPTER, --anywhere ADAPTER
                        Sequence of an adapter that may be ligated to the 5' or
                        3' end (paired data: of the first read). Both types of
                        matches as described under -a and -g are allowed. If the
                        first base of the read is part of the match, the
                        behavior is as with -g, otherwise as with -a. This
                        option is mostly for rescuing failed library
                        preparations - do not use if you know which end your
                        adapter was ligated to!
  -e E, --error-rate E, --errors E
                        Maximum allowed error rate (if 0 <= E < 1), or absolute
                        number of errors for full-length adapter match (if E is
                        an integer >= 1). Error rate = no. of errors divided by
                        length of matching region. Default: 0.1 (10%)
  --no-indels           Allow only mismatches in alignments. Default: allow both
                        mismatches and indels
  -n COUNT, --times COUNT
                        Remove up to COUNT adapters from each read. Default: 1
  -O MINLENGTH, --overlap MINLENGTH
                        Require MINLENGTH overlap between read and adapter for
                        an adapter to be found. Default: 3
  --match-read-wildcards
                        Interpret IUPAC wildcards in reads. Default: False
  -N, --no-match-adapter-wildcards
                        Do not interpret IUPAC wildcards in adapters.
  --action {trim,retain,mask,lowercase,crop,none}
                        What to do if a match was found. trim: trim adapter and
                        up- or downstream sequence; retain: trim, but retain
                        adapter; mask: replace with 'N' characters; lowercase:
                        convert to lowercase; crop: trim up and downstream
                        sequence; none: leave unchanged. Default: trim
  --rc, --revcomp       Check both the read and its reverse complement for
                        adapter matches. If match is on reverse-complemented
                        version, output that one. Default: check only read

Additional read modifications:
  -u LEN, --cut LEN     Remove LEN bases from each read (or R1 if paired; use -U
                        option for R2). If LEN is positive, remove bases from
                        the beginning. If LEN is negative, remove bases from the
                        end. Can be used twice if LENs have different signs.
                        Applied *before* adapter trimming.
  --nextseq-trim 3'CUTOFF
                        NextSeq-specific quality trimming (each read). Trims
                        also dark cycles appearing as high-quality G bases.
  -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF
                        Trim low-quality bases from 5' and/or 3' ends of each
                        read before adapter removal. Applied to both reads if
                        data is paired. If one value is given, only the 3' end
                        is trimmed. If two comma-separated cutoffs are given,
                        the 5' end is trimmed with the first cutoff, the 3' end
                        with the second.
  --quality-base N      Assume that quality values in FASTQ are encoded as
                        ascii(quality + N). This needs to be set to 64 for some
                        old Illumina FASTQ files. Default: 33
  --poly-a              Trim poly-A tails
  --length LENGTH, -l LENGTH
                        Shorten reads to LENGTH. Positive values remove bases at
                        the end while negative ones remove bases at the
                        beginning. This and the following modifications are
                        applied after adapter trimming.
  --trim-n              Trim N's on ends of reads.
  --length-tag TAG      Search for TAG followed by a decimal number in the
                        description field of the read. Replace the decimal
                        number with the correct length of the trimmed read. For
                        example, use --length-tag 'length=' to correct fields
                        like 'length=123'.
  --strip-suffix STRIP_SUFFIX
                        Remove this suffix from read names if present. Can be
                        given multiple times.
  -x PREFIX, --prefix PREFIX
                        Add this prefix to read names. Use {name} to insert the
                        name of the matching adapter.
  -y SUFFIX, --suffix SUFFIX
                        Add this suffix to read names; can also include {name}
  --rename TEMPLATE     Rename reads using TEMPLATE containing variables such as
                        {id}, {adapter_name} etc. (see documentation)
  --zero-cap, -z        Change negative quality values to zero.

Filtering of processed reads:
  Filters are applied after above read modifications. Paired-end reads are
  always discarded pairwise (see also --pair-filter).

  -m LEN[:LEN2], --minimum-length LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0
  -M LEN[:LEN2], --maximum-length LEN[:LEN2]
                        Discard reads longer than LEN. Default: no limit
  --max-n COUNT         Discard reads with more than COUNT 'N' bases. If COUNT
                        is a number between 0 and 1, it is interpreted as a
                        fraction of the read length.
  --max-expected-errors ERRORS, --max-ee ERRORS
                        Discard reads whose expected number of errors (computed
                        from quality values) exceeds ERRORS.
  --max-average-error-rate ERROR_RATE, --max-aer ERROR_RATE
                        as --max-expected-errors (see above), but divided by
                        length to account for reads of varying length.
  --discard-trimmed, --discard
                        Discard reads that contain an adapter. Use also -O to
                        avoid discarding too many randomly matching reads.
  --discard-untrimmed, --trimmed-only
                        Discard reads that do not contain an adapter.
  --discard-casava      Discard reads that did not pass CASAVA filtering (header
                        has :Y:).

Output:
  --quiet               Print only error messages.
  --report {full,minimal}
                        Which type of report to print: 'full' or 'minimal'.
                        Default: full
  --json FILE           Dump report in JSON format to FILE
  -o FILE, --output FILE
                        Write trimmed reads to FILE. FASTQ or FASTA format is
                        chosen depending on input. Summary report is sent to
                        standard output. Use '{name}' for demultiplexing (see
                        docs). Default: write to standard output
  --fasta               Output FASTA to standard output even on FASTQ input.
  --compression-level N
                        Compression level for compressed output files. Default:
                        1
  -Z                    DEPRECATED because compression level 1 is now the
                        default.
  --info-file FILE      Write information about each read and its adapter
                        matches into FILE. See the documentation for the file
                        format.
  -r FILE, --rest-file FILE
                        When the adapter matches in the middle of a read, write
                        the rest (after the adapter) to FILE.
  --wildcard-file FILE  When the adapter has N wildcard bases, write adapter
                        bases matching wildcard positions to FILE. (Inaccurate
                        with indels.)
  --too-short-output FILE
                        Write reads that are too short (according to length
                        specified by -m) to FILE. Default: discard reads
  --too-long-output FILE
                        Write reads that are too long (according to length
                        specified by -M) to FILE. Default: discard reads
  --untrimmed-output FILE
                        Write reads that do not contain any adapter to FILE.
                        Default: output to same file as trimmed reads

Paired-end options:
  The -A/-G/-B/-U/-Q options work like their lowercase counterparts, but are
  applied to R2 (second read in pair)

  -A ADAPTER            3' adapter to be removed from R2
  -G ADAPTER            5' adapter to be removed from R2
  -B ADAPTER            5'/3 adapter to be removed from R2
  -U LENGTH             Remove LENGTH bases from R2
  -Q [5'CUTOFF,]3'CUTOFF
                        Quality-trimming cutoff for R2. Default: same as for R1
  -L LENGTH             Shorten R2 to LENGTH. Default: same as for R1
  -p FILE, --paired-output FILE
                        Write R2 to FILE.
  --pair-adapters       Treat adapters given with -a/-A etc. as pairs. Either
                        both or none are removed from each read pair.
  --pair-filter {any,both,first}
                        Which of the reads in a paired-end read have to match
                        the filtering criterion in order for the pair to be
                        filtered. Default: any
  --interleaved         Read and/or write interleaved paired-end reads.
  --untrimmed-paired-output FILE
                        Write second read in a pair to this FILE when no adapter
                        was found. Use with --untrimmed-output. Default: output
                        to same file as trimmed reads
  --too-short-paired-output FILE
                        Write second read in a pair to this file if pair is too
                        short.
  --too-long-paired-output FILE
                        Write second read in a pair to this file if pair is too
                        long.
```
