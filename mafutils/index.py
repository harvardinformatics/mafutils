#############################################################################
# Create block and scaffold indexes for a MAF file.
#
# Gregg Thomas, December 2023
# Refactored into mafutils package form, April 2026
#############################################################################

import hashlib
import os
import shutil
import tempfile
from typing import Annotated, Optional

import typer

from mafutils.lib import common as COMMON


def process_maf_block(block):
    """
    Returns (row, ref_len, aln_cols, num_seqs) for one block: the index row as
    strings, plus the three integers the header aggregates need.

    The integers are returned rather than re-parsed from `row` on purpose --
    this loop runs once per block (233M times on a whole-genome MAF) and is the
    hot path the v0.6.0 indexing speedup optimised, so int()-ing its own output
    back again would be a needless per-block cost.
    """
    ref_seq = block[1].split()
    ref_scaff = ref_seq[1].split(".", 1)[1]
    line_len = str(len(block[1]))
    num_seqs = len(block) - 1          # every s line, including the reference
    seq_len = len(ref_seq[6])          # block width, gap columns included
    row = [ref_scaff, ref_seq[2], ref_seq[3], str(seq_len), line_len, str(num_seqs)]
    return row, int(ref_seq[3]), seq_len, num_seqs


def run_index(maf_file, block_index_path, scaffold_index_path):
    maf_compression = COMMON.detectCompression(maf_file)
    size = os.path.getsize(maf_file)
    mtime = os.path.getmtime(maf_file)
    hash_obj = hashlib.md5()

    out_dir = os.path.dirname(block_index_path) or "."
    with tempfile.TemporaryDirectory(prefix="maf_index_tmp_", dir=out_dir) as tmp_dir:
        block_tmp_path = os.path.join(tmp_dir, "block.idx")
        scaffold_tmp_path = os.path.join(tmp_dir, "scaffold.idx")

        # Rows are written to temp files first because the header (line 1 of
        # each real output file) needs size/mtime/hash, and the hash isn't
        # final until the whole file has been read via openMafHashing below.
        with COMMON.openMafHashing(maf_file, maf_compression, hash_obj) as maf_stream, \
                open(block_tmp_path, "w") as block_stream, \
                open(scaffold_tmp_path, "w") as scaffold_stream:

            current_scaffold = None
            region_start_byte = None
            region_end_byte = None

            # none/gz get the fast binary+byte-counting path; bgzip must use
            # text mode + tell(), since its virtual offsets aren't additive
            # (see iterMafBlocks/openMafHashing).
            use_binary = maf_compression in ("none", "gz")

            # Header aggregates, accumulated as we go so `mafutils info` can
            # report them later from the header alone rather than rescanning.
            n_blocks = 0
            total_ref_bases = 0
            total_aln_cols = 0
            total_seq_lines = 0
            max_seqs = 0
            distinct_scaffolds = set()

            for block, block_start, block_end in COMMON.iterMafBlocks(maf_stream, binary=use_binary):
                block_info, ref_len, aln_cols, num_seqs = process_maf_block(block)
                ref_scaffold = block_info[0]

                n_blocks += 1
                total_ref_bases += ref_len
                total_aln_cols += aln_cols
                total_seq_lines += num_seqs
                if num_seqs > max_seqs:
                    max_seqs = num_seqs
                distinct_scaffolds.add(ref_scaffold)

                mdx_line = block_info + [str(block_start), str(block_end)]
                block_stream.write("\t".join(mdx_line) + "\n")

                if current_scaffold is None:
                    current_scaffold = ref_scaffold
                    region_start_byte = block_start
                    region_end_byte = block_end
                elif ref_scaffold != current_scaffold:
                    scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")
                    current_scaffold = ref_scaffold
                    region_start_byte = block_start
                    region_end_byte = block_end
                else:
                    region_end_byte = block_end

            if current_scaffold is not None:
                scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")

        content_hash = f"md5:{hash_obj.hexdigest()}"

        # Identical in both headers on purpose -- validate cross-checks the two
        # for exact equality, and these describe the MAF, not the index file.
        aggregates = {
            "blocks": n_blocks,
            "scaffolds": len(distinct_scaffolds),
            "ref_bases": total_ref_bases,
            "aln_cols": total_aln_cols,
            "seq_lines": total_seq_lines,
            "max_seqs": max_seqs,
        }

        with open(block_index_path, "w") as block_stream:
            COMMON.writeIndexHeader(block_stream, maf_file, maf_compression, size, mtime, content_hash, aggregates)
            with open(block_tmp_path, "r") as tmp_fp:
                shutil.copyfileobj(tmp_fp, block_stream)

        with open(scaffold_index_path, "w") as scaffold_stream:
            COMMON.writeIndexHeader(scaffold_stream, maf_file, maf_compression, size, mtime, content_hash, aggregates)
            with open(scaffold_tmp_path, "r") as tmp_fp:
                shutil.copyfileobj(tmp_fp, scaffold_stream)


def index_command(
    maf_file: Annotated[str, typer.Argument(help="Input MAF file (.maf, .maf.gz, or bgzip-compressed .maf)")],
    block_index: Annotated[Optional[str], typer.Argument(help="Output block index path (default: <MAF_FILE>.block.idx)")] = None,
    scaffold_index: Annotated[Optional[str], typer.Argument(help="Output scaffold index path (default: <MAF_FILE>.scaffold.idx)")] = None,
) -> None:
    if block_index is None:
        block_index = COMMON.deriveBlockIndexPath(maf_file)
    if scaffold_index is None:
        scaffold_index = COMMON.deriveScaffoldIndexPath(maf_file)
    run_index(maf_file, block_index, scaffold_index)
