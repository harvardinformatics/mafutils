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
    ref_seq = block[1].split()
    ref_scaff = ref_seq[1].split(".", 1)[1]
    ref_len = block[1].split()[3]
    line_len = str(len(block[1]))
    num_seqs = str(len(block) - 1)
    seq_len = str(len(ref_seq[6]))
    return [ref_scaff, ref_seq[2], ref_len, seq_len, line_len, num_seqs]


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

            for block_text, block_start, block_end in COMMON.iterMafBlocks(maf_stream):
                block = block_text.split("\n")
                block_info = process_maf_block(block)
                ref_scaffold = block_info[0]

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

        with open(block_index_path, "w") as block_stream:
            COMMON.writeIndexHeader(block_stream, maf_file, maf_compression, size, mtime, content_hash)
            with open(block_tmp_path, "r") as tmp_fp:
                shutil.copyfileobj(tmp_fp, block_stream)

        with open(scaffold_index_path, "w") as scaffold_stream:
            COMMON.writeIndexHeader(scaffold_stream, maf_file, maf_compression, size, mtime, content_hash)
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
