#############################################################################
# Create block and scaffold indexes for a MAF file.
#
# Gregg Thomas, December 2023
# Refactored into mafutils package form, April 2026
#############################################################################

import gzip
from typing import Annotated

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

    if maf_compression == "gz":
        maf_stream = gzip.open(maf_file, "rt")
    elif maf_compression == "none":
        maf_stream = open(maf_file, "r")
    else:
        raise ValueError(f"Unsupported compression type: {maf_compression}")

    with maf_stream, open(block_index_path, "w") as block_stream, open(scaffold_index_path, "w") as scaffold_stream:
        line = "init"
        block = []

        current_scaffold = None
        region_start_byte = None
        region_end_byte = None

        while line != "":
            line = maf_stream.readline()
            if line.startswith("#") or line.strip() == "":
                continue

            if line.startswith("a"):
                header_line_len = len(line)

                if block:
                    block_end = maf_stream.tell() - header_line_len
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

                block_start = maf_stream.tell() - header_line_len
                block = [line.strip()]
            else:
                block.append(line.strip())

        if block:
            block_end = maf_stream.tell()
            block_info = process_maf_block(block)
            ref_scaffold = block_info[0]

            mdx_line = block_info + [str(block_start), str(block_end)]
            block_stream.write("\t".join(mdx_line) + "\n")

            if current_scaffold != ref_scaffold:
                if current_scaffold is not None:
                    scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")
                current_scaffold = ref_scaffold
                region_start_byte = block_start
                region_end_byte = block_end
            else:
                region_end_byte = block_end

            scaffold_stream.write(f"{current_scaffold}\t{region_start_byte}\t{region_end_byte}\n")


def index_command(
    maf_file: Annotated[str, typer.Argument(help="Input MAF file (.maf or .maf.gz)")],
    block_index: Annotated[str, typer.Argument(help="Output block index path")],
    scaffold_index: Annotated[str, typer.Argument(help="Output scaffold index path")],
) -> None:
    run_index(maf_file, block_index, scaffold_index)
