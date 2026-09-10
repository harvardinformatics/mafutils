#############################################################################
# Common functions to support the utility scripts for the pipeline.
#
# Gregg Thomas, April 2025
#############################################################################

import gzip
import hashlib
import io
import os
import sys

from mafutils.lib import bgzf

#############################################################################

def _isBgzip(filename):
# BGZF blocks are gzip members with a self-identifying 'BC' extra subfield.
# Header layout (see the BAM/SAM spec and bgzf._bgzf_header):
#   bytes 0-3:   1f 8b 08 04  (gzip magic + FEXTRA flag set)
#   bytes 10-11: XLEN (extra field length, little-endian uint16)
#   bytes 12-13: subfield id, must be b"BC" for BGZF

    with open(filename, "rb") as fp:
        header = fp.read(18)

    if len(header) < 14:
        return False
    if header[0:4] != b"\x1f\x8b\x08\x04":
        return False
    extra_len = int.from_bytes(header[10:12], "little")
    if extra_len < 6:
        return False
    return header[12:14] == b"BC"

#############################################################################

def detectCompression(filename):
# Detect compression of a file by examining the first lines in the file

    compression_type = "none";

    magic_dict = {
            b"\x1f\x8b\x08": "gz",
            # b"\x1f\x8b\x08\x08": "gz",
            b"\x42\x5a\x68": "bz2",
            b"\x50\x4b\x03\x04": "zip"
        }
    # An encoded set of possible "magic strings" that start different types of compressed files
    # From: https://www.garykessler.net/library/file_sigs.html
    # \x is the escape code for hex values
    # b converts strings to bytes

    max_len = max(len(x) for x in magic_dict)
    # The number of characters to read from the beginning of the file should be the length of
    # the longest magic string

    file_start = open(filename, "rb").read(max_len);
    # Read the beginning of the file up to the length of the longest magic string

    for magic_string in magic_dict:
        if file_start.startswith(magic_string):
            compression_type = magic_dict[magic_string];
    # Check each magic string against the start of the file

    if compression_type == "gz" and _isBgzip(filename):
        compression_type = "bgzip"
    # Plain gzip and bgzip share the same leading magic bytes; bgzip is only
    # distinguishable by its 'BC' extra subfield.

    return compression_type;

#############################################################################

def openMaf(filename, compression, mode="rt"):
    """
    Open a MAF file with the reader appropriate for its compression type.
    "none" -> builtin open, "gz" -> gzip.open, "bgzip" -> bgzf.open (this
    package's vendored copy of Biopython's Bio.bgzf, see mafutils/lib/bgzf.py)
    -- real random access via virtual offsets, unlike plain gzip.
    """
    if compression == "none":
        return open(filename, mode)
    elif compression == "gz":
        return gzip.open(filename, mode)
    elif compression == "bgzip":
        return bgzf.open(filename, mode)
    else:
        raise ValueError(f"Unsupported MAF compression: {compression}")

#############################################################################

class _HashingRawIO(io.RawIOBase):
    """
    Wraps a raw binary file object; every byte physically read from it is
    fed into a hash object before being returned. Used only at index-build
    time so mafutils can compute a whole-file content hash in the same pass
    already used to find block offsets, for all compression types -- the
    hash reflects the exact on-disk bytes (compressed representation for
    gz/bgzip), since that's what the index's offsets are tied to, and it's
    read before any decompression happens.
    """

    def __init__(self, raw, hash_obj):
        self._raw = raw
        self._hash = hash_obj

    def readable(self):
        return True

    def seekable(self):
        return self._raw.seekable()

    def tell(self):
        return self._raw.tell()

    def seek(self, offset, whence=0):
        # bgzf.BgzfReader issues "seeks" as part of purely sequential reading
        # (e.g. re-asserting the current position when loading each new
        # block) -- those are genuine no-ops here. A real jump would mean
        # bytes get hashed out of the true sequential order, silently
        # corrupting the hash, so only allow seeking to the position we're
        # already at; anything else fails loudly instead of doing that.
        if whence != 0:
            raise io.UnsupportedOperation("_HashingRawIO only supports absolute seeks")
        current = self._raw.tell()
        if offset != current:
            raise io.UnsupportedOperation(
                f"_HashingRawIO only supports a single sequential pass; refusing to "
                f"seek from {current} to {offset}."
            )
        return current

    def readinto(self, b):
        n = self._raw.readinto(b)
        if n:
            self._hash.update(bytes(b[:n]))
        return n

    def close(self):
        try:
            self._raw.close()
        finally:
            super().close()


def openMafHashing(filename, compression, hash_obj):
    """
    Like openMaf(), but threads a _HashingRawIO underneath the reader so
    every byte read from disk is fed into hash_obj -- for gz/bgzip this
    means hashing the raw compressed bytes, before decompression. Used only
    by mafutils index.

    For "none"/"gz" this returns a BINARY stream, to be consumed by
    iterMafBlocks(..., binary=True): those two compression types have plain
    additive stream positions, so byte offsets can be tracked with a simple
    counter instead of calling tell(). That matters enormously at scale --
    on a TextIOWrapper, tell() must snapshot the incremental UTF-8 decoder's
    state to build a seekable cookie, making it ~4.6x more expensive than
    the readline() it accompanies, and iterMafBlocks needs a position before
    every single line. Measured on a 1.5GB real MAF: 11.89s (text+tell) vs
    1.38s (binary+counting) for the same, byte-identical offsets.

    For "bgzip" this still returns a TEXT stream, because BGZF virtual
    offsets pack a compressed-block offset and an in-block offset into one
    int -- they are NOT additive, so a byte counter cannot reproduce them
    and tell() is genuinely required. See iterMafBlocks for the full
    explanation.
    """
    hashing_raw = _HashingRawIO(open(filename, "rb"), hash_obj)

    if compression == "none":
        return io.BufferedReader(hashing_raw)
    elif compression == "gz":
        return io.BufferedReader(gzip.GzipFile(fileobj=hashing_raw, mode="rb"))
    elif compression == "bgzip":
        return bgzf.BgzfReader(fileobj=hashing_raw, mode="rt")
    else:
        raise ValueError(f"Unsupported MAF compression: {compression}")


def computeFileHash(maf_file, algo="md5"):
    """
    Reads maf_file's raw bytes in chunks and returns "<algo>:<hexdigest>".
    Used to check a file against a previously-stored hash (mafutils
    validate, and --verify-hash on fetch/stats/gc) -- unlike
    openMafHashing(), this is a dedicated read purely for hashing, since
    there's no other pass to piggyback on at validation time.
    """
    hash_obj = hashlib.new(algo)
    with open(maf_file, "rb") as fp:
        for chunk in iter(lambda: fp.read(1024 * 1024), b""):
            hash_obj.update(chunk)
    return f"{algo}:{hash_obj.hexdigest()}"

#############################################################################

COPY_CHUNK_SIZE = 8 * 1024 * 1024


def copyMafRangeToStream(handle, compression, offset_start, offset_end, out_stream, chunk_size=COPY_CHUNK_SIZE):
    """
    Copies the byte range [offset_start, offset_end) from an open binary-mode
    MAF handle straight into out_stream, in bounded chunks. Returns the number
    of bytes written.

    Use this instead of readMafBlockBytes() whenever the range can be large --
    notably whole-scaffold extraction. readMafBlockBytes() materializes the
    entire range as one bytes object, which is fine for a single alignment
    block (KBs) but fatal for a scaffold: on a real 7.4TB 241-species MAF,
    chr1 spans bytes 16 -> 597,495,640,164, so that became a single ~597GB
    read() and died with a MemoryError whose str() is the empty string
    (hence a completely blank error message in the log).

    Memory here is O(chunk_size) regardless of range size, for every
    compression type -- including bgzip, whose line-by-line path must also
    write as it goes rather than accumulating lines to join at the end.
    """
    handle.seek(offset_start)
    written = 0

    if compression == "bgzip":
        # Virtual offsets are comparable but not subtractable, so read
        # forward until the handle's position reaches offset_end. Write each
        # line straight out instead of accumulating them.
        while handle.tell() < offset_end:
            line = handle.readline()
            if not line:
                break
            out_stream.write(line)
            written += len(line)
        return written

    remaining = offset_end - offset_start
    while remaining > 0:
        chunk = handle.read(min(chunk_size, remaining))
        if not chunk:
            break
        out_stream.write(chunk)
        written += len(chunk)
        remaining -= len(chunk)
    return written


def readMafBlockBytes(handle, compression, offset_start, offset_end):
    """
    Seeks an open binary-mode MAF handle to offset_start and reads exactly
    the bytes of one block, returning them.

    Only appropriate for genuinely block-sized ranges, since it holds the
    whole range in memory -- see copyMafRangeToStream() for large ranges
    such as whole scaffolds.

    For "none"/"gz", offset_start/offset_end are raw byte positions /
    decompressed-stream positions respectively, so a simple subtraction
    gives the correct read length.

    For "bgzip", offsets are BGZF virtual offsets (compressed block offset
    packed with an in-block offset) -- subtracting two virtual offsets does
    NOT give a meaningful byte count. Virtual offsets are only valid to
    *compare* (for ordering), not to arithmetic on, so instead we read
    forward line-by-line until the handle's position reaches offset_end.
    """
    handle.seek(offset_start)
    if compression == "bgzip":
        chunks = []
        while handle.tell() < offset_end:
            line = handle.readline()
            if not line:
                break
            chunks.append(line)
        return b"".join(chunks)
    else:
        return handle.read(offset_end - offset_start)

#############################################################################

def deriveBlockIndexPath(maf_file):
    return f"{maf_file}.block.idx"


def deriveScaffoldIndexPath(maf_file):
    return f"{maf_file}.scaffold.idx"

#############################################################################

INDEX_HEADER_PREFIX = "# mafutils-index"


def writeIndexHeader(stream, maf_file, compression, size, mtime, content_hash):
    stream.write(
        f"{INDEX_HEADER_PREFIX} format=2 maf={os.path.basename(maf_file)} "
        f"compression={compression} size={size} mtime={mtime} hash={content_hash}\n"
    )


def readIndexHeader(index_file):
    """
    Peeks at the first line of an index file. Returns a dict of the header's
    key=value fields, or None if the index predates this header (an older
    mafutils version built it, or the file is otherwise headerless).
    """
    with open(index_file, "r", encoding="utf-8") as fp:
        first_line = fp.readline()

    if not first_line.startswith(INDEX_HEADER_PREFIX):
        return None

    fields = {}
    for token in first_line.strip().split()[2:]:
        if "=" in token:
            key, _, value = token.partition("=")
            fields[key] = value
    return fields or None


MTIME_TOLERANCE_SECONDS = 2.0


def compareIndexHeader(header, maf_file, detected_compression, check_hash=False):
    """
    Compares a parsed index header (from readIndexHeader) against the actual
    file at maf_file. Does not log or raise -- returns a dict with one entry
    per checked field (compression/size/mtime/hash), each either None (the
    header didn't record that field, e.g. a format=1 index predating it) or
    {"match": bool, "expected": ..., "actual": ...}.

    check_hash controls whether the (expensive -- reads the whole file)
    hash comparison is performed at all; callers that only want the cheap
    size/mtime check should leave it False.
    """
    result = {"compression": None, "size": None, "mtime": None, "hash": None}
    if header is None:
        return result

    header_compression = header.get("compression")
    if header_compression is not None:
        result["compression"] = {
            "match": header_compression == detected_compression,
            "expected": header_compression,
            "actual": detected_compression,
        }

    header_size = header.get("size")
    if header_size is not None:
        actual_size = os.path.getsize(maf_file)
        result["size"] = {
            "match": int(header_size) == actual_size,
            "expected": header_size,
            "actual": str(actual_size),
        }

    header_mtime = header.get("mtime")
    if header_mtime is not None:
        actual_mtime = os.path.getmtime(maf_file)
        result["mtime"] = {
            "match": abs(float(header_mtime) - actual_mtime) <= MTIME_TOLERANCE_SECONDS,
            "expected": header_mtime,
            "actual": str(actual_mtime),
        }

    header_hash = header.get("hash")
    if header_hash is not None and check_hash:
        algo = header_hash.split(":", 1)[0]
        actual_hash = computeFileHash(maf_file, algo=algo)
        result["hash"] = {
            "match": header_hash == actual_hash,
            "expected": header_hash,
            "actual": actual_hash,
        }

    return result


def validateIndexHeader(header, maf_file, detected_compression, LOG, strict=False):
    """
    Validates a parsed index header (from readIndexHeader) against the
    actual MAF file being processed now, and fails the run (via sys.exit)
    on any conclusive problem.

    strict=False (default): cheap check -- compression mismatch is an
    error, size mismatch is an error (files differing in size are
    definitely not the same content the index was built from), mtime
    mismatch is only a warning (weak signal -- files get touched/copied
    without content changing). A missing header only warns.

    strict=True (--verify-hash): compares against the stored content hash
    instead -- authoritative, but requires reading the whole file. A
    missing header, or a header with no stored hash, is an error in this
    mode: an explicit request for certainty can't be silently downgraded.
    """
    if header is None:
        if strict:
            LOG.error(
                f"--verify-hash requested but index for {maf_file} has no mafutils "
                f"header to verify against (built by an older mafutils version)."
            )
            sys.exit(1)
        LOG.warning(
            f"Index for {maf_file} has no mafutils header (built by an older "
            f"mafutils version); cannot verify it matches this file's "
            f"compression ({detected_compression})."
        )
        return

    result = compareIndexHeader(header, maf_file, detected_compression, check_hash=strict)

    compression_check = result["compression"]
    if compression_check is not None and not compression_check["match"]:
        LOG.error(
            f"Index/MAF compression mismatch for {maf_file}: index was built "
            f"with compression={compression_check['expected']}, but this file was "
            f"detected as compression={compression_check['actual']}. Rebuild the "
            f"index with `mafutils index` against this exact file."
        )
        sys.exit(1)

    if strict:
        hash_check = result["hash"]
        if hash_check is None:
            LOG.error(
                f"--verify-hash requested but index for {maf_file} has no stored "
                f"hash (built by an older mafutils version). Rebuild the index to "
                f"enable hash verification."
            )
            sys.exit(1)
        if not hash_check["match"]:
            LOG.error(
                f"Index/MAF content mismatch for {maf_file}: index was built from "
                f"a file with hash {hash_check['expected']}, but this file currently "
                f"hashes to {hash_check['actual']}. Rebuild the index with "
                f"`mafutils index` against this exact file."
            )
            sys.exit(1)
        LOG.info(f"Verified {maf_file} content hash matches index.")
        return

    size_check = result["size"]
    if size_check is not None and not size_check["match"]:
        LOG.error(
            f"Index/MAF size mismatch for {maf_file}: index was built from a "
            f"{size_check['expected']}-byte file, but this file is currently "
            f"{size_check['actual']} bytes. Rebuild the index with `mafutils index` "
            f"against this exact file."
        )
        sys.exit(1)

    mtime_check = result["mtime"]
    if mtime_check is not None and not mtime_check["match"]:
        LOG.warning(
            f"Index/MAF modification time differs for {maf_file} (index built at "
            f"mtime={mtime_check['expected']}, file currently has "
            f"mtime={mtime_check['actual']}). This alone doesn't necessarily mean "
            f"the content changed -- rerun with --verify-hash for a definitive check "
            f"if unsure."
        )

#############################################################################

def iterMafBlocks(stream, binary=False):
    """
    Splits an open MAF stream into raw per-block lines, split on 'a'-prefixed
    header lines (skipping blank/comment lines). Yields
    (block_lines, block_start_offset, block_end_offset) tuples, where the
    offsets bracket each block and are suitable for later random-access
    seeking back into the same file.

    block_lines[0] (the 'a' header) and block_lines[1] (the reference/first
    sequence line) are always returned as stripped `str` -- they are the
    only two lines whose *content* mafutils index reads. The remaining
    sequence lines are returned exactly as readline() produced them, and in
    binary mode are still `bytes`: only their *count* (len(block_lines)) is
    ever used downstream, so decoding/stripping them would be pure overhead.
    That matters enormously at scale -- a block can carry hundreds of
    sequence lines (one per species) and a whole-genome MAF hundreds of
    millions of blocks.

    Two position-tracking modes, because the two cases are fundamentally
    different:

    binary=True (for "none"/"gz", via openMafHashing): the stream is binary
    and positions are tracked with a plain additive byte counter. Valid
    because both have simple additive stream positions (raw file bytes, and
    decompressed-stream bytes, respectively). This exists for speed: on a
    TextIOWrapper, tell() must snapshot the incremental UTF-8 decoder state
    to build a seekable cookie, costing ~4.6x more than the readline() it
    accompanies -- and a position is needed before every single line.
    Measured on a 1.5GB real MAF: 11.89s (text+tell) vs 1.38s
    (binary+counting), producing byte-identical offsets over 269,308 blocks.

    binary=False (for "bgzip"): positions come from stream.tell() and are
    NEVER derived by arithmetic (e.g. tell() - len(line)). BGZF virtual
    offsets pack a compressed-block offset and an in-block offset into one
    int, so they are not additive -- subtracting or adding a plain byte
    count can "borrow" across the packing whenever an 'a' line straddles a
    real BGZF block boundary, landing on a byte position that was never a
    valid block start (confirmed against a real ~114k-block file, where it
    corrupted the offset by exactly the size of that borrow). This is also
    why the binary counting path above must not be used for bgzip.
    """
    if binary:
        empty, hash_prefix, a_prefix = b"", b"#", b"a"
    else:
        empty, hash_prefix, a_prefix = "", "#", "a"

    block_lines = []
    block_start = None
    pos = 0 if binary else stream.tell()
    next_line_is_reference = False

    while True:
        line = stream.readline()
        if line == empty:
            break
        next_pos = pos + len(line) if binary else stream.tell()

        if line.startswith(hash_prefix) or not line.strip():
            pos = next_pos
            continue

        if line.startswith(a_prefix):
            if block_lines:
                yield block_lines, block_start, pos
            block_start = pos
            block_lines = [line.decode().strip() if binary else line.strip()]
            next_line_is_reference = True
        elif next_line_is_reference:
            block_lines.append(line.decode().strip() if binary else line.strip())
            next_line_is_reference = False
        else:
            block_lines.append(line)

        pos = next_pos

    if block_lines:
        yield block_lines, block_start, pos if binary else stream.tell()

#############################################################################
