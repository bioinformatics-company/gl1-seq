#!/usr/bin/env python3
"""
seqnft_pipeline.py — SeqNFT + GL1 dual-mode pipeline (on-chain/off-chain identical semantics)
===========================================================================================

UPDATED FOR GL1ENCv2 (gl1enc_v2.py)
-----------------------------------
This version adapts the pipeline to use **gl1enc_v2.py** as the single codec engine
for compact 2/4/6-bit encodings + ASCII fallback:

  - DNA2   (2-bit)  : A/C/G/T only
  - DNA4   (4-bit)  : IUPAC nucleotides + '-' (gap)
  - SIXBIT (6-bit)  : 64-symbol alphabet for proteins / general symbols
  - ASCII  (8-bit)  : fallback

The pipeline is a **local simulation** of the on-chain SeqNFT system described in
`SEQNFT_PIPELINE_SPEC_v2.md`.

It provides:
- Content-addressed object store (git-like)
- Commit DAG (history, branching-ready)
- SeqNFT hierarchy: Assembly (root) -> Chromosome/Record (parent) -> Chunk (child)
- Two deployment modes:
    * fully-onchain: store bytes in the "chain simulator"
    * anchor/ipfs: store only CID-like pointers (simulated) + sha256
- Export/import to file sidecars:
    * .gl1  (GL1F v1)  main payload container (chunked)
    * .gl1x (GL1X v2)  random access + masks/exceptions + chunk hashes + Merkle root
    * .gl1q (GL1Q v2)  FASTQ qualities (optional) with TOC
    * .gl1i (GL1I v1)  minimizer index (optional) with TOC/buckets
- Region queries ("view")
- Evolution-style mutation commits ("mutate")

Notes
-----
- This is designed to "roughly mimic on-chain behaviour locally" with deterministic hashing.
- Hashing uses SHA-256 to match the GL1 family utilities in this reference pipeline.
- If you later port to Solidity, you can switch the hash function to keccak256, but you MUST
  keep canonical encodings stable and commit to one hash family on-chain.

License
-------
MIT-style permission is implied for reference code; adjust as needed for your project.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as _dt
import gzip
import hashlib
import io
import json
import os
import random
import shutil
import struct
import tempfile
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union


# --- Optional dependency: gl1enc_v2.py in same directory (or installed) ---
# You can also rename gl1enc_v2.py -> gl1enc.py and it will be found.
try:
    import gl1enc_v2 as gl1enc  # type: ignore
except Exception:
    try:
        import gl1enc as gl1enc  # type: ignore
    except Exception as e:  # pragma: no cover
        gl1enc = None
        _GL1ENC_IMPORT_ERROR = e


# =============================================================================
# GL1F/GL1X/GL1Q/GL1I constants (pipeline-side reference)
# =============================================================================

# --- Kinds (record semantic type) ---
KIND_NUC = 1
KIND_PROT = 2

# --- Codecs (payload encoding). These match gl1enc.SeqEncoding numeric IDs ---
#   0 DNA2
#   1 ASCII
#   2 DNA4
#   3 SIXBIT
CODEC_DNA2 = 0
CODEC_ASCII = 1
CODEC_DNA4 = 2
CODEC_SIXBIT = 3

# --- Record flags (bitfield) ---
FLAG_HAS_MASK = 1 << 0   # record has lowercase mask intervals
FLAG_HAS_EXCP = 1 << 1   # record has exceptions (used only with DNA2 payloads in this pipeline)
FLAG_HAS_QUAL = 1 << 2   # record has qualities (GL1Q)

# --- GL1F (main) ---
MAGIC_FILE = b"GL1F"
FILE_VERSION = 1
FILE_HDR_FMT = "<4sHHII"   # magic, version(u16), flags(u16), record_count(u32), reserved(u32)
FILE_HDR_SIZE = struct.calcsize(FILE_HDR_FMT)

MAGIC_REC = b"GL1R"
REC_VERSION = 1
REC_HDR_FMT = "<4sHBBHIIHHI"
#  magic(4s)
#  ver(u16)
#  kind(u8)
#  codec(u8)
#  flags(u16)
#  length(u32)
#  chunk_symbols(u32)
#  name_len(u16)
#  reserved(u16)
#  chunk_count(u32)
REC_HDR_SIZE = struct.calcsize(REC_HDR_FMT)

CHUNK_HDR_FMT = "<II"  # sym_count(u32), payload_len(u32)
CHUNK_HDR_SIZE = struct.calcsize(CHUNK_HDR_FMT)

# --- GL1X (sidecar) ---
MAGIC_X = b"GL1X"
X_VERSION = 2
X_HDR_FMT = "<4sHII"  # magic, version(u16), record_count(u32), reserved(u32)
X_HDR_SIZE = struct.calcsize(X_HDR_FMT)

# --- GL1Q (qualities) ---
MAGIC_Q = b"GL1Q"
Q_VERSION = 2
# magic, version(u16), reserved(u16), toc_off(u64), record_count(u32), reserved2(u32)
Q_HDR_FMT = "<4sHHQII"
Q_HDR_SIZE = struct.calcsize(Q_HDR_FMT)

Q_REC_META_FMT = "<III"  # total_len(u32), chunk_symbols(u32), chunk_count(u32)
Q_REC_META_SIZE = struct.calcsize(Q_REC_META_FMT)

Q_CHUNK_HDR_FMT = "<II"  # sym_count(u32), zlen(u32)
Q_CHUNK_HDR_SIZE = struct.calcsize(Q_CHUNK_HDR_FMT)

MAGIC_QTOC = b"QTOC"
QTOC_VERSION = 1
QTOC_HDR_FMT = "<4sHI"  # magic, version(u16), record_count(u32)
QTOC_HDR_SIZE = struct.calcsize(QTOC_HDR_FMT)

# --- GL1I (minimizer index) ---
MAGIC_I = b"GL1I"
I_VERSION = 1
I_FLAG_CANON = 1 << 0
I_FLAG_NUC = 1 << 1      # index built for nucleotide k-mers (2-bit)
I_FLAG_PROT = 1 << 2     # index built for protein k-mers (6-bit)
# magic, version(u16), k(u8), w(u8), flags(u16), bucket_count(u32), record_count(u32), reserved(u32)
I_HDR_FMT = "<4sHBBHIII"
I_HDR_SIZE = struct.calcsize(I_HDR_FMT)

# per-bucket directory entry: offset(u64), key_count(u32), reserved(u32)
I_BUCKET_DIR_FMT = "<QII"
I_BUCKET_DIR_SIZE = struct.calcsize(I_BUCKET_DIR_FMT)


# =============================================================================
# Canonical encoding / hashing utilities
# =============================================================================

def canonical_json(obj: Any) -> str:
    """
    Deterministic JSON encoding used for object hashing.

    - sorted keys
    - compact separators
    - UTF-8
    """
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_bytes(data: bytes) -> bytes:
    return hashlib.sha256(data).digest()


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def obj_id(obj_type: str, core: Dict[str, Any]) -> str:
    """
    Content-addressed object id (git-like), stable across on/off-chain simulation.

    object_id = sha256( b"seqnft:" + type + b"\\0" + canonical_json(core) )
    """
    raw = canonical_json(core).encode("utf-8")
    return sha256_hex(b"seqnft:" + obj_type.encode("utf-8") + b"\0" + raw)


def merkle_root_sha256(leaves: List[bytes]) -> bytes:
    """
    Merkle root using sha256(concat(left,right)) with duplicate-last rule.
    Leaves are 32-byte digests.
    """
    if not leaves:
        return hashlib.sha256(b"").digest()
    level = leaves[:]
    while len(level) > 1:
        nxt: List[bytes] = []
        i = 0
        while i < len(level):
            a = level[i]
            b = level[i + 1] if i + 1 < len(level) else a
            nxt.append(hashlib.sha256(a + b).digest())
            i += 2
        level = nxt
    return level[0]


def merkle_root_sha256_hex(leaves_hex: List[str]) -> str:
    """
    Merkle root returning hex.
    Leaves are provided as hex digests (32 bytes each).
    """
    leaves = [bytes.fromhex(x) for x in leaves_hex]
    return merkle_root_sha256(leaves).hex()


# =============================================================================
# Repo layout
# =============================================================================

def _subdir_for(hex_id: str) -> str:
    return hex_id[:2]


@dataclass
class RepoPaths:
    root: Path
    dot: Path
    objects: Path
    blobs: Path
    refs: Path
    heads: Path
    chain: Path
    chain_tokens: Path
    chain_blobs: Path
    chain_ipfs: Path

    @staticmethod
    def from_root(root: Path) -> "RepoPaths":
        dot = root / ".seqnft"
        return RepoPaths(
            root=root,
            dot=dot,
            objects=dot / "objects",
            blobs=dot / "blobs",
            refs=dot / "refs",
            heads=dot / "refs" / "heads",
            chain=dot / "chain_sim",
            chain_tokens=dot / "chain_sim" / "tokens",
            chain_blobs=dot / "chain_sim" / "onchain_blobs",
            chain_ipfs=dot / "chain_sim" / "ipfs_blobs",
        )


# =============================================================================
# Object store (CAS)
# =============================================================================

class ObjectStore:
    def __init__(self, paths: RepoPaths):
        self.p = paths

    def ensure_dirs(self) -> None:
        for d in [
            self.p.dot, self.p.objects, self.p.blobs, self.p.refs, self.p.heads,
            self.p.chain, self.p.chain_tokens, self.p.chain_blobs, self.p.chain_ipfs
        ]:
            d.mkdir(parents=True, exist_ok=True)

    def put_blob(self, data: bytes) -> str:
        bid = sha256_hex(data)
        sub = self.p.blobs / _subdir_for(bid)
        sub.mkdir(parents=True, exist_ok=True)
        fp = sub / f"{bid}.bin"
        if not fp.exists():
            fp.write_bytes(data)
        return bid

    def get_blob(self, bid: str) -> bytes:
        fp = self.p.blobs / _subdir_for(bid) / f"{bid}.bin"
        return fp.read_bytes()

    def has_blob(self, bid: str) -> bool:
        return (self.p.blobs / _subdir_for(bid) / f"{bid}.bin").exists()

    def put_object(self, obj_type: str, core: Dict[str, Any]) -> str:
        oid = obj_id(obj_type, core)
        sub = self.p.objects / _subdir_for(oid)
        sub.mkdir(parents=True, exist_ok=True)
        fp = sub / f"{oid}.json"
        if not fp.exists():
            fp.write_text(canonical_json({"type": obj_type, "core": core}) + "\n", encoding="utf-8")
        return oid

    def get_object(self, oid: str) -> Tuple[str, Dict[str, Any]]:
        fp = self.p.objects / _subdir_for(oid) / f"{oid}.json"
        data = json.loads(fp.read_text(encoding="utf-8"))
        return data["type"], data["core"]

    def has_object(self, oid: str) -> bool:
        return (self.p.objects / _subdir_for(oid) / f"{oid}.json").exists()


# =============================================================================
# Refs / HEAD (git-like)
# =============================================================================

class Refs:
    def __init__(self, paths: RepoPaths):
        self.p = paths

    def head_path(self) -> Path:
        return self.p.refs / "HEAD"

    def main_path(self) -> Path:
        return self.p.heads / "main"

    def get_head(self) -> Optional[str]:
        hp = self.head_path()
        if hp.exists():
            return hp.read_text(encoding="utf-8").strip() or None
        return None

    def set_head(self, commit_id: str) -> None:
        self.head_path().write_text(commit_id + "\n", encoding="utf-8")
        self.main_path().write_text(commit_id + "\n", encoding="utf-8")

    def get_main(self) -> Optional[str]:
        mp = self.main_path()
        if mp.exists():
            return mp.read_text(encoding="utf-8").strip() or None
        return None


# =============================================================================
# Chain simulator (deployments)
# =============================================================================

class ChainSim:
    """
    Local "chain" simulation:
      - tokens are JSON files keyed by tokenId (object id hex)
      - fully-onchain chunks copy blob bytes into chain_sim/onchain_blobs/
      - anchor/ipfs chunks copy blob bytes into chain_sim/ipfs_blobs/ and record a CID-like string

    This mimics the on-chain/off-chain dual mode without requiring Solidity.
    """

    def __init__(self, paths: RepoPaths, store: ObjectStore):
        self.p = paths
        self.store = store

    def _token_file(self, token_id: str) -> Path:
        sub = self.p.chain_tokens / _subdir_for(token_id)
        sub.mkdir(parents=True, exist_ok=True)
        return sub / f"{token_id}.json"

    def mint_token(self, token_id: str, payload: Dict[str, Any]) -> None:
        fp = self._token_file(token_id)
        if fp.exists():
            return
        fp.write_text(canonical_json(payload) + "\n", encoding="utf-8")

    def get_token(self, token_id: str) -> Dict[str, Any]:
        fp = self._token_file(token_id)
        return json.loads(fp.read_text(encoding="utf-8"))

    def put_onchain_blob(self, blob_id: str) -> str:
        data = self.store.get_blob(blob_id)
        sub = self.p.chain_blobs / _subdir_for(blob_id)
        sub.mkdir(parents=True, exist_ok=True)
        fp = sub / f"{blob_id}.bin"
        if not fp.exists():
            fp.write_bytes(data)
        return f"chain://{blob_id}"

    def put_ipfs_blob(self, blob_id: str) -> str:
        data = self.store.get_blob(blob_id)
        sub = self.p.chain_ipfs / _subdir_for(blob_id)
        sub.mkdir(parents=True, exist_ok=True)
        fp = sub / f"{blob_id}.bin"
        if not fp.exists():
            fp.write_bytes(data)
        # pseudo CID: stable, content-addressed
        return f"ipfs://sha256-{blob_id}"


# =============================================================================
# GL1X / GL1Q / GL1I helper structures
# =============================================================================

@dataclass
class XRec:
    name: str
    record_off: int
    kind: int
    codec: int
    flags: int
    length: int
    chunk_symbols: int
    chunk_count: int
    chunk_offs: List[int]                 # absolute offsets to CHUNK_HDR in .gl1
    mask_intervals: List[Tuple[int, int]] # absolute (start, len)
    excp: List[Tuple[int, int]]           # absolute (pos, code) where code is DNA4 nibble (0..15)
    chunk_hashes: List[bytes]             # sha256(payload) per chunk (32 bytes)
    merkle_root: bytes                    # merkle_root_sha256(chunk_hashes)


@dataclass
class QTocEntry:
    name: str
    offset: int


# =============================================================================
# File path helpers
# =============================================================================

def gl1x_path(gl1_path: str) -> str:
    if gl1_path.endswith(".gl1"):
        return gl1_path[:-4] + ".gl1x"
    return gl1_path + ".gl1x"


def gl1q_path(gl1_path: str) -> str:
    if gl1_path.endswith(".gl1"):
        return gl1_path[:-4] + ".gl1q"
    return gl1_path + ".gl1q"


def gl1i_path(gl1_path: str) -> str:
    if gl1_path.endswith(".gl1"):
        return gl1_path[:-4] + ".gl1i"
    return gl1_path + ".gl1i"


# =============================================================================
# GL1X (sidecar) read/write
# =============================================================================

def write_gl1x(path: str, xrecs: List[XRec]) -> None:
    """
    Write GL1X v2 sidecar.

    The format is defined in SEQNFT_PIPELINE_SPEC_v2.md.
    """
    with open(path, "wb") as f:
        f.write(struct.pack(X_HDR_FMT, MAGIC_X, X_VERSION, len(xrecs), 0))
        for xr in xrecs:
            nb = xr.name.encode("utf-8")
            f.write(struct.pack("<H", len(nb)))
            f.write(nb)
            f.write(struct.pack("<QBBHIII", int(xr.record_off), int(xr.kind), int(xr.codec), int(xr.flags),
                                int(xr.length), int(xr.chunk_symbols), int(xr.chunk_count)))
            f.write(struct.pack("<II", len(xr.mask_intervals), len(xr.excp)))

            # chunk offsets
            for off in xr.chunk_offs:
                f.write(struct.pack("<Q", int(off)))

            # mask intervals (abs)
            for s, ln in xr.mask_intervals:
                f.write(struct.pack("<II", int(s), int(ln)))

            # exceptions (abs): (pos u32, code u8)
            for p, code in xr.excp:
                f.write(struct.pack("<IB", int(p), int(code) & 0xFF))

            # chunk hashes
            for h in xr.chunk_hashes:
                if len(h) != 32:
                    raise ValueError("chunk_hash must be 32 bytes")
                f.write(h)

            if len(xr.merkle_root) != 32:
                raise ValueError("merkle_root must be 32 bytes")
            f.write(xr.merkle_root)


def read_gl1x(path: str) -> List[XRec]:
    with open(path, "rb") as f:
        hdr = f.read(X_HDR_SIZE)
        if len(hdr) != X_HDR_SIZE:
            raise ValueError("Bad GL1X header")
        magic, ver, rec_count, _ = struct.unpack(X_HDR_FMT, hdr)
        if magic != MAGIC_X:
            raise ValueError("Not a GL1X file")
        if ver != X_VERSION:
            raise ValueError(f"Unsupported GL1X version: {ver}")

        out: List[XRec] = []
        for _ in range(rec_count):
            (nlen,) = struct.unpack("<H", f.read(2))
            name = f.read(nlen).decode("utf-8", "replace")
            record_off, kind, codec, flags, length, chunk_symbols, chunk_count = struct.unpack("<QBBHIII", f.read(struct.calcsize("<QBBHIII")))
            mask_count, excp_count = struct.unpack("<II", f.read(8))

            chunk_offs = [struct.unpack("<Q", f.read(8))[0] for _ in range(chunk_count)]
            mask = [struct.unpack("<II", f.read(8)) for _ in range(mask_count)]
            excp: List[Tuple[int, int]] = []
            for _j in range(excp_count):
                p, code = struct.unpack("<IB", f.read(5))
                excp.append((p, code))

            chunk_hashes = [f.read(32) for _ in range(chunk_count)]
            merkle_root = f.read(32)

            out.append(XRec(
                name=name,
                record_off=record_off,
                kind=kind,
                codec=codec,
                flags=flags,
                length=length,
                chunk_symbols=chunk_symbols,
                chunk_count=chunk_count,
                chunk_offs=chunk_offs,
                mask_intervals=mask,
                excp=excp,
                chunk_hashes=chunk_hashes,
                merkle_root=merkle_root,
            ))
        return out


# =============================================================================
# GL1Q (qualities) read/write
# =============================================================================

def write_gl1q_init(fq: io.BufferedRandom) -> None:
    """
    Initialize a GL1Q file with a placeholder header; finalize_gl1q() will fix it.
    """
    fq.seek(0)
    fq.write(struct.pack(Q_HDR_FMT, MAGIC_Q, Q_VERSION, 0, 0, 0, 0))


def finalize_gl1q(fq: io.BufferedRandom, toc: List[QTocEntry], record_count: int) -> None:
    """
    Append TOC and patch header with toc offset and record_count.
    """
    toc_off = fq.tell()
    fq.write(struct.pack(QTOC_HDR_FMT, MAGIC_QTOC, QTOC_VERSION, record_count))
    for e in toc:
        nb = e.name.encode("utf-8")
        fq.write(struct.pack("<H", len(nb)))
        fq.write(nb)
        fq.write(struct.pack("<Q", int(e.offset)))

    # patch header
    fq.seek(0)
    fq.write(struct.pack(Q_HDR_FMT, MAGIC_Q, Q_VERSION, 0, int(toc_off), int(record_count), 0))
    fq.flush()


def read_gl1q_toc(path: str) -> Dict[str, int]:
    """
    Return {record_name: offset} for GL1Q.
    """
    with open(path, "rb") as f:
        hdr = f.read(Q_HDR_SIZE)
        if len(hdr) != Q_HDR_SIZE:
            raise ValueError("Bad GL1Q header")
        magic, ver, _rsv, toc_off, rec_count, _rsv2 = struct.unpack(Q_HDR_FMT, hdr)
        if magic != MAGIC_Q:
            raise ValueError("Not a GL1Q file")
        if ver != Q_VERSION:
            raise ValueError(f"Unsupported GL1Q version: {ver}")

        if toc_off == 0:
            return {}

        f.seek(toc_off)
        t_hdr = f.read(QTOC_HDR_SIZE)
        t_magic, t_ver, t_count = struct.unpack(QTOC_HDR_FMT, t_hdr)
        if t_magic != MAGIC_QTOC:
            raise ValueError("Bad GL1Q TOC magic")
        if t_ver != QTOC_VERSION:
            raise ValueError(f"Unsupported GL1Q TOC version: {t_ver}")
        if t_count != rec_count:
            # allow but warn by ignoring header count mismatch
            pass

        out: Dict[str, int] = {}
        for _ in range(t_count):
            (nlen,) = struct.unpack("<H", f.read(2))
            name = f.read(nlen).decode("utf-8", "replace")
            (off,) = struct.unpack("<Q", f.read(8))
            out[name] = off
        return out


def read_gl1q_record_comp_chunks(gl1q_path_: str, toc: Dict[str, int], name: str) -> List[bytes]:
    """
    Read a GL1Q record block and return per-chunk compressed bytes:
      (Q_CHUNK_HDR + zbytes) for each chunk.
    """
    off = toc.get(name)
    if off is None:
        raise KeyError(name)
    out: List[bytes] = []
    with open(gl1q_path_, "rb") as fq:
        fq.seek(off)
        (nlen,) = struct.unpack("<H", fq.read(2))
        nm = fq.read(nlen).decode("utf-8", "replace")
        if nm != name:
            raise ValueError("GL1Q TOC mismatch")
        total_len, chunk_symbols, chunk_count = struct.unpack(Q_REC_META_FMT, fq.read(Q_REC_META_SIZE))
        for _ in range(chunk_count):
            hdr = fq.read(Q_CHUNK_HDR_SIZE)
            if len(hdr) != Q_CHUNK_HDR_SIZE:
                raise ValueError("Truncated GL1Q chunk header")
            sym, zlen = struct.unpack(Q_CHUNK_HDR_FMT, hdr)
            z = fq.read(zlen)
            if len(z) != zlen:
                raise ValueError("Truncated GL1Q chunk payload")
            out.append(hdr + z)
    return out


# =============================================================================
# GL1I (minimizer index) read/write + utilities
# =============================================================================

def _kmer_key_nuc2(kmer: str, canonical: bool) -> Optional[int]:
    """
    Encode an A/C/G/T kmer into a uint64 using 2-bit encoding.
    Returns None if kmer contains non-ACGT.
    """
    km = kmer.upper()
    v = 0
    for ch in km:
        v <<= 2
        if ch == "A":
            pass
        elif ch == "C":
            v |= 1
        elif ch == "G":
            v |= 2
        elif ch == "T":
            v |= 3
        else:
            return None

    if canonical:
        # compute reverse complement key as well and take min
        # Reverse complement of A/C/G/T
        rc_v = 0
        for ch in reversed(km):
            rc_v <<= 2
            if ch == "A":
                rc_v |= 3
            elif ch == "C":
                rc_v |= 2
            elif ch == "G":
                rc_v |= 1
            elif ch == "T":
                rc_v |= 0
        v = min(v, rc_v)
    return v


def _kmer_key_sixbit(kmer: str, canonical: bool) -> Optional[int]:
    """
    Encode a SIXBIT kmer into a uint64 (6 bits per symbol).
    Returns None if a symbol isn't in SIXBIT alphabet or if it overflows 64 bits.
    """
    if gl1enc is None:
        raise RuntimeError("gl1enc not available")
    km = kmer.upper()
    bits = 6 * len(km)
    if bits > 64:
        return None

    v = 0
    for ch in km:
        if ch not in gl1enc._SIXBIT_CHAR_TO_CODE:  # type: ignore
            return None
        v = (v << 6) | int(gl1enc._SIXBIT_CHAR_TO_CODE[ch])  # type: ignore

    if canonical:
        # canonical for proteins is not well-defined; keep as-is.
        # We keep this flag for API symmetry.
        pass
    return v


def write_gl1i(path: str, *, k: int, w: int, canonical: bool, bucket_count: int,
              records: List[Tuple[str, int]], buckets: List[Dict[int, List[Tuple[int, int]]]], flags: int) -> None:
    """
    Write GL1I v1 to `path`.

    buckets[b][key] = list of (rec_id, pos) postings
    """
    if bucket_count != len(buckets):
        raise ValueError("bucket_count mismatch")
    if k <= 0 or w <= 0:
        raise ValueError("k and w must be positive")
    if bucket_count <= 0:
        raise ValueError("bucket_count must be positive")

    with open(path, "wb") as f:
        f.write(struct.pack(I_HDR_FMT, MAGIC_I, I_VERSION, int(k) & 0xFF, int(w) & 0xFF,
                            int(flags), int(bucket_count), int(len(records)), 0))

        # records table
        for name, length in records:
            nb = name.encode("utf-8")
            if len(nb) > 65535:
                raise ValueError("Record name too long")
            f.write(struct.pack("<H", len(nb)))
            f.write(nb)
            f.write(struct.pack("<I", int(length)))

        # placeholder bucket directory
        bucket_dir_off = f.tell()
        for _ in range(bucket_count):
            f.write(struct.pack(I_BUCKET_DIR_FMT, 0, 0, 0))

        # bucket payloads
        bucket_offsets: List[int] = []
        bucket_key_counts: List[int] = []

        for b in range(bucket_count):
            bucket_offsets.append(f.tell())
            key_map = buckets[b]
            keys_sorted = sorted(key_map.keys())
            bucket_key_counts.append(len(keys_sorted))
            f.write(struct.pack("<I", len(keys_sorted)))
            for key in keys_sorted:
                postings = key_map[key]
                f.write(struct.pack("<Q", int(key)))
                f.write(struct.pack("<I", len(postings)))
                for rec_id, pos in postings:
                    f.write(struct.pack("<II", int(rec_id), int(pos)))

        # patch bucket directory
        end_off = f.tell()
        f.seek(bucket_dir_off)
        for b in range(bucket_count):
            f.write(struct.pack(I_BUCKET_DIR_FMT, int(bucket_offsets[b]), int(bucket_key_counts[b]), 0))
        f.seek(end_off)


def read_gl1i_header(path: str) -> Tuple[int, int, int, int, List[Tuple[str, int]], int]:
    """
    Read GL1I header and record table and bucket directory offset.

    Returns: (k, w, flags, bucket_count, records, bucket_dir_off)
    """
    with open(path, "rb") as f:
        hdr = f.read(I_HDR_SIZE)
        if len(hdr) != I_HDR_SIZE:
            raise ValueError("Bad GL1I header")
        magic, ver, k, w, flags, bucket_count, record_count, _rsv = struct.unpack(I_HDR_FMT, hdr)
        if magic != MAGIC_I:
            raise ValueError("Not a GL1I file")
        if ver != I_VERSION:
            raise ValueError(f"Unsupported GL1I version: {ver}")

        records: List[Tuple[str, int]] = []
        for _ in range(record_count):
            (nlen,) = struct.unpack("<H", f.read(2))
            name = f.read(nlen).decode("utf-8", "replace")
            (length,) = struct.unpack("<I", f.read(4))
            records.append((name, length))

        bucket_dir_off = f.tell()
        # skip directory here; read in query
        return int(k), int(w), int(flags), int(bucket_count), records, bucket_dir_off


def query_gl1i(path: str, key: int) -> List[Tuple[int, int]]:
    """
    Query GL1I for a key and return postings [(rec_id, pos), ...].
    """
    with open(path, "rb") as f:
        hdr = f.read(I_HDR_SIZE)
        magic, ver, k, w, flags, bucket_count, record_count, _rsv = struct.unpack(I_HDR_FMT, hdr)
        if magic != MAGIC_I or ver != I_VERSION:
            raise ValueError("Not a supported GL1I file")

        # read record table
        for _ in range(record_count):
            (nlen,) = struct.unpack("<H", f.read(2))
            f.seek(nlen, io.SEEK_CUR)
            f.seek(4, io.SEEK_CUR)

        # read bucket directory
        bucket_dir: List[Tuple[int, int]] = []  # (offset, key_count)
        for _ in range(bucket_count):
            off, key_count, _r = struct.unpack(I_BUCKET_DIR_FMT, f.read(I_BUCKET_DIR_SIZE))
            bucket_dir.append((off, key_count))

        b = int(key) % int(bucket_count)
        off, key_count = bucket_dir[b]
        if off == 0:
            return []
        f.seek(off)
        (kc2,) = struct.unpack("<I", f.read(4))
        # kc2 should match key_count
        for _ in range(kc2):
            (kkey,) = struct.unpack("<Q", f.read(8))
            (pcnt,) = struct.unpack("<I", f.read(4))
            if kkey == int(key):
                hits = [struct.unpack("<II", f.read(8)) for _ in range(pcnt)]
                return [(int(r), int(p)) for (r, p) in hits]
            else:
                f.seek(pcnt * 8, io.SEEK_CUR)
        return []


# =============================================================================
# SeqNFT Repo facade
# =============================================================================

class SeqNFTRepo:
    def __init__(self, root: Path):
        self.paths = RepoPaths.from_root(root)
        self.store = ObjectStore(self.paths)
        self.refs = Refs(self.paths)
        self.chain = ChainSim(self.paths, self.store)

    @staticmethod
    def discover(start: Path) -> "SeqNFTRepo":
        cur = start.resolve()
        for p in [cur] + list(cur.parents):
            if (p / ".seqnft").exists():
                return SeqNFTRepo(p)
        raise SystemExit("Not inside a SeqNFT repo (missing .seqnft). Run `seqnft init` first.")

    def init(self) -> None:
        self.store.ensure_dirs()
        cfg = self.paths.dot / "config.json"
        if not cfg.exists():
            cfg.write_text(canonical_json({
                "spec": "seqnft-gl1-v2",
                "hash": "sha256",
                "gl1enc": getattr(gl1enc, "GL1ENC_FORMAT", None) if gl1enc else None,
                "created": int(time.time()),
            }) + "\n", encoding="utf-8")

    # -------------------------------------------------------------------------
    # Object builders
    # -------------------------------------------------------------------------

    def _mk_chunk_object(
        self,
        chunk_blob: bytes,
        kind: int,
        codec: int,
        sym_count: int,
        payload_hash_hex: str,
        mask_rel: List[Tuple[int, int]],
        excp_rel: List[Tuple[int, int]],
        qual_chunk_bytes: Optional[bytes],
    ) -> str:
        """
        Store:
          - chunk_blob (CHUNK_HDR + payload)
          - optional quality chunk bytes (Q_CHUNK_HDR + zlib bytes)
          - chunk-local metadata (mask/excp) in canonical JSON
        """
        blob_id = self.store.put_blob(chunk_blob)
        qual_blob_id = self.store.put_blob(qual_chunk_bytes) if qual_chunk_bytes is not None else None

        core = {
            "v": 2,
            "kind": int(kind),
            "codec": int(codec),
            "sym_count": int(sym_count),
            "chunk_blob": blob_id,
            "payload_sha256": payload_hash_hex,
            "mask": [[int(a), int(b)] for a, b in mask_rel],
            "excp": [[int(p), int(c)] for p, c in excp_rel],
            "qual_chunk": qual_blob_id,
        }
        return self.store.put_object("chunk", core)

    def _mk_chrom_object(
        self,
        name: str,
        kind: int,
        codec: int,
        flags: int,
        length: int,
        chunk_symbols: int,
        chunk_ids: List[str],
        payload_hashes_hex: List[str],
        mask_leaves_hex: List[str],
        excp_leaves_hex: List[str],
        qual_leaves_hex: List[str],
    ) -> str:
        seq_payload_root = merkle_root_sha256_hex(payload_hashes_hex)

        mask_root = merkle_root_sha256_hex(mask_leaves_hex) if mask_leaves_hex else ""
        excp_root = merkle_root_sha256_hex(excp_leaves_hex) if excp_leaves_hex else ""
        qual_root = merkle_root_sha256_hex(qual_leaves_hex) if qual_leaves_hex else ""

        core = {
            "v": 2,
            "name": name,
            "kind": int(kind),
            "codec": int(codec),
            "flags": int(flags),
            "length": int(length),
            "chunk_symbols": int(chunk_symbols),
            "chunks": chunk_ids,
            "seq_payload_root": seq_payload_root,
            "mask_root": mask_root,
            "excp_root": excp_root,
            "qual_root": qual_root,
        }
        return self.store.put_object("chrom", core)

    def _mk_assembly_object(self, assembly_name: str, chrom_ids: List[str]) -> str:
        # assembly root is merkle over sha256(chrom_id_bytes) leaves (pipeline-level)
        chrom_leaves = [sha256_hex(bytes.fromhex(cid)) for cid in chrom_ids]
        assembly_root = merkle_root_sha256_hex(chrom_leaves) if chrom_leaves else merkle_root_sha256_hex([])
        core = {
            "v": 2,
            "assembly": assembly_name,
            "created": int(time.time()),
            "chromosomes": chrom_ids,
            "assembly_root": assembly_root,
        }
        return self.store.put_object("assembly", core)

    def _mk_commit_object(
        self,
        tree_id: str,
        parents: List[str],
        author: str,
        message: str,
        tags: Optional[Dict[str, Any]] = None,
    ) -> str:
        core = {
            "v": 2,
            "tree": tree_id,
            "parents": parents,
            "author": author,
            "message": message,
            "time": int(time.time()),
            "tags": tags or {},
        }
        return self.store.put_object("commit", core)

    # -------------------------------------------------------------------------
    # Helpers: input parsing and normalization
    # -------------------------------------------------------------------------

    @staticmethod
    def _open_text_auto(path: str):
        if path.endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return open(path, "rt", encoding="utf-8", errors="replace")

    @staticmethod
    def _detect_format(path: str) -> str:
        with SeqNFTRepo._open_text_auto(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("@"):
                    return "fastq"
                break
        raise ValueError("Cannot detect input format; expected FASTA (>) or FASTQ (@).")

    @staticmethod
    def _remove_ws_keep_case(s: str) -> str:
        # remove ASCII whitespace used by GL1ENC normalization, but keep case
        return s.replace(" ", "").replace("\n", "").replace("\r", "").replace("\t", "")

    @staticmethod
    def _mask_intervals_from_raw(raw_no_ws: str) -> List[Tuple[int, int]]:
        """
        Compute lowercase intervals from a sequence that has no whitespace and keeps original case.
        """
        out: List[Tuple[int, int]] = []
        run: Optional[int] = None
        for i, ch in enumerate(raw_no_ws):
            is_low = "a" <= ch <= "z"
            if is_low and run is None:
                run = i
            if (not is_low) and run is not None:
                out.append((run, i - run))
                run = None
        if run is not None:
            out.append((run, len(raw_no_ws) - run))
        return out

    @staticmethod
    def _guess_kind(normalized_upper: str, *, assume_fastq: bool) -> int:
        """
        Heuristic kind inference. FASTQ is assumed nucleotide.
        """
        if assume_fastq:
            return KIND_NUC
        if not normalized_upper:
            return KIND_NUC
        nuc_alphabet = set("ACGTUNRYSMWKVHDBN-")
        nuc_hits = sum(1 for ch in normalized_upper if ch in nuc_alphabet)
        frac = nuc_hits / max(1, len(normalized_upper))
        return KIND_NUC if frac >= 0.90 else KIND_PROT

    @staticmethod
    def _sanitize_unknowns(seq: str, kind: int, unknown_to_x: bool) -> str:
        """
        Replace unknown characters with 'N' (nuc) or 'X' (prot) if unknown_to_x.
        Otherwise return unchanged.
        """
        if not unknown_to_x:
            return seq
        if kind == KIND_NUC:
            allowed = set(gl1enc._DNA4_CHAR_TO_NIBBLE.keys()) if gl1enc is not None else set("ACGTN-")
            repl = "N"
        else:
            allowed = set(gl1enc._SIXBIT_CHAR_TO_CODE.keys()) if gl1enc is not None else set("ACDEFGHIKLMNPQRSTVWYBXZ*")
            repl = "X"
        return "".join((ch if ch in allowed else repl) for ch in seq)

    @staticmethod
    def _select_codec(
        kind: int,
        norm_upper: str,
        mode: str,
        nuc2_excp_threshold: float,
        strict: bool,
        unknown_to_x: bool,
    ) -> Tuple[int, bool]:
        """
        Decide record codec and whether to use exceptions with DNA2 payloads.

        Returns: (codec_id, use_exceptions)
        """
        if gl1enc is None:
            raise RuntimeError(f"gl1enc_v2 not available: {_GL1ENC_IMPORT_ERROR}")

        s = norm_upper

        # Optionally sanitize unknowns to keep compact encodings viable.
        s2 = SeqNFTRepo._sanitize_unknowns(s, kind, unknown_to_x)

        if kind == KIND_NUC:
            if mode == "auto":
                if gl1enc.is_pure_acgt(s2, assume_normalized=True):
                    return CODEC_DNA2, False
                if gl1enc.is_dna4_eligible(s2, assume_normalized=True):
                    non_acgt = sum(1 for ch in s2 if ch not in "ACGT")
                    frac = non_acgt / max(1, len(s2))
                    if frac <= nuc2_excp_threshold:
                        return CODEC_DNA2, True
                    return CODEC_DNA4, False
                # not DNA4 eligible => ASCII fallback
                return CODEC_ASCII, False

            if mode == "2bit":
                # DNA2; if not pure ACGT we store exceptions where possible
                if gl1enc.is_pure_acgt(s2, assume_normalized=True):
                    return CODEC_DNA2, False
                # only safe if exceptions can represent symbols (DNA4)
                if gl1enc.is_dna4_eligible(s2, assume_normalized=True):
                    return CODEC_DNA2, True
                if strict:
                    raise ValueError("2bit mode requested but sequence has non-DNA4 symbols (cannot store exceptions).")
                return CODEC_ASCII, False

            if mode == "4bit":
                if gl1enc.is_dna4_eligible(s2, assume_normalized=True):
                    return CODEC_DNA4, False
                if strict:
                    raise ValueError("4bit mode requested but sequence is not DNA4-eligible.")
                return CODEC_ASCII, False

            if mode == "6bit":
                # not meaningful for nucleotide; treat as auto
                return SeqNFTRepo._select_codec(kind, s2, "auto", nuc2_excp_threshold, strict, unknown_to_x)

            raise ValueError(f"Unknown mode: {mode}")

        # KIND_PROT
        if mode == "auto":
            if gl1enc.is_sixbit_eligible(s2, assume_normalized=True):
                return CODEC_SIXBIT, False
            return CODEC_ASCII, False

        if mode == "6bit":
            if gl1enc.is_sixbit_eligible(s2, assume_normalized=True):
                return CODEC_SIXBIT, False
            if strict:
                raise ValueError("6bit mode requested but sequence is not SIXBIT-eligible.")
            return CODEC_ASCII, False

        # 2bit/4bit requested for protein => fallback or error
        if strict:
            raise ValueError(f"{mode} mode requested but record looks like protein; use --mode auto or --mode 6bit.")
        return CODEC_ASCII, False

    # -------------------------------------------------------------------------
    # High-level operations
    # -------------------------------------------------------------------------

    def mint_from_fasta_fastq(
        self,
        inp: str,
        assembly_name: str,
        chunk_symbols: int,
        mode: str,
        strict: bool,
        unknown_to_x: bool,
        store_qual: bool,
        nuc2_excp_threshold: float,
        author: str,
        message: str,
    ) -> str:
        """
        Ingest FASTA/FASTQ directly into the SeqNFT object graph using GL1ENCv2.

        - Splits each record into chunks of size `chunk_symbols` (except last).
        - Encodes each chunk payload using codec selected by mode/auto.
        - Stores mask intervals to preserve lowercase when requested during `view`.
        - Stores exceptions list only when codec is DNA2 and record contains non-ACGT symbols.
        - Stores per-chunk qualities (zlib) if input is FASTQ and store_qual is enabled.

        Returns: commit object id (hex)
        """
        if gl1enc is None:
            raise SystemExit(f"gl1enc_v2.py is required. Import error: {_GL1ENC_IMPORT_ERROR}")

        self.init()

        parent = self.refs.get_head()
        parents = [parent] if parent else []

        fmt = self._detect_format(inp)

        # Parse records
        records: List[Tuple[str, str, Optional[str], List[Tuple[int, int]]]] = []
        # tuple: (name, normalized_upper_seq, qual_or_none, mask_intervals_abs)

        if fmt == "fasta":
            with self._open_text_auto(inp) as f:
                name: Optional[str] = None
                seq_parts: List[str] = []
                for line in f:
                    if not line:
                        continue
                    if line.startswith(">"):
                        # flush previous
                        if name is not None:
                            raw = self._remove_ws_keep_case("".join(seq_parts))
                            mask = self._mask_intervals_from_raw(raw)
                            norm = raw.upper()
                            records.append((name, norm, None, mask))
                        header = line[1:].strip()
                        name = header.split()[0] if header else f"rec{len(records)}"
                        seq_parts = []
                    else:
                        seq_parts.append(line.strip())
                # flush last
                if name is not None:
                    raw = self._remove_ws_keep_case("".join(seq_parts))
                    mask = self._mask_intervals_from_raw(raw)
                    norm = raw.upper()
                    records.append((name, norm, None, mask))

        elif fmt == "fastq":
            with self._open_text_auto(inp) as f:
                while True:
                    h = f.readline()
                    if not h:
                        break
                    h = h.rstrip("\n\r")
                    if not h:
                        continue
                    if not h.startswith("@"):
                        raise ValueError("FASTQ parse error: expected '@' header")
                    name = h[1:].strip().split()[0] or f"rec{len(records)}"

                    # sequence lines until '+' line
                    seq_lines: List[str] = []
                    while True:
                        line = f.readline()
                        if not line:
                            raise ValueError("FASTQ parse error: unexpected EOF in sequence block")
                        line = line.rstrip("\n\r")
                        if line.startswith("+"):
                            break
                        seq_lines.append(line.strip())

                    raw_seq = self._remove_ws_keep_case("".join(seq_lines))
                    mask = self._mask_intervals_from_raw(raw_seq)
                    norm_seq = raw_seq.upper()

                    # quality lines until length matches seq length
                    q_lines: List[str] = []
                    q_len = 0
                    while q_len < len(raw_seq):
                        ql = f.readline()
                        if not ql:
                            raise ValueError("FASTQ parse error: unexpected EOF in quality block")
                        ql = ql.rstrip("\n\r")
                        q_lines.append(ql)
                        q_len += len(ql)
                    qual = "".join(q_lines)
                    if len(qual) != len(raw_seq):
                        raise ValueError("FASTQ parse error: quality length != sequence length")

                    records.append((name, norm_seq, qual if store_qual else None, mask))
        else:
            raise ValueError(f"Unsupported input format: {fmt}")

        if not records:
            raise ValueError("No records found in input")

        chrom_ids: List[str] = []

        for name, norm_seq, qual, mask_abs in records:
            kind = self._guess_kind(norm_seq, assume_fastq=(fmt == "fastq"))
            # sanitize unknowns only if requested; this may change content, so it's opt-in
            norm_seq2 = self._sanitize_unknowns(norm_seq, kind, unknown_to_x)

            codec, use_excp = self._select_codec(kind, norm_seq2, mode, nuc2_excp_threshold, strict, unknown_to_x)

            length = len(norm_seq2)
            cs = chunk_symbols if chunk_symbols > 0 else length
            chunk_count = (length + cs - 1) // cs if cs > 0 else 0

            flags = 0
            if mask_abs:
                flags |= FLAG_HAS_MASK
            if use_excp:
                flags |= FLAG_HAS_EXCP
            if qual is not None:
                flags |= FLAG_HAS_QUAL

            payload_hashes_hex: List[str] = []
            mask_leaf_hex: List[str] = []
            excp_leaf_hex: List[str] = []
            qual_leaf_hex: List[str] = []
            chunk_oids: List[str] = []

            for ci in range(chunk_count):
                chunk_start = ci * cs
                chunk_end = min(length, chunk_start + cs)
                chunk_seq = norm_seq2[chunk_start:chunk_end]
                sym = len(chunk_seq)

                # mask relative to chunk
                mask_rel = _slice_intervals(mask_abs, chunk_start, chunk_end)

                # quality chunk (optional)
                qc_bytes: Optional[bytes] = None
                if qual is not None:
                    q_chunk = qual[chunk_start:chunk_end].encode("ascii", "replace")
                    z = zlib.compress(q_chunk, 6)
                    qc_bytes = struct.pack(Q_CHUNK_HDR_FMT, len(q_chunk), len(z)) + z

                excp_rel: List[Tuple[int, int]] = []

                # encode payload according to record codec
                if codec == CODEC_DNA2:
                    base_chars: List[str] = []
                    for j, ch in enumerate(chunk_seq):
                        if ch in "ACGT":
                            base_chars.append(ch)
                        else:
                            base_chars.append("A")
                            if use_excp:
                                # must be DNA4 representable by construction (or sanitized via unknown_to_x)
                                if ch in gl1enc._DNA4_CHAR_TO_NIBBLE:  # type: ignore
                                    code = int(gl1enc._DNA4_CHAR_TO_NIBBLE[ch])  # type: ignore
                                elif unknown_to_x:
                                    code = int(gl1enc._DNA4_CHAR_TO_NIBBLE["N"])  # type: ignore
                                else:
                                    raise ValueError(f"Symbol not representable for exceptions: {ch!r}")
                                excp_rel.append((j, code))
                    payload = gl1enc.pack_dna2("".join(base_chars))  # type: ignore

                elif codec == CODEC_DNA4:
                    if not gl1enc.is_dna4_eligible(chunk_seq, assume_normalized=True):  # type: ignore
                        raise ValueError("Chunk not DNA4-eligible under DNA4 record codec")
                    payload = gl1enc.pack_dna4(chunk_seq)  # type: ignore

                elif codec == CODEC_SIXBIT:
                    if not gl1enc.is_sixbit_eligible(chunk_seq, assume_normalized=True):  # type: ignore
                        raise ValueError("Chunk not SIXBIT-eligible under SIXBIT record codec")
                    payload = gl1enc.pack_sixbit(chunk_seq)  # type: ignore

                else:
                    payload = chunk_seq.encode("ascii", "replace")

                # build chunk blob: CHUNK_HDR + payload
                chunk_blob = struct.pack(CHUNK_HDR_FMT, sym, len(payload)) + payload

                ph = sha256_hex(payload)
                payload_hashes_hex.append(ph)

                # per-chunk meta leaves for chrom roots
                mask_leaf_hex.append(sha256_hex(canonical_json([[a, b] for a, b in mask_rel]).encode("utf-8")))
                excp_leaf_hex.append(sha256_hex(canonical_json([[p, c] for p, c in excp_rel]).encode("utf-8")))
                qual_leaf_hex.append(sha256_hex(qc_bytes) if qc_bytes is not None else sha256_hex(b""))

                chunk_oid = self._mk_chunk_object(
                    chunk_blob=chunk_blob,
                    kind=kind,
                    codec=codec,
                    sym_count=sym,
                    payload_hash_hex=ph,
                    mask_rel=mask_rel,
                    excp_rel=excp_rel,
                    qual_chunk_bytes=qc_bytes,
                )
                chunk_oids.append(chunk_oid)

            chrom_oid = self._mk_chrom_object(
                name=name,
                kind=kind,
                codec=codec,
                flags=flags,
                length=length,
                chunk_symbols=cs,
                chunk_ids=chunk_oids,
                payload_hashes_hex=payload_hashes_hex,
                mask_leaves_hex=mask_leaf_hex,
                excp_leaves_hex=excp_leaf_hex,
                qual_leaves_hex=qual_leaf_hex if qual is not None else [],
            )
            chrom_ids.append(chrom_oid)

        assembly_oid = self._mk_assembly_object(assembly_name=assembly_name, chrom_ids=chrom_ids)
        commit_oid = self._mk_commit_object(
            tree_id=assembly_oid,
            parents=parents,
            author=author,
            message=message,
            tags={"source": Path(inp).name, "chunk_symbols": chunk_symbols, "mode": mode, "gl1enc": gl1enc.GL1ENC_FORMAT},  # type: ignore
        )
        self.refs.set_head(commit_oid)
        return commit_oid

    # -------------------------------------------------------------------------
    # Commit/tree accessors
    # -------------------------------------------------------------------------

    def get_commit(self, commit_id: Optional[str]) -> Tuple[str, Dict[str, Any]]:
        if commit_id is None:
            commit_id = self.refs.get_head()
        if not commit_id:
            raise SystemExit("No commits yet. Mint/import something first.")
        typ, core = self.store.get_object(commit_id)
        if typ != "commit":
            raise SystemExit("Not a commit id")
        return commit_id, core

    def get_assembly(self, assembly_id: str) -> Dict[str, Any]:
        typ, core = self.store.get_object(assembly_id)
        if typ != "assembly":
            raise SystemExit("Not an assembly id")
        return core

    def iter_chromosomes(self, assembly_id: str) -> Iterator[Tuple[str, Dict[str, Any]]]:
        assembly = self.get_assembly(assembly_id)
        for cid in assembly["chromosomes"]:
            typ, core = self.store.get_object(cid)
            if typ != "chrom":
                raise SystemExit("Corrupt assembly: non-chrom id")
            yield cid, core

    def find_chrom_by_name(self, assembly_id: str, name: str) -> Tuple[str, Dict[str, Any]]:
        for cid, core in self.iter_chromosomes(assembly_id):
            if core["name"] == name:
                return cid, core
        raise SystemExit(f"Chromosome/record not found: {name}")

    # -------------------------------------------------------------------------
    # Decode / View
    # -------------------------------------------------------------------------

    def view(self, commit_id: Optional[str], chrom_name: str, start: int, end: int, preserve_case: bool) -> str:
        if gl1enc is None:
            raise SystemExit(f"gl1enc_v2.py is required. Import error: {_GL1ENC_IMPORT_ERROR}")

        _, ccore = self.get_commit(commit_id)
        assembly_id = ccore["tree"]
        _chrom_id, chrom = self.find_chrom_by_name(assembly_id=assembly_id, name=chrom_name)

        length = int(chrom["length"])
        start = max(0, start)
        end = min(length, end)
        if end <= start:
            return ""

        cs = int(chrom["chunk_symbols"]) or length
        c0 = start // cs
        c1 = (end - 1) // cs

        out = bytearray()
        for ci in range(c0, c1 + 1):
            chunk_oid = chrom["chunks"][ci]
            _, chunk = self.store.get_object(chunk_oid)
            blob = self.store.get_blob(chunk["chunk_blob"])

            sym, blen = struct.unpack(CHUNK_HDR_FMT, blob[:CHUNK_HDR_SIZE])
            payload = blob[CHUNK_HDR_SIZE:CHUNK_HDR_SIZE + blen]

            chunk_start = ci * cs
            chunk_end = chunk_start + sym
            a = max(start, chunk_start)
            b = min(end, chunk_end)
            if b <= a:
                continue

            seq = _decode_chunk_payload(
                kind=int(chrom["kind"]),
                codec=int(chrom["codec"]),
                payload=payload,
                sym=sym,
                excp_rel=chunk.get("excp", []),
                mask_rel=chunk.get("mask", []),
                preserve_case=preserve_case,
            )

            out.extend(seq[(a - chunk_start):(b - chunk_start)])

        return out.decode("ascii")

    # -------------------------------------------------------------------------
    # Export to GL1 files
    # -------------------------------------------------------------------------

    def export_gl1(self, commit_id: Optional[str], out_prefix: str) -> Tuple[str, Optional[str]]:
        """
        Export commit tree to:
          - <out_prefix>.gl1
          - <out_prefix>.gl1x
          - optionally <out_prefix>.gl1q if qualities exist

        Returns (gl1_path, gl1q_path_or_none)
        """
        if gl1enc is None:
            raise SystemExit(f"gl1enc_v2.py is required. Import error: {_GL1ENC_IMPORT_ERROR}")

        commit_id, ccore = self.get_commit(commit_id)
        assembly_id = ccore["tree"]
        chroms = list(self.iter_chromosomes(assembly_id))

        gl1_path = out_prefix + ".gl1"
        x_path = gl1x_path(gl1_path)
        q_path = gl1q_path(gl1_path)

        xrecs: List[XRec] = []
        q_present_any = False

        with open(gl1_path, "wb") as out:
            out.write(struct.pack(FILE_HDR_FMT, MAGIC_FILE, FILE_VERSION, 0, len(chroms), 0))

            # prepare GL1Q if needed
            fq: Optional[io.BufferedRandom] = open(q_path, "w+b")
            if fq is not None:
                write_gl1q_init(fq)
            qtoc: List[QTocEntry] = []
            q_record_count = 0

            try:
                for _chrom_id, chrom in chroms:
                    name = chrom["name"]
                    nb = name.encode("utf-8")
                    record_off = out.tell()

                    chunk_symbols = int(chrom["chunk_symbols"])
                    chunk_count = len(chrom["chunks"])
                    flags = int(chrom["flags"])

                    # detect qualities exist
                    has_q = False
                    for ch_oid in chrom["chunks"][:1]:
                        _, ch = self.store.get_object(ch_oid)
                        if ch.get("qual_chunk"):
                            has_q = True
                            break
                    if has_q:
                        q_present_any = True

                    out.write(struct.pack(
                        REC_HDR_FMT, MAGIC_REC, REC_VERSION,
                        int(chrom["kind"]), int(chrom["codec"]), flags,
                        int(chrom["length"]),
                        int(chunk_symbols),
                        len(nb),
                        0,
                        int(chunk_count),
                    ))
                    out.write(nb)

                    # GL1Q record start
                    if has_q and fq is not None:
                        qtoc.append(QTocEntry(name, fq.tell()))
                        q_record_count += 1
                        qnb = name.encode("utf-8")
                        fq.write(struct.pack("<H", len(qnb)))
                        fq.write(qnb)
                        fq.write(struct.pack(Q_REC_META_FMT, int(chrom["length"]), int(chunk_symbols), int(chunk_count)))

                    # write chunks
                    chunk_offs: List[int] = []
                    chunk_hashes: List[bytes] = []
                    mask_abs: List[Tuple[int, int]] = []
                    excp_abs: List[Tuple[int, int]] = []

                    cs = chunk_symbols or int(chrom["length"])

                    for ci, chunk_oid in enumerate(chrom["chunks"]):
                        _, chunk = self.store.get_object(chunk_oid)
                        blob = self.store.get_blob(chunk["chunk_blob"])

                        # record chunk offset (points to CHUNK_HDR)
                        chunk_offs.append(out.tell())

                        # write blob bytes directly
                        out.write(blob)

                        # parse payload hash leaf
                        sym, blen = struct.unpack(CHUNK_HDR_FMT, blob[:CHUNK_HDR_SIZE])
                        payload = blob[CHUNK_HDR_SIZE:CHUNK_HDR_SIZE + blen]
                        chunk_hashes.append(sha256_bytes(payload))

                        # flatten mask/excp to absolute coords (for GL1X)
                        chunk_start = ci * cs
                        if chunk.get("mask"):
                            for s, ln in chunk["mask"]:
                                mask_abs.append((chunk_start + int(s), int(ln)))
                        if chunk.get("excp"):
                            for p, code in chunk["excp"]:
                                excp_abs.append((chunk_start + int(p), int(code)))

                        # write qualities if present
                        if has_q and fq is not None:
                            qbid = chunk.get("qual_chunk")
                            if qbid is None:
                                # fill with '!' if missing: should not happen if built from FASTQ
                                q = b"!" * int(sym)
                                z = zlib.compress(q, 6)
                                fq.write(struct.pack(Q_CHUNK_HDR_FMT, len(q), len(z)))
                                fq.write(z)
                            else:
                                qchunk_bytes = self.store.get_blob(qbid)
                                fq.write(qchunk_bytes)

                    mr_bytes = merkle_root_sha256(chunk_hashes)
                    xrecs.append(XRec(
                        name=name,
                        record_off=record_off,
                        kind=int(chrom["kind"]),
                        codec=int(chrom["codec"]),
                        flags=int(chrom["flags"]),
                        length=int(chrom["length"]),
                        chunk_symbols=int(chunk_symbols),
                        chunk_count=int(chunk_count),
                        chunk_offs=chunk_offs,
                        mask_intervals=sorted(mask_abs),
                        excp=sorted(excp_abs),
                        chunk_hashes=chunk_hashes,
                        merkle_root=mr_bytes,
                    ))

                # finalize GL1Q
                if fq is not None:
                    if q_present_any:
                        finalize_gl1q(fq, qtoc, q_record_count)
                    else:
                        fq.close()
                        os.unlink(q_path)
                        fq = None

            finally:
                try:
                    if fq is not None:
                        fq.close()
                except Exception:
                    pass

        write_gl1x(x_path, xrecs)
        return gl1_path, (q_path if q_present_any else None)

    # -------------------------------------------------------------------------
    # Import from GL1 files
    # -------------------------------------------------------------------------

    def import_gl1(self, gl1_path_: str, assembly_name: str, author: str, message: str) -> str:
        """
        Ingest an existing .gl1 + optional .gl1x and .gl1q into the object graph.

        - If .gl1x exists, uses it for mask/exceptions slicing and chunk offsets.
        - Otherwise imports payload-only (no mask/exceptions).
        - If .gl1q exists, imports qualities.

        Returns: commit id.
        """
        self.init()

        x_path = gl1x_path(gl1_path_)
        q_path = gl1q_path(gl1_path_)

        xrecs: Optional[List[XRec]] = None
        if os.path.exists(x_path):
            xrecs = read_gl1x(x_path)

        qtoc: Optional[Dict[str, int]] = None
        if os.path.exists(q_path):
            qtoc = read_gl1q_toc(q_path)

        chrom_ids: List[str] = []

        with open(gl1_path_, "rb") as f:
            hdr = f.read(FILE_HDR_SIZE)
            if len(hdr) != FILE_HDR_SIZE:
                raise ValueError("Bad GL1F header")
            magic, ver, _flags, record_count, _rsv = struct.unpack(FILE_HDR_FMT, hdr)
            if magic != MAGIC_FILE or ver != FILE_VERSION:
                raise ValueError("Not a supported GL1F file")

            for ri in range(record_count):
                record_off = f.tell()
                rh = f.read(REC_HDR_SIZE)
                if len(rh) != REC_HDR_SIZE:
                    raise ValueError("Truncated record header")
                (rmagic, rver, kind, codec, flags, length, chunk_symbols, name_len, _rsv2, chunk_count) = struct.unpack(REC_HDR_FMT, rh)
                if rmagic != MAGIC_REC or rver != REC_VERSION:
                    raise ValueError("Bad record header magic/version")
                name = f.read(name_len).decode("utf-8", "replace")

                # determine chunk offsets / metadata
                if xrecs is not None:
                    # match by name
                    xr = next((x for x in xrecs if x.name == name), None)
                    if xr is None:
                        raise ValueError(f"Missing GL1X record: {name}")
                    chunk_offs = xr.chunk_offs
                    mask_abs = xr.mask_intervals
                    excp_abs = xr.excp
                    chunk_hashes_bytes = xr.chunk_hashes
                else:
                    # compute offsets by scanning
                    chunk_offs = []
                    pos0 = f.tell()
                    for _ in range(chunk_count):
                        chunk_offs.append(f.tell())
                        chdr = f.read(CHUNK_HDR_SIZE)
                        if len(chdr) != CHUNK_HDR_SIZE:
                            raise ValueError("Truncated chunk header")
                        sym, blen = struct.unpack(CHUNK_HDR_FMT, chdr)
                        payload = f.read(blen)
                        if len(payload) != blen:
                            raise ValueError("Truncated chunk payload")
                    # reset to start of chunks for actual import
                    f.seek(pos0)
                    mask_abs = []
                    excp_abs = []
                    chunk_hashes_bytes = []

                # read qualities for this record if present
                qual_chunks: List[Optional[bytes]] = []
                if qtoc is not None and name in qtoc:
                    qual_chunks = read_gl1q_record_comp_chunks(q_path, qtoc, name)
                    if len(qual_chunks) != chunk_count:
                        raise ValueError("GL1Q chunk count mismatch")

                # build chunks
                payload_hashes_hex: List[str] = []
                mask_leaf_hex: List[str] = []
                excp_leaf_hex: List[str] = []
                qual_leaf_hex: List[str] = []
                chunk_oids: List[str] = []

                cs = int(chunk_symbols) or int(length)

                # pre-sort global lists
                mask_global = sorted(mask_abs)
                excp_global = sorted(excp_abs)

                for ci in range(chunk_count):
                    # ensure file pointer is at the chunk offset (if xrecs exist)
                    if xrecs is not None:
                        f.seek(chunk_offs[ci])
                    sym, blen = struct.unpack(CHUNK_HDR_FMT, f.read(CHUNK_HDR_SIZE))
                    payload = f.read(blen)
                    chunk_blob = struct.pack(CHUNK_HDR_FMT, sym, blen) + payload

                    ph = sha256_hex(payload)
                    payload_hashes_hex.append(ph)

                    chunk_start = ci * cs
                    chunk_end = chunk_start + sym

                    mask_rel = _slice_intervals(mask_global, chunk_start, chunk_end)
                    excp_rel = _slice_points(excp_global, chunk_start, chunk_end)

                    mask_leaf_hex.append(sha256_hex(canonical_json([[a, b] for a, b in mask_rel]).encode("utf-8")))
                    excp_leaf_hex.append(sha256_hex(canonical_json([[p, c] for p, c in excp_rel]).encode("utf-8")))

                    qc_bytes = qual_chunks[ci] if qual_chunks else None
                    qual_leaf_hex.append(sha256_hex(qc_bytes) if qc_bytes is not None else sha256_hex(b""))

                    chunk_oid = self._mk_chunk_object(
                        chunk_blob=chunk_blob,
                        kind=int(kind),
                        codec=int(codec),
                        sym_count=int(sym),
                        payload_hash_hex=ph,
                        mask_rel=mask_rel,
                        excp_rel=excp_rel,
                        qual_chunk_bytes=qc_bytes,
                    )
                    chunk_oids.append(chunk_oid)

                chrom_oid = self._mk_chrom_object(
                    name=name,
                    kind=int(kind),
                    codec=int(codec),
                    flags=int(flags),
                    length=int(length),
                    chunk_symbols=int(chunk_symbols),
                    chunk_ids=chunk_oids,
                    payload_hashes_hex=payload_hashes_hex,
                    mask_leaves_hex=mask_leaf_hex,
                    excp_leaves_hex=excp_leaf_hex,
                    qual_leaves_hex=qual_leaf_hex if qual_chunks else [],
                )
                chrom_ids.append(chrom_oid)

        assembly_oid = self._mk_assembly_object(assembly_name=assembly_name, chrom_ids=chrom_ids)
        parent = self.refs.get_head()
        parents = [parent] if parent else []
        commit_oid = self._mk_commit_object(tree_id=assembly_oid, parents=parents, author=author, message=message,
                                            tags={"source": Path(gl1_path_).name, "import": True})
        self.refs.set_head(commit_oid)
        return commit_oid

    # -------------------------------------------------------------------------
    # Mutation / evolution
    # -------------------------------------------------------------------------

    def mutate_random_snps(
        self,
        commit_id: Optional[str],
        chrom_name: str,
        n_snps: int,
        seed: int,
        author: str,
        message: str,
        preserve_case: bool,
    ) -> str:
        """
        Create a new commit by applying random A/C/G/T substitutions on a single chromosome.
        Only affected chunks are rewritten; unchanged chunks are reused (git-like efficiency).

        Notes:
        - This mutation op is defined for nucleotide records. If the record is protein, it errors.
        """
        if gl1enc is None:
            raise SystemExit(f"gl1enc_v2.py is required. Import error: {_GL1ENC_IMPORT_ERROR}")
        random.seed(seed)

        commit_id, ccore = self.get_commit(commit_id)
        assembly_id = ccore["tree"]
        assembly = self.get_assembly(assembly_id)

        chrom_ids = list(assembly["chromosomes"])
        chrom_id, chrom = self.find_chrom_by_name(assembly_id, chrom_name)

        if int(chrom["kind"]) != KIND_NUC:
            raise SystemExit("mutate_random_snps only supports nucleotide records in this reference implementation.")

        length = int(chrom["length"])
        if length <= 0:
            raise SystemExit("Chromosome length is 0")

        cs = int(chrom["chunk_symbols"]) or length

        # pick SNP positions
        positions = sorted({random.randrange(0, length) for _ in range(n_snps)})
        if not positions:
            raise SystemExit("No SNP positions selected")

        # map chunk -> list of positions within that chunk
        edits: Dict[int, List[int]] = {}
        for p in positions:
            ci = p // cs
            edits.setdefault(ci, []).append(p)

        new_chunk_ids = list(chrom["chunks"])

        for ci, pos_list in edits.items():
            chunk_oid = chrom["chunks"][ci]
            _, chunk = self.store.get_object(chunk_oid)
            blob = self.store.get_blob(chunk["chunk_blob"])
            sym, blen = struct.unpack(CHUNK_HDR_FMT, blob[:CHUNK_HDR_SIZE])
            payload = blob[CHUNK_HDR_SIZE:CHUNK_HDR_SIZE + blen]

            # decode semantic sequence for this chunk (apply excp and mask optionally)
            seq = _decode_chunk_payload(
                kind=int(chrom["kind"]),
                codec=int(chrom["codec"]),
                payload=payload,
                sym=sym,
                excp_rel=chunk.get("excp", []),
                mask_rel=chunk.get("mask", []),
                preserve_case=preserve_case,
            )

            chunk_start = ci * cs
            # mutate bases
            for abs_p in pos_list:
                rel = abs_p - chunk_start
                if rel < 0 or rel >= len(seq):
                    continue
                old = chr(seq[rel])
                # pick new base different from old (case-insensitive)
                old_u = old.upper()
                if old_u not in "ACGT":
                    old_u = "A"
                newb = random.choice([b for b in "ACGT" if b != old_u])
                if preserve_case and ("a" <= old <= "z"):
                    newb = newb.lower()
                seq[rel] = ord(newb)

            # re-encode chunk using the SAME record codec
            # (payload always stored uppercase; mask carries lowercase info)
            seq_str_case = seq.decode("ascii")
            # recompute mask (if preserve_case) from seq_str_case before uppercasing
            mask_rel_new = _mask_intervals_from_seq(seq_str_case.encode("ascii"))

            # payload uses uppercase
            seq_norm = seq_str_case.upper()

            codec = int(chrom["codec"])
            use_excp = bool(int(chrom["flags"]) & FLAG_HAS_EXCP)

            excp_rel_new: List[Tuple[int, int]] = []
            payload2: bytes

            if codec == CODEC_DNA2:
                base_chars: List[str] = []
                for j, ch in enumerate(seq_norm):
                    if ch in "ACGT":
                        base_chars.append(ch)
                    else:
                        base_chars.append("A")
                        if use_excp:
                            if ch in gl1enc._DNA4_CHAR_TO_NIBBLE:  # type: ignore
                                code = int(gl1enc._DNA4_CHAR_TO_NIBBLE[ch])  # type: ignore
                            else:
                                code = int(gl1enc._DNA4_CHAR_TO_NIBBLE["N"])  # type: ignore
                            excp_rel_new.append((j, code))
                payload2 = gl1enc.pack_dna2("".join(base_chars))  # type: ignore

            elif codec == CODEC_DNA4:
                if not gl1enc.is_dna4_eligible(seq_norm, assume_normalized=True):  # type: ignore
                    # promote to ASCII
                    codec = CODEC_ASCII
                    payload2 = seq_norm.encode("ascii", "replace")
                else:
                    payload2 = gl1enc.pack_dna4(seq_norm)  # type: ignore

            elif codec == CODEC_ASCII:
                payload2 = seq_norm.encode("ascii", "replace")

            else:
                # shouldn't happen for nucleotides
                codec = CODEC_ASCII
                payload2 = seq_norm.encode("ascii", "replace")

            new_blob = struct.pack(CHUNK_HDR_FMT, int(sym), len(payload2)) + payload2
            ph = sha256_hex(payload2)

            new_chunk_oid = self._mk_chunk_object(
                chunk_blob=new_blob,
                kind=int(chrom["kind"]),
                codec=codec,
                sym_count=int(sym),
                payload_hash_hex=ph,
                mask_rel=mask_rel_new,
                excp_rel=excp_rel_new,
                qual_chunk_bytes=self.store.get_blob(chunk["qual_chunk"]) if chunk.get("qual_chunk") else None,
            )
            new_chunk_ids[ci] = new_chunk_oid

        # rebuild chrom object (roots will change)
        payload_hashes_hex: List[str] = []
        mask_leaves_hex: List[str] = []
        excp_leaves_hex: List[str] = []
        qual_leaves_hex: List[str] = []

        for ch_oid in new_chunk_ids:
            _, ch = self.store.get_object(ch_oid)
            # parse payload hash from stored field
            payload_hashes_hex.append(ch["payload_sha256"])
            mask_leaves_hex.append(sha256_hex(canonical_json(ch.get("mask", [])).encode("utf-8")))
            excp_leaves_hex.append(sha256_hex(canonical_json(ch.get("excp", [])).encode("utf-8")))
            if ch.get("qual_chunk"):
                qual_leaves_hex.append(sha256_hex(self.store.get_blob(ch["qual_chunk"])))

        new_chrom_oid = self._mk_chrom_object(
            name=chrom["name"],
            kind=int(chrom["kind"]),
            codec=int(chrom["codec"]),
            flags=int(chrom["flags"]),
            length=int(chrom["length"]),
            chunk_symbols=int(chrom["chunk_symbols"]),
            chunk_ids=new_chunk_ids,
            payload_hashes_hex=payload_hashes_hex,
            mask_leaves_hex=mask_leaves_hex,
            excp_leaves_hex=excp_leaves_hex,
            qual_leaves_hex=qual_leaves_hex if qual_leaves_hex else [],
        )

        # new assembly: replace just this chrom id
        new_chroms = [new_chrom_oid if cid == chrom_id else cid for cid in chrom_ids]
        new_assembly_oid = self._mk_assembly_object(assembly_name=assembly["assembly"], chrom_ids=new_chroms)

        new_commit_oid = self._mk_commit_object(
            tree_id=new_assembly_oid,
            parents=[commit_id],
            author=author,
            message=message,
            tags={"mutation": "random_snps", "chrom": chrom_name, "n": n_snps, "seed": seed},
        )
        self.refs.set_head(new_commit_oid)
        return new_commit_oid

    # -------------------------------------------------------------------------
    # Deploy (chain simulation)
    # -------------------------------------------------------------------------

    def deploy(self, commit_id: Optional[str], mode: str, out_json: str) -> str:
        """
        Create a chain-sim deployment:
          - mode=onchain: copy all referenced blobs to chain_sim/onchain_blobs/
          - mode=ipfs: copy all referenced blobs to chain_sim/ipfs_blobs/ and record cid
        Mint token JSON files for assembly/chrom/chunk objects.
        """
        commit_id, ccore = self.get_commit(commit_id)
        assembly_id = ccore["tree"]

        # mint commit token
        self.chain.mint_token(commit_id, {"type": "CommitNFT", "commit": commit_id, "tree": assembly_id, "parents": ccore["parents"], "message": ccore["message"]})

        assembly = self.get_assembly(assembly_id)
        self.chain.mint_token(assembly_id, {"type": "AssemblyNFT", "assembly_id": assembly_id, "assembly": assembly["assembly"], "chromosomes": assembly["chromosomes"], "assembly_root": assembly["assembly_root"]})

        deployed_blobs: Dict[str, str] = {}  # blob_id -> uri
        deployed_chunks: List[str] = []
        deployed_chroms: List[str] = []

        for chrom_id, chrom in self.iter_chromosomes(assembly_id):
            deployed_chroms.append(chrom_id)
            self.chain.mint_token(chrom_id, {"type": "ChromosomeNFT", "chrom_id": chrom_id, "name": chrom["name"], "length": chrom["length"], "chunk_symbols": chrom["chunk_symbols"], "chunks": chrom["chunks"],
                                             "seq_payload_root": chrom["seq_payload_root"], "mask_root": chrom["mask_root"], "excp_root": chrom["excp_root"], "qual_root": chrom["qual_root"]})

            for chunk_id in chrom["chunks"]:
                deployed_chunks.append(chunk_id)
                _, ch = self.store.get_object(chunk_id)
                # deploy chunk blob
                bid = ch["chunk_blob"]
                if bid not in deployed_blobs:
                    deployed_blobs[bid] = self.chain.put_onchain_blob(bid) if mode == "onchain" else self.chain.put_ipfs_blob(bid)
                # deploy quality blob if present
                if ch.get("qual_chunk"):
                    qbid = ch["qual_chunk"]
                    if qbid not in deployed_blobs:
                        deployed_blobs[qbid] = self.chain.put_onchain_blob(qbid) if mode == "onchain" else self.chain.put_ipfs_blob(qbid)

                self.chain.mint_token(chunk_id, {"type": "ChunkNFT", "chunk_id": chunk_id, "kind": ch["kind"], "codec": ch["codec"], "sym_count": ch["sym_count"],
                                                 "chunk_blob": {"blob_id": bid, "uri": deployed_blobs[bid], "payload_sha256": ch["payload_sha256"]},
                                                 "mask": ch.get("mask", []), "excp": ch.get("excp", []), "qual_chunk": ({"blob_id": ch["qual_chunk"], "uri": deployed_blobs[ch["qual_chunk"]]} if ch.get("qual_chunk") else None)})

        deployment = {
            "deployment_v": 2,
            "mode": mode,
            "commit": commit_id,
            "assembly": assembly_id,
            "chromosomes": deployed_chroms,
            "chunks": deployed_chunks,
            "blobs": deployed_blobs,
        }
        if out_json == "-":
            return json.dumps(deployment, indent=2)
        Path(out_json).write_text(json.dumps(deployment, indent=2) + "\n", encoding="utf-8")
        return out_json

    # -------------------------------------------------------------------------
    # Minimizer index (GL1I)
    # -------------------------------------------------------------------------

    def build_minimizer_index(self, commit_id: Optional[str], k: int, w: int, canonical: bool, bucket_count: int) -> str:
        """
        Build GL1I minimizer index for the assembly *directly* from stored chunks
        (no temp export needed).

        For nucleotide records:
          - only kmers with A/C/G/T are indexed (ambiguous symbols are skipped)

        Stores the resulting .gl1i bytes as a blob and returns an Index object id.
        """
        if gl1enc is None:
            raise SystemExit(f"gl1enc_v2.py is required. Import error: {_GL1ENC_IMPORT_ERROR}")
        if k <= 0 or w <= 0:
            raise ValueError("k and w must be positive")
        if bucket_count <= 0:
            raise ValueError("bucket_count must be positive")
        if k > 31:
            # 2-bit packing into uint64 supports up to 31 bases (62 bits)
            raise ValueError("k too large for nuc2 uint64 key (max 31 in this reference)")

        commit_id, ccore = self.get_commit(commit_id)
        assembly_id = ccore["tree"]

        # Build records list and posting buckets
        records_meta: List[Tuple[str, int]] = []
        # buckets[bucket] -> dict(key -> postings[(rec_id,pos)])
        buckets: List[Dict[int, List[Tuple[int, int]]]] = [dict() for _ in range(bucket_count)]

        flags = I_FLAG_CANON if canonical else 0
        flags |= I_FLAG_NUC

        rec_id = 0
        for _cid, chrom in self.iter_chromosomes(assembly_id):
            name = chrom["name"]
            kind = int(chrom["kind"])
            length = int(chrom["length"])
            records_meta.append((name, length))

            if kind != KIND_NUC:
                rec_id += 1
                continue

            # stream decode full record in uppercase (ignore preserve_case)
            seq = self.view(commit_id, name, 0, length, preserve_case=False)
            if len(seq) != length:
                # tolerate but warn by using decoded length
                length = len(seq)

            # compute kmer keys for each position; None if invalid
            keys: List[Optional[int]] = [None] * max(0, length - k + 1)
            for i in range(0, length - k + 1):
                key = _kmer_key_nuc2(seq[i:i+k], canonical)
                keys[i] = key

            # minimizer over window w (over k-mers)
            last_min: Optional[Tuple[int, int]] = None  # (key, pos)
            for i in range(0, len(keys) - w + 1):
                # find min in window
                min_key: Optional[int] = None
                min_pos: Optional[int] = None
                for j in range(i, i + w):
                    kj = keys[j]
                    if kj is None:
                        continue
                    if min_key is None or kj < min_key:
                        min_key = kj
                        min_pos = j
                if min_key is None or min_pos is None:
                    continue
                cur = (min_key, min_pos)
                if last_min == cur:
                    continue
                last_min = cur
                b = min_key % bucket_count
                buckets[b].setdefault(min_key, []).append((rec_id, int(min_pos)))

            rec_id += 1

        # write to temp gl1i, then store as blob
        tmpdir = Path(tempfile.mkdtemp(prefix="seqnft_gl1i_"))
        try:
            out_i = str(tmpdir / "index.gl1i")
            write_gl1i(out_i, k=k, w=w, canonical=canonical, bucket_count=bucket_count, records=records_meta, buckets=buckets, flags=flags)
            data = Path(out_i).read_bytes()
            bid = self.store.put_blob(data)
            core = {
                "v": 2,
                "kind": "minimizer_gl1i",
                "source_commit": commit_id,
                "k": int(k),
                "w": int(w),
                "canonical": bool(canonical),
                "bucket_count": int(bucket_count),
                "gl1i_blob": bid,
                "sha256": sha256_hex(data),
            }
            return self.store.put_object("index", core)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def query_minimizer_index(self, index_oid: str, kmer: Optional[str], key_hex: Optional[str]) -> List[str]:
        """
        Query stored GL1I index and return lines of "record\\tpos".
        """
        typ, core = self.store.get_object(index_oid)
        if typ != "index" or core.get("kind") != "minimizer_gl1i":
            raise SystemExit("Not a minimizer index object")
        bid = core["gl1i_blob"]
        data = self.store.get_blob(bid)

        tmp = tempfile.NamedTemporaryFile("wb", delete=False, suffix=".gl1i")
        try:
            tmp.write(data)
            tmp.flush()
            tmp.close()

            k, w, flags, bucket_count, recs, bucket_dir_off = read_gl1i_header(tmp.name)
            canonical = bool(flags & I_FLAG_CANON)

            if kmer is not None:
                if bool(flags & I_FLAG_NUC):
                    key = _kmer_key_nuc2(kmer, canonical)
                else:
                    key = _kmer_key_sixbit(kmer, canonical)
                if key is None:
                    return []
            elif key_hex is not None:
                key = int(key_hex, 16)
            else:
                raise SystemExit("Provide --kmer or --key")

            hits = query_gl1i(tmp.name, int(key))
            out: List[str] = []
            for rec_id, pos in hits:
                name = recs[rec_id][0] if 0 <= rec_id < len(recs) else f"rec{rec_id}"
                out.append(f"{name}\t{pos}")
            return out
        finally:
            try:
                os.unlink(tmp.name)
            except OSError:
                pass


# =============================================================================
# Helpers (slicing, decoding, mask/exceptions)
# =============================================================================

def _slice_intervals(global_intervals: List[Tuple[int, int]], abs_start: int, abs_end: int) -> List[Tuple[int, int]]:
    """
    Slice global (start,len) intervals into [abs_start,abs_end), returning relative (start,len).
    """
    out: List[Tuple[int, int]] = []
    for s, ln in global_intervals:
        a = max(int(s), abs_start)
        b = min(int(s) + int(ln), abs_end)
        if b <= a:
            continue
        out.append((a - abs_start, b - a))
    return out


def _slice_points(global_points: List[Tuple[int, int]], abs_start: int, abs_end: int) -> List[Tuple[int, int]]:
    """
    Slice global (pos,code) points into [abs_start,abs_end), returning relative (pos,code).
    """
    out: List[Tuple[int, int]] = []
    for p, code in global_points:
        p = int(p)
        if abs_start <= p < abs_end:
            out.append((p - abs_start, int(code)))
    return out


def _mask_intervals_from_seq(seq_bytes: bytes) -> List[Tuple[int, int]]:
    """
    Recompute lowercase intervals within a chunk from raw (possibly mixed-case) bytes.
    """
    out: List[Tuple[int, int]] = []
    run: Optional[int] = None
    for i, b in enumerate(seq_bytes):
        is_low = 97 <= b <= 122
        if is_low and run is None:
            run = i
        if (not is_low) and run is not None:
            out.append((run, i - run))
            run = None
    if run is not None:
        out.append((run, len(seq_bytes) - run))
    return out


def _apply_mask_inplace(seq: bytearray, mask_rel: List[Sequence[int]]) -> None:
    """
    Apply lowercase mask intervals (relative) to seq in-place.
    """
    for s, ln in mask_rel:
        s = int(s); ln = int(ln)
        if ln <= 0:
            continue
        end = min(len(seq), s + ln)
        for i in range(max(0, s), end):
            b = seq[i]
            if 65 <= b <= 90:  # A-Z
                seq[i] = b + 32  # to lowercase


def _apply_excp_inplace(seq: bytearray, excp_rel: List[Sequence[int]]) -> None:
    """
    Apply exceptions (relative) to seq in-place. Exception code is a DNA4 nibble (0..15).
    """
    if gl1enc is None:
        raise RuntimeError("gl1enc not available")
    for p, code in excp_rel:
        p = int(p); code = int(code) & 0xFF
        if 0 <= p < len(seq):
            # map nibble -> canonical IUPAC/gap char
            ch = gl1enc._DNA4_NIBBLE_TO_CHAR.get(code, "N")  # type: ignore
            seq[p] = ord(ch)


def _decode_chunk_payload(
    *,
    kind: int,
    codec: int,
    payload: bytes,
    sym: int,
    excp_rel: List[Sequence[int]],
    mask_rel: List[Sequence[int]],
    preserve_case: bool,
) -> bytearray:
    """
    Decode payload bytes into uppercase ASCII letters, then optionally apply
    exceptions and mask to reconstruct semantic sequence.

    - For codec DNA2: payload stores packed A/C/G/T only. If excp_rel is non-empty,
      exceptions are applied to restore ambiguous symbols.
    - For codec DNA4: payload stores packed IUPAC nucleotides (and '-' gap).
    - For codec SIXBIT: payload stores packed 64-symbol alphabet.
    - For codec ASCII: payload stores raw uppercase ASCII bytes.

    Returns: bytearray of length `sym`.
    """
    if gl1enc is None:
        raise RuntimeError(f"gl1enc_v2 not available: {_GL1ENC_IMPORT_ERROR}")

    if codec == CODEC_DNA2:
        s = gl1enc.unpack_dna2(payload, int(sym))  # type: ignore
        seq = bytearray(s.encode("ascii"))
        if excp_rel:
            _apply_excp_inplace(seq, excp_rel)
    elif codec == CODEC_DNA4:
        s = gl1enc.unpack_dna4(payload, int(sym))  # type: ignore
        seq = bytearray(s.encode("ascii"))
    elif codec == CODEC_SIXBIT:
        s = gl1enc.unpack_sixbit(payload, int(sym))  # type: ignore
        seq = bytearray(s.encode("ascii"))
    else:
        # ASCII fallback
        if len(payload) != int(sym):
            # allow but truncate/pad
            s = payload.decode("ascii", "replace")
            s = (s + ("N" * int(sym)))[: int(sym)]
            seq = bytearray(s.encode("ascii"))
        else:
            seq = bytearray(payload)

    if preserve_case and mask_rel:
        _apply_mask_inplace(seq, mask_rel)
    return seq


# =============================================================================
# CLI
# =============================================================================

def cmd_init(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo(Path(args.repo))
    repo.init()
    print(f"Initialized SeqNFT repo at {repo.paths.dot}")


def cmd_mint(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    commit = repo.mint_from_fasta_fastq(
        inp=args.inp,
        assembly_name=args.assembly,
        chunk_symbols=args.chunk,
        mode=args.mode,
        strict=args.strict,
        unknown_to_x=args.unknown_to_x,
        store_qual=args.store_qual,
        nuc2_excp_threshold=args.nuc2_excp_threshold,
        author=args.author,
        message=args.message,
    )
    print(commit)


def cmd_import(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    commit = repo.import_gl1(args.gl1, assembly_name=args.assembly, author=args.author, message=args.message)
    print(commit)


def cmd_log(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    head = args.commit or repo.refs.get_head()
    if not head:
        print("(no commits)")
        return
    cur = head
    n = 0
    while cur and n < args.n:
        typ, c = repo.store.get_object(cur)
        if typ != "commit":
            break
        ts = _dt.datetime.fromtimestamp(int(c["time"])).isoformat()
        print(f"commit {cur}")
        print(f"Author: {c['author']}")
        print(f"Date:   {ts}")
        print(f"\n    {c['message']}\n")
        parents = c.get("parents") or []
        cur = parents[0] if parents else None
        n += 1


def cmd_view(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    s = repo.view(args.commit, args.chrom, args.start, args.end, preserve_case=args.preserve_case)
    if args.to == "raw":
        print(s, end="")
    else:
        print(f">{args.chrom}:{args.start}-{args.end}")
        w = args.wrap
        for i in range(0, len(s), w):
            print(s[i:i+w])


def cmd_export(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    gl1_path_out, q_path_out = repo.export_gl1(args.commit, args.out_prefix)
    print(gl1_path_out)
    print(gl1x_path(gl1_path_out))
    if q_path_out:
        print(q_path_out)


def cmd_mutate(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    newc = repo.mutate_random_snps(
        commit_id=args.commit,
        chrom_name=args.chrom,
        n_snps=args.n,
        seed=args.seed,
        author=args.author,
        message=args.message,
        preserve_case=args.preserve_case,
    )
    print(newc)


def cmd_deploy(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    out = repo.deploy(args.commit, mode=args.mode, out_json=args.out)
    print(out)


def cmd_index_build(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    oid = repo.build_minimizer_index(args.commit, k=args.k, w=args.w, canonical=(not args.no_canonical), bucket_count=args.bucket_count)
    print(oid)


def cmd_index_query(args: argparse.Namespace) -> None:
    repo = SeqNFTRepo.discover(Path(args.repo) if args.repo else Path.cwd())
    lines = repo.query_minimizer_index(args.index, kmer=args.kmer, key_hex=args.key)
    for ln in lines:
        print(ln)


def build_cli() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="seqnft", description="SeqNFT+GL1 pipeline (local chain simulator) — GL1ENCv2")
    ap.add_argument("--repo", default=None, help="Repo root (default: discover from cwd). For init only: where to create repo.")

    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("init", help="Initialize a SeqNFT repo (.seqnft).")
    p.add_argument("--repo", default=".", help="Repo directory to initialize (default: current directory).")
    p.set_defaults(fn=cmd_init)

    p = sub.add_parser("mint", help="Mint a new commit from FASTA/FASTQ using GL1ENCv2 encodings.")
    p.add_argument("-i", "--in", dest="inp", required=True)
    p.add_argument("--assembly", default="assembly")
    p.add_argument("--chunk", type=int, default=1_000_000)
    p.add_argument("--mode", choices=["auto", "2bit", "4bit", "6bit"], default="auto",
                   help="Encoding mode. auto chooses best; 2bit forces DNA2(+excp); 4bit forces DNA4; 6bit forces SIXBIT for proteins.")
    p.add_argument("--strict", action="store_true", help="Fail if chosen mode cannot represent input symbols.")
    p.add_argument("--unknown-to-x", action="store_true", help="Replace unknown symbols with 'N' (nuc) or 'X' (prot) to keep compact encoding.")
    p.add_argument("--store-qual", action="store_true", help="If input is FASTQ, store qualities in GL1Q sidecar blobs.")
    p.add_argument("--nuc2-excp-threshold", type=float, default=0.001,
                   help="In auto mode (nuc): if fraction of non-ACGT <= threshold, store DNA2 payload + sparse exceptions, else use DNA4.")
    p.add_argument("--author", default="anon")
    p.add_argument("-m", "--message", default="mint")
    p.set_defaults(fn=cmd_mint)

    p = sub.add_parser("import", help="Import an existing .gl1 (+ optional .gl1x/.gl1q) into the repo as a new commit.")
    p.add_argument("--gl1", required=True, help="Path to .gl1 file")
    p.add_argument("--assembly", default="assembly")
    p.add_argument("--author", default="anon")
    p.add_argument("-m", "--message", default="import")
    p.set_defaults(fn=cmd_import)

    p = sub.add_parser("log", help="Show commit history (linear parents[0]).")
    p.add_argument("--commit", default=None, help="Start commit (default HEAD).")
    p.add_argument("-n", type=int, default=20)
    p.set_defaults(fn=cmd_log)

    p = sub.add_parser("view", help="Extract a region from a chromosome in a commit.")
    p.add_argument("--commit", default=None)
    p.add_argument("--chrom", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--preserve-case", action="store_true")
    p.add_argument("--to", choices=["fasta", "raw"], default="fasta")
    p.add_argument("--wrap", type=int, default=60)
    p.set_defaults(fn=cmd_view)

    p = sub.add_parser("export", help="Export commit to .gl1/.gl1x/(.gl1q).")
    p.add_argument("--commit", default=None)
    p.add_argument("-o", "--out-prefix", required=True, help="Output prefix (writes <prefix>.gl1 etc).")
    p.set_defaults(fn=cmd_export)

    p = sub.add_parser("mutate", help="Create a new commit by applying random SNPs to a nucleotide chromosome.")
    p.add_argument("--commit", default=None)
    p.add_argument("--chrom", required=True)
    p.add_argument("-n", type=int, default=1000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--author", default="anon")
    p.add_argument("-m", "--message", default="mutate")
    p.add_argument("--preserve-case", action="store_true")
    p.set_defaults(fn=cmd_mutate)

    p = sub.add_parser("deploy", help="Deploy a commit into the chain simulator (mint token JSON + store blobs).")
    p.add_argument("--commit", default=None)
    p.add_argument("--mode", choices=["onchain", "ipfs"], default="onchain")
    p.add_argument("-o", "--out", default="-", help="deployment.json output path or '-' for stdout")
    p.set_defaults(fn=cmd_deploy)

    p = sub.add_parser("index-build", help="Build minimizer index GL1I for commit (no export needed).")
    p.add_argument("--commit", default=None)
    p.add_argument("--k", type=int, default=15)
    p.add_argument("--w", type=int, default=10)
    p.add_argument("--no-canonical", action="store_true")
    p.add_argument("--bucket-count", type=int, default=256)
    p.set_defaults(fn=cmd_index_build)

    p = sub.add_parser("index-query", help="Query a stored minimizer GL1I index.")
    p.add_argument("--index", required=True, help="Index object id returned by index-build")
    p.add_argument("--kmer", default=None)
    p.add_argument("--key", default=None, help="hex key like 0x1234... or deadbeef...")
    p.set_defaults(fn=cmd_index_query)

    return ap


def main(argv: Optional[List[str]] = None) -> int:
    ap = build_cli()
    args = ap.parse_args(argv)
    args.fn(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
