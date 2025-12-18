# SeqNFT + GL1 Unified Pipeline Specification (v1)

This spec defines a **single pipeline** that runs with identical semantics in:
- **off-chain mode** (local files, IPFS/object storage, cloud),
- **on-chain mode** (smart contracts storing bytes or anchoring CIDs).

It is built on top of the existing **GL1** format family:

- `.gl1` (GL1F v1) payload container (bitpacked chunks)
- `.gl1x` (GL1X v2) sidecar: random access + masking + exceptions + chunk hashes + Merkle root
- `.gl1q` (GL1Q v2) sidecar: FASTQ qualities with TOC
- `.gl1i` (GL1I v1) minimizer postings index

The key additional concept is **SeqNFT objects + commits**, which make GL1 **history-aware** (git-like) and **NFT-native**.

---

## 1. Goals

1) **Identical semantics on-chain and off-chain**  
   Same hashes, same roots, same record/chunk layout, same region extraction output.

2) **Future-proof “git-like” history for sequences**  
   Commits form a DAG; new commits reuse unchanged chunks; only modified chunks are new objects.

3) **NFT-native hierarchy**  
   - Chromosome (record) = parent NFT  
   - Chunk = child NFT  
   - Assembly/Commit = root NFTs (versioned dataset)

4) Two deployment modes:
   - **Fully on-chain**: chunk bytes are stored directly in contracts.
   - **Anchor/IPFS**: contracts store only `CID + sha256`, clients fetch bytes from IPFS and verify.

5) Ready for:
   - scalable region queries (“view”)
   - scalable search (“minimizer index”)
   - deterministic derivations (gene→protein, variants→new commits)
   - alignment-ready expansion (index + block format hooks)

---

## 2. Object Model (Content Addressed)

All objects are immutable and addressed by `sha256` of a **canonical encoding** of the object core.

### 2.1 Blob
Raw bytes addressed by `sha256(bytes)`.
Used for:
- GL1 chunk blobs (CHUNK_HDR + payload)
- compressed quality-chunk bytes (GL1Q chunk hdr + zlib bytes)
- `.gl1i` bytes (minimizer index)

### 2.2 Chunk object (`type="chunk"`)
Represents a **single GL1 chunk** plus chunk-local meta.

Core fields (conceptual):
- `chunk_blob`: blob id holding exact GL1 bytes `CHUNK_HDR + payload`
- `payload_sha256`: sha256(payload) for GL1X-compatible Merkle leaves
- `mask`: list of `(start,len)` **within chunk**
- `excp`: list of `(pos, code)` **within chunk** (NUC4 NT16 code)
- `qual_chunk` (optional): blob id holding the compressed quality chunk bytes

This preserves core GL1 features:
- masking
- exceptions
- chunk hashing and Merkle

### 2.3 Chromosome object (`type="chrom"`)
A record manifest = **parent NFT**.

Core fields:
- `name`, `kind`, `codec`, `flags`, `length`, `chunk_symbols`
- `chunks`: ordered list of chunk object ids
- roots:
  - `seq_payload_root`: Merkle root over `sha256(payload)` leaves (GL1X-compatible)
  - `mask_root`, `excp_root`, `qual_root`: Merkle roots over per-chunk meta leaves (pipeline-level integrity)

### 2.4 Assembly object (`type="assembly"`)
Dataset manifest of chromosomes.

Core fields:
- `assembly`: name/version string
- `chromosomes`: list of chromosome ids
- `assembly_root`: Merkle root over chromosome ids (pipeline-level)

### 2.5 Commit object (`type="commit"`)
Git-like commit node; can be minted as a “Version NFT”.

Core fields:
- `tree`: assembly id
- `parents`: list of parent commit ids
- `author`, `message`, `time`
- `tags` (optional): structured metadata (aligner version, pipeline params, etc.)

---

## 3. Canonical Hashing

### 3.1 Canonical JSON
For v1, canonical object encoding is JSON with:
- sorted keys
- separators `(",", ":")`
- UTF-8 bytes

Then:
```
object_id = sha256( b"seqnft:" + type + b"\\0" + canonical_json(core) )
```

This makes hashing consistent off-chain and on-chain (Solidity can reproduce by deterministic encoding).

### 3.2 Merkle
Merkle root uses:
- leaf hashes are 32-byte sha256 digests
- parent = sha256(left || right)
- if odd count, duplicate last (“dup-last”)

This matches GL1’s Merkle utilities.

---

## 4. NFT Mapping (On-chain)

### 4.1 Token IDs
Token IDs SHOULD be the 256-bit value of the object id:
- `tokenId = uint256(sha256_hex)`
or stored as `bytes32`.

### 4.2 Hierarchy
- ChromosomeNFT token id == chromosome object id
- ChunkNFT token id == chunk object id
- ChromosomeNFT “children” are the `chunks[]` list (implicit hierarchy)

### 4.3 Deployment Modes
- **fully-onchain**: ChunkNFT stores the blob bytes directly or via a data contract.
- **anchor/ipfs**: ChunkNFT stores only `(cid, sha256, byte_length)`.
In both cases, `sha256` is authoritative.

---

## 5. GL1 Compatibility (Import / Export)

### 5.1 Export
An Assembly object can always be exported to:
- `.gl1`: concatenating stored chunk blobs in record order
- `.gl1x`: built from chromosome/chunk metadata and computed offsets
- `.gl1q`: if qual chunks exist, reconstructed from stored quality-chunk blobs
- `.gl1i`: stored as an index blob if built

### 5.2 Import
A `.gl1` + `.gl1x` (+ optional `.gl1q`) is ingested by:
- storing each chunk’s exact bytes as a blob
- partitioning mask/exceptions into chunk-local lists
- building chromosome + assembly objects
- creating a commit

---

## 6. Search and Alignment Readiness

### 6.1 Minimizer Index (GL1I)
Index objects can store a `.gl1i` blob produced with fixed `(k,w,canonical,bucket_count)` parameters.
This enables fast seed hits for aligners and k-mer/minimizer search.

### 6.2 Alignment-Ready Extension (future)
Add `alnset` objects referencing:
- reference commit
- read-set commit
- alignment blocks (compressed, columnar)
- region index (CSI-like bins -> block pointers)
These can also be Merkle-committed and deployed as IndexNFTs.

---

## 7. Evolution via Commits
Because chromosome manifests are lists of chunk ids:
- a mutation affecting only a region rewrites only the overlapping chunks
- all other chunk ids are reused
- the new chromosome id and assembly id differ, but storage is incremental

This makes evolutionary simulation natural:
- each evolutionary step = commit
- history is preserved and replayable
- diffs are chunk-level by default; base-level diffs can be computed by region decode.

---

## 8. Summary
This spec makes GL1 into a **history-aware, NFT-native sequence system**:
- GL1 remains the payload codec and file compatibility layer.
- SeqNFT objects add Git-like commits, on-chain deployment, and explicit integrity roots.
- The same object graph runs both on-chain and off-chain.
