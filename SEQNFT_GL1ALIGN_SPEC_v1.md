# GL1 Align Stack + SeqNFT: Superscalable, Verifiable Alignment Architecture (v1)

**Spec ID:** `GL1ALIGNv1`  
**Targets:** GL1ENCv2 + SeqNFT pipeline (commit DAG + chunk NFTs)  
**Primary goal:** make alignments and downstream analysis **(a) scalable across thousands of genomes**, **(b) content-addressed & reproducible**, and **(c) deployable with identical semantics on-chain and off-chain** (bytes stored on-chain or anchored via IPFS).

---

## 0. Abstract

We define an alignment stack for **SeqNFT + GL1** where:

1. **Genomes are stored as immutable chunk objects** (child NFTs), grouped into chromosomes (parent NFTs), grouped into assemblies and commits (version NFTs).
2. **Indexes are built at the chunk level** so they naturally deduplicate across commits/genomes and are incrementally updateable.
3. **Alignment results are stored as content-addressed blocks** with a Merkle root, creating an “alignment receipt” that is reproducible and verifiable on-chain or off-chain.
4. The *same object graph* (chunks, commits, indexes, alignments) supports two deployment modes:
   - fully on-chain (bytes stored directly),
   - anchored/IPFS (store CID + sha256, fetch & verify client-side).

The result is a pipeline that behaves like **git for genomes + verifiable “BAM on-chain semantics”**, while retaining the high-performance codec benefits of GL1.

---

## 1. Why this scales across thousands of genomes

### 1.1 Chunk-level deduplication is the unit of scale
Most genomics tooling treats each genome as a monolithic FASTA. In SeqNFT, the unit is a **chunk NFT**:

- identical sequence chunks across genomes map to the same `chunk_id` (content addressed),
- indices and alignments reference chunk_ids, not raw files.

So with a cohort (10k genomes), storage and indexing scale closer to:

> **unique chunks** rather than **genomes × full length**.

This is the same conceptual win that git has over copying whole files: unchanged blocks are reused.

### 1.2 Sharded seed directories are IPFS-native
A global seed index across thousands of genomes is too large to store as one blob. Instead, we define **GL1D** as an array of bucket blobs, each retrievable independently:

- query a seed → compute bucket id → fetch one small bucket blob → verify sha256 → decode postings.

This works naturally with IPFS and on-chain anchoring (store only roots and CIDs).

### 1.3 Incremental updates are native
A new commit changes only the affected chunks. So you:

- build indices only for **new chunks**,
- update seed directory buckets only for affected seeds,
- append alignment blocks only for new read sets.

No re-indexing of the full reference is required.

---

## 2. Concepts and invariants

### 2.1 Determinism invariants
All objects and blobs are content-addressed by `sha256` in the reference pipeline:

- **blob id** = `sha256(blob_bytes)`
- **object id** = `sha256( b"seqnft:" + type + b"\0" + canonical_json(core) )`

Canonical JSON means:
- sorted keys
- separators `(",", ":")`
- UTF-8 encoding

### 2.2 Semantic sequence invariants
SeqNFT chunks store:
- a GL1 chunk payload (`DNA2`, `DNA4`, `SIXBIT`, `ASCII`) and
- optional exceptions list (`DNA2` case) to restore IUPAC symbols.

**Semantic sequence** is *payload decoded + exceptions applied*.  
All indices and alignment computations are defined on the semantic sequence.

---

## 3. Alignment stack overview

The GL1 Align Stack introduces these artifacts:

### 3.1 GL1I2: per-chunk seed index
**Purpose:** index seeds within one chunk, deterministically.

- Input: semantic chunk sequence
- Output: mapping `seed_key -> [positions]` (positions are within chunk)

Stored as:
- `chunk_index` object (metadata)
- `gl1i2_blob` (binary bytes)

### 3.2 GL1D: sharded seed directory
**Purpose:** global mapping from seeds to chunk_ids.

- Input: many GL1I2 chunk indices
- Output: bucket blobs mapping `seed_key -> [chunk_id]`

Stored as:
- `seed_dir` object (metadata, bucket references, Merkle root)
- `bucket blobs` (binary bytes per bucket)

### 3.3 GL1A: alignment blocks + Merkle root
**Purpose:** store alignments as immutable, verifiable blocks.

- Input: mapping results (read → {ref chrom, pos, cigar, score,…})
- Output: blocks of alignments partitioned by reference coordinate

Stored as:
- `alnset` object (metadata + ordered block list)
- `ALN1` block blobs (compressed JSONL)
- Merkle root over all block hashes (`aln_root`)

---

## 4. Object model additions (SeqNFT CAS objects)

All objects below are content addressed.

### 4.1 `chunk_index` object (`type="chunk_index"`)
Represents a GL1I2 index for one chunk under fixed seeding params.

Core fields (v1):
- `spec`: `"GL1ALIGNv1"`
- `chunk`: chunk object id (the chunk NFT id)
- `chunk_kind`, `chunk_codec`
- `seed`: `{scheme,k,w,canonical,kind,max_occ}`
- `gl1i2_blob`: blob id (sha256 of GL1I2 bytes)
- `sha256`: sha256 of GL1I2 bytes
- `key_count`, `sym_count`
- `built` timestamp

### 4.2 `seed_dir` object (`type="seed_dir"`)
Represents a sharded directory built for **one seeding configuration**.

Core fields (v1):
- `spec`: `"GL1ALIGNv1"`
- `seed`: `{scheme,k,w,canonical,kind,max_occ}`
- `bucket_count`
- `buckets`: list of bucket blob ids (length = bucket_count)
- `bucket_sha256`: list of sha256 per bucket blob
- `dir_root`: Merkle root over `bucket_sha256` digests
- `built`

### 4.3 `alnset` object (`type="alnset"`)
Represents the full alignment output of reads vs a reference commit.

Core fields:
- `spec`: `"GL1ALIGNv1"`
- `ref_commit`: commit id mapped against
- `reads`: `{path, sha256}` (source identity)
- `seed_dir`: seed_dir id used
- `params`: seed params and aligner params
- `block_bp`: block size used (e.g., 1,000,000)
- `blocks`: list of block descriptors
- `aln_root`: Merkle root over ordered block digests
- `stats`: `{reads_total, reads_mapped}`
- `built`

Block descriptor fields:
- `chrom`, `start`, `end`
- `block_blob` (blob id)
- `sha256` (sha256 of block bytes)
- `count` record count

---

## 5. Binary formats

This section defines compact binary blobs to keep IPFS storage and on-chain anchors small.

### 5.1 GL1I2 chunk index (`gl1i2_blob`)
**Magic:** `b"GLI2"`  
**Version:** `1`

**Header (little endian):**
```
I2_HDR_FMT = "<4sHBBBBHIIQQQ"
magic      : 4s  ("GLI2")
ver        : u16 (1)
scheme     : u8  (1=minimizer)
kind       : u8  (1=nuc, 2=prot)
k          : u8
w          : u8
flags      : u16 (bit0=canonical)
sym_count  : u32  (#symbols in chunk)
key_count  : u32  (#distinct seed keys)
dir_off    : u64  (offset of directory; fixed = header size)
data_off   : u64  (offset of postings)
reserved   : u64
```

**Directory entries (sorted by key):**
```
I2_DIR_FMT = "<QQII"
key        : u64
pos_off    : u64  (absolute offset to postings bytes)
byte_len   : u32  (postings byte length)
pos_count  : u32  (#positions for this key)
```

**Postings bytes:**
- positions are delta-varint encoded:
  - first position as varint absolute
  - remaining positions as varint deltas

### 5.2 GL1D bucket blob (`buckets[i]`)
**Magic:** `b"G1DB"`  
**Version:** `1`

**Header:**
```
DB_HDR_FMT = "<4sHIIQQ"
magic      : 4s ("G1DB")
ver        : u16 (1)
bucket_id  : u32
key_count  : u32
dir_off    : u64 (header size)
data_off   : u64 (directory end)
```

**Directory entries (sorted by key):**
```
DB_DIR_FMT = "<QQII"
key         : u64
chunk_off   : u64  (absolute offset to chunk id list)
chunk_count : u32
reserved    : u32
```

**Chunk id lists:**
- each chunk id is 32 bytes (bytes32)
- duplicates removed in builder

### 5.3 ALN1 alignment block (`block_blob`)
**Magic:** `b"ALN1"`  
**Version:** `1`

Header:
```
ALN_HDR_FMT = "<4sHIII"
magic     : 4s ("ALN1")
ver       : u16 (1)
rec_count : u32
jsonl_len : u32 (uncompressed)
zlen      : u32 (compressed)
```

Payload:
- zlib-compressed JSONL (one alignment record per line)

---

## 6. Algorithms

### 6.1 Index building (chunk → GL1I2)
For each chunk:
1. decode semantic sequence bytes
2. generate seed list for `scheme=minimizer(k,w,canonical)`:
   - build A/C/G/T-only k-mers
   - select minimizers using O(n) deque
3. accumulate postings: `key -> [positions]`
4. optional filter: drop keys with too many occurrences (`max_occ`)
5. write GL1I2 (header + sorted directory + postings)

### 6.2 Directory building (GL1I2 set → GL1D buckets)
For each GL1I2:
- read only keys from directory
- emit records `(key, chunk_id)` into temp bucket file `key % bucket_count`

Per bucket:
- external sort records by (key, chunk_id)
- group by key, deduplicate chunk_ids
- write G1DB bucket blob
- compute sha256 and store as blob id

Compute `dir_root = Merkle(bucket_sha256[])`.

### 6.3 Mapping (read → candidate regions → alignment)
The reference implementation does:
1. seed the read
2. query seed_dir for each key → accumulate votes per chunk_id
3. take top chunks
4. for each chunk, compute best offsets between read and chunk using GL1I2 positions
5. translate chunk offsets into reference coordinates using the commit’s chunk layout
6. extend using Smith-Waterman on a small reference window

Production engines can replace step (6) with WFA/affine DP and step (4) with full minimap-style chaining while preserving the same index artifacts.

---

## 7. On-chain / IPFS deployment

### 7.1 Anchor-only mode
Store on-chain:
- `dir_root` and `bucket_sha256[i]`
- per bucket: `cid` (IPFS), `sha256`, `byte_len`
- `aln_root` and alignment block hashes + CIDs

Clients fetch bucket/blob by CID and verify sha256 before using.

### 7.2 Fully on-chain mode
Store bucket blobs and alignment blocks directly in data contracts, referenced by `bytes32 blob_id`.

This can be expensive but yields maximum on-chain availability and “no external dependency” proofs.

---

## 8. Why this beats the *system-level* limitations of existing tooling

This stack is not “a new DP algorithm” in isolation; it’s an **end-to-end verifiable system**:

- **Git-native evolution:** alignments can be reproduced at any commit in history.
- **Deduplicated indices:** build once per chunk; reused across genomes/commits.
- **Proof-carrying alignment receipts:** output root commits the full result set.
- **IPFS-friendly sharding:** only fetch what you query (buckets / blocks).
- **Composable analysis:** alignment blocks are immutable building blocks for variant calling, annotations, and pangenome analytics.

Throughput depends on implementation language and hardware, but the *architecture* scales natively.

---

## 9. Reference Python implementation

The reference code is provided as `seqnft_align_pipeline.py` and runs next to `seqnft_pipeline.py`.

Typical workflow:

```bash
# 1) Create repo and mint a reference genome
python seqnft_pipeline.py init --repo .
python seqnft_pipeline.py mint -i ref.fa --assembly hg --chunk 1000000 --mode auto -m "hg mint"

# capture commit id from stdout as $REF_COMMIT

# 2) Build a sharded seed directory (this builds per-chunk indices too)
python seqnft_align_pipeline.py seed-dir-build --commit $REF_COMMIT --k 15 --w 10 --bucket-count 4096

# capture seed_dir id as $SEED_DIR

# 3) Align reads
python seqnft_align_pipeline.py align --ref-commit $REF_COMMIT --seed-dir $SEED_DIR --reads reads.fq.gz

# capture alnset id as $ALNSET

# 4) Query alignments in a region
python seqnft_align_pipeline.py aln-view --alnset $ALNSET --chrom chr1 --start 100000 --end 101000 --to table
```

---

## 10. Roadmap (production-grade extensions)

The spec is intentionally modular. Next steps that preserve compatibility:

1. **More seed schemes:** syncmers and strobemers (new `scheme_id` values).
2. **Chaining:** minimap-style chaining of seed hits into long-read candidates.
3. **Affine WFA extension:** fast long-read alignment while still emitting ALN blocks.
4. **Cohort membership shards:** chunk_id → bitset(genome_ids) for multi-genome pangenome mapping.
5. **Variant objects:** store VCF-like diffs as commits and align directly to variant graphs.

---

## 11. Summary

GL1 + SeqNFT already provides:
- compact storage (2/4/6-bit),
- masking and exceptions,
- deterministic hashes and Merkle roots,
- commit history and chunk reuse.

GL1ALIGNv1 adds the missing layer:
- **chunk indices**
- **sharded directories**
- **verifiable alignment outputs**

…so you get a scalable genomics substrate that is naturally deployable on-chain or IPFS.



---

## 12. Cohort mode (thousands of genomes)

The core “thousands of genomes” scaling feature is **cohort mode**, which adds:

1. `cohort` object: the list of commit ids representing genomes.
2. `membership` object: sharded mapping `chunk_id -> genome_ids` (sparse, varint-coded).
3. cohort-wide `seed_dir`: built across all unique chunks in the cohort.

### 12.1 `cohort` object
`type="cohort"`

Core fields:
- `name`
- `commits`: ordered list of commit ids

The order defines `genome_id` indices (`0..N-1`) for membership bitsets/lists.

### 12.2 `membership` object
`type="membership"`, `kind="cohort_membership"`

Core fields:
- `cohort`: cohort object id
- `bucket_count`
- `buckets`: list of bucket blob ids
- `bucket_sha256`: list of sha256 per bucket
- `root`: Merkle root over bucket hashes

### 12.3 CMB1 membership bucket format
**Magic:** `b"CMB1"`  
**Version:** `1`

Header:
```
CMB_HDR_FMT = "<4sHIIQQ"
magic      : 4s ("CMB1")
ver        : u16 (1)
bucket_id  : u32
entry_count: u32
dir_off    : u64
data_off   : u64
```

Directory entry:
```
CMB_DIR_FMT = "<32sQII"
chunk_id    : 32 bytes
off         : u64  (absolute offset to data)
byte_len    : u32
gid_count   : u32
```

Data:
- varint-delta-coded list of genome ids:
  - first gid absolute
  - remaining as delta vs previous
- sorted, unique

### 12.4 Cohort alignment behavior
Given reads:
1. seed the read
2. use cohort `seed_dir` to retrieve candidate chunks
3. use `membership` to translate chunks into candidate genomes (commit ids)
4. align read to the top-N genomes (reference implementation uses N=3)

For higher performance, production engines can:
- reuse chaining and long-read alignment,
- cache buckets and chunk indices aggressively,
- incorporate cohort graph/haplotype models.

---

## 13. Deployment artifacts in the reference code

The reference `seqnft_align_pipeline.py` includes deployment helpers that work with the existing chain simulator (`ChainSim`) from `seqnft_pipeline.py`:

- `deploy-seed-dir` stores bucket blobs (on-chain or IPFS) and mints a SeedDir token JSON.
- `deploy-alnset` stores ALN1 blocks (on-chain or IPFS) and mints an AlignmentSet token JSON.
- `deploy-membership` stores membership buckets and mints a CohortMembership token JSON.

These JSON files are emitted into the repo’s `chain/` folder exactly like commit/chunk/assembly deploys.
