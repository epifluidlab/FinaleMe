# Add UXM-Compatible Output Mode (pat.gz + beta) to FinaleMe Decode

## Context

FinaleMe already predicts per-CpG methylation state for each cfDNA fragment during Viterbi decoding, but currently only outputs **aggregated** per-CpG summary statistics. The user wants a new output mode that preserves **per-fragment haplotype-level** methylation patterns in the `.pat.gz` and `.beta` formats used by [UXM_deconv](https://github.com/nloyfer/UXM_deconv), enabling direct tissues-of-origin deconvolution without format conversion.

---

## UXM File Format Reference

### `.pat.gz` (bgzip-compressed, tab-separated, no header)
| Column | Description |
|--------|-------------|
| chr | Chromosome |
| CpG index | 1-based genome-wide CpG site index |
| Pattern | `C` = methylated, `T` = unmethylated (one char per consecutive CpG) |
| Count | Number of identical fragments with same chr, start CpG index, and pattern |

Sorted by CpG index. Identical (chr, startIndex, pattern) records merged with count.

### `.beta` (binary)
- `NR_SITES` rows x 2 columns of `uint8`
- Column 0: methylated count, Column 1: total coverage (both capped at 255, proportionally scaled)
- Row index = CpG site index - 1 (0-based array for 1-based CpG index)

---

## Implementation Plan

### Step 1: Add CLI options and validation

**File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) ~line 144 (after `-t` option)

Add two new options:
```java
@Option(name="-patOutput", usage="output UXM-compatible .pat.gz and .beta files for deconvolution. Requires -cpgIndexFile. default: false")
public boolean patOutput = false;

@Option(name="-cpgIndexFile", usage="CG_motif bedgraph/bed file listing all CpG positions genome-wide, for building CpG index. Required with -patOutput.")
public String cpgIndexFile = null;
```

In `doMain()` (~line 200), add validation before `decodeHmm()`:
```java
if (patOutput && cpgIndexFile == null) {
    throw new IllegalArgumentException("-patOutput requires -cpgIndexFile");
}
```

### Step 2: Add CpgIndex inner class and loader

**File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) (after MatrixObj inner class, ~line 1207)

Add `CpgIndex` static inner class with:
- `LinkedHashMap<String, int[]> chrPositions` — chr → sorted CpG start positions (preserves chromosome order from file)
- `HashMap<String, Integer> chrOffsets` — chr → cumulative offset (total CpGs on preceding chromosomes)
- `int totalSites` — total CpG count genome-wide
- `int getGlobalIndex(String chr, int start)` — binary search in `chrPositions[chr]`, returns `chrOffsets[chr] + localIndex + 1` (1-based), or -1 if not found

Memory: ~112 MB for ~28M CpG positions as `int[]` arrays. Efficient.

Add `private CpgIndex loadCpgIndex(String cpgIndexFile)` method:
- Read bedgraph/bed file (support `.gz` via GZIPInputStream, same pattern as line 242-249)
- Parse chr (col 0) and start (col 1) from each line
- Group into per-chromosome ArrayList, then convert to sorted `int[]`
- Compute cumulative offsets in file order (chr1, chr2, ..., chrX, chrY)

### Step 3: Modify decodeHmm() to collect per-fragment patterns

**File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) `decodeHmm()` (lines 713-861)

Changes:
1. Add `CpgIndex cpgIndex` parameter to method signature (line 713)
2. Update call site in `doMain()` (line 225) to pass the loaded `cpgIndex` (or `null`)
3. Before the futures loop (~line 739), when `patOutput` is true:
   ```java
   HashMap<String, Integer> patRecords = new HashMap<>(); // key="chr\tCpGIdx\tpattern", value=count
   long skippedFragments = 0;
   ```
4. Inside the sequential aggregation loop (lines 766-821), after the existing `methySummary` logic, add pattern collection block guarded by `if (patOutput && cpgIndex != null)`:
   - For the current fragment (index `j`), use `locRow` to get all CpG locations
   - Parse chr and start from `locRow.get(0)` to get the fragment's first CpG
   - Look up global CpG index for each CpG in the fragment
   - Verify contiguity (indices form consecutive sequence); if not, increment `skippedFragments` and skip
   - Build pattern string: `'C'` if `hiddenState[i] == methylatedState`, else `'T'`
   - Merge into `patRecords` HashMap: key = `chr + "\t" + startCpGIndex + "\t" + pattern`, value += 1

### Step 4: Write `.pat.gz` output

**File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) (after existing output writing, ~line 849)

When `patOutput` is true:
1. Derive output filename: replace `.bed.gz` with `.pat.gz` in `outputFile`, or append `.pat.gz`
2. Convert `patRecords` entries to a list, sort by CpG index (numerically)
3. Write using `htsjdk.samtools.util.BlockCompressedOutputStream` (already available in htsjdk bundled in fat jar):
   ```java
   BlockCompressedOutputStream patOut = new BlockCompressedOutputStream(patFile);
   for each sorted record:
       patOut.write((chr + "\t" + cpgIndex + "\t" + pattern + "\t" + count + "\n").getBytes());
   patOut.close();
   ```
4. Log skipped fragment count if > 0

### Step 5: Write `.beta` output

**File**: [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) (immediately after `.pat.gz` writing)

When `patOutput` is true:
1. Derive output filename: replace `.pat.gz` with `.beta`
2. Allocate `int[] methyCounts = new int[cpgIndex.totalSites]` and `int[] totalCounts = new int[cpgIndex.totalSites]` (~224 MB for hg19)
3. Iterate sorted pat records; for each record, for each position `p` in pattern:
   - Array index = startCpGIndex + p - 1 (convert 1-based to 0-based)
   - If `pattern.charAt(p) == 'C'`: `methyCounts[idx] += count`
   - `totalCounts[idx] += count`
4. Write binary file: for each CpG site, write 2 bytes (uint8):
   - If `totalCounts[i] > 255`: scale `methyCounts[i] = round(methyCounts[i] * 255.0 / totalCounts[i])`, `totalCounts[i] = 255`
   - Write `(byte) methyCounts[i]`, `(byte) totalCounts[i]`
5. Use `BufferedOutputStream` wrapping `FileOutputStream` for efficiency

### Step 6: Add imports

```java
import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.util.LinkedHashMap;
```

---

## File Changes Summary

| File | Changes |
|------|---------|
| [FinaleMe.java](src/main/java/edu/northwestern/epifluidlab/finaleme/hmm/FinaleMe.java) | CLI options, CpgIndex class, loadCpgIndex(), modify decodeHmm() signature + pattern collection, pat.gz + beta writers |

All changes are in a single file, following the existing pattern (MatrixObj is already an inner class).

---

## Verification

1. Build: `mvn clean compile assembly:single`
2. Run existing tests: `mvn test` (should still pass — no behavior change without `-patOutput`)
3. Run decode with new flags on test data:
   ```
   java -cp "target/FinaleMe-0.60-jar-with-dependencies.jar:lib/jahmm-0.6.2.jar" \
     edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe \
     test.model input_matrix.bed.gz prediction.bed.gz \
     -decodeModeOnly -patOutput -cpgIndexFile CG_motif.bedgraph
   ```
4. Verify `.pat.gz`: `zcat prediction.pat.gz | head` — check 4 columns, C/T patterns, sorted by CpG index
5. Verify `.beta`: `python3 -c "import numpy as np; d=np.fromfile('prediction.beta',dtype=np.uint8).reshape(-1,2); print(d.shape, d[d[:,1]>0][:5])"` — shape should be (NR_SITES, 2)
6. Validate with UXM_deconv: load the generated files and run deconvolution
