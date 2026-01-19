# Trinity Version Comparison: DC1L1 vs DC1R3

## Summary Comparison

### DC1L1 (Successful) ✓
- **Trinity Version:** v2.9.1
- **Trinity.fasta:** 267 MByte ✓
- **Runtime:** 14.7 hours (52,769 seconds)
- **Chrysalis:** 14.1 hours (50,879 seconds)
- **Phase 2:** 21 minutes (1,278 seconds)
- **Reads:** 29,792,544
- **KMERs:** 326,673,959

### DC1R3 (Problematic) ✗
- **Trinity Version:** v2.15.2 (newer!)
- **Trinity.fasta:** 0 MByte ✗
- **Runtime:** 79 hours / 3.3 days (284,798 seconds)
- **Chrysalis:** 78.6 hours (282,850 seconds) ⚠️
- **Phase 2:** 23.5 minutes (1,413 seconds)
- **Reads:** 26,548,730 (similar)
- **KMERs:** 255,027,636 (fewer)

---

## Key Findings

### 1. **Different Trinity Versions**

**DC1L1:** Trinity v2.9.1 (older, working version)  
**DC1R3:** Trinity v2.15.2 (newer version, but failed)

**This is likely the root cause!** Different Trinity versions may have:
- Different output formats
- Different file locations
- Different assembly algorithms
- Bugs in newer versions

---

### 2. **Massive Runtime Difference**

| Phase | DC1L1 | DC1R3 | Difference |
|-------|-------|-------|------------|
| **Total** | 14.7 hrs | 79 hrs | **5.4x longer** |
| **Chrysalis** | 14.1 hrs | 78.6 hrs | **5.6x longer** |
| **Phase 2** | 21 min | 23.5 min | Similar |

**DC1R3 took 5.4x longer** - this suggests:
- The newer version may be less efficient
- Or there were issues during processing that caused slowdowns
- Or the algorithm changed significantly

---

### 3. **Chrysalis Phase Took Much Longer**

Chrysalis (read clustering) took **78.6 hours** in DC1R3 vs **14.1 hours** in DC1L1.

This could indicate:
- Algorithm changes in v2.15.2 that are slower
- Memory issues causing swapping
- I/O bottlenecks
- Or actual problems during clustering

---

### 4. **Phase 2 Completed But No Output**

Both show Phase 2 completed:
- DC1L1: 1,278 seconds → **Trinity.fasta created (267 MB)**
- DC1R3: 1,413 seconds → **Trinity.fasta shows 0 MB**

This suggests:
- Phase 2 ran but the final merge/creation of Trinity.fasta failed
- Or the file was created in a different location
- Or there's a bug in v2.15.2's final output step

---

## What to Check

### 1. Verify Actual File Location

The newer Trinity version might output to a different location:

```bash
cd /projects/tholl_lab_1/daisy_analysis/01_processed/DC1R3_trinity_assembly

# Check for Trinity.fasta in various locations
find . -name "Trinity.fasta" -type f
find . -name "*.fasta" -type f | grep -i trinity

# Check the main directory
ls -lh Trinity.fasta 2>/dev/null || echo "Not in main directory"

# Check if there's a trinity_out_dir subdirectory (newer versions use this)
ls -lh trinity_out_dir/Trinity.fasta 2>/dev/null || echo "Not in trinity_out_dir"
```

### 2. Check Trinity Version Behavior

Different versions may have different output structures:

```bash
# Compare directory structures
echo "=== DC1L1 structure ==="
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/DC1L1_trinity_assembly/ | head -20

echo "=== DC1R3 structure ==="
ls -lh /projects/tholl_lab_1/daisy_analysis/01_processed/DC1R3_trinity_assembly/ | head -20
```

### 3. Check for Alternative Output Files

Newer Trinity might create files with different names:

```bash
cd DC1R3_trinity_assembly
find . -name "*.fasta" -type f -exec ls -lh {} \;
```

### 4. Check Logs for Version-Specific Issues

```bash
# Look for version-specific errors or warnings
grep -i "version\|v2.15\|output\|fasta" ../trinity_all_*.out | grep -i "DC1R3"
```

---

## Solutions

### Solution 1: Use Same Trinity Version as DC1L1

Since DC1L1 worked with v2.9.1, use that version for consistency:

```bash
# Check current Trinity version
Trinity --version

# If you need to switch versions, you may need to:
# 1. Install Trinity v2.9.1 in your conda environment
# 2. Or use a different conda environment with v2.9.1
```

### Solution 2: Check for File in Different Location

Newer Trinity versions might output to `trinity_out_dir/`:

```bash
cd DC1R3_trinity_assembly
if [ -d "trinity_out_dir" ]; then
    ls -lh trinity_out_dir/Trinity.fasta
fi
```

### Solution 3: Re-run with Same Version as DC1L1

For consistency across all samples, re-run DC1R3 with Trinity v2.9.1:

```bash
# Make sure you're using the same version
conda activate rnaseq
Trinity --version  # Should show v2.9.1

# Then re-run (after backing up/removing incomplete run)
```

### Solution 4: Check if File Was Created Elsewhere

Sometimes Trinity creates files in unexpected locations:

```bash
# Search the entire analysis directory
find /projects/tholl_lab_1/daisy_analysis -name "Trinity.fasta" -type f -newer DC1R3_trinity_assembly/Trinity.timing
```

---

## Recommendations

1. **Use consistent Trinity versions** across all samples
   - DC1L1 used v2.9.1 and worked
   - DC1R3 used v2.15.2 and failed
   - Use v2.9.1 for all samples

2. **Check actual file location** - newer versions may output differently

3. **If file truly missing**, re-run with v2.9.1 to match DC1L1

4. **Consider version compatibility** - mixing versions can cause issues downstream

---

## Expected Behavior

Based on DC1L1 (working example):
- Trinity.fasta should be **~267 MB** (or similar, depending on assembly size)
- File should be in the main output directory
- Runtime should be **~15 hours**, not 79 hours

DC1R3 deviates significantly from this, suggesting version-related issues.


