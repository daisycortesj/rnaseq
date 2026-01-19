# Trinity Troubleshooting: Trinity.fasta 0 MByte

## Problem

Your Trinity output shows:
```
Trinity.fasta 0 MByte
```

This is **NOT expected** - a successful Trinity run should produce a non-zero Trinity.fasta file.

---

## What to Check

### 1. Verify if Trinity.fasta Actually Exists

```bash
cd /projects/tholl_lab_1/daisy_analysis/01_processed/DC1R3_trinity_assembly
ls -lh Trinity.fasta
```

**Expected:** File should exist and be > 0 bytes (typically hundreds of MB for your data size)

**If file doesn't exist or is 0 bytes:** Trinity may have failed during final assembly

---

### 2. Check for Error Messages

```bash
# Check the Trinity log output
grep -i "error\|fail\|warn" trinity_all_*.out | tail -50

# Check for any error files
ls -lh *.err *.log 2>/dev/null
```

---

### 3. Check Phase 2 Completion Status

```bash
cd DC1R3_trinity_assembly

# Check if all commands completed
wc -l recursive_trinity.cmds
wc -l recursive_trinity.cmds.completed

# Compare sizes (should match if complete)
ls -lh recursive_trinity.cmds*
```

**Expected:** Both files should have the same number of lines when complete

---

### 4. Check for Butterfly Output

Phase 2 (Butterfly) creates the final assembly. Check if Butterfly completed:

```bash
# Look for Butterfly output in read_partitions
find read_partitions -name "*.allProbPaths.fasta" | head -5

# Check if any Butterfly processes failed
grep -i "butterfly\|error" recursive_trinity.cmds.completed | tail -20
```

---

### 5. Check Disk Space

Trinity needs significant disk space. Check if you ran out:

```bash
df -h /projects/tholl_lab_1/daisy_analysis/01_processed/DC1R3_trinity_assembly
```

---

## Common Causes

### Cause 1: Phase 2 Didn't Fully Complete

Even though the summary says it completed, some Butterfly processes may have failed silently.

**Solution:**
```bash
# Check completion percentage
completed=$(wc -l < recursive_trinity.cmds.completed)
total=$(wc -l < recursive_trinity.cmds)
percent=$(echo "scale=2; $completed / $total * 100" | bc)
echo "Completion: $percent%"
```

If < 100%, Phase 2 didn't complete.

---

### Cause 2: Butterfly Failed to Merge Assemblies

Butterfly needs to merge all the individual cluster assemblies into Trinity.fasta.

**Check:**
```bash
# Look for Butterfly output files
find read_partitions -name "*.allProbPaths.fasta" | wc -l

# Should have many files (one per cluster)
```

---

### Cause 3: Insufficient Memory During Final Merge

The final merge step can require a lot of memory.

**Check:**
```bash
# Look for out-of-memory errors
grep -i "killed\|oom\|memory" trinity_all_*.out
```

---

### Cause 4: File System Issues

Network filesystems can cause issues with large file creation.

**Check:**
```bash
# Try to create a test file
touch test_write.txt && rm test_write.txt
echo "File system is writable"
```

---

## Solutions

### Solution 1: Re-run Trinity (If Phase 2 Incomplete)

If Phase 2 didn't complete, you may need to re-run:

```bash
# Remove incomplete output (CAREFUL - backup first!)
# Only do this if you're sure it failed
# rm -rf /projects/tholl_lab_1/daisy_analysis/01_processed/DC1R3_trinity_assembly

# Re-run Trinity
sbatch scripts/run_trinity.sbatch 00_1_DC DC1R3
```

---

### Solution 2: Manually Complete Phase 2

If most of Phase 2 completed, you might be able to manually trigger the final merge:

```bash
cd DC1R3_trinity_assembly

# Check if you can manually run Butterfly merge
# This is complex and may require Trinity expertise
```

**Note:** This is advanced and may not work. Usually better to re-run.

---

### Solution 3: Check Trinity Logs for Specific Errors

```bash
# Look for specific error patterns
grep -i "butterfly.*fail\|cannot.*merge\|error.*assembly" trinity_all_*.out
```

---

## Comparison with Working Sample (DC1L1)

DC1L1 has:
- `Trinity.fasta` (268M) ✓
- `Trinity.fasta.gene_trans_map` (17M) ✓
- Complete `recursive_trinity.cmds.completed` ✓

DC1R3 should have the same files. If it doesn't, something failed.

---

## Quick Diagnostic Commands

Run these to get a full picture:

```bash
cd /projects/tholl_lab_1/daisy_analysis/01_processed/DC1R3_trinity_assembly

echo "=== File Check ==="
ls -lh Trinity.fasta 2>/dev/null || echo "Trinity.fasta NOT FOUND"

echo "=== Completion Status ==="
if [ -f recursive_trinity.cmds ] && [ -f recursive_trinity.cmds.completed ]; then
    total=$(wc -l < recursive_trinity.cmds)
    done=$(wc -l < recursive_trinity.cmds.completed)
    echo "Commands: $done / $total"
    if [ "$total" -gt 0 ]; then
        percent=$(echo "scale=1; $done * 100 / $total" | bc)
        echo "Completion: ${percent}%"
    fi
fi

echo "=== Butterfly Output ==="
find read_partitions -name "*.allProbPaths.fasta" 2>/dev/null | wc -l | xargs echo "Butterfly output files:"

echo "=== Disk Space ==="
df -h . | tail -1

echo "=== Recent Errors ==="
grep -i "error\|fail" ../trinity_all_*.out 2>/dev/null | tail -10
```

---

## Expected vs Actual

**Expected Output:**
```
Trinity.fasta 268 MByte  (or similar, > 0)
```

**Your Output:**
```
Trinity.fasta 0 MByte
```

This indicates a problem that needs investigation.

---

## Next Steps

1. **First:** Check if Trinity.fasta actually exists and its real size
2. **Second:** Check completion status of Phase 2
3. **Third:** Look for error messages in logs
4. **Fourth:** If incomplete, consider re-running Trinity

The fact that DC1L1 worked suggests your setup is correct, so this might be a transient issue (disk space, memory, network filesystem hiccup, etc.).


