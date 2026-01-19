# Trinity Version Output Location Differences

## The Issue: Different File Locations by Version

Trinity v2.9.1 and v2.15.2 output files in **different locations** with **different names**.

---

## Version Comparison

### Trinity v2.9.1 (Older - DC1L1, DC1L2, etc.)

**Output location:** Inside the assembly directory

```
01_processed/
└── DC1L1_trinity_assembly/
    ├── Trinity.fasta                    ← HERE (267 MB)
    ├── Trinity.fasta.gene_trans_map     ← HERE
    └── [other files]
```

**File path:**
```bash
01_processed/DC1L1_trinity_assembly/Trinity.fasta
```

---

### Trinity v2.15.2 (Newer - DC1R3, DC2L1, etc.)

**Output location:** In the parent directory with prefix

```
01_processed/
├── DC1R3_trinity_assembly/               ← Directory
│   └── [intermediate files]
├── DC1R3_trinity_assembly.Trinity.fasta  ← HERE (164 MB)
└── DC1R3_trinity_assembly.Trinity.fasta.gene_trans_map  ← HERE
```

**File path:**
```bash
01_processed/DC1R3_trinity_assembly.Trinity.fasta
```

---

## How to Find Your Files

### For v2.9.1 assemblies (older):
```bash
cd 01_processed/DC1L1_trinity_assembly
ls -lh Trinity.fasta
```

### For v2.15.2 assemblies (newer):
```bash
cd 01_processed
ls -lh DC1R3_trinity_assembly.Trinity.fasta
```

---

## Quick Check Script

To find all Trinity.fasta files regardless of version:

```bash
# Find all Trinity.fasta files
find 01_processed -name "*Trinity.fasta" -type f

# Or more specifically:
find 01_processed -name "Trinity.fasta" -o -name "*_trinity_assembly.Trinity.fasta"
```

---

## For Your RSEM Quantification Script

Your `run_trinity_rsem_all.sbatch` script looks for:
```bash
transcripts="${asm_dir}/Trinity.fasta"
```

This will work for v2.9.1 but **NOT for v2.15.2**.

**Fix needed:** Update the script to check both locations:

```bash
# Check for Trinity.fasta in both possible locations
if [[ -f "${asm_dir}/Trinity.fasta" ]]; then
    transcripts="${asm_dir}/Trinity.fasta"
elif [[ -f "${asm_dir}.Trinity.fasta" ]]; then
    transcripts="${asm_dir}.Trinity.fasta"
else
    echo "No Trinity.fasta found for ${asm_dir}, skipping"
    continue
fi
```

---

## Summary

| Version | File Location | File Name Pattern |
|---------|---------------|-------------------|
| v2.9.1 | Inside assembly dir | `{assembly_dir}/Trinity.fasta` |
| v2.15.2 | Parent directory | `{assembly_dir}.Trinity.fasta` |

**DC1R3 is complete and working!** The file is just in a different location than expected.

---

## Next Steps

1. ✅ **DC1R3 is complete** - File exists at `DC1R3_trinity_assembly.Trinity.fasta` (164 MB)

2. **Update RSEM script** to handle both versions

3. **For consistency**, you could:
   - Copy/move v2.15.2 files to match v2.9.1 structure
   - Or update all scripts to check both locations
   - Or standardize on one version

4. **Continue with RSEM quantification** - The assembly is ready!


