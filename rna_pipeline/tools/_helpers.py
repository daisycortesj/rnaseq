def sjdb_overhang_from_readlen(readlen: int) -> int:
    """
    STAR recommendation: sjdbOverhang = read_len - 1, minimum = 1.
    This small helper ensures STAR is called with a valid overhang value
    even if readlen is not provided or is invalid.
    """
    try:
        rl = int(readlen)
    except (TypeError, ValueError):
        rl = 150  # fallback default read length
    return max(1, rl - 1)
