"""Basic DNA sequence helpers used by the PCR engine."""

VALID_BASES = set("ACGTN")

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def clean_sequence(seq: str) -> str:
    """Return an uppercase sequence with spaces, numbers and line breaks removed."""
    cleaned = "".join(ch for ch in seq.upper() if ch.isalpha())

    bad = set(cleaned) - VALID_BASES
    if bad:
        raise ValueError(f"Invalid characters in sequence: {sorted(bad)}")

    return cleaned


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]



def parse_fasta(text: str) -> dict[str, str]:
    """Read FASTA text into {name: sequence}. Name = first word after '>'."""
    records = {}
    name = None
    chunks = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = clean_sequence("".join(chunks))
            parts = line[1:].split()
            name = parts[0] if parts else f"seq{len(records) + 1}"
            chunks = []
        else:
            if name is None:
                raise ValueError("FASTA must start with a '>' header line.")
            chunks.append(line)

    if name is not None:
        records[name] = clean_sequence("".join(chunks))
    if not records:
        raise ValueError("No sequences found.")
    return records