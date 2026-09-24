"""Core logic for the short-read assembler.

Reads may originate from either strand, so every read is used exactly once
in either its forward ("+") or reverse-complement ("-") orientation.  An
assembly is a permutation of the reads together with one orientation per
read; adjacent reads are joined by the maximal exact suffix/prefix overlap
between the previous read and the next one.

The optimal assembly minimises the final assembled length.  Ties are broken
first by the lexicographically smallest assembled string, then by the
lexicographically smallest id/orientation path (compared element by element
as ``(id, orientation)`` pairs along the assembly order).
"""
from __future__ import annotations

from dataclasses import dataclass

MIN_READS = 2
MAX_READS = 10
MIN_READ_LEN = 4
MAX_READ_LEN = 30

ORIENTATIONS = ("+", "-")
_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class AssemblyError(ValueError):
    """Raised when input data or a read set fails validation."""


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an A/C/G/T string."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Read:
    id: str
    seq: str


@dataclass(frozen=True)
class Placement:
    """One read's placement: orientation and half-open [start, end) span."""

    id: str
    orientation: str
    start: int
    end: int


@dataclass(frozen=True)
class Assembly:
    sequence: str
    placements: tuple[Placement, ...]

    @property
    def length(self) -> int:
        return len(self.sequence)


def max_overlap(left: str, right: str) -> int:
    """Largest k such that the length-k suffix of ``left`` equals the
    length-k prefix of ``right``."""
    for k in range(min(len(left), len(right)), 0, -1):
        if left.endswith(right[:k]):
            return k
    return 0


def validate_reads(reads: list[Read]) -> None:
    """Reject an entire read set that violates any input constraint."""
    count = len(reads)
    if not MIN_READS <= count <= MAX_READS:
        raise AssemblyError(
            f"expected between {MIN_READS} and {MAX_READS} reads, got {count}"
        )
    seen: set[str] = set()
    for read in reads:
        if not read.id or any(ord(c) < 0x21 or ord(c) > 0x7E for c in read.id):
            raise AssemblyError(
                f"invalid read id {read.id!r}: need non-empty printable ASCII"
            )
        if read.id in seen:
            raise AssemblyError(f"duplicate read id {read.id!r}")
        seen.add(read.id)
        if not MIN_READ_LEN <= len(read.seq) <= MAX_READ_LEN:
            raise AssemblyError(
                f"read {read.id}: length {len(read.seq)} outside "
                f"{MIN_READ_LEN}..{MAX_READ_LEN}"
            )
        bad = sorted(set(read.seq) - _BASES)
        if bad:
            raise AssemblyError(
                f"read {read.id}: invalid bases {''.join(bad)!r}, expected A/C/G/T"
            )
    # Reverse-complement pairs (necessarily equal length) are rejected...
    for i in range(count):
        for j in range(i + 1, count):
            if reads[i].seq == reverse_complement(reads[j].seq):
                raise AssemblyError(
                    f"reads {reads[i].id} and {reads[j].id} are reverse complements"
                )
    # ...as is any read fully contained in another, in either orientation.
    for container in reads:
        for other in reads:
            if other is container or len(other.seq) > len(container.seq):
                continue
            if (
                other.seq in container.seq
                or reverse_complement(other.seq) in container.seq
            ):
                raise AssemblyError(
                    f"read {other.id} is contained in read {container.id}"
                )


def assemble(reads: list[Read]) -> Assembly:
    """Return the optimal assembly of ``reads`` (see module docstring)."""
    validate_reads(reads)
    n = len(reads)
    ids = [read.id for read in reads]
    index_of = {rid: i for i, rid in enumerate(ids)}
    # oriented[i][0] is the forward sequence, oriented[i][1] the reverse
    # complement.
    oriented = [(read.seq, reverse_complement(read.seq)) for read in reads]

    # overlap[(i, oi, j, oj)]: max suffix/prefix overlap between read i and
    # read j in the given orientations.
    overlap = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            for oi in range(2):
                for oj in range(2):
                    overlap[(i, oi, j, oj)] = max_overlap(
                        oriented[i][oi], oriented[j][oj]
                    )

    # Held-Karp DP over subsets.  best[(mask, j, oj)] holds the optimal
    # (length, string, path) assembly of exactly the reads in ``mask`` that
    # ends with read j in orientation oj.  Candidates for the same state are
    # compared by the same (length, string, path) ordering used for the
    # final answer: equal-length candidates have equal-length strings, so
    # their relative order is preserved when identical suffixes are appended
    # by later transitions.
    best: dict[tuple[int, int, int], tuple[int, str, tuple]] = {}
    for j in range(n):
        for oj in range(2):
            seq = oriented[j][oj]
            best[(1 << j, j, oj)] = (
                len(seq),
                seq,
                ((ids[j], ORIENTATIONS[oj]),),
            )

    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for j in range(n):
            if not mask & (1 << j):
                continue
            for oj in range(2):
                state = best.get((mask, j, oj))
                if state is None:
                    continue
                cur_len, cur_str, cur_path = state
                for k in range(n):
                    if mask & (1 << k):
                        continue
                    for ok in range(2):
                        ov = overlap[(j, oj, k, ok)]
                        nxt = oriented[k][ok]
                        cand = (
                            cur_len + len(nxt) - ov,
                            cur_str + nxt[ov:],
                            cur_path + ((ids[k], ORIENTATIONS[ok]),),
                        )
                        key = (mask | (1 << k), k, ok)
                        if key not in best or cand < best[key]:
                            best[key] = cand

    _, sequence, path = min(
        best[(full, j, oj)] for j in range(n) for oj in range(2)
    )

    # Reconstruct per-read coordinates along the winning path.
    placements = []
    offset = 0
    prev_seq = None
    for rid, sign in path:
        seq = oriented[index_of[rid]][0 if sign == "+" else 1]
        if prev_seq is not None:
            offset -= max_overlap(prev_seq, seq)
        placements.append(
            Placement(id=rid, orientation=sign, start=offset, end=offset + len(seq))
        )
        offset += len(seq)
        prev_seq = seq
    return Assembly(sequence=sequence, placements=tuple(placements))


def reads_from_data(data) -> list[Read]:
    """Build and validate reads from decoded JSON input."""
    if not isinstance(data, dict) or not isinstance(data.get("reads"), list):
        raise AssemblyError('input must be a JSON object with a "reads" array')
    reads = []
    for pos, item in enumerate(data["reads"]):
        if not isinstance(item, dict):
            raise AssemblyError(f"reads[{pos}] must be an object")
        rid, seq = item.get("id"), item.get("seq")
        if not isinstance(rid, str) or not isinstance(seq, str):
            raise AssemblyError(f'reads[{pos}] must have string "id" and "seq"')
        reads.append(Read(id=rid, seq=seq))
    validate_reads(reads)
    return reads


def assembly_to_data(assembly: Assembly) -> dict:
    """Serialise an assembly to the JSON output document."""
    return {
        "assembly": assembly.sequence,
        "length": assembly.length,
        "placements": [
            {
                "id": p.id,
                "orientation": p.orientation,
                "start": p.start,
                "end": p.end,
            }
            for p in assembly.placements
        ],
    }
