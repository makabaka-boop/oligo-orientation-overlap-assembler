"""Tests for the assembler.

The optimiser is cross-checked against an exhaustive enumeration of every
permutation and every orientation assignment for inputs of up to five
reads, alongside directed cases for zero-overlap reads, reverse-complement
palindromes and tied optima.
"""
from __future__ import annotations

import itertools
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

from assembler.core import (
    AssemblyError,
    Read,
    assemble,
    assembly_to_data,
    max_overlap,
    reads_from_data,
    reverse_complement,
    validate_reads,
)

ROOT = Path(__file__).resolve().parent.parent
SEED = 20260924


def rc(seq: str) -> str:
    return reverse_complement(seq)


def oriented(seq: str, sign: str) -> str:
    return seq if sign == "+" else rc(seq)


def brute_force(reads):
    """Optimal (sequence, path) by exhaustive permutation x orientation search."""
    best = None
    for perm in itertools.permutations(range(len(reads))):
        for signs in itertools.product("+-", repeat=len(reads)):
            seqs = [oriented(reads[i].seq, s) for i, s in zip(perm, signs)]
            assembled = seqs[0]
            for left, right in zip(seqs, seqs[1:]):
                assembled += right[max_overlap(left, right):]
            path = tuple((reads[i].id, s) for i, s in zip(perm, signs))
            key = (len(assembled), assembled, path)
            if best is None or key < best[0]:
                best = (key, assembled, path)
    return best[1], best[2]


def assert_matches_brute_force(reads):
    expected_seq, expected_path = brute_force(reads)
    result = assemble(reads)
    assert result.sequence == expected_seq
    assert [(p.id, p.orientation) for p in result.placements] == list(expected_path)
    return result


def assert_consistent(reads, asm):
    """Placements must reconstruct the assembled sequence exactly."""
    by_id = {r.id: r for r in reads}
    assert len(asm.placements) == len(reads)
    assert len({p.id for p in asm.placements}) == len(reads)
    for p in asm.placements:
        seq = oriented(by_id[p.id].seq, p.orientation)
        assert p.orientation in ("+", "-")
        assert 0 <= p.start < p.end <= len(asm.sequence)
        assert p.end - p.start == len(seq)
        assert asm.sequence[p.start : p.end] == seq
    starts = [p.start for p in asm.placements]
    assert starts[0] == 0
    assert all(a < b for a, b in zip(starts, starts[1:]))
    assert asm.placements[-1].end == len(asm.sequence)


def random_valid_reads(rng, n, min_len=4, max_len=10, palindrome_chance=0.0):
    while True:
        reads = []
        for i in range(n):
            length = rng.randrange(min_len, max_len + 1)
            if palindrome_chance and rng.random() < palindrome_chance:
                half = "".join(rng.choice("ACGT") for _ in range(length // 2))
                seq = half + rc(half)
            else:
                seq = "".join(rng.choice("ACGT") for _ in range(length))
            reads.append(Read(id=f"r{i}", seq=seq))
        try:
            validate_reads(reads)
            return reads
        except AssemblyError:
            continue


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------


def test_reverse_complement():
    assert rc("AACCGG") == "CCGGTT"
    assert rc("ACGT") == "ACGT"  # reverse-complement palindrome
    assert rc("AATT") == "AATT"  # reverse-complement palindrome
    rng = random.Random(SEED)
    for _ in range(200):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randrange(4, 31)))
        assert rc(rc(seq)) == seq


def test_max_overlap():
    assert max_overlap("AACCGG", "CGGTTA") == 3
    assert max_overlap("ACGT", "CACG") == 0
    assert max_overlap("AABC", "ABCC") == 3
    assert max_overlap("AAAC", "AAAC") == 4  # defensive; validation forbids it
    assert max_overlap("ACGTACGT", "GT") == 2


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


def test_validate_accepts_palindromic_read():
    validate_reads([Read("p", "ACGT"), Read("q", "CGTAC")])


@pytest.mark.parametrize(
    "reads, match",
    [
        ([Read("a", "ACGT")], "between 2 and 10"),
        ([Read(f"r{i}", "ACGT") for i in range(11)], "between 2 and 10"),
        ([Read("a", "ACGT"), Read("a", "TGCA")], "duplicate read id"),
        ([Read("", "ACGT"), Read("b", "TGCA")], "invalid read id"),
        ([Read("a b", "ACGT"), Read("b", "TGCA")], "invalid read id"),
        ([Read("réad", "ACGT"), Read("b", "TGCA")], "invalid read id"),
        ([Read("a", "ACG"), Read("b", "TGCA")], "length 3 outside"),
        ([Read("a", "A" * 31), Read("b", "TGCA")], "length 31 outside"),
        ([Read("a", "ACGX"), Read("b", "TGCA")], "invalid bases"),
        ([Read("a", "acgt"), Read("b", "TGCA")], "invalid bases"),
        ([Read("a", "AAAACCCC"), Read("b", "GGGGTTTT")], "reverse complements"),
        ([Read("a", "AACCGGTT"), Read("b", "CCGG")], "contained"),
        ([Read("a", "AACCGGTT"), Read("b", "AACC")], "contained"),  # rc(b)=GGTT
        ([Read("a", "ACGTAC"), Read("b", "ACGTAC")], "contained"),
    ],
)
def test_validate_rejects(reads, match):
    with pytest.raises(AssemblyError, match=match):
        validate_reads(reads)


def test_reads_from_data_errors():
    with pytest.raises(AssemblyError, match='"reads" array'):
        reads_from_data([1, 2])
    with pytest.raises(AssemblyError, match='"reads" array'):
        reads_from_data({"reads": "nope"})
    with pytest.raises(AssemblyError, match="must be an object"):
        reads_from_data({"reads": [42]})
    with pytest.raises(AssemblyError, match='string "id" and "seq"'):
        reads_from_data({"reads": [{"id": 1, "seq": "ACGT"}]})
    with pytest.raises(AssemblyError, match="reverse complements"):
        reads_from_data(
            {"reads": [{"id": "a", "seq": "AAAACCCC"}, {"id": "b", "seq": "GGGGTTTT"}]}
        )


def test_reads_from_data_ignores_extra_keys():
    reads = reads_from_data(
        {
            "reads": [
                {"id": "a", "seq": "ACGT", "quality": "!!!!"},
                {"id": "b", "seq": "CGTA"},
            ],
            "meta": {"source": "test"},
        }
    )
    assert reads == [Read("a", "ACGT"), Read("b", "CGTA")]


# ---------------------------------------------------------------------------
# directed assembly cases
# ---------------------------------------------------------------------------


def test_simple_overlap():
    # The 5-base overlap only appears once r0 is reverse complemented
    # (CCGGTT -> CGGTTA).  The mirror chain (r1 reversed, then r0 forward)
    # also has length 7 but assembles to the larger string "TAACCGG".
    reads = [Read("r0", "AACCGG"), Read("r1", "CGGTTA")]
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    assert asm.sequence == "CCGGTTA"
    assert [(p.id, p.orientation, p.start, p.end) for p in asm.placements] == [
        ("r0", "-", 0, 6),
        ("r1", "+", 1, 7),
    ]


def test_optimum_uses_reverse_complement():
    # Only by reversing r1 (AAAAGGGG -> CCCCTTTT) does the 4-base overlap
    # with r0 appear.
    reads = [Read("r0", "AAAACCCC"), Read("r1", "AAAAGGGG")]
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    assert asm.sequence == "AAAACCCCTTTT"
    assert [(p.id, p.orientation, p.start, p.end) for p in asm.placements] == [
        ("r0", "+", 0, 8),
        ("r1", "-", 4, 12),
    ]


def test_three_read_chain():
    reads = [Read("r0", "AACA"), Read("r1", "ACAT"), Read("r2", "CATT")]
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    assert asm.sequence == "AACATT"
    assert [(p.id, p.orientation, p.start, p.end) for p in asm.placements] == [
        ("r0", "+", 0, 4),
        ("r1", "+", 1, 5),
        ("r2", "+", 2, 6),
    ]


def test_no_overlap_reads_are_concatenated():
    reads = [Read("r0", "ATGC"), Read("r1", "GAAAC")]
    for a, b in itertools.permutations(reads):
        for sa in "+-":
            for sb in "+-":
                assert max_overlap(oriented(a.seq, sa), oriented(b.seq, sb)) == 0
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    assert asm.length == 4 + 5
    # Every candidate ties on length; the lexicographically smallest string
    # wins, which here keeps both reads forward.
    assert asm.sequence == "ATGCGAAAC"
    assert [(p.id, p.orientation, p.start, p.end) for p in asm.placements] == [
        ("r0", "+", 0, 4),
        ("r1", "+", 4, 9),
    ]


def test_palindromic_read_prefers_forward_orientation():
    reads = [Read("r0", "ACGT"), Read("r1", "CGTAC")]  # rc(ACGT) == ACGT
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    assert asm.sequence == "ACGTAC"
    # r0 is a reverse-complement palindrome: both orientations give the same
    # assembled string, so the id/orientation path decides ("+" < "-").
    assert [(p.id, p.orientation, p.start, p.end) for p in asm.placements] == [
        ("r0", "+", 0, 4),
        ("r1", "+", 1, 6),
    ]


def test_two_palindromes_tie_broken_by_path():
    reads = [Read("r0", "ACGT"), Read("r1", "AATT")]  # both rc-palindromes
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    # No overlaps exist; all orientations are equivalent, so the smallest
    # string (AATTACGT) and then the smallest path settle it.
    assert asm.sequence == "AATTACGT"
    assert [(p.id, p.orientation) for p in asm.placements] == [
        ("r1", "+"),
        ("r0", "+"),
    ]


def test_tied_optima_pick_smallest_string():
    # Two chains achieve the minimal length 8: r0 reversed -> r1 -> r2
    # reversed gives "CCGGTTAA", the mirror chain gives "TTAACCGG".
    reads = [Read("r0", "AACCGG"), Read("r1", "CGGTTA"), Read("r2", "TTAAC")]
    asm = assert_matches_brute_force(reads)
    assert_consistent(reads, asm)
    assert asm.length == 8
    assert asm.sequence == "CCGGTTAA"
    assert [(p.id, p.orientation, p.start, p.end) for p in asm.placements] == [
        ("r0", "-", 0, 6),
        ("r1", "+", 1, 7),
        ("r2", "-", 3, 8),
    ]


# ---------------------------------------------------------------------------
# exhaustive cross-checks (all permutations x orientations, <= 5 reads)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_random_sets_match_brute_force(n):
    rng = random.Random(SEED + n)
    for _ in range(25):
        reads = random_valid_reads(rng, n)
        asm = assert_matches_brute_force(reads)
        assert_consistent(reads, asm)


def test_random_palindrome_heavy_sets_match_brute_force():
    rng = random.Random(SEED + 1000)
    for _ in range(25):
        n = rng.randrange(2, 6)
        reads = random_valid_reads(rng, n, palindrome_chance=0.6)
        asm = assert_matches_brute_force(reads)
        assert_consistent(reads, asm)


def test_all_four_mer_pairs_match_brute_force():
    mers = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
    checked = 0
    for i, first in enumerate(mers):
        for second in mers[i + 1 :]:
            reads = [Read("r0", first), Read("r1", second)]
            try:
                validate_reads(reads)
            except AssemblyError:
                continue  # reverse-complement pair
            assert_matches_brute_force(reads)
            checked += 1
    assert checked == len(mers) * (len(mers) - 1) // 2 - 120  # 120 rc pairs


# ---------------------------------------------------------------------------
# larger inputs
# ---------------------------------------------------------------------------


def test_ten_reads():
    rng = random.Random(SEED + 10)
    reads = random_valid_reads(rng, 10, min_len=8, max_len=20)
    asm = assemble(reads)
    assert_consistent(reads, asm)
    # The optimum is no worse than assembling in input order, all forward.
    naive = reads[0].seq
    for left, right in zip(reads, reads[1:]):
        naive += right.seq[max_overlap(left.seq, right.seq):]
    assert asm.length <= len(naive)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def run_cli(args=(), stdin_text=None):
    return subprocess.run(
        [sys.executable, "-m", "assembler", *args],
        input=stdin_text,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


def test_cli_stdin_three_reads():
    payload = {
        "reads": [
            {"id": "read1", "seq": "AACCGG"},
            {"id": "read2", "seq": "CGGTTA"},
            {"id": "read3", "seq": "TTAAC"},
        ]
    }
    proc = run_cli(stdin_text=json.dumps(payload))
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    # The output carries the assembled string and every placement, not just
    # the optimal length.
    assert out == {
        "assembly": "CCGGTTAA",
        "length": 8,
        "placements": [
            {"id": "read1", "orientation": "-", "start": 0, "end": 6},
            {"id": "read2", "orientation": "+", "start": 1, "end": 7},
            {"id": "read3", "orientation": "-", "start": 3, "end": 8},
        ],
    }


def test_cli_file_input(tmp_path):
    payload = {"reads": [{"id": "a", "seq": "AAAACCCC"}, {"id": "b", "seq": "AAAAGGGG"}]}
    path = tmp_path / "reads.json"
    path.write_text(json.dumps(payload))
    proc = run_cli(args=[str(path)])
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["assembly"] == "AAAACCCCTTTT"
    assert [p["orientation"] for p in out["placements"]] == ["+", "-"]


def test_cli_pretty(tmp_path):
    payload = {"reads": [{"id": "a", "seq": "AAAACCCC"}, {"id": "b", "seq": "AAAAGGGG"}]}
    pretty = run_cli(args=["--pretty"], stdin_text=json.dumps(payload))
    plain = run_cli(stdin_text=json.dumps(payload))
    assert pretty.returncode == plain.returncode == 0
    assert "\n" in pretty.stdout.strip()
    assert json.loads(pretty.stdout) == json.loads(plain.stdout)


def test_cli_invalid_json():
    proc = run_cli(stdin_text="{not json")
    assert proc.returncode == 1
    assert "invalid JSON" in json.loads(proc.stderr)["error"]


def test_cli_rejected_read_set():
    payload = {"reads": [{"id": "a", "seq": "AAAACCCC"}, {"id": "b", "seq": "GGGGTTTT"}]}
    proc = run_cli(stdin_text=json.dumps(payload))
    assert proc.returncode == 1
    assert "reverse complements" in json.loads(proc.stderr)["error"]


def test_cli_missing_file():
    proc = run_cli(args=["does-not-exist.json"])
    assert proc.returncode == 2
    assert "error" in json.loads(proc.stderr)


def test_cli_output_round_trips_through_core():
    payload = {
        "reads": [{"id": f"r{i}", "seq": s} for i, s in enumerate(
            ["AACCGG", "CGGTTA", "TTAAC", "ACGTAC", "GTACGA"]
        )]
    }
    proc = run_cli(stdin_text=json.dumps(payload))
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    reads = reads_from_data(payload)
    asm = assemble(reads)
    assert out == assembly_to_data(asm)
    assert_consistent(reads, asm)
