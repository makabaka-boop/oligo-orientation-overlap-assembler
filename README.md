# read-assembler

Short DNA reads may come from either strand. This assembler chooses, for
every read, exactly one orientation — forward (`+`) or reverse complement
(`-`) — and an order in which to lay the reads down, joining adjacent reads
by the **maximal exact suffix/prefix overlap** between the previous read and
the next one. It searches every permutation and orientation combination
(with a Held–Karp dynamic program, so 10 reads are no problem) and returns
the assembly with the **shortest final string**.

Ties are broken deterministically:

1. lexicographically smallest assembled string;
2. lexicographically smallest id/orientation path, compared element by
   element as `(id, orientation)` pairs along the assembly order
   (`"+"` sorts before `"-"`).

The output is the assembled string plus every read's orientation and
half-open `[start, end)` coordinates — never just the optimal length.

## Input

One JSON object with a `reads` array on stdin (or in a file argument):

```json
{
  "reads": [
    { "id": "read1", "seq": "AACCGG" },
    { "id": "read2", "seq": "CGGTTA" },
    { "id": "read3", "seq": "TTAAC" }
  ]
}
```

Constraints — violating any of them rejects the **entire** data set:

- 2 to 10 reads;
- each `id` unique, non-empty printable ASCII;
- each `seq` 4 to 30 bases of `A`/`C`/`G`/`T` only;
- no two reads are reverse complements of each other;
- no read is fully contained in another read, in either orientation.

## Output

```json
{
  "assembly": "CCGGTTAA",
  "length": 8,
  "placements": [
    { "id": "read1", "orientation": "-", "start": 0, "end": 6 },
    { "id": "read2", "orientation": "+", "start": 1, "end": 7 },
    { "id": "read3", "orientation": "-", "start": 3, "end": 8 }
  ]
}
```

`placements` is listed in assembly order; overlapping neighbours share
coordinates, and `assembly[start:end]` always equals the read in its chosen
orientation.

Exit codes: `0` success, `1` invalid JSON or rejected read set (a JSON
`{"error": ...}` document goes to stderr), `2` input/usage errors.

## Usage

Local (Python ≥ 3.10, no dependencies):

```sh
python -m assembler < reads.json
python -m assembler reads.json --pretty
```

As the Compose `assembler` service:

```sh
docker compose build
docker compose run --rm -T assembler < reads.json
```

(`-T` disables the pseudo-TTY so stdin pipes through cleanly.)

## Development

```sh
python -m pytest            # run the test suite
python -m pytest -q
```

The tests cross-check the optimiser against an exhaustive enumeration of
all permutations and orientation assignments for every input of up to five
reads (including every valid pair of 4-mers), and cover zero-overlap reads,
reverse-complement palindromes and tied optima.

## Layout

```
assembler/
  core.py   # validation, reverse complement, overlap, Held-Karp optimiser
  cli.py    # JSON stdin/file -> JSON stdout command line interface
tests/
  test_assembler.py
Dockerfile, docker-compose.yml
```
