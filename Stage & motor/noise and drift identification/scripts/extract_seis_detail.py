#!/usr/bin/env python3
'''One pass over a SEIS CSV, collecting two things the summary caches lack.

  raw snippets   a few seconds of the actual waveform, for the supporting
                 time-domain figure
  block PSDs     the recording split into equal blocks, each given its own
                 Welch estimate, so the spread between blocks shows how
                 stationary the record is

The second is the honest substitute for repeated recordings. It measures
stationarity *within* one record, which is weaker than repeating the
measurement but is what these data can support, and it is the assumption the
table-over-floor ratio rests on.
'''

import argparse
from pathlib import Path

import numpy as np

from compute_seis_welch import WelchAccumulator, active_exclusions
from render_seis_sequence import csv_blocks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--result', type=Path, required=True)
    parser.add_argument('--exclusions', type=Path,
                        default=ROOT / 'analysis' / 'exclusion_windows.csv')
    parser.add_argument('--snippet-start', type=float, default=150.0)
    parser.add_argument('--snippet-seconds', type=float, default=5.0)
    parser.add_argument('--blocks', type=int, default=6)
    parser.add_argument('--nperseg', type=int, default=200000)
    args = parser.parse_args()

    with np.load(args.summary) as summary:
        fs = 1.0 / float(summary['dt'])
        samples = int(summary['samples'])

    excluded = active_exclusions(args.exclusions, args.source.name)
    masked = [(int(a * fs), int(b * fs)) for a, b, _ in excluded]

    edges = np.linspace(0, samples, args.blocks + 1).astype(int)
    accumulators = [WelchAccumulator(fs, args.nperseg)
                    for _ in range(args.blocks)]
    active = 0

    snip_a = int(args.snippet_start * fs)
    snip_b = snip_a + int(args.snippet_seconds * fs)
    snippet = []

    cursor = 0
    for table in csv_blocks(args.source):
        # Column 0 is the timestamp; the accumulator wants the three axes
        # only, while the snippet keeps the time column for plotting.
        table = np.atleast_2d(table)
        values = table[:, 1:4]
        start, end = cursor, cursor + len(values)
        cursor = end

        if start < snip_b and end > snip_a:
            snippet.append(table[max(0, snip_a - start):snip_b - start])

        # Walk the chunk forward, handing each accumulator only the samples
        # inside its own block and never anything inside an exclusion.
        position = start
        while position < end:
            while active + 1 < args.blocks and position >= edges[active + 1]:
                active += 1
            stop = min(end, edges[active + 1]) if active + 1 < len(edges) else end
            piece = values[position - start:stop - start]
            offset = position
            for ma, mb in masked:
                if offset < mb and offset + len(piece) > ma:
                    head = piece[:max(0, ma - offset)]
                    if len(head):
                        accumulators[active].feed(head)
                    accumulators[active].reset_gap()
                    piece = piece[max(0, mb - offset):]
                    offset = max(offset, mb)
            if len(piece):
                accumulators[active].feed(piece)
            position = stop
        print(f'  {cursor:,}/{samples:,} samples', flush=True)

    frequency = None
    spectra, counts = [], []
    for index, accumulator in enumerate(accumulators):
        if accumulator.count == 0:
            continue
        frequency, psd = accumulator.finish()
        spectra.append(psd)
        counts.append(accumulator.count)

    raw = (np.vstack(snippet) if snippet
           else np.zeros((0, 4), dtype=np.float32))
    np.savez_compressed(
        args.result, frequency_hz=frequency,
        block_psd=np.array(spectra), block_counts=np.array(counts),
        block_edges_s=edges / fs, snippet=raw,
        snippet_start_s=args.snippet_start, fs_hz=fs, nperseg=args.nperseg)
    print(f'wrote {args.result}  '
          f'({len(spectra)} blocks, {counts}; snippet {raw.shape})')


if __name__ == '__main__':
    main()
