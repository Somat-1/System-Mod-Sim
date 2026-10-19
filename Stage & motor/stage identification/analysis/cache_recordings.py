#!/usr/bin/env python3
'''Parse the raw hardware exports into npz caches the plot scripts load.

The three raw CSVs total 1.3 GB and are local-only (see the repository
.gitignore). Everything downstream reads the caches written here, so the raw
files are touched exactly once per campaign.

    IDdata.csv        TwinCAT YT Scope export, EL5101 counter, 1 kHz, 1 channel
    tempTvalm.csv     Dewesoft export, 3 accelerometer axes, 20 kHz
    tempTvalm2.csv    same format, longer recording
'''

import argparse
import io
import re
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE_DIR = HERE.parent
DEFAULT_CACHE = HERE / 'cache'

READ_BLOCK = 64 * 1024 * 1024
DATA_ROW = re.compile(rb'^\d+\t-?\d+')


def _blocks(stream, first_line_handler=None):
    '''Yield whole-line byte blocks, so a row is never split across parses.'''
    header = stream.readline()
    if first_line_handler is not None:
        first_line_handler(header)
    tail = b''
    while True:
        chunk = stream.read(READ_BLOCK)
        if not chunk:
            break
        chunk = tail + chunk
        cut = chunk.rfind(b'\n') + 1
        tail = chunk[cut:]
        if cut:
            yield chunk[:cut]
    if tail.strip():
        yield tail if tail.endswith(b'\n') else tail + b'\n'


def parse_scope(path, cache_dir):
    '''EL5101 counter export -> counts (uint32) plus the scope header metadata.'''
    meta, counts = {}, []
    t0 = time.time()
    with path.open('rb') as stream:
        # The header is a free-form key/tab/value block; data starts at the
        # first line that is "<sample index><tab><counter>".
        while True:
            pos = stream.tell()
            line = stream.readline()
            if not line:
                raise ValueError(f'no sample rows found in {path}')
            if DATA_ROW.match(line):
                stream.seek(pos)
                break
            key, _, value = line.decode('latin-1').strip().partition('\t')
            if key:
                meta.setdefault(key, value.strip())
        tail = b''
        while True:
            chunk = stream.read(READ_BLOCK)
            if not chunk:
                break
            chunk = tail + chunk
            cut = chunk.rfind(b'\n') + 1
            tail, chunk = chunk[cut:], chunk[:cut]
            if not chunk:
                continue
            # The export ends with a literal "EOF" line.
            chunk = chunk.replace(b'EOF', b'')
            block = np.loadtxt(io.BytesIO(chunk), delimiter='\t', dtype=np.int64)
            counts.append(np.atleast_2d(block))

    table = np.vstack(counts)
    index, counter = table[:, 0], table[:, 1]
    if not np.array_equal(index, np.arange(index.size)):
        raise ValueError('sample index is not contiguous; the export has gaps')

    sample_ms = float(meta.get('SampleTime[ms]', 1.0))
    out = cache_dir / f'{path.stem}.npz'
    np.savez_compressed(out, counter=counter.astype(np.int64),
                        sample_ms=sample_ms,
                        channel=meta.get('SymbolName', ''),
                        starttime=meta.get('Starttime of export', ''),
                        endtime=meta.get('Endtime of export', ''))
    print(f'{path.name}: {counter.size} samples @ {1000 / sample_ms:.0f} Hz '
          f'= {counter.size * sample_ms / 1000:.1f} s  ({time.time() - t0:.0f} s) '
          f'-> {out.name}')
    return out


def parse_accel(path, cache_dir):
    '''Dewesoft 3-axis export -> time base plus float32 acceleration columns.'''
    units, rows = [], []
    t0 = time.time()
    with path.open('rb') as stream:
        for block in _blocks(stream, lambda h: units.extend(
                h.decode('latin-1').strip().split(','))):
            rows.append(np.loadtxt(io.BytesIO(block), delimiter=',',
                                   dtype=np.float32))

    table = np.vstack(rows)
    t, accel = table[:, 0].astype(np.float64), table[:, 1:]
    step = np.diff(t[:1000])
    dt = float(np.median(step))
    if not np.allclose(step, dt, rtol=0, atol=dt * 1e-3):
        raise ValueError('non-uniform sample spacing in the first 1000 rows')

    out = cache_dir / f'{path.stem}.npz'
    np.savez_compressed(out, accel=accel.astype(np.float32), dt=dt,
                        duration=float(t[-1]), units=np.array(units[1:]))
    print(f'{path.name}: {accel.shape[0]} samples x {accel.shape[1]} axes @ '
          f'{1 / dt / 1000:.1f} kHz = {t[-1]:.2f} s  ({time.time() - t0:.0f} s) '
          f'-> {out.name}')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    ap.add_argument('--scope', type=Path, default=STAGE_DIR / 'IDdata.csv')
    ap.add_argument('--accel', type=Path, nargs='*', default=[
        STAGE_DIR / 'chirp' / 'tempTvalm.csv',
        STAGE_DIR / 'chirp' / 'tempTvalm2.csv'])
    ap.add_argument('--force', action='store_true',
                    help='re-parse even if the cache file already exists')
    args = ap.parse_args()

    args.cache.mkdir(parents=True, exist_ok=True)
    jobs = [(args.scope, parse_scope)] + [(p, parse_accel) for p in args.accel]
    for path, parser in jobs:
        if not path.exists():
            print(f'{path.name}: missing, skipped')
            continue
        if (args.cache / f'{path.stem}.npz').exists() and not args.force:
            print(f'{path.name}: cached already, skipped')
            continue
        parser(path, args.cache)


if __name__ == '__main__':
    main()
