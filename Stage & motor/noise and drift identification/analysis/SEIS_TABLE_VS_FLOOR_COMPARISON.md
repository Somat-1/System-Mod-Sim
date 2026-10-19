# SEIS table and floor measurement comparison

## Scope

This note compares the two three-channel SEIS voltage recordings made with the
same acquisition device:

- **Table:** `data/2026-10-15_seis/TvalSeis.csv`
- **Floor:** `data/2026-10-15_seis_floor/TvalSeisfloor.csv`

The comparison is made in recorded voltage units. No voltage-to-acceleration
sensitivity was supplied, so the results must not be interpreted as calibrated
acceleration or as a direct transmissibility measurement.

## Processing held constant

| Setting | Value |
|---|---:|
| Sample rate | approximately 20 kHz |
| Welch window | periodic Hann |
| Window length | 1,048,576 samples (52.43 s) |
| Overlap | 50% (524,288 samples) |
| Frequency resolution | approximately 0.0191 Hz |
| Detrending | mean removed separately from every window |
| Reported comparison band | 0.02–400 Hz |
| Full rendered band | approximately 0.019 Hz–10 kHz |

The table estimate is gap-aware: no Welch window crosses the active 26–30 min
exclusion. The floor measurement has no exclusion.

## Measurement coverage

| Quantity | Table | Floor |
|---|---:|---:|
| Raw duration | 2027.290 s (33.79 min) | 972.249 s (16.20 min) |
| Duration admitted to Welch estimate | 1787.290 s | 972.249 s |
| Excluded duration | 240.000 s | 0 s |
| Welch windows averaged | 65 | 36 |

The floor PSD therefore has fewer independent averages and should be expected
to look somewhat less statistically smooth than the table PSD.

## Time-domain comparison

The quantities below use the three-channel magnitude relative to each
recording's initial baseline. Table statistics omit the active 26–30 min mask.

| Statistic | Table [mV] | Floor [mV] | Floor / table |
|---|---:|---:|---:|
| Median 50 ms RMS magnitude | 0.511 | 2.962 | 5.80× |
| 95th-percentile 50 ms RMS magnitude | 0.811 | 3.715 | 4.58× |
| Maximum 50 ms peak magnitude | 14.952 | 26.606 | 1.78× |

The floor measurement has a markedly higher continuous vibration/noise level,
not merely a larger isolated transient.

## Band-integrated spectral RMS

Each value is the square root of the integrated one-sided PSD in that band and
is therefore expressed as RMS voltage. Channels are not combined.

| Band [Hz] | Table Ch. 1 [mV] | Table Ch. 2 [mV] | Table Ch. 3 [mV] | Floor Ch. 1 [mV] | Floor Ch. 2 [mV] | Floor Ch. 3 [mV] |
|---|---:|---:|---:|---:|---:|---:|
| 0.02–1 | 0.0051 | 0.0075 | 0.0056 | 0.0069 | 0.0044 | 0.0051 |
| 1–10 | 0.1990 | 0.2218 | 0.2053 | 0.0451 | 0.0544 | 0.1750 |
| 10–100 | 0.0821 | 0.1043 | 0.2516 | 1.1102 | 1.0875 | 1.2567 |
| 100–400 | 0.0532 | 0.1339 | 0.2542 | 0.7850 | 0.7819 | 1.9630 |
| **0.02–400 total** | **0.2218** | **0.2794** | **0.4125** | **1.3604** | **1.3405** | **2.3374** |

### Floor-to-table RMS ratio

| Band [Hz] | Channel 1 | Channel 2 | Channel 3 |
|---|---:|---:|---:|
| 0.02–1 | 1.36× | 0.58× | 0.91× |
| 1–10 | 0.23× | 0.25× | 0.85× |
| 10–100 | 13.53× | 10.43× | 4.99× |
| 100–400 | 14.75× | 5.84× | 7.72× |
| **0.02–400 total** | **6.13×** | **4.80×** | **5.67×** |

## Dominant individual spectral lines below 400 Hz

| Recording | Channel 1 | Channel 2 | Channel 3 |
|---|---:|---:|---:|
| Table | 2.098 Hz | 9.232 Hz | 48.981 Hz |
| Floor | 49.019 Hz | 49.019 Hz | 49.019 Hz |

The floor recording contains a coherent line at approximately 49.019 Hz on all
three channels. Its peak PSD is 17.37, 15.40, and 14.38 mV²/Hz on channels
1–3, respectively. The table recording also contains energy near 49 Hz, but it
is dominant only on channel 3 and is much weaker there (0.237 mV²/Hz).

The near-50 Hz coherence is consistent with either mains-related pickup or a
mechanically coupled 50 Hz source. The present recordings alone cannot
distinguish those mechanisms; a controlled electrical isolation or mounting
test would be required.

## Main findings

1. **Very-low-frequency behaviour is broadly comparable.** From 0.02–1 Hz,
   neither placement is uniformly higher across all three channels.
2. **The table is stronger on channels 1–2 from 1–10 Hz.** Floor RMS in this
   band is only about one quarter of the table value on those channels.
3. **The floor is substantially stronger from 10–400 Hz.** Depending on the
   channel and band, its integrated RMS is approximately 5–15 times higher.
4. **The floor result is dominated by a common 49.019 Hz line.** Its appearance
   on all channels suggests a shared excitation or shared electrical coupling.
5. **Across the complete 0.02–400 Hz focus band, the floor is 4.8–6.1 times
   higher in RMS voltage.** This conclusion is consistent with the higher
   median and 95th-percentile time-domain magnitude.

## Figures

### Frequency-domain montage

The large montage keeps the measurements and recorded axes separate. The first
two columns are the SEIS/table recording (`TvalSeis.csv`); the last two are the
floor recording (`TvalSeisfloor.csv`). Rows correspond to channels 1, 2, and 3.
Every panel therefore contains exactly one measurement and one recorded axis.
The faint coloured curve is that axis's estimate, and the heavy black curve is
a smoothed trace derived from the same axis only. Direct FFT amplitude spectral
density and Welch PSD are shown side by side. The FFT columns are not true FRFs
because no independent input/reference channel was recorded.

![Table and floor FFT/PSD montage](plots/seis_frequency_montage.png)

### Table raw-data context

![Table raw-data crop with 26–30 min excluded](plots/seis_04_raw_exclusion_crop.png)

### Floor measurement

![Floor raw sequence](plots/seis_floor_01_recorded_sequence.png)

## Interpretation limits

- The sensor mounting and environmental conditions differ between recordings.
- Voltage has not been converted to acceleration.
- The table estimate omits an operator-confirmed disturbance, whereas all floor
  samples are included.
- The recordings have different durations and numbers of Welch averages.
- These results compare recorded output levels; they do not by themselves
  identify whether a line is mechanical vibration, electrical pickup, or both.
