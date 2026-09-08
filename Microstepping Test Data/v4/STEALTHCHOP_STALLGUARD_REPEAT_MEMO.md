# Memo — repeat v4 with StealthChop and StallGuard

## Scope

The current v4 MRES/trajectory campaign is the baseline motion measurement. It
must not change chopper mode or enable StallGuard. After the baseline dataset is
complete, repeat the same measurement sequence with **StealthChop and
StallGuard enabled together** as a separate experiment and dataset.

Do not retrofit those settings into `run_mres_trajectory_campaign.py`. Clone or
wrap the baseline sequence so the trajectories, order, markers, MRES values,
current, and block names remain comparable while the driver-mode configuration
is explicit and separately versioned.

## Requirements for the repeat

1. Preserve the baseline CSV and IDS recording unchanged.
2. Use the verified ESP/TMC UART path, or first confirm that the dedicated
   controller exposes the required TMC register access. The EVO motion command
   protocol used by the baseline runner does not itself prove that direct TMC
   register reads are available.
3. Before acquisition, explicitly configure StealthChop and read back the
   relevant driver configuration rather than inferring the mode from motor
   sound.
4. Read and log the raw StallGuard result/status continuously with timestamps
   aligned to the motion-event CSV.
5. Recalibrate the StallGuard threshold after the mechanical assembly and motor
   current are fixed. Store the unloaded baseline, chosen threshold, filter or
   debounce length, and ignored acceleration/end-zone windows in the run
   metadata.
6. Keep operation running after a detection. The indicator should be red only
   while the live signal is below the threshold and should return to green when
   load recovers; detections must also be logged as edges/events for analysis.
7. Repeat exactly MRES `1, 4, 16, 32`, the 15-cycle/30 s oscillations, and both
   25 mm command-generation approaches at 27.5, 70, and 200 full steps/s.
8. Retain the same unique separation markers and trigger boundaries. Add a
   dataset-level mode tag such as `STEALTHCHOP_STALLGUARD` so it cannot be
   confused with the baseline run.

## Comparison notes

Compare the baseline and repeat only when mechanics, supply voltage, current,
origin, preload, and acquisition settings are unchanged. Report missed-step or
stall indications per block and also retain the continuous raw StallGuard trace;
a thresholded LED state alone is not sufficient for later analysis.

Status: memo only. No StealthChop or StallGuard configuration is added to the
current v4 runner.
