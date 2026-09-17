(*
    Unnotched high-band chirp inspired by ESP chirp v3.

    Sequence: 5 s settle -> 250 to 1000 Hz linear sweep (96 s) ->
    5 s settle -> 1000 to 250 Hz linear sweep (96 s) -> 5 s settle.

    This block streams direct position setpoints on the EL70x1 internal
    microstep grid. It deliberately contains no notch and never commands the
    known problematic band below 250 Hz.
*)
FUNCTION_BLOCK FB_Chirp250To1000
VAR_INPUT
    bStart                          : BOOL;
    bAbort                          : BOOL;
    bDriveReady                     : BOOL;
    bDriveError                     : BOOL;
    bAcknowledgeMarginalSampling    : BOOL;
    diActualPosition                : DINT;
    udiTaskCycleUs                  : UDINT := 250;
    uiTerminalMicrostepsPerFullStep : UINT := 64;
    lrAmplitudeFullSteps            : LREAL := 3.0;
    lrStartFrequencyHz              : LREAL := 250.0;
    lrEndFrequencyHz                : LREAL := 1000.0;
    lrSweepDurationS                : LREAL := 96.0;
    lrSettleDurationS               : LREAL := 5.0;
    lrSpeedRangeFullStepsPerSecond  : LREAL := 32000.0;
END_VAR
VAR_OUTPUT
    bBusy                         : BOOL;
    bDone                         : BOOL;
    bAborted                      : BOOL;
    bError                        : BOOL;
    bMarginalTiming               : BOOL;
    udiErrorId                    : UDINT;
    eState                        : E_BeckhoffSequenceState := SEQ_IDLE;
    uiSegment                     : UINT; (* 0 idle, 1 pre, 2 up, 3 mid, 4 down, 5 tail *)
    diPositionCommand             : DINT;
    lrCurrentFrequencyHz          : LREAL;
    lrStateElapsedS               : LREAL;
    lrTotalElapsedS               : LREAL;
    lrSamplesPerCycleAtEnd        : LREAL;
    lrRequiredPeakVelocityFullStepS: LREAL;
    uiActiveMicrostepsPerFullStep   : UINT;
END_VAR
VAR
    fbStartEdge      : R_TRIG;
    diOrigin         : DINT;
    uiInternalState  : UINT;
    lrDtS            : LREAL;
    lrSweepRateHzPerS: LREAL;
    lrPhaseRad       : LREAL;
    lrOffsetCounts   : LREAL;
    lrTwoPi          : LREAL := 6.2831853071795864769;
END_VAR

fbStartEdge(CLK := bStart);
bDone := FALSE;

IF udiTaskCycleUs > 0 THEN
    lrDtS := UDINT_TO_LREAL(udiTaskCycleUs) * 0.000001;
ELSE
    lrDtS := 0.0;
END_IF

IF (udiTaskCycleUs > 0) AND (lrEndFrequencyHz > 0.0) THEN
    lrSamplesPerCycleAtEnd := 1000000.0 /
        (UDINT_TO_LREAL(udiTaskCycleUs) * lrEndFrequencyHz);
ELSE
    lrSamplesPerCycleAtEnd := 0.0;
END_IF

lrRequiredPeakVelocityFullStepS :=
    lrTwoPi * lrEndFrequencyHz * lrAmplitudeFullSteps;
bMarginalTiming := (lrSamplesPerCycleAtEnd < 10.0);

IF bAbort AND bBusy THEN
    (* An abort freezes the commanded position at the terminal's reported
       internal position. Normal completion, unlike abort, returns to origin. *)
    diPositionCommand := diActualPosition;
    bBusy := FALSE;
    bAborted := TRUE;
    eState := SEQ_ABORTED;
    uiSegment := 0;
    uiInternalState := 0;
ELSIF bBusy AND (bDriveError OR NOT bDriveReady) THEN
    diPositionCommand := diActualPosition;
    bBusy := FALSE;
    bError := TRUE;
    udiErrorId := 16#0001;
    eState := SEQ_ERROR;
    uiSegment := 0;
    uiInternalState := 0;
ELSE
    CASE uiInternalState OF
        0: (* Idle and validation *)
            IF fbStartEdge.Q THEN
                bAborted := FALSE;
                bError := FALSE;
                udiErrorId := 0;
                lrStateElapsedS := 0.0;
                lrTotalElapsedS := 0.0;
                lrCurrentFrequencyHz := 0.0;

                IF NOT bDriveReady OR bDriveError THEN
                    bError := TRUE;
                    udiErrorId := 16#0201;
                    eState := SEQ_ERROR;
                ELSIF (udiTaskCycleUs <> 250) THEN
                    bError := TRUE;
                    udiErrorId := 16#0202; (* supplied EL70x1 requires the documented 250 us task *)
                    eState := SEQ_ERROR;
                ELSIF (uiTerminalMicrostepsPerFullStep = 0) OR
                      (lrAmplitudeFullSteps <= 0.0) OR
                      (lrSweepDurationS <= 0.0) OR
                      (lrSettleDurationS < 0.0) THEN
                    bError := TRUE;
                    udiErrorId := 16#0203;
                    eState := SEQ_ERROR;
                ELSIF (lrStartFrequencyHz < 250.0) OR
                      (lrEndFrequencyHz > 1000.0) OR
                      (lrEndFrequencyHz <= lrStartFrequencyHz) THEN
                    bError := TRUE;
                    udiErrorId := 16#0204; (* stay inside the requested safe band *)
                    eState := SEQ_ERROR;
                ELSIF lrSamplesPerCycleAtEnd < 4.0 THEN
                    bError := TRUE;
                    udiErrorId := 16#0205; (* too few setpoint samples per period *)
                    eState := SEQ_ERROR;
                ELSIF bMarginalTiming AND NOT bAcknowledgeMarginalSampling THEN
                    bError := TRUE;
                    udiErrorId := 16#0206; (* explicit acknowledgement required *)
                    eState := SEQ_ERROR;
                ELSIF lrSpeedRangeFullStepsPerSecond <= 0.0 THEN
                    bError := TRUE;
                    udiErrorId := 16#0207;
                    eState := SEQ_ERROR;
                ELSIF lrRequiredPeakVelocityFullStepS >
                      (0.9 * lrSpeedRangeFullStepsPerSecond) THEN
                    bError := TRUE;
                    udiErrorId := 16#0208; (* requested sine exceeds speed margin *)
                    eState := SEQ_ERROR;
                ELSE
                    uiActiveMicrostepsPerFullStep := uiTerminalMicrostepsPerFullStep;
                    diOrigin := diActualPosition;
                    diPositionCommand := diOrigin;
                    lrSweepRateHzPerS :=
                        (lrEndFrequencyHz - lrStartFrequencyHz) / lrSweepDurationS;
                    bBusy := TRUE;
                    eState := SEQ_SETTLING;
                    uiSegment := 1;
                    uiInternalState := 10;
                END_IF
            END_IF

        10: (* Pre-sweep settling *)
            diPositionCommand := diOrigin;
            lrCurrentFrequencyHz := 0.0;
            lrStateElapsedS := lrStateElapsedS + lrDtS;
            lrTotalElapsedS := lrTotalElapsedS + lrDtS;
            IF lrStateElapsedS >= lrSettleDurationS THEN
                lrStateElapsedS := 0.0;
                eState := SEQ_RUNNING;
                uiSegment := 2;
                uiInternalState := 20;
            END_IF

        20: (* Linear up-chirp, 250 -> 1000 Hz, no notch *)
            lrCurrentFrequencyHz :=
                lrStartFrequencyHz + lrSweepRateHzPerS * lrStateElapsedS;
            IF lrCurrentFrequencyHz > lrEndFrequencyHz THEN
                lrCurrentFrequencyHz := lrEndFrequencyHz;
            END_IF
            lrPhaseRad := lrTwoPi *
                (lrStartFrequencyHz * lrStateElapsedS +
                 0.5 * lrSweepRateHzPerS * lrStateElapsedS * lrStateElapsedS);
            lrOffsetCounts := lrAmplitudeFullSteps *
                UINT_TO_LREAL(uiActiveMicrostepsPerFullStep) * SIN(lrPhaseRad);
            diPositionCommand := diOrigin + LREAL_TO_DINT(lrOffsetCounts);
            lrStateElapsedS := lrStateElapsedS + lrDtS;
            lrTotalElapsedS := lrTotalElapsedS + lrDtS;
            IF lrStateElapsedS >= lrSweepDurationS THEN
                diPositionCommand := diOrigin;
                lrCurrentFrequencyHz := 0.0;
                lrStateElapsedS := 0.0;
                eState := SEQ_SETTLING;
                uiSegment := 3;
                uiInternalState := 30;
            END_IF

        30: (* Mid-sweep settling at the original position *)
            diPositionCommand := diOrigin;
            lrStateElapsedS := lrStateElapsedS + lrDtS;
            lrTotalElapsedS := lrTotalElapsedS + lrDtS;
            IF lrStateElapsedS >= lrSettleDurationS THEN
                lrStateElapsedS := 0.0;
                eState := SEQ_RUNNING;
                uiSegment := 4;
                uiInternalState := 40;
            END_IF

        40: (* Linear down-chirp, 1000 -> 250 Hz, no notch *)
            lrCurrentFrequencyHz :=
                lrEndFrequencyHz - lrSweepRateHzPerS * lrStateElapsedS;
            IF lrCurrentFrequencyHz < lrStartFrequencyHz THEN
                lrCurrentFrequencyHz := lrStartFrequencyHz;
            END_IF
            lrPhaseRad := lrTwoPi *
                (lrEndFrequencyHz * lrStateElapsedS -
                 0.5 * lrSweepRateHzPerS * lrStateElapsedS * lrStateElapsedS);
            lrOffsetCounts := lrAmplitudeFullSteps *
                UINT_TO_LREAL(uiActiveMicrostepsPerFullStep) * SIN(lrPhaseRad);
            diPositionCommand := diOrigin + LREAL_TO_DINT(lrOffsetCounts);
            lrStateElapsedS := lrStateElapsedS + lrDtS;
            lrTotalElapsedS := lrTotalElapsedS + lrDtS;
            IF lrStateElapsedS >= lrSweepDurationS THEN
                diPositionCommand := diOrigin;
                lrCurrentFrequencyHz := 0.0;
                lrStateElapsedS := 0.0;
                eState := SEQ_SETTLING;
                uiSegment := 5;
                uiInternalState := 50;
            END_IF

        50: (* Tail settling at origin *)
            diPositionCommand := diOrigin;
            lrStateElapsedS := lrStateElapsedS + lrDtS;
            lrTotalElapsedS := lrTotalElapsedS + lrDtS;
            IF lrStateElapsedS >= lrSettleDurationS THEN
                diPositionCommand := diOrigin;
                bBusy := FALSE;
                bDone := TRUE;
                eState := SEQ_COMPLETE;
                uiSegment := 0;
                uiInternalState := 60;
            END_IF

        60: (* One-scan complete state *)
            eState := SEQ_IDLE;
            uiInternalState := 0;

        ELSE
            diPositionCommand := diActualPosition;
            bBusy := FALSE;
            bError := TRUE;
            udiErrorId := 16#02FF;
            eState := SEQ_ERROR;
            uiSegment := 0;
            uiInternalState := 0;
    END_CASE
END_IF
END_FUNCTION_BLOCK
