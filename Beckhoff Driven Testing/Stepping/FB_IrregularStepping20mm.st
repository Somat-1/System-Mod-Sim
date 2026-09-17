(*
    Reproducible irregular stepping over a centered 20 mm total span.

    Default limits are -10 mm to +10 mm relative to the internal position
    captured at start. Pseudo-random magnitudes are generated from a seed,
    while target signs alternate to force back-and-forth travel. Each move is
    logged in output arrays, followed by a 5 s dwell after the terminal's
    internal position reaches the target. The final logged move is to origin.
*)
FUNCTION_BLOCK FB_IrregularStepping20mm
VAR_INPUT
    bStart                          : BOOL;
    bAbort                          : BOOL;
    bDriveReady                     : BOOL;
    bDriveError                     : BOOL;
    diActualPosition                : DINT;
    udiTaskCycleUs                  : UDINT := 250;
    uiTerminalMicrostepsPerFullStep : UINT := 64;
    lrFullStepMm                    : LREAL := 0.010;
    lrTotalRangeMm                  : LREAL := 20.0;
    lrMinimumAbsTargetMm            : LREAL := 0.5;
    uiRandomMoves                   : UINT := 16;
    udiSeed                         : UDINT := 16#6D2B79F5;
    tInitialSettle                  : TIME := T#5S;
    tDwellBetweenMoves              : TIME := T#5S;
    tMoveTimeout                    : TIME := T#10S;
    udiPositionToleranceCounts      : UDINT := 2;
END_VAR
VAR_OUTPUT
    bBusy                    : BOOL;
    bDone                    : BOOL;
    bAborted                 : BOOL;
    bError                   : BOOL;
    udiErrorId               : UDINT;
    eState                   : E_BeckhoffSequenceState := SEQ_IDLE;
    diPositionCommand        : DINT;
    diFollowingErrorCounts   : DINT;
    uiMoveIndex              : UINT;
    uiTotalCommands          : UINT;
    uiLogCount               : UINT;
    lrElapsedS               : LREAL;
    lrLastTargetOffsetMm     : LREAL;
    lrLastMoveDistanceMm     : LREAL;
    aCommandedTargetOffsetMm : ARRAY[1..31] OF LREAL;
    aCommandedMoveDistanceMm : ARRAY[1..31] OF LREAL;
    aCommandTimeS            : ARRAY[1..31] OF LREAL;
    aCommandedPositionCounts : ARRAY[1..31] OF DINT;
END_VAR
VAR
    fbStartEdge             : R_TRIG;
    fbDwell                 : TON;
    fbMoveTimeout           : TON;
    uiInternalState         : UINT;
    uiRunRandomMoves         : UINT;
    uiRunMicrosteps          : UINT;
    lrRunFullStepMm          : LREAL;
    lrRunMinimumAbsTargetMm  : LREAL;
    uiClearIndex            : UINT;
    udiRngState             : UDINT;
    diOrigin                : DINT;
    diRelativeTargetCounts  : DINT;
    diToleranceCounts       : DINT;
    lrDtS                   : LREAL;
    lrHalfRangeMm           : LREAL;
    lrRandomFraction        : LREAL;
    lrRequestedMagnitudeMm  : LREAL;
    lrRequestedTargetMm     : LREAL;
    lrPreviousTargetOffsetMm: LREAL;
END_VAR

fbStartEdge(CLK := bStart);
bDone := FALSE;
diFollowingErrorCounts := diPositionCommand - diActualPosition;

IF udiTaskCycleUs > 0 THEN
    lrDtS := UDINT_TO_LREAL(udiTaskCycleUs) * 0.000001;
ELSE
    lrDtS := 0.0;
END_IF

IF bBusy THEN
    lrElapsedS := lrElapsedS + lrDtS;
END_IF

IF bAbort AND bBusy THEN
    diPositionCommand := diActualPosition;
    bBusy := FALSE;
    bAborted := TRUE;
    eState := SEQ_ABORTED;
    uiInternalState := 0;
    fbDwell(IN := FALSE);
    fbMoveTimeout(IN := FALSE);
ELSIF bBusy AND (bDriveError OR NOT bDriveReady) THEN
    diPositionCommand := diActualPosition;
    bBusy := FALSE;
    bError := TRUE;
    udiErrorId := 16#0001;
    eState := SEQ_ERROR;
    uiInternalState := 0;
    fbDwell(IN := FALSE);
    fbMoveTimeout(IN := FALSE);
ELSE
    CASE uiInternalState OF
        0: (* Idle and validation *)
            fbDwell(IN := FALSE);
            fbMoveTimeout(IN := FALSE);
            IF fbStartEdge.Q THEN
                bAborted := FALSE;
                bError := FALSE;
                udiErrorId := 0;
                uiLogCount := 0;
                lrElapsedS := 0.0;
                lrLastTargetOffsetMm := 0.0;
                lrLastMoveDistanceMm := 0.0;

                IF NOT bDriveReady OR bDriveError THEN
                    bError := TRUE;
                    udiErrorId := 16#0301;
                    eState := SEQ_ERROR;
                ELSIF (udiTaskCycleUs = 0) OR
                      (uiTerminalMicrostepsPerFullStep = 0) OR
                      (lrFullStepMm <= 0.0) THEN
                    bError := TRUE;
                    udiErrorId := 16#0302;
                    eState := SEQ_ERROR;
                ELSIF (lrTotalRangeMm <= 0.0) OR
                      (lrMinimumAbsTargetMm <= 0.0) OR
                      (lrMinimumAbsTargetMm >= (0.5 * lrTotalRangeMm)) THEN
                    bError := TRUE;
                    udiErrorId := 16#0303;
                    eState := SEQ_ERROR;
                ELSIF (uiRandomMoves = 0) OR (uiRandomMoves > 30) THEN
                    bError := TRUE;
                    udiErrorId := 16#0304;
                    eState := SEQ_ERROR;
                ELSIF udiPositionToleranceCounts > 1000 THEN
                    bError := TRUE;
                    udiErrorId := 16#0305;
                    eState := SEQ_ERROR;
                ELSE
                    FOR uiClearIndex := 1 TO 31 DO
                        aCommandedTargetOffsetMm[uiClearIndex] := 0.0;
                        aCommandedMoveDistanceMm[uiClearIndex] := 0.0;
                        aCommandTimeS[uiClearIndex] := 0.0;
                        aCommandedPositionCounts[uiClearIndex] := 0;
                    END_FOR

                    diOrigin := diActualPosition;
                    diPositionCommand := diOrigin;
                    lrHalfRangeMm := 0.5 * lrTotalRangeMm;
                    uiRunRandomMoves := uiRandomMoves;
                    uiRunMicrosteps := uiTerminalMicrostepsPerFullStep;
                    lrRunFullStepMm := lrFullStepMm;
                    lrRunMinimumAbsTargetMm := lrMinimumAbsTargetMm;
                    lrPreviousTargetOffsetMm := 0.0;
                    uiMoveIndex := 1;
                    uiTotalCommands := uiRunRandomMoves + 1; (* final origin return *)
                    diToleranceCounts := UDINT_TO_DINT(udiPositionToleranceCounts);
                    udiRngState := udiSeed;
                    IF udiRngState = 0 THEN
                        udiRngState := 16#6D2B79F5;
                    END_IF
                    bBusy := TRUE;
                    eState := SEQ_SETTLING;
                    uiInternalState := 10;
                END_IF
            END_IF

        10: (* Initial settling at captured origin *)
            fbDwell(IN := TRUE, PT := tInitialSettle);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                uiInternalState := 20;
            END_IF

        20: (* Generate, quantize, command, and log the next target *)
            fbDwell(IN := FALSE);
            fbMoveTimeout(IN := FALSE);

            IF uiMoveIndex <= uiRunRandomMoves THEN
                (* xorshift32: deterministic, non-zero-seed pseudo-random source *)
                udiRngState := udiRngState XOR SHL(udiRngState, 13);
                udiRngState := udiRngState XOR SHR(udiRngState, 17);
                udiRngState := udiRngState XOR SHL(udiRngState, 5);
                lrRandomFraction :=
                    UDINT_TO_LREAL(udiRngState MOD UDINT#1000000) / 999999.0;
                lrRequestedMagnitudeMm := lrRunMinimumAbsTargetMm +
                    lrRandomFraction * (lrHalfRangeMm - lrRunMinimumAbsTargetMm);

                IF (uiMoveIndex MOD 2) = 1 THEN
                    lrRequestedTargetMm := lrRequestedMagnitudeMm;
                ELSE
                    lrRequestedTargetMm := -lrRequestedMagnitudeMm;
                END_IF
            ELSE
                lrRequestedTargetMm := 0.0; (* final return to origin *)
            END_IF

            diRelativeTargetCounts := LREAL_TO_DINT(
                (lrRequestedTargetMm / lrRunFullStepMm) *`r`n                UINT_TO_LREAL(uiRunMicrosteps));
            diPositionCommand := diOrigin + diRelativeTargetCounts;

            (* Log the actual quantized command, not the unrepresentable request. *)
            lrLastTargetOffsetMm :=
                DINT_TO_LREAL(diRelativeTargetCounts) * lrRunFullStepMm /`r`n                UINT_TO_LREAL(uiRunMicrosteps);
            lrLastMoveDistanceMm :=
                lrLastTargetOffsetMm - lrPreviousTargetOffsetMm;
            uiLogCount := uiMoveIndex;
            aCommandedTargetOffsetMm[uiMoveIndex] := lrLastTargetOffsetMm;
            aCommandedMoveDistanceMm[uiMoveIndex] := lrLastMoveDistanceMm;
            aCommandTimeS[uiMoveIndex] := lrElapsedS;
            aCommandedPositionCounts[uiMoveIndex] := diPositionCommand;

            eState := SEQ_RUNNING;
            uiInternalState := 30;

        30: (* Wait for internal target arrival before starting the dwell *)
            fbMoveTimeout(IN := TRUE, PT := tMoveTimeout);
            IF (diFollowingErrorCounts <= diToleranceCounts) AND
               (diFollowingErrorCounts >= -diToleranceCounts) THEN
                fbMoveTimeout(IN := FALSE);
                eState := SEQ_SETTLING;
                uiInternalState := 40;
            ELSIF fbMoveTimeout.Q THEN
                fbMoveTimeout(IN := FALSE);
                diPositionCommand := diActualPosition;
                bBusy := FALSE;
                bError := TRUE;
                udiErrorId := 16#0306; (* internal target not reached in time *)
                eState := SEQ_ERROR;
                uiInternalState := 0;
            END_IF

        40: (* Five-second dwell after reaching each commanded target *)
            fbDwell(IN := TRUE, PT := tDwellBetweenMoves);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                lrPreviousTargetOffsetMm := lrLastTargetOffsetMm;
                IF uiMoveIndex < uiTotalCommands THEN
                    uiMoveIndex := uiMoveIndex + 1;
                    eState := SEQ_RUNNING;
                    uiInternalState := 20;
                ELSE
                    diPositionCommand := diOrigin;
                    bBusy := FALSE;
                    bDone := TRUE;
                    eState := SEQ_COMPLETE;
                    uiInternalState := 50;
                END_IF
            END_IF

        50: (* One-scan complete state, then idle at origin *)
            eState := SEQ_IDLE;
            uiInternalState := 0;

        ELSE
            diPositionCommand := diActualPosition;
            bBusy := FALSE;
            bError := TRUE;
            udiErrorId := 16#03FF;
            eState := SEQ_ERROR;
            uiInternalState := 0;
    END_CASE
END_IF
END_FUNCTION_BLOCK