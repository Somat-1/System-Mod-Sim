(*
    Oscillating one-microstep equilibrium test derived from the v5/v6 ESP32
    microstepping-only block. There are no markers and no long trajectories.

    The EL70x1 remains on its fine internal grid (normally 1/64 step). The
    requested v5 MRES ladder is represented by position increments of
    terminal_resolution / requested_MRES counts: 64, 32, 16, 8, 4, 2, 1
    counts for MRES 1, 2, 4, 8, 16, 32, 64 on a 1/64 terminal grid.
*)
FUNCTION_BLOCK FB_OscillatingStepping
VAR_INPUT
    bStart                         : BOOL;
    bAbort                         : BOOL;
    bDriveReady                    : BOOL;
    bDriveError                    : BOOL;
    diActualPosition               : DINT;
    uiTerminalMicrostepsPerFullStep: UINT := 64;
    uiCyclesPerLevel               : UINT := 15;
    tEndpointDwell                 : TIME := T#5S;
    tLeadSettle                    : TIME := T#5S;
    tInterLevelSettle              : TIME := T#5S;
    tTailSettle                    : TIME := T#5S;
END_VAR
VAR_OUTPUT
    bBusy                  : BOOL;
    bDone                  : BOOL;
    bAborted               : BOOL;
    bError                 : BOOL;
    udiErrorId             : UDINT;
    eState                 : E_BeckhoffSequenceState := SEQ_IDLE;
    uiCurrentMres          : UINT;
    uiCurrentCycle         : UINT;
    udiIncrementCounts     : UDINT;
    diPositionCommand      : DINT;
    diFollowingErrorCounts : DINT;
END_VAR
VAR
    fbStartEdge    : R_TRIG;
    fbDwell        : TON;
    usiMresIndex   : USINT;
    diOrigin       : DINT;
    uiInternalState: UINT;
END_VAR

fbStartEdge(CLK := bStart);
bDone := FALSE;
diFollowingErrorCounts := diPositionCommand - diActualPosition;

IF bAbort AND bBusy THEN
    (* Freeze at the reported internal position; do not make an automatic
       recovery move after an emergency/abort request. *)
    diPositionCommand := diActualPosition;
    bBusy := FALSE;
    bAborted := TRUE;
    eState := SEQ_ABORTED;
    uiInternalState := 0;
    fbDwell(IN := FALSE);
ELSIF bBusy AND (bDriveError OR NOT bDriveReady) THEN
    diPositionCommand := diActualPosition;
    bBusy := FALSE;
    bError := TRUE;
    udiErrorId := 16#0001; (* drive lost ready or reported an error *)
    eState := SEQ_ERROR;
    uiInternalState := 0;
    fbDwell(IN := FALSE);
ELSE
    CASE uiInternalState OF
        0: (* Idle *)
            fbDwell(IN := FALSE);
            IF fbStartEdge.Q THEN
                bAborted := FALSE;
                bError := FALSE;
                udiErrorId := 0;

                IF NOT bDriveReady OR bDriveError THEN
                    bError := TRUE;
                    udiErrorId := 16#0101;
                    eState := SEQ_ERROR;
                ELSIF uiTerminalMicrostepsPerFullStep = 0 THEN
                    bError := TRUE;
                    udiErrorId := 16#0102;
                    eState := SEQ_ERROR;
                ELSIF uiCyclesPerLevel = 0 THEN
                    bError := TRUE;
                    udiErrorId := 16#0103;
                    eState := SEQ_ERROR;
                ELSE
                    diOrigin := diActualPosition;
                    diPositionCommand := diOrigin;
                    usiMresIndex := 0;
                    uiCurrentMres := 1;
                    udiIncrementCounts := UINT_TO_UDINT(uiTerminalMicrostepsPerFullStep);
                    uiCurrentCycle := 1;
                    bBusy := TRUE;
                    eState := SEQ_SETTLING;
                    uiInternalState := 10;
                END_IF
            END_IF

        10: (* Initial 5 s settling at the captured origin *)
            fbDwell(IN := TRUE, PT := tLeadSettle);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                uiInternalState := 20;
            END_IF

        20: (* Command one requested-MRES microstep forward *)
            fbDwell(IN := FALSE);
            diPositionCommand := diOrigin + UDINT_TO_DINT(udiIncrementCounts);
            eState := SEQ_RUNNING;
            uiInternalState := 30;

        30: (* Dwell/settle at the forward endpoint *)
            fbDwell(IN := TRUE, PT := tEndpointDwell);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                uiInternalState := 40;
            END_IF

        40: (* Return to the same origin; no accumulating position drift *)
            fbDwell(IN := FALSE);
            diPositionCommand := diOrigin;
            eState := SEQ_RETURNING;
            uiInternalState := 50;

        50: (* Dwell/settle at origin *)
            fbDwell(IN := TRUE, PT := tEndpointDwell);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                IF uiCurrentCycle < uiCyclesPerLevel THEN
                    uiCurrentCycle := uiCurrentCycle + 1;
                    uiInternalState := 20;
                ELSIF usiMresIndex < 6 THEN
                    usiMresIndex := usiMresIndex + 1;
                    CASE usiMresIndex OF
                        1: uiCurrentMres := 2;
                        2: uiCurrentMres := 4;
                        3: uiCurrentMres := 8;
                        4: uiCurrentMres := 16;
                        5: uiCurrentMres := 32;
                        6: uiCurrentMres := 64;
                    END_CASE

                    IF (uiTerminalMicrostepsPerFullStep MOD uiCurrentMres) <> 0 THEN
                        diPositionCommand := diActualPosition;
                        bBusy := FALSE;
                        bError := TRUE;
                        udiErrorId := 16#0104; (* terminal grid cannot represent this level *)
                        eState := SEQ_ERROR;
                        uiInternalState := 0;
                    ELSE
                        udiIncrementCounts := UINT_TO_UDINT(
                            uiTerminalMicrostepsPerFullStep / uiCurrentMres);
                        uiCurrentCycle := 1;
                        eState := SEQ_SETTLING;
                        uiInternalState := 60;
                    END_IF
                ELSE
                    eState := SEQ_SETTLING;
                    uiInternalState := 70;
                END_IF
            END_IF

        60: (* Additional settling between requested MRES levels *)
            fbDwell(IN := TRUE, PT := tInterLevelSettle);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                uiInternalState := 20;
            END_IF

        70: (* Tail settling at the original position *)
            fbDwell(IN := TRUE, PT := tTailSettle);
            IF fbDwell.Q THEN
                fbDwell(IN := FALSE);
                diPositionCommand := diOrigin;
                bBusy := FALSE;
                bDone := TRUE;
                eState := SEQ_COMPLETE;
                uiInternalState := 80;
            END_IF

        80: (* One-scan complete state, then idle while holding origin *)
            eState := SEQ_IDLE;
            uiInternalState := 0;

        ELSE
            diPositionCommand := diActualPosition;
            bBusy := FALSE;
            bError := TRUE;
            udiErrorId := 16#01FF;
            eState := SEQ_ERROR;
            uiInternalState := 0;
    END_CASE
END_IF
END_FUNCTION_BLOCK
