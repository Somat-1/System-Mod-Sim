(*
    Top-level test harness. Create a TwinCAT PROGRAM named MAIN_Tester and
    paste this declaration/implementation. Assign it to the same task that
    updates the EL70x1 PDOs.

    Operator sequence:
      1. Commission/link the PDOs described in GVL_BeckhoffIO.st.
      2. Set bArmDrive TRUE and verify bDriveReady.
      3. Select TEST_STEPPING, TEST_IRREGULAR_STEPPING, or TEST_CHIRP.
      4. Pulse bStartTest for one PLC scan.
      5. Pulse bAbortAll at any time to freeze the position command.
*)
PROGRAM MAIN_Tester
VAR
    (* Operator controls: expose these in a visualization/watch list. *)
    eSelectedTest                    : E_BeckhoffTestSelection := TEST_IDLE;
    bArmDrive                        : BOOL;
    bStartTest                       : BOOL;
    bAbortAll                        : BOOL;
    bResetDrive                      : BOOL;
    bReduceHoldingTorque             : BOOL;
    bAcknowledgeMarginalChirpSampling: BOOL;
    udiIrregularSeed                 : UDINT := 16#6D2B79F5;
    uiIrregularRandomMoves           : UINT := 16;

    (* Must equal the actual cyclic task period. Use 250 us for chirp. *)
    udiTaskCycleUs                   : UDINT := 250;
    uiTerminalMicrostepsPerFullStep  : UINT := 64;

    fbStepping  : FB_OscillatingStepping;
    fbIrregular : FB_IrregularStepping20mm;
    fbChirp     : FB_Chirp250To1000;

    bDriveHealthy    : BOOL;
    bAnyTestBusy     : BOOL;
    diActualPosition : DINT;
    diSelectedCommand: DINT;

    (* Scope/log channels. Record these together with the external sensor. *)
    uiActiveSegment             : UINT;
    uiSteppingMres              : UINT;
    uiSteppingCycle             : UINT;
    uiIrregularCommandIndex     : UINT;
    lrIrregularTargetOffsetMm   : LREAL;
    lrIrregularMoveDistanceMm   : LREAL;
    lrCommandFrequencyHz        : LREAL;
    lrTestElapsedS              : LREAL;
    diFollowingErrorCounts      : DINT;
    udiActiveErrorId            : UDINT;
END_VAR

(* Convert the terminal's unsigned PDO bit pattern to a signed position. *)
diActualPosition := UDINT_TO_DINT(GVL_BeckhoffIO.udiStmInternalPosition);
bDriveHealthy := GVL_BeckhoffIO.bStmReady AND
                 NOT GVL_BeckhoffIO.bStmError AND
                 NOT GVL_BeckhoffIO.bStmSyncError;

fbStepping(
    bStart := bStartTest AND NOT bAbortAll AND
              (eSelectedTest = TEST_STEPPING) AND bArmDrive,
    bAbort := bAbortAll OR NOT bArmDrive OR
              (eSelectedTest <> TEST_STEPPING),
    bDriveReady := bDriveHealthy,
    bDriveError := GVL_BeckhoffIO.bStmError OR GVL_BeckhoffIO.bStmSyncError,
    diActualPosition := diActualPosition,
    uiTerminalMicrostepsPerFullStep := uiTerminalMicrostepsPerFullStep,
    uiCyclesPerLevel := 15,
    tEndpointDwell := T#5S,
    tLeadSettle := T#5S,
    tInterLevelSettle := T#5S,
    tTailSettle := T#5S);

fbIrregular(
    bStart := bStartTest AND NOT bAbortAll AND
              (eSelectedTest = TEST_IRREGULAR_STEPPING) AND bArmDrive,
    bAbort := bAbortAll OR NOT bArmDrive OR
              (eSelectedTest <> TEST_IRREGULAR_STEPPING),
    bDriveReady := bDriveHealthy,
    bDriveError := GVL_BeckhoffIO.bStmError OR GVL_BeckhoffIO.bStmSyncError,
    diActualPosition := diActualPosition,
    udiTaskCycleUs := udiTaskCycleUs,
    uiTerminalMicrostepsPerFullStep := uiTerminalMicrostepsPerFullStep,
    lrFullStepMm := 0.010,
    lrTotalRangeMm := 20.0,
    lrMinimumAbsTargetMm := 0.5,
    uiRandomMoves := uiIrregularRandomMoves,
    udiSeed := udiIrregularSeed,
    tInitialSettle := T#5S,
    tDwellBetweenMoves := T#5S,
    tMoveTimeout := T#10S,
    udiPositionToleranceCounts := 2);

fbChirp(
    bStart := bStartTest AND NOT bAbortAll AND
              (eSelectedTest = TEST_CHIRP) AND bArmDrive,
    bAbort := bAbortAll OR NOT bArmDrive OR
              (eSelectedTest <> TEST_CHIRP),
    bDriveReady := bDriveHealthy,
    bDriveError := GVL_BeckhoffIO.bStmError OR GVL_BeckhoffIO.bStmSyncError,
    bAcknowledgeMarginalSampling := bAcknowledgeMarginalChirpSampling,
    diActualPosition := diActualPosition,
    udiTaskCycleUs := udiTaskCycleUs,
    uiTerminalMicrostepsPerFullStep := uiTerminalMicrostepsPerFullStep,
    lrAmplitudeFullSteps := 3.0,
    lrStartFrequencyHz := 250.0,
    lrEndFrequencyHz := 1000.0,
    lrSweepDurationS := 96.0,
    lrSettleDurationS := 5.0,
    lrSpeedRangeFullStepsPerSecond := 32000.0);

bAnyTestBusy := fbStepping.bBusy OR fbIrregular.bBusy OR fbChirp.bBusy;

(* Only a running block owns the setpoint. Otherwise hold the terminal's
   reported internal position so selecting a test cannot cause motion. *)
IF fbStepping.bBusy THEN
    diSelectedCommand := fbStepping.diPositionCommand;
ELSIF fbIrregular.bBusy THEN
    diSelectedCommand := fbIrregular.diPositionCommand;
ELSIF fbChirp.bBusy THEN
    diSelectedCommand := fbChirp.diPositionCommand;
ELSE
    diSelectedCommand := diActualPosition;
END_IF

(* Physical output ownership is centralized here. *)
GVL_BeckhoffIO.bStmEnable := bArmDrive;
GVL_BeckhoffIO.bStmReset := bResetDrive;
GVL_BeckhoffIO.bStmReduceTorque := bReduceHoldingTorque AND NOT bAnyTestBusy;
GVL_BeckhoffIO.udiStmPosition := DINT_TO_UDINT(diSelectedCommand);

(* Stable channels for TwinCAT Scope or a logger. *)
IF fbStepping.bBusy THEN
    uiActiveSegment := UINT#1;
    uiSteppingMres := fbStepping.uiCurrentMres;
    uiSteppingCycle := fbStepping.uiCurrentCycle;
    uiIrregularCommandIndex := 0;
    lrIrregularTargetOffsetMm := 0.0;
    lrIrregularMoveDistanceMm := 0.0;
    lrCommandFrequencyHz := 0.0;
    lrTestElapsedS := 0.0;
    diFollowingErrorCounts := fbStepping.diFollowingErrorCounts;
    udiActiveErrorId := fbStepping.udiErrorId;
ELSIF fbIrregular.bBusy THEN
    uiActiveSegment := UINT#2;
    uiSteppingMres := 0;
    uiSteppingCycle := 0;
    uiIrregularCommandIndex := fbIrregular.uiMoveIndex;
    lrIrregularTargetOffsetMm := fbIrregular.lrLastTargetOffsetMm;
    lrIrregularMoveDistanceMm := fbIrregular.lrLastMoveDistanceMm;
    lrCommandFrequencyHz := 0.0;
    lrTestElapsedS := fbIrregular.lrElapsedS;
    diFollowingErrorCounts := fbIrregular.diFollowingErrorCounts;
    udiActiveErrorId := fbIrregular.udiErrorId;
ELSIF fbChirp.bBusy THEN
    uiActiveSegment := fbChirp.uiSegment;
    uiSteppingMres := 0;
    uiSteppingCycle := 0;
    uiIrregularCommandIndex := 0;
    lrIrregularTargetOffsetMm := 0.0;
    lrIrregularMoveDistanceMm := 0.0;
    lrCommandFrequencyHz := fbChirp.lrCurrentFrequencyHz;
    lrTestElapsedS := fbChirp.lrTotalElapsedS;
    diFollowingErrorCounts := diSelectedCommand - diActualPosition;
    udiActiveErrorId := fbChirp.udiErrorId;
ELSE
    uiActiveSegment := 0;
    uiSteppingMres := 0;
    uiSteppingCycle := 0;
    uiIrregularCommandIndex := fbIrregular.uiLogCount;
    lrIrregularTargetOffsetMm := fbIrregular.lrLastTargetOffsetMm;
    lrIrregularMoveDistanceMm := fbIrregular.lrLastMoveDistanceMm;
    lrCommandFrequencyHz := 0.0;
    lrTestElapsedS := 0.0;
    diFollowingErrorCounts := diSelectedCommand - diActualPosition;
    IF fbStepping.bError THEN
        udiActiveErrorId := fbStepping.udiErrorId;
    ELSIF fbIrregular.bError THEN
        udiActiveErrorId := fbIrregular.udiErrorId;
    ELSE
        udiActiveErrorId := fbChirp.udiErrorId;
    END_IF
END_IF
END_PROGRAM