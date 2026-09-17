(*
    Create a TwinCAT GVL named GVL_BeckhoffIO and paste this declaration.
    Link every variable individually to the EL70x1 PDO named in the comment.

    Required PDO assignment:
      SM2: 0x1602 STM Control + 0x1603 STM Position
      SM3: 0x1A03 STM Status + 0x1A07 STM Internal position
*)
VAR_GLOBAL
    (* EL70x1 -> PLC, 0x1A03 STM Status *)
    bStmReadyToEnable      AT %I* : BOOL;  (* 0x6010:01 *)
    bStmReady              AT %I* : BOOL;  (* 0x6010:02 *)
    bStmWarning            AT %I* : BOOL;  (* 0x6010:03 *)
    bStmError              AT %I* : BOOL;  (* 0x6010:04 *)
    bStmMovingPositive     AT %I* : BOOL;  (* 0x6010:05 *)
    bStmMovingNegative     AT %I* : BOOL;  (* 0x6010:06 *)
    bStmTorqueReduced      AT %I* : BOOL;  (* 0x6010:07 *)
    bStmSyncError          AT %I* : BOOL;  (* 0x6010:0E *)

    (* EL70x1 -> PLC, 0x1A07 STM Internal position *)
    udiStmInternalPosition AT %I* : UDINT; (* 0x6010:14, internal microstep position *)

    (* PLC -> EL70x1, 0x1602 STM Control *)
    bStmEnable             AT %Q* : BOOL;  (* 0x7010:01 *)
    bStmReset              AT %Q* : BOOL;  (* 0x7010:02, use as a pulse *)
    bStmReduceTorque       AT %Q* : BOOL;  (* 0x7010:03 *)

    (* PLC -> EL70x1, 0x1603 STM Position *)
    udiStmPosition         AT %Q* : UDINT; (* 0x7010:11, position setpoint *)
END_VAR