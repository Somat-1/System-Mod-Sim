(* Import or paste this DUT before the function blocks and MAIN_Tester. *)
TYPE E_BeckhoffTestSelection :
(
    TEST_IDLE := 0,
    TEST_STEPPING := 1,
    TEST_CHIRP := 2,
    TEST_IRREGULAR_STEPPING := 3
);
END_TYPE

TYPE E_BeckhoffSequenceState :
(
    SEQ_IDLE := 0,
    SEQ_SETTLING := 10,
    SEQ_RUNNING := 20,
    SEQ_RETURNING := 30,
    SEQ_COMPLETE := 40,
    SEQ_ABORTED := 50,
    SEQ_ERROR := 60
);
END_TYPE