# Firmware for the MRES test — requirements

No sketch here yet. The stage-identification sketch
([`../../firmware/stage_identification`](../../firmware/stage_identification/stage_identification.ino))
is the base to extend: same pins, same PCNT encoder, same plan-upload model.
It cannot run `../sequence/mres_plan_S*.csv` as it stands, for four reasons.

| Needed | Current sketch | Why |
|---|---|---|
| targets in u256 (1/256 FS) | whole full steps | the B3 fine tier ends at 9/16 FS |
| staircase increments of k microsteps | 1 microstep per increment | B3b is 1/16 FS = MRES/16 microsteps per increment |
| MSCNT in telemetry | not read | θ_e for the electrical-cycle harmonics (§6) |
| configurable IHOLD / TPOWERDOWN / intpol | fixed in `configureDriver()` | §2: operating values, recorded |

## Plan row

One `PLAN ADD` per CSV row, positions in u256:

```
PLAN ADD <kind> <target_u256> <speed_fs> <dwell_ms> <mres> <n_inc> <inc_usteps> <inc_dwell_ms> <tag>
```

- `kind` is `move`, `stair` or `set_mres`.
- `set_mres` refuses unless the position is a whole full step **and** MSCNT
  read back is on the full-step grid. It writes CHOPCONF.MRES, reads it back,
  and emits `# MRES,<mres>,mscnt=<n>` before its dwell.
- `stair` checks `n_inc * inc_usteps * (256 / mres) == |target - start|`.

## Part 1

Homing to the two end switches is not specified yet. When it is, it belongs
in this sketch too (`HOME`, reporting both switch positions and the centre),
and it has to finish at x_start on a full step.

## Identity and replies `run_mres.py` waits for

```
PING                 -> # PONG,fw=mres_test,ready=1,driver=1,...
HOLD <m> <tpd> <hd>  -> # HOLD_OK,...
INTPOL <0|1>         -> # INTPOL_OK,intpol=<0|1>
DRIVER               -> # DRIVER,spreadcycle=1,intpol=..,stealth=0,ihold_mult=..,tpowerdown=..
ZERO                 -> # ZERO_OK
PLAN CLEAR / ADD     -> # PLAN_CLEARED / # PLAN_ADD,index=<i>,count=<n>,tag=<tag>
STATUS               -> # STATUS,...,moves=<n>,...
PLAN RUN             -> # RUN_START ... # RUN_COMPLETE,pos=<u256>
```

Any `# ERR` stops the host. The simulator in `run_mres.py` (`SimulatedBoard`)
is the executable form of this list.

## Telemetry

```
# S,<t_ms>,<pos_u256>,<encoder>,<mscnt>,<row>
# MOVE,<row>,<kind>,<tag>,target=<u256>,mres=<m>
# ARRIVED,<row>,pos=<u256>,enc=<counts>,mscnt=<n>
```

## Driver

`HOLD <ihold_multiplier> <tpowerdown> <iholddelay>` and `INTPOL <0|1>`, each
read back and echoed in `# DRIVER`. The host writes them from
`mres_design.DRIVER` before each session and stores the echo in every block
file header.
