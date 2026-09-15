/*
 * TEMPORARY tool sketch -- not part of the v4 campaign.
 *
 * Manual range-finding plus automatic back-and-forth sweep, for a stage with
 * no end-of-travel sensors. All motion is operator-commanded over USB CDC
 * serial: small jog/nudge commands to explore toward each hard mechanical
 * stop, SETMIN/SETMAX to record what was found, then SWEEP to run
 * automatically between them. The run can be paused and resumed at any time,
 * and optional waypoints make it pause automatically at chosen in-between
 * positions.
 *
 * STEP GPIO5, DIR GPIO6, UART ESP_RX GPIO18 / ESP_TX GPIO17, RGB LED GPIO48,
 * EN/ENN externally grounded -- same verified pin map as the v4 campaign.
 *
 * Chopper mode is forced to SpreadCycle with interpolation off and verified
 * by readback, exactly like the campaign firmware: this project found that
 * the TMC2209's power-on defaults (StealthChop, MicroPlyer interpolation on)
 * make commanded microsteps land at the wrong place, which would corrupt this
 * sketch's own position bookkeeping since there is no independent position
 * feedback here (no encoder, no limit switches) -- the firmware's pulse count
 * IS the only position reference.
 *
 * R_SENSE_OHM is 0.11, the sense resistor value actually fitted on this
 * board. Earlier sketches in this project used 0.03 (a Watterott/eval-board
 * value, not this board's), which under-delivers current by roughly 3x
 * (~122 mA actual against a 360 mA request). That is corrected here: this
 * sketch targets 360 mA RMS and actually delivers close to it, so it has
 * more holding/running torque available than any earlier run -- relevant
 * because this tool is specifically used to jog up against unknown hard
 * stops.
 *
 * MRES is fixed for the whole session (set once at boot, not changeable at
 * runtime) so that every position -- including SETMIN/SETMAX and MARK
 * waypoints captured mid-NUDGE -- stays exactly reachable without rounding.
 * Edit MRES_FIXED below and reflash if a different resolution is wanted.
 *
 * Commands over USB CDC (see HELP for the live list):
 *   J+ [n] / J- [n]   Jog n full steps (default 5), unclamped until a range
 *                      is set -- this is how you explore toward each stop.
 *   NUDGE+ / NUDGE-    Single microstep, for the last careful approach.
 *   SETMIN / SETMAX    Record the current position as a soft limit.
 *   CLEARRANGE         Forget the recorded limits.
 *   MARK               Add the current position as a pause waypoint
 *                      (requires a range, and the position must be inside it).
 *   CLEARMARKS         Forget all waypoints.
 *   GOTO <full_steps>  Move directly to an absolute position (needs a range).
 *   GOTOMM <mm>        Same, specified in millimetres.
 *   GOTOMIN / GOTOMAX  Move directly to the recorded MIN/MAX limit.
 *   SPEED <full_steps/s>  Set the rate used by jog/GOTO/SWEEP (not NUDGE).
 *   SWEEP [cycles]     Run back-and-forth between MIN and MAX. 0 or omitted
 *                      means run until ABORT. Pauses (holds for RESUME) at
 *                      every waypoint it passes, and dwells briefly at each
 *                      end before reversing.
 *   PAUSE / RESUME     Pause any in-progress motion or sweep; resume it.
 *   ABORT              Stop everything now (not resumable -- SWEEP again to
 *                      restart). Also cancels a PAUSE/waypoint hold.
 *   STATUS / HELP
 *
 * Only ABORT, PAUSE, RESUME, STATUS, and HELP are accepted while a move or
 * sweep is in progress; anything else is rejected with "# BUSY" rather than
 * queued, so a jog cannot be issued mid-SWEEP by accident. Pause the sweep
 * (or ABORT it) first if manual jogging is needed.
 */

#include <Arduino.h>
#include <TMCStepper.h>

constexpr uint8_t STEP_PIN = 5;
constexpr uint8_t DIR_PIN = 6;
constexpr uint8_t UART_RX_PIN = 18;
constexpr uint8_t UART_TX_PIN = 17;
constexpr uint8_t LED_PIN = 48;
constexpr uint8_t DRIVER_ADDRESS = 0;
constexpr float R_SENSE_OHM = 0.11F;  // Actual resistor on this board.
constexpr uint32_t DRIVER_BAUD = 115200;
constexpr uint16_t CURRENT_RMS_MA = 360;

constexpr uint16_t MRES_FIXED = 16;  // Fixed for the whole session; see header.
constexpr uint16_t FULL_STEPS_PER_MM = 100;  // 2 mm/rev lead, 200 full steps/rev.

constexpr int32_t JOG_DEFAULT_FULL_STEPS = 5;
constexpr uint32_t NUDGE_RATE_FULL_STEPS_S = 8;
constexpr uint32_t DEFAULT_TRAVEL_RATE_FULL_STEPS_S = 40;
constexpr uint32_t MIN_TRAVEL_RATE_FULL_STEPS_S = 5;
constexpr uint32_t MAX_TRAVEL_RATE_FULL_STEPS_S = 300;
constexpr uint32_t ENDPOINT_DWELL_MS = 1000;
constexpr uint32_t STEP_HIGH_US = 3;
constexpr uint32_t DIR_SETUP_US = 20;
constexpr uint8_t MAX_WAYPOINTS = 8;

HardwareSerial TmcSerial(1);
TMC2209Stepper driver(&TmcSerial, R_SENSE_OHM, DRIVER_ADDRESS);

bool driverReady = false;
bool busy = false;
volatile bool abortRequested = false;
volatile bool pauseRequested = false;

int64_t positionMicrosteps = 0;  // Absolute, 1 unit = 1 microstep at MRES_FIXED.
bool haveMin = false, haveMax = false, rangeSet = false;
int64_t rangeMinCandidate = 0, rangeMaxCandidate = 0;
int64_t rangeMinMicrosteps = 0, rangeMaxMicrosteps = 0;

int64_t waypoints[MAX_WAYPOINTS];
uint8_t waypointCount = 0;

uint32_t activeRateMilliHz = 0;
uint32_t travelRateMilliHz = DEFAULT_TRAVEL_RATE_FULL_STEPS_S * 1000;

void colour(uint8_t red, uint8_t green, uint8_t blue) {
  rgbLedWrite(LED_PIN, red, green, blue);
}

void printHelp() {
  Serial.println("# HELP,J+ [n] / J- [n]        jog n full steps (default 5)");
  Serial.println("# HELP,NUDGE+ / NUDGE-        jog one microstep");
  Serial.println("# HELP,SETMIN / SETMAX        record current position as a limit");
  Serial.println("# HELP,CLEARRANGE             forget recorded limits");
  Serial.println("# HELP,MARK / CLEARMARKS      add/clear a pause waypoint at current position");
  Serial.println("# HELP,GOTO <full_steps>      absolute move (needs a range)");
  Serial.println("# HELP,GOTOMM <mm>            absolute move in millimetres");
  Serial.println("# HELP,GOTOMIN / GOTOMAX      move directly to the recorded limit");
  Serial.println("# HELP,SPEED <full_steps/s>   set jog/GOTO/SWEEP rate");
  Serial.println("# HELP,SWEEP [cycles]         auto back-and-forth, 0=until ABORT");
  Serial.println("# HELP,PAUSE / RESUME         pause/resume in-progress motion");
  Serial.println("# HELP,ABORT                  stop now, not resumable");
  Serial.println("# HELP,STATUS / HELP");
}

void printStatus() {
  Serial.printf(
      "# STATUS,position_full_steps=%lld,position_mm=%.4f,off_grid_microsteps=%lld,"
      "range_set=%u,min_full_steps=%lld,max_full_steps=%lld,waypoints=%u,"
      "speed_full_steps_s=%lu,mres=%u,busy=%u,paused=%u,driver=%u\n",
      static_cast<long long>(positionMicrosteps / MRES_FIXED),
      static_cast<double>(positionMicrosteps) / MRES_FIXED / FULL_STEPS_PER_MM,
      static_cast<long long>(positionMicrosteps % MRES_FIXED),
      static_cast<unsigned>(rangeSet),
      rangeSet ? static_cast<long long>(rangeMinMicrosteps / MRES_FIXED) : 0,
      rangeSet ? static_cast<long long>(rangeMaxMicrosteps / MRES_FIXED) : 0,
      waypointCount, static_cast<unsigned long>(travelRateMilliHz / 1000),
      MRES_FIXED, static_cast<unsigned>(busy), static_cast<unsigned>(pauseRequested),
      static_cast<unsigned>(driverReady));
}

// Called frequently (idle loop, and every ~32 pulses / each dwell tick during
// motion) so ABORT/PAUSE/RESUME/STATUS/HELP always land promptly regardless
// of what else is running. Anything else is rejected while busy rather than
// queued, so a jog can't be issued mid-SWEEP by accident.
void pollSerial() {
  if (!Serial.available()) return;
  String raw = Serial.readStringUntil('\n');
  raw.trim();
  if (raw.length() == 0) return;
  String upper = raw;
  upper.toUpperCase();
  int sp = upper.indexOf(' ');
  String cmd = sp < 0 ? upper : upper.substring(0, sp);
  String args = sp < 0 ? String("") : raw.substring(sp + 1);
  args.trim();

  if (cmd == "ABORT") {
    abortRequested = true;
    pauseRequested = false;
    Serial.println("# ABORT_REQUESTED");
    return;
  }
  if (cmd == "PAUSE") {
    pauseRequested = true;
    Serial.println("# PAUSE_REQUESTED");
    return;
  }
  if (cmd == "RESUME") {
    pauseRequested = false;
    Serial.println("# RESUME_REQUESTED");
    return;
  }
  if (cmd == "STATUS") {
    printStatus();
    return;
  }
  if (cmd == "HELP") {
    printHelp();
    return;
  }

  if (busy) {
    Serial.printf("# BUSY,command_ignored=%s\n", cmd.c_str());
    return;
  }

  if (cmd == "J+" || cmd == "J-") {
    doJog(cmd == "J+" ? 1 : -1, args);
  } else if (cmd == "NUDGE+" || cmd == "NUDGE-") {
    doNudge(cmd == "NUDGE+" ? 1 : -1);
  } else if (cmd == "SETMIN") {
    setLimit(true);
  } else if (cmd == "SETMAX") {
    setLimit(false);
  } else if (cmd == "CLEARRANGE") {
    rangeSet = false;
    haveMin = false;
    haveMax = false;
    Serial.println("# RANGE_CLEARED");
  } else if (cmd == "MARK") {
    addWaypoint();
  } else if (cmd == "CLEARMARKS") {
    waypointCount = 0;
    Serial.println("# WAYPOINTS_CLEARED");
  } else if (cmd == "GOTO") {
    doGoto(args, false);
  } else if (cmd == "GOTOMM") {
    doGoto(args, true);
  } else if (cmd == "GOTOMIN" || cmd == "GOTOMAX") {
    doGotoLimit(cmd == "GOTOMIN");
  } else if (cmd == "SPEED") {
    setSpeed(args);
  } else if (cmd == "SWEEP") {
    uint32_t cycles = args.length() ? static_cast<uint32_t>(args.toInt()) : 0;
    runSweep(cycles);
  } else if (cmd.length() > 0) {
    Serial.printf("# UNKNOWN_COMMAND=%s\n", cmd.c_str());
  }
}

bool waitWhilePaused() {
  while (pauseRequested && !abortRequested) {
    pollSerial();
    delay(5);
  }
  return !abortRequested;
}

// Not compensated for time spent paused inside it: a PAUSE during an
// endpoint dwell shortens the remaining dwell by the real time elapsed.
// Cosmetic only, not safety-relevant.
bool cancellableDwellWithPause(uint32_t durationMs) {
  const uint32_t start = millis();
  while (millis() - start < durationMs) {
    pollSerial();
    if (abortRequested) return false;
    if (pauseRequested && !waitWhilePaused()) return false;
    delay(2);
  }
  return true;
}

bool configureMotionRate(uint32_t fullStepRateMilliHz) {
  const uint64_t pulseRateMilliHz =
      static_cast<uint64_t>(fullStepRateMilliHz) * MRES_FIXED;
  if (pulseRateMilliHz == 0 || pulseRateMilliHz > UINT32_MAX) return false;
  const uint64_t periodUs = 1000000000ULL / pulseRateMilliHz;
  if (periodUs <= STEP_HIGH_US + 2) return false;
  activeRateMilliHz = fullStepRateMilliHz;
  return true;
}

// Updates positionMicrosteps per pulse actually sent, not on completion, so
// position stays correct even when cut short by ABORT (the expected way this
// tool gets used: PAUSE or ABORT right when the stage looks close to a hard
// stop).
bool moveMicrosteps(int32_t signedMicrosteps) {
  if (signedMicrosteps == 0) return true;
  if (activeRateMilliHz == 0) return false;
  const bool positive = signedMicrosteps > 0;
  const uint32_t count = positive
      ? static_cast<uint32_t>(signedMicrosteps)
      : static_cast<uint32_t>(-static_cast<int64_t>(signedMicrosteps));
  const uint64_t pulseRateMilliHz =
      static_cast<uint64_t>(activeRateMilliHz) * MRES_FIXED;
  const uint32_t periodUs = static_cast<uint32_t>(
      (1000000000ULL + pulseRateMilliHz / 2) / pulseRateMilliHz);
  digitalWrite(DIR_PIN, positive ? HIGH : LOW);
  delayMicroseconds(DIR_SETUP_US);
  for (uint32_t pulse = 0; pulse < count; ++pulse) {
    const uint64_t startedUs = esp_timer_get_time();
    digitalWrite(STEP_PIN, HIGH);
    delayMicroseconds(STEP_HIGH_US);
    digitalWrite(STEP_PIN, LOW);
    positionMicrosteps += positive ? 1 : -1;
    while (esp_timer_get_time() - startedUs < periodUs) {
      const uint64_t remaining = periodUs - (esp_timer_get_time() - startedUs);
      if (remaining > 1500) delayMicroseconds(500);
    }
    if ((pulse & 31U) == 31U) {
      pollSerial();
      if (pauseRequested && !waitWhilePaused()) return false;
      if (abortRequested) return false;
    }
  }
  return !abortRequested;
}

bool moveFullSteps(int32_t signedFullSteps, uint32_t rateMilliHz) {
  if (!configureMotionRate(rateMilliHz)) return false;
  const int64_t pulses = static_cast<int64_t>(signedFullSteps) * MRES_FIXED;
  if (pulses > INT32_MAX || pulses < INT32_MIN) return false;
  return moveMicrosteps(static_cast<int32_t>(pulses));
}

// Travels to targetMicrosteps, stopping and holding (waiting for RESUME) at
// every waypoint strictly between the start and the target, in travel order.
bool travelWithWaypoints(int64_t targetMicrosteps) {
  while (true) {
    const int64_t from = positionMicrosteps;
    const bool ascending = targetMicrosteps > from;
    int64_t nextStop = targetMicrosteps;
    bool stopIsWaypoint = false;
    for (uint8_t i = 0; i < waypointCount; ++i) {
      const int64_t w = waypoints[i];
      const bool between =
          ascending ? (w > from && w < targetMicrosteps) : (w < from && w > targetMicrosteps);
      if (!between) continue;
      const bool closer = !stopIsWaypoint || (ascending ? (w < nextStop) : (w > nextStop));
      if (closer) {
        nextStop = w;
        stopIsWaypoint = true;
      }
    }
    const int32_t deltaFullSteps = static_cast<int32_t>((nextStop - from) / MRES_FIXED);
    if (deltaFullSteps != 0) {
      colour(0, 8, 28);
      if (!moveFullSteps(deltaFullSteps, travelRateMilliHz)) return false;
    }
    if (stopIsWaypoint) {
      Serial.printf(
          "# WAYPOINT_HOLD,position_full_steps=%lld (send RESUME to continue)\n",
          static_cast<long long>(positionMicrosteps / MRES_FIXED));
      colour(20, 16, 0);
      pauseRequested = true;
      if (!waitWhilePaused()) return false;
    }
    if (nextStop == targetMicrosteps) break;
  }
  return true;
}

void doJog(int direction, const String &args) {
  long n = args.length() ? args.toInt() : JOG_DEFAULT_FULL_STEPS;
  if (n <= 0) n = JOG_DEFAULT_FULL_STEPS;
  int64_t target = positionMicrosteps + static_cast<int64_t>(direction) * n * MRES_FIXED;
  if (rangeSet) {
    if (target < rangeMinMicrosteps) target = rangeMinMicrosteps;
    if (target > rangeMaxMicrosteps) target = rangeMaxMicrosteps;
  }
  const int32_t deltaFullSteps =
      static_cast<int32_t>((target - positionMicrosteps) / MRES_FIXED);
  if (deltaFullSteps == 0) {
    Serial.println("# JOG_NOOP,at_limit_or_zero");
    return;
  }
  busy = true;
  abortRequested = false;
  colour(0, 20, 0);
  const bool ok = moveFullSteps(deltaFullSteps, travelRateMilliHz);
  busy = false;
  colour(0, 18, 18);
  Serial.printf("# JOG_%s,requested_full_steps=%ld,applied_full_steps=%ld,ok=%u\n",
                direction > 0 ? "FWD" : "REV", n, static_cast<long>(deltaFullSteps),
                static_cast<unsigned>(ok));
  printStatus();
}

void doNudge(int direction) {
  const int64_t target = positionMicrosteps + direction;
  if (rangeSet) {
    if (target < rangeMinMicrosteps) {
      Serial.println("# NUDGE_BLOCKED,at_min");
      return;
    }
    if (target > rangeMaxMicrosteps) {
      Serial.println("# NUDGE_BLOCKED,at_max");
      return;
    }
  }
  busy = true;
  abortRequested = false;
  colour(0, 20, 0);
  bool ok = false;
  // configureMotionRate expects a milliHz-scaled rate (Hz x 1000), matching
  // travelRateMilliHz elsewhere in this file -- NUDGE_RATE_FULL_STEPS_S is a
  // plain Hz constant, so it must be scaled here too or every NUDGE runs
  // 1000x slower than intended (~7.8s instead of ~7.8ms per microstep).
  if (configureMotionRate(NUDGE_RATE_FULL_STEPS_S * 1000)) ok = moveMicrosteps(direction);
  busy = false;
  colour(0, 18, 18);
  Serial.printf("# NUDGE_%s,ok=%u\n", direction > 0 ? "FWD" : "REV", static_cast<unsigned>(ok));
  printStatus();
}

void setLimit(bool isMin) {
  if (isMin) {
    rangeMinCandidate = positionMicrosteps;
    haveMin = true;
  } else {
    rangeMaxCandidate = positionMicrosteps;
    haveMax = true;
  }
  if (haveMin && haveMax) {
    rangeMinMicrosteps = min(rangeMinCandidate, rangeMaxCandidate);
    rangeMaxMicrosteps = max(rangeMinCandidate, rangeMaxCandidate);
    rangeSet = true;
    const double travelMm = static_cast<double>(rangeMaxMicrosteps - rangeMinMicrosteps) /
                            MRES_FIXED / FULL_STEPS_PER_MM;
    Serial.printf("# RANGE_SET,min_full_steps=%lld,max_full_steps=%lld,travel_mm=%.3f\n",
                  static_cast<long long>(rangeMinMicrosteps / MRES_FIXED),
                  static_cast<long long>(rangeMaxMicrosteps / MRES_FIXED), travelMm);
  } else {
    Serial.printf("# %s_CANDIDATE_SET,position_full_steps=%lld (send %s to complete range)\n",
                  isMin ? "MIN" : "MAX",
                  static_cast<long long>(positionMicrosteps / MRES_FIXED),
                  isMin ? "SETMAX" : "SETMIN");
  }
}

void addWaypoint() {
  if (!rangeSet) {
    Serial.println("# MARK_REJECTED,reason=range_not_set");
    return;
  }
  if (positionMicrosteps < rangeMinMicrosteps || positionMicrosteps > rangeMaxMicrosteps) {
    Serial.println("# MARK_REJECTED,reason=outside_range");
    return;
  }
  if (waypointCount >= MAX_WAYPOINTS) {
    Serial.println("# MARK_REJECTED,reason=list_full");
    return;
  }
  waypoints[waypointCount] = positionMicrosteps;
  Serial.printf("# MARK_ADDED,index=%u,position_full_steps=%lld\n", waypointCount,
                static_cast<long long>(positionMicrosteps / MRES_FIXED));
  ++waypointCount;
}

void doGoto(const String &args, bool millimetres) {
  if (!rangeSet) {
    Serial.println("# GOTO_REJECTED,reason=range_not_set");
    return;
  }
  if (args.length() == 0) {
    Serial.println("# GOTO_REJECTED,reason=missing_argument");
    return;
  }
  const long targetFullSteps =
      millimetres ? lround(args.toFloat() * FULL_STEPS_PER_MM) : args.toInt();
  const int64_t targetMicrosteps = static_cast<int64_t>(targetFullSteps) * MRES_FIXED;
  if (targetMicrosteps < rangeMinMicrosteps || targetMicrosteps > rangeMaxMicrosteps) {
    Serial.println("# GOTO_REJECTED,reason=outside_range");
    return;
  }
  busy = true;
  abortRequested = false;
  const bool ok = travelWithWaypoints(targetMicrosteps);
  busy = false;
  colour(0, 18, 18);
  Serial.printf("# GOTO_COMPLETE,ok=%u\n", static_cast<unsigned>(ok));
  printStatus();
}

void doGotoLimit(bool toMin) {
  if (!rangeSet) {
    Serial.println("# GOTO_REJECTED,reason=range_not_set");
    return;
  }
  busy = true;
  abortRequested = false;
  colour(0, 8, 28);
  const bool ok =
      travelWithWaypoints(toMin ? rangeMinMicrosteps : rangeMaxMicrosteps);
  busy = false;
  colour(0, 18, 18);
  Serial.printf("# GOTO_COMPLETE,ok=%u\n", static_cast<unsigned>(ok));
  printStatus();
}

void setSpeed(const String &args) {
  if (args.length() == 0) {
    Serial.printf("# SPEED,full_steps_s=%lu\n",
                  static_cast<unsigned long>(travelRateMilliHz / 1000));
    return;
  }
  long value = args.toInt();
  if (value < static_cast<long>(MIN_TRAVEL_RATE_FULL_STEPS_S))
    value = MIN_TRAVEL_RATE_FULL_STEPS_S;
  if (value > static_cast<long>(MAX_TRAVEL_RATE_FULL_STEPS_S))
    value = MAX_TRAVEL_RATE_FULL_STEPS_S;
  travelRateMilliHz = static_cast<uint32_t>(value) * 1000;
  Serial.printf("# SPEED_SET,full_steps_s=%ld\n", value);
}

bool runSweep(uint32_t cycles) {
  if (!rangeSet) {
    Serial.println("# SWEEP_REJECTED,reason=range_not_set");
    return false;
  }
  busy = true;
  abortRequested = false;
  pauseRequested = false;
  const uint32_t legsTarget = cycles == 0 ? 0 : cycles * 2;
  uint32_t legIndex = 0;
  bool headingToMax = true;
  Serial.printf(
      "# SWEEP_START,cycles=%lu (0=infinite),min_full_steps=%lld,max_full_steps=%lld\n",
      static_cast<unsigned long>(cycles),
      static_cast<long long>(rangeMinMicrosteps / MRES_FIXED),
      static_cast<long long>(rangeMaxMicrosteps / MRES_FIXED));
  bool ok = true;
  while (legsTarget == 0 || legIndex < legsTarget) {
    const int64_t target = headingToMax ? rangeMaxMicrosteps : rangeMinMicrosteps;
    Serial.printf("# SWEEP_LEG_START,leg=%lu,heading=%s\n",
                  static_cast<unsigned long>(legIndex + 1), headingToMax ? "MAX" : "MIN");
    colour(0, 8, 28);
    if (!travelWithWaypoints(target)) {
      ok = false;
      break;
    }
    Serial.printf("# SWEEP_LEG_END,leg=%lu,position_full_steps=%lld\n",
                  static_cast<unsigned long>(legIndex + 1),
                  static_cast<long long>(positionMicrosteps / MRES_FIXED));
    ++legIndex;
    colour(20, 16, 0);
    if (!cancellableDwellWithPause(ENDPOINT_DWELL_MS)) {
      ok = false;
      break;
    }
    headingToMax = !headingToMax;
  }
  busy = false;
  colour(0, 18, 18);
  Serial.printf(ok ? "# SWEEP_COMPLETE,legs=%lu\n" : "# SWEEP_ABORTED,legs=%lu\n",
                static_cast<unsigned long>(legIndex));
  return ok;
}

// See esp32_v4_mres_trajectory_campaign.ino for why a bare CHOPCONF read is
// ambiguous: TMCStepper's UART read gives up after 2 retries and returns 0
// with CRCerror set, which looks exactly like a genuine MRES-code-0 readback.
bool readChopconfChecked(uint32_t &value) {
  for (uint8_t attempt = 0; attempt < 5; ++attempt) {
    const uint32_t candidate = driver.CHOPCONF();
    if (!driver.CRCerror && (candidate & 0x0FUL) != 0UL) {
      value = candidate;
      return true;
    }
    delay(5);
  }
  return false;
}

bool setMresFixed(uint16_t mres) {
  uint8_t code = 8;
  for (uint16_t value = mres; value > 1; value >>= 1) --code;
  constexpr uint32_t MRES_MASK = 0x0F000000UL;
  for (uint8_t attempt = 1; attempt <= 5; ++attempt) {
    uint32_t chopconf = 0;
    if (!readChopconfChecked(chopconf)) {
      delay(10);
      continue;
    }
    driver.CHOPCONF((chopconf & ~MRES_MASK) | (static_cast<uint32_t>(code) << 24));
    uint32_t verify = 0;
    if (!readChopconfChecked(verify)) {
      delay(10);
      continue;
    }
    if (((verify & MRES_MASK) >> 24) == code) {
      Serial.printf("# MRES_OK,mres=%u,code=%u,attempt=%u\n", mres, code, attempt);
      return true;
    }
    delay(10);
  }
  Serial.printf("# MRES_READBACK_FAILED,wrote=%u\n", code);
  return false;
}

bool configureDriver() {
  TmcSerial.begin(DRIVER_BAUD, SERIAL_8N1, UART_RX_PIN, UART_TX_PIN);
  driver.begin();
  driver.toff(5);
  driver.pdn_disable(true);
  driver.mstep_reg_select(true);
  driver.I_scale_analog(false);
  driver.rms_current(CURRENT_RMS_MA);
  driver.ihold(driver.irun());
  driver.iholddelay(0);
  driver.TPOWERDOWN(0);

  driver.en_spreadCycle(true);
  driver.TPWMTHRS(0);
  driver.intpol(false);

  const uint32_t gconf = driver.GCONF();
  const uint32_t chopconf = driver.CHOPCONF();
  const bool spreadCycleOn = ((gconf >> 2) & 1U) != 0U;
  const bool interpolationOff = ((chopconf >> 28) & 1U) == 0U;
  Serial.printf("# CHOPPER_CONFIG,gconf=0x%08lX,chopconf=0x%08lX,en_spreadcycle=%u,intpol=%u\n",
                static_cast<unsigned long>(gconf), static_cast<unsigned long>(chopconf),
                static_cast<unsigned>(spreadCycleOn), static_cast<unsigned>(!interpolationOff));
  if (!spreadCycleOn || !interpolationOff) {
    Serial.println("# CHOPPER_CONFIG_FAILED,refusing_to_run");
    return false;
  }
  const uint32_t drvStatus = driver.DRV_STATUS();
  const bool stealthActive = ((drvStatus >> 30) & 1U) != 0U;
  if (stealthActive) {
    Serial.println("# STEALTHCHOP_STILL_ACTIVE,refusing_to_run");
    return false;
  }
  bool uartOk = false;
  for (uint8_t attempt = 0; attempt < 10; ++attempt) {
    if (driver.test_connection() == 0 && driver.version() == 0x21) {
      uartOk = true;
      break;
    }
    delay(150);
  }
  if (!uartOk) {
    Serial.println("# TMC_UART_FAILED");
    return false;
  }
  Serial.println("# TMC_UART_OK,version=0x21,rx=18,tx=17");
  Serial.printf("# CURRENT,requested_rms_ma=%u,irun_code=%u,vsense=%u,actual_rms_ma=%u\n",
                CURRENT_RMS_MA, driver.irun(), driver.vsense(), driver.rms_current());
  return setMresFixed(MRES_FIXED);
}

void failWithCode(uint8_t code) {
  Serial.printf("# FAILURE_CODE=%u\n", code);
  while (true) {
    for (uint8_t blink = 0; blink < code; ++blink) {
      colour(32, 0, 0);
      delay(250);
      colour(0, 0, 0);
      delay(250);
    }
    delay(1500);
  }
}

void setup() {
  pinMode(STEP_PIN, OUTPUT);
  pinMode(DIR_PIN, OUTPUT);
  digitalWrite(STEP_PIN, LOW);
  digitalWrite(DIR_PIN, LOW);

  colour(24, 24, 24);
  delay(1000);
  colour(0, 0, 0);

  Serial.begin(115200);
  Serial.setTimeout(20);
  const uint32_t serialStart = millis();
  while (!Serial && millis() - serialStart < 3000) delay(10);
  delay(4000);  // Let the USB-CDC capture attach before anything is printed.

  Serial.println("# TMP_V4_MANUAL_RANGE_SWEEP (temporary tool; no motion until commanded)");
  Serial.println("# STEP=GPIO5 DIR=GPIO6 UART_RX=GPIO18 UART_TX=GPIO17 EN_EXTERNALLY_GROUNDED");
  printHelp();

  driverReady = configureDriver();
  if (!driverReady) failWithCode(1);

  colour(0, 18, 18);
  Serial.println("# READY");
  printStatus();
}

void loop() {
  pollSerial();
  delay(5);
}
