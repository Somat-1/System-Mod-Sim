# Stepper Isolated — v1 (2026-10-02): rotary encoder bring-up

This version sets up a rotary encoder read by an ESP32-S3. It contains the
wiring record, the ESP32 firmware that streams the encoder state over USB, and
a desktop monitor that plots the signals live, records them and exports CSV.

```
v1_2026-10-02/
├── README.md                          ← this document
├── encoder_monitor.py                 ← live monitor / recorder / exporter (host)
├── requirements.txt
├── firmware/esp32s3_encoder_stream/   ← Arduino sketch for the ESP32-S3
└── recordings/                        ← default export location (created on first save)
```

---

## 1. Connections

| Encoder signal | ESP32-S3 pin | Function in firmware |
|---|---|---|
| A | **GPIO 4** | PCNT quadrature input (edge on A, level on B) |
| B | **GPIO 5** | PCNT quadrature input (edge on B, level on A) |
| Z (index) | **GPIO 6** | GPIO interrupt on the rising edge; latches the count |
| GND | GND | Common ground: required |
| V+ | 3V3 or 5V | Depends on the encoder (see below) |

```
   Rotary encoder                 ESP32-S3
  ┌────────────┐               ┌──────────────┐
  │          A ├───────────────┤ GPIO4        │
  │          B ├───────────────┤ GPIO5        │
  │          Z ├───────────────┤ GPIO6        │
  │        GND ├───────────────┤ GND          │
  │         V+ ├───────────────┤ 3V3 / 5V     │
  └────────────┘               │          USB ├──── host (USB CDC, /dev/cu.usbmodem…)
                               └──────────────┘
```

**Pin choice.** GPIO 4, 5 and 6 are general-purpose pins on the ESP32-S3. They
are not strapping pins (those are 0, 3, 45 and 46), not used for USB (19 and
20), and not used for flash or PSRAM. They are safe to use for this.

**Signal levels: check before powering.** ESP32-S3 inputs are **3.3 V only,
not 5 V tolerant**.
- If the encoder runs from 3.3 V, or has open-collector outputs with no
  external pull-ups, wiring it directly is fine. The firmware turns on the
  internal pull-ups (about 45 kΩ).
- If the encoder runs from 5 V **and has push-pull outputs**, put a level
  shifter or a resistor divider (for example 10 kΩ / 20 kΩ) on A, B and Z.
- If the encoder has open-collector outputs with pull-ups to 5 V on its board,
  remove those pull-ups or move them to 3.3 V.
- If the encoder has differential outputs (A+/A−…), use an RS-422 receiver
  such as an AM26LS32 powered from 3.3 V.

**Encoder details. Still to fill in:**

| Item | Value |
|---|---|
| Model / part number | _TBD_ |
| Lines per rev (PPR) | _TBD_. Counts per rev = 4 × PPR (x4 decoding) |
| Supply voltage | _TBD_ |
| Output type | _TBD_ (open-collector / push-pull / line driver) |
| Level shifting | Voltage divider fitted on the encoder outputs (resistor values _TBD_) |
| Mechanical mounting | _TBD_ (shaft coupling to the stepper, etc.) |

If you don't know the PPR, turn the shaft through more than two revolutions
in the monitor. The monitor measures the counts between Z pulses and shows
them as `CPR≈…`; **Use Z** copies that value into the CPR field.

---

## 2. How the firmware reads the encoder

- **A/B → PCNT in hardware, x4 decoding.** The ESP32-S3 pulse-counter
  peripheral counts every edge on both channels. Counting does not depend on
  the CPU or the USB link, so no counts are lost even at high shaft speed. The
  driver extends the 16-bit hardware counter to 32 bits through watch points
  at ±30000. A 1 µs glitch filter rejects short spikes. This limits the
  maximum edge rate to about 500 kHz, far above what this setup needs.
- **Z → GPIO interrupt.** On each rising edge of Z, the interrupt stores the
  current count (`z_pos`) and increments `z_count`. It also sets a latch
  flag, so an index pulse shorter than one sample period still appears in the
  stream.
- **Sampling.** The firmware samples at a fixed rate: 1000 Hz by default,
  adjustable from 10 to 5000 Hz. Each sample is time-stamped with
  `esp_timer` (µs) and carries a sequence number, so the host can detect
  dropped samples.

> The **count** is exact: it comes from the hardware counter. The **A/B/Z
> levels** are just the pin states at each sample instant. Once the edge rate
> is above the sample rate, the A/B traces alias. They are there to check
> wiring and phasing, not to measure position.

### Serial protocol (USB CDC, any baud)

Data line, one per sample:

```
D,<seq>,<t_us>,<count>,<levels>,<z_count>,<z_pos>
```

| Field | Meaning |
|---|---|
| `seq` | sample counter (gaps = dropped samples) |
| `t_us` | ESP32 time at sampling [µs] |
| `count` | accumulated x4 quadrature count (signed) |
| `levels` | bit0 A, bit1 B, bit2 Z (instantaneous), bit3 Z edge seen since previous sample |
| `z_count` | number of index pulses since boot/zero |
| `z_pos` | count captured at the most recent index pulse |

Lines that start with `#` are info messages and acknowledgements (`#INFO`,
`#ACK`, `#ERR`).

Commands, each ending with a newline:

| Command | Effect |
|---|---|
| `R<hz>` | set the sample rate (10–5000), e.g. `R2000` |
| `Z` | zero the count and the index statistics |
| `S0` / `S1` | stop / start streaming |
| `I` | print the configuration |

Which direction counts as positive depends on the A/B phasing. To reverse it,
swap the A and B wires, or swap the edge actions in `setupPcnt()`.

---

## 3. Build and flash

The board and options are the same as in the rest of the repository:

```sh
cd "Stepper Isolated/v1_2026-10-02/firmware"
arduino-cli compile --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' esp32s3_encoder_stream
arduino-cli upload  --fqbn 'esp32:esp32:esp32s3:USBMode=hwcdc,CDCOnBoot=cdc' \
  --port /dev/cu.usbmodem3101 esp32s3_encoder_stream
```

**Use the board's native USB port.** On this Mac it enumerates as
`/dev/cu.usbmodem3101`. Both flashing and the data stream go through this
port.

**Avoid the UART port.** It goes through a WCH CH343 chip (USB ID `1a86:55d3`)
and enumerates as `/dev/cu.usbmodem5970…`. While the WCH driver
extension is not approved, macOS's generic driver corrupts esptool's data
packets. Every upload then fails with `Failed to write to target RAM (result
was 0107: Checksum error)` or `0105: The format of the received message is
invalid`.

**If the chip stays in the bootloader after an upload** (no `#INFO` line, and
`esptool … chip-id` reports `Staying in bootloader`), either tap **RESET** or
run:

```sh
~/Library/Arduino15/packages/esp32/tools/esptool_py/5.3.0/esptool --chip esp32s3 \
  --port /dev/cu.usbmodem3101 --before no-reset --after watchdog-reset chip-id
```

---

## 4. Monitor application

Requirements: Python 3 with `numpy`, `matplotlib`, `tkinter` and `pyserial`.
On this Mac they are provided by `Stepper Isolated/.venv`, which is
Homebrew Python 3.14 with system site-packages and pyserial added.

**Just run it.** `python encoder_monitor.py` works from any terminal,
including the one built into VS Code. If the script finds it is not running
in that venv, or that `pyexpat` can't load (see §5), it re-launches itself
under `.venv` with `DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib`.

**From VS Code:** in Run and Debug (⇧⌘D), pick **Encoder monitor (ESP32-S3)**
or **Encoder monitor (simulated)** and press F5. Both are defined in
`.vscode/launch.json` at the repository root.

To rebuild the venv on another machine:

```sh
cd "Stepper Isolated"
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib python3 -m venv --system-site-packages .venv
DYLD_LIBRARY_PATH=/opt/homebrew/opt/expat/lib .venv/bin/python -m pip install -r v1_2026-10-02/requirements.txt
```

On a machine whose Python isn't affected, drop the `DYLD_LIBRARY_PATH=`
prefix.

```sh
python encoder_monitor.py                          # choose the port in the window
python encoder_monitor.py --port /dev/cu.usbmodem3101 --cpr 4000
python encoder_monitor.py --simulate               # synthetic encoder, no hardware
```

The window has three plots that share one time axis:
1. **Position**: in counts, or in degrees once a CPR is set.
2. **Speed**: in counts/s, or in rpm once a CPR is set, differenced over 10 ms.
3. **Digital traces**: A, B and Z. Z includes the edge latch, so short index
   pulses stay visible.

The status line shows the count, revolutions, angle within the revolution,
speed, the current A/B/Z levels, the number of index pulses with the measured
CPR, the actual sample rate, and the number of dropped samples.

| Control | Action |
|---|---|
| Port / ↻ / Connect | choose the USB port; `SIMULATE` is always listed |
| Rate [Hz] + Set | set the ESP32 sample rate |
| CPR (x4) / Use Z | counts per revolution for unit conversion; **Use Z** takes it from the index spacing |
| Window [s] | length of the visible time window (up to 60 s) |
| Zero | reset the count on the ESP32 |
| Freeze view | stop redrawing the plots so you can inspect or zoom them; data and recording continue |
| ● Record / ■ Stop | capture every sample from now until Stop |
| Export CSV… | save the last recording (default `recordings/encoder_<timestamp>.csv`) |
| Save plot PNG… | save the current figure |

When you start a new recording or close the window, the monitor asks before
discarding a recording that hasn't been exported.

### Headless recording

To record from a script or over SSH, with no window:

```sh
python encoder_monitor.py --headless --port /dev/cu.usbmodem3101 \
  --rate 2000 --duration 10 --zero --out recordings/spin_test.csv
```

### CSV format

The file begins with `#` metadata lines (source, rate, CPR, start time, sample
count, dropped samples, pins), followed by these columns:

```
host_time_s,seq,t_us,count,A,B,Z,Z_edge,z_count,z_pos
```

For analysis, use `t_us` (ESP32 clock) as the time base. `host_time_s` is
assigned when the USB packet arrives and arrives in bursts.

---

## 5. Status and known issues (2026-10-02)

- The firmware **compiles** with esp32 core 3.3.10. Flash use is 321 kB (24 %).
- The monitor was tested in `--simulate` mode, both headless and in the CSV
  export path.
- **Flashed and running on the board** through the native USB port
  `/dev/cu.usbmodem3101`. The stream shows the expected `#INFO` line and
  about 1 kHz of `D,` lines with consecutive sequence numbers. The encoder has
  not yet been checked while turning: in the first capture the shaft was
  still, with count 0 and A = B = Z = 0.
- **Homebrew Python on this Mac: worked around, not fixed.** Homebrew
  Python 3.13 and 3.14 link `pyexpat` against `/usr/lib/libexpat.1.dylib`,
  which lacks `_XML_SetAllocTrackerActivationThreshold`. As a result,
  matplotlib, `plistlib` and `pip` all fail to import. Homebrew's own `expat`
  2.8.3 has the symbol, so the monitor loads it through `DYLD_LIBRARY_PATH`
  (see §4). Python outside this project is unchanged. A permanent fix would
  need a Homebrew rebuild, or a macOS update that brings a newer system
  libexpat.
- **Verified 2026-10-02:** plain `python3 encoder_monitor.py` re-launches
  correctly. Recording from the board without a window gave 3049 samples in
  3 s with 0 dropped. The GUI opened on `/dev/cu.usbmodem3101` without
  errors.
