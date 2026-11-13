"""Automated Frequency Response Analyzer (Bode Plotter) for Hantek DSO2D10.

========================================================================
Key Features:
1. Sweep Settings Dialog:
   - Number of points, frequency range, generator amplitude, Auto-V/div on/off,
     estimated sweep duration.
2. Auto-V/div (Auto-Ranging):
   - CH1 set once for the stimulus, CH2 re-ranged per point (1-2-5 grid incl. probe
     factor, read-back, stale first acquisition discarded). V/div stored per point;
     original settings restored after the sweep. ADC clipping is flagged.
3. Trend-Based Parasitic Detection (both sides):
   - Passband anchor = highest flat region; stopbands followed outwards on both sides.
   - Parasitic takeover detected by magnitude rise, phase reversal or ADC floor;
     valid points must lie 10 dB above the extrapolated parasitic level.
   - -3 dB cutoff evaluated within the reliable range only.
4. Sub-Sample Interpolated Zero-Crossing Dimensioning (Delta-t):
   - Computes real rising zero-crossings for both CH1 and CH2 via linear interpolation.
   - Correctly anchors phase lead (Delta-t < 0 / CH2 ahead) and phase lag (Delta-t > 0).
5. Dynamic Axis Scaling:
   - Voltage axis follows the largest measured peak; phase axis follows the reliable range.
6. Bode Plot Export:
   - "Save Bode Plot" button saves only the two Bode diagrams (PNG/PDF/SVG).
"""

from __future__ import annotations

import re
import sys
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.widgets import Button, Slider
import numpy as np
import pyvisa

# Global UI configuration
plt.rcParams["toolbar"] = "None"

TIMEBASE_STEPS: Tuple[float, ...] = (
    2e-9, 5e-9, 10e-9, 20e-9, 50e-9,
    100e-9, 200e-9, 500e-9,
    1e-6, 2e-6, 5e-6, 10e-6, 20e-6, 50e-6,
    100e-6, 200e-6, 500e-6,
    1e-3, 2e-3, 5e-3, 10e-3, 20e-3, 50e-3,
    100e-3, 200e-3, 500e-3,
    1.0,
)

VISA_RESOURCE_ID: str = "USB0::1183::20574::CN5606029103430::0::INSTR"

# Default sweep settings (pre-filled in the start dialog)
DEFAULT_FREQ_START_HZ: float = 10.0
DEFAULT_FREQ_STOP_HZ: float = 25_000_000.0
DEFAULT_NUM_POINTS: int = 100
DEFAULT_DDS_VPP: float = 2.0

# Permitted ranges for the start dialog
FREQ_MIN_HZ: float = 1.0
FREQ_MAX_HZ: float = 25_000_000.0      # DDS upper limit of the DSO2D10
NUM_POINTS_MIN: int = 2
NUM_POINTS_MAX: int = 2000
DDS_VPP_MIN: float = 0.1
DDS_VPP_MAX: float = 7.0

# ADC geometry: 8-bit samples, 25 counts per vertical division (8 div = 200 counts)
ADC_COUNTS_PER_DIV: float = 25.0
ADC_CLIP_MARGIN: int = 1               # samples within this many counts of 0/255 count as clipped
ADC_CLIP_MIN_SAMPLES: int = 3          # tolerate single glitch samples

# Empirical timing model for the duration estimate (approx. 2 min for 100 points)
SWEEP_OVERHEAD_PER_POINT_S: float = 1.15

# Oscilloscope plot: y-limit = factor * largest measured peak (1.72 -> +-1.6 V at 0.93 V peak)
PLOT_Y_HEADROOM: float = 1.72

# ------------------------------------------------------------------------------
# Auto-V/div (auto-ranging)
# ------------------------------------------------------------------------------
DEFAULT_AUTO_VDIV: bool = True
# Device range at 1X probe: 2 mV/div ... 10 V/div (verified at 10X: 20 mV ... 100 V).
# The scope silently ignores values outside the range or off the 1-2-5 grid.
VDIV_MIN_AT_1X: float = 2e-3
VDIV_MAX_AT_1X: float = 10.0
CH1_TARGET_MAX_DIV: float = 6.0        # stimulus may span at most ~6 div on CH1
AUTORANGE_UPPER_COUNTS: int = 175      # > 7 div peak-peak -> one step coarser
AUTORANGE_LOWER_COUNTS: int = 60       # < 2.4 div peak-peak -> one step finer
AUTORANGE_MAX_SWITCHES: int = 4        # per measurement point
# Verified on the DSO2D10: the first acquisition after :CHANnel<n>:SCALe may still
# contain (partly) old-scale data, even after 200 ms. The second one is always valid.
# So after every switch: short wait, discard one acquisition, then measure.
AUTORANGE_POST_SWITCH_MIN_WAIT_S: float = 0.05
AUTORANGE_POST_SWITCH_MAX_WAIT_S: float = 0.6
# Duration estimate surcharge: typical number of range switches x cost per switch
# (set + read-back ~0.7 s, discarded acquisition, re-measurement)
AUTORANGE_ESTIMATE_SWITCHES: int = 8
AUTORANGE_ESTIMATE_COST_PER_SWITCH_S: float = 2.0


def build_vdiv_steps(probe_factor: float) -> List[float]:
    """1-2-5 V/div steps within the device range, including the probe factor."""
    steps = []
    for exponent in range(-4, 4):
        for mantissa in (1.0, 2.0, 5.0):
            value = float(f"{mantissa * 10.0**exponent * probe_factor:.6g}")
            if VDIV_MIN_AT_1X * probe_factor * 0.999 <= value <= VDIV_MAX_AT_1X * probe_factor * 1.001:
                steps.append(value)
    return steps


def nearest_step_index(steps: List[float], vdiv: float) -> int:
    return int(np.argmin([abs(np.log(step / vdiv)) for step in steps]))


def format_vdiv(vdiv: float) -> str:
    return f"{vdiv*1000:g} mV" if vdiv < 1.0 else f"{vdiv:g} V"


def get_optimal_timebase(freq_hz: float) -> float:
    target_timebase = 3.5 / (14.0 * freq_hz)
    return min(TIMEBASE_STEPS, key=lambda tb: abs(tb - target_timebase))


def get_settle_time(freq_hz: float) -> float:
    return max(0.015, min(0.12, 2.5 / freq_hz))


def format_frequency(freq_val: float) -> str:
    if freq_val < 1e3:
        return f"{freq_val:.1f} Hz"
    if freq_val < 1e6:
        return f"{freq_val/1e3:.2f} kHz"
    return f"{freq_val/1e6:.2f} MHz"


_SI_PREFIXES = {"": 1.0, "k": 1e3, "m": 1e6, "meg": 1e6}


def parse_frequency(text: str) -> float:
    """Parse '10', '1,5k', '25M', '2.5 MHz', '100 Hz' into Hz (case-insensitive prefixes)."""
    cleaned = text.strip().replace(",", ".").replace(" ", "")
    cleaned = re.sub(r"(?i)hz$", "", cleaned)
    match = re.fullmatch(r"(?i)([0-9]*\.?[0-9]+(?:e[+-]?[0-9]+)?)(meg|k|m)?", cleaned)
    if not match:
        raise ValueError(f"Invalid frequency: '{text}'")
    return float(match.group(1)) * _SI_PREFIXES[(match.group(2) or "").lower()]


def parse_number(text: str) -> float:
    return float(text.strip().replace(",", "."))


def estimate_sweep_duration(f_start: float, f_stop: float, num_points: int, auto_vdiv: bool = False) -> float:
    freqs = np.geomspace(f_start, f_stop, num_points)
    duration = float(sum(
        SWEEP_OVERHEAD_PER_POINT_S + get_settle_time(f) + 14.0 * get_optimal_timebase(f)
        for f in freqs
    ))
    if auto_vdiv:
        duration += AUTORANGE_ESTIMATE_SWITCHES * AUTORANGE_ESTIMATE_COST_PER_SWITCH_S
    return duration


def ask_sweep_settings() -> Optional[Tuple[float, float, int, float, bool]]:
    """Show the start dialog. Returns (f_start, f_stop, num_points, dds_vpp, auto_vdiv) or None if cancelled."""
    import tkinter as tk
    from tkinter import ttk

    result: Dict[str, Any] = {}

    root = tk.Tk()
    root.title("Bode Analyzer – Sweep Settings")
    root.resizable(False, False)

    frame = ttk.Frame(root, padding=14)
    frame.grid()

    var_points = tk.StringVar(value=str(DEFAULT_NUM_POINTS))
    var_start = tk.StringVar(value=format_frequency(DEFAULT_FREQ_START_HZ))
    var_stop = tk.StringVar(value=format_frequency(DEFAULT_FREQ_STOP_HZ))
    var_vpp = tk.StringVar(value=f"{DEFAULT_DDS_VPP:.1f}")
    var_info = tk.StringVar()
    var_auto = tk.BooleanVar(value=DEFAULT_AUTO_VDIV)
    var_vdiv_hint = tk.StringVar()

    rows = (
        ("Measurement points:", var_points, f"{NUM_POINTS_MIN} … {NUM_POINTS_MAX}"),
        ("Start frequency:", var_start, f"from {format_frequency(FREQ_MIN_HZ)}  (e.g. 10, 1.5k, 2M)"),
        ("Stop frequency:", var_stop, f"up to {format_frequency(FREQ_MAX_HZ)}"),
        ("Generator amplitude [Vpp]:", var_vpp, f"{DDS_VPP_MIN} … {DDS_VPP_MAX}"),
    )
    entries = []
    for row_idx, (label, var, hint) in enumerate(rows):
        ttk.Label(frame, text=label).grid(row=row_idx, column=0, sticky="w", pady=3)
        entry = ttk.Entry(frame, textvariable=var, width=14)
        entry.grid(row=row_idx, column=1, padx=8, pady=3)
        entries.append(entry)
        ttk.Label(frame, text=hint, foreground="gray").grid(row=row_idx, column=2, sticky="w")

    ttk.Checkbutton(
        frame, text="Auto-V/div (adjust vertical scale automatically)", variable=var_auto,
    ).grid(row=len(rows), column=0, columnspan=3, sticky="w", pady=(8, 0))

    ttk.Label(frame, textvariable=var_vdiv_hint, foreground="gray").grid(
        row=len(rows) + 1, column=0, columnspan=3, sticky="w", pady=(2, 2)
    )

    info_label = ttk.Label(frame, textvariable=var_info)
    info_label.grid(row=len(rows) + 2, column=0, columnspan=3, sticky="w", pady=(4, 8))

    def validate() -> Optional[Tuple[float, float, int, float, bool]]:
        auto_vdiv = bool(var_auto.get())
        var_vdiv_hint.set(
            "CH1 is set once, CH2 is re-ranged at every measurement point.\n"
            "The original V/div settings are restored afterwards."
            if auto_vdiv else
            "V/div is taken from the oscilloscope – set it beforehand so that\n"
            "CH1 and CH2 stay within the screen."
        )
        try:
            points_f = parse_number(var_points.get())
            if not points_f.is_integer():
                raise ValueError("Measurement points must be an integer")
            points = int(points_f)
            f_start = parse_frequency(var_start.get())
            f_stop = parse_frequency(var_stop.get())
            vpp = parse_number(var_vpp.get())
        except ValueError as exc:
            var_info.set(f"⚠ {exc}")
            info_label.configure(foreground="red")
            return None

        error = None
        if not NUM_POINTS_MIN <= points <= NUM_POINTS_MAX:
            error = f"Measurement points must be between {NUM_POINTS_MIN} and {NUM_POINTS_MAX}"
        elif not FREQ_MIN_HZ <= f_start <= FREQ_MAX_HZ or not FREQ_MIN_HZ <= f_stop <= FREQ_MAX_HZ:
            error = f"Frequencies must be between {format_frequency(FREQ_MIN_HZ)} and {format_frequency(FREQ_MAX_HZ)}"
        elif f_start >= f_stop:
            error = "Start frequency must be lower than stop frequency"
        elif not DDS_VPP_MIN <= vpp <= DDS_VPP_MAX:
            error = f"Amplitude must be between {DDS_VPP_MIN} and {DDS_VPP_MAX} Vpp"

        if error:
            var_info.set(f"⚠ {error}")
            info_label.configure(foreground="red")
            return None

        duration = estimate_sweep_duration(f_start, f_stop, points, auto_vdiv)
        minutes, seconds = divmod(int(round(duration)), 60)
        var_info.set(f"Estimated sweep duration: approx. {minutes} min {seconds:02d} s")
        info_label.configure(foreground="black")
        return f_start, f_stop, points, vpp, auto_vdiv

    def on_start(_event: Any = None) -> None:
        settings = validate()
        if settings is not None:
            result["settings"] = settings
            root.destroy()

    def on_reset() -> None:
        var_points.set(str(DEFAULT_NUM_POINTS))
        var_start.set(format_frequency(DEFAULT_FREQ_START_HZ))
        var_stop.set(format_frequency(DEFAULT_FREQ_STOP_HZ))
        var_vpp.set(f"{DEFAULT_DDS_VPP:.1f}")
        var_auto.set(DEFAULT_AUTO_VDIV)

    for var in (var_points, var_start, var_stop, var_vpp, var_auto):
        var.trace_add("write", lambda *_: validate())

    buttons = ttk.Frame(frame)
    buttons.grid(row=len(rows) + 3, column=0, columnspan=3, sticky="e")
    ttk.Button(buttons, text="Defaults", command=on_reset).grid(row=0, column=0, padx=4)
    ttk.Button(buttons, text="Cancel", command=root.destroy).grid(row=0, column=1, padx=4)
    ttk.Button(buttons, text="Start sweep", command=on_start).grid(row=0, column=2, padx=4)

    root.bind("<Return>", on_start)
    root.bind("<Escape>", lambda _e: root.destroy())
    validate()
    entries[0].focus_set()
    root.mainloop()

    return result.get("settings")


def detect_clipping(raw_offset_binary: np.ndarray) -> bool:
    """True if the waveform hits the ADC rails (signal exceeds the V/div range)."""
    at_rail = (raw_offset_binary <= ADC_CLIP_MARGIN) | (raw_offset_binary >= 255 - ADC_CLIP_MARGIN)
    return int(np.count_nonzero(at_rail)) >= ADC_CLIP_MIN_SAMPLES


def fetch_waveform_raw(scope_inst: pyvisa.resources.MessageBasedResource) -> bytes:
    samples_total = -1
    samples_data: Optional[bytearray] = None
    samples_received = 0

    for _ in range(6):
        try:
            scope_inst.write("PRIVate:WAVeform:DATA:ALL?")
            raw_payload = scope_inst.read_raw()
        except Exception:
            continue

        if not raw_payload or chr(raw_payload[0]) != "#" or len(raw_payload) < 29:
            time.sleep(0.015)
            continue

        total_samples = int(raw_payload[11:20].decode())
        current_offset = int(raw_payload[20:29].decode())

        if samples_total == -1:
            samples_total = total_samples
            samples_data = bytearray(samples_total)

        chunk_length = len(raw_payload) - 128
        if chunk_length > 0 and samples_data is not None:
            samples_data[current_offset:current_offset + chunk_length] = raw_payload[128:128 + chunk_length]
            samples_received += chunk_length

        if samples_total > 0 and samples_received >= samples_total:
            return bytes(samples_data)

    return b""


RAW_BLOCK_LEN: int = 2000
RAW_MIN_BYTES: int = 4 * RAW_BLOCK_LEN


def split_channels(raw_bytes: bytes) -> Tuple[np.ndarray, np.ndarray]:
    """De-interleave the raw buffer into CH1/CH2 offset-binary samples (0..255)."""
    raw_buffer = np.frombuffer(raw_bytes, dtype=np.uint8)
    blk = RAW_BLOCK_LEN
    ch1 = np.concatenate([raw_buffer[0:blk], raw_buffer[2*blk:3*blk]]) ^ 0x80
    ch2 = np.concatenate([raw_buffer[blk:2*blk], raw_buffer[3*blk:4*blk]]) ^ 0x80
    return ch1, ch2


def acquire_channels(
    scope_inst: pyvisa.resources.MessageBasedResource, discard_first: bool = False
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Fetch one acquisition. With discard_first, one acquisition is thrown away first
    (needed after a V/div switch, see AUTORANGE_POST_SWITCH_* notes)."""
    if discard_first:
        fetch_waveform_raw(scope_inst)
    raw_bytes = fetch_waveform_raw(scope_inst)
    if len(raw_bytes) < RAW_MIN_BYTES:
        time.sleep(0.02)
        raw_bytes = fetch_waveform_raw(scope_inst)
        if len(raw_bytes) < RAW_MIN_BYTES:
            return None
    return split_channels(raw_bytes)


def set_vdiv(scope_inst: pyvisa.resources.MessageBasedResource, channel: int, vdiv: float) -> float:
    """Set V/div and return the value actually applied by the scope (read back)."""
    scope_inst.write(f":CHANnel{channel}:SCALe {vdiv:g}")
    return float(scope_inst.query(f":CHANnel{channel}:SCALe?").strip())


def wait_after_vdiv_switch(timebase_s: float) -> None:
    time.sleep(float(np.clip(2.0 * 14.0 * timebase_s,
                             AUTORANGE_POST_SWITCH_MIN_WAIT_S, AUTORANGE_POST_SWITCH_MAX_WAIT_S)))


def autorange_direction(samples_offset_binary: np.ndarray) -> int:
    """+1 = one step coarser, -1 = one step finer, 0 = range is fine."""
    peak_to_peak_counts = int(np.ptp(samples_offset_binary))
    if detect_clipping(samples_offset_binary) or peak_to_peak_counts > AUTORANGE_UPPER_COUNTS:
        return +1
    if peak_to_peak_counts < AUTORANGE_LOWER_COUNTS:
        return -1
    return 0


def compute_spectral_point(
    ch1_voltage: np.ndarray,
    ch2_voltage: np.ndarray,
    sample_interval_s: float,
    frequency_hz: float,
) -> Tuple[float, float, float]:
    num_samples = len(ch1_voltage)
    window = np.hanning(num_samples)
    w_ch1 = ch1_voltage * window
    w_ch2 = ch2_voltage * window

    time_vector = np.arange(num_samples) * sample_interval_s
    carrier_reference = np.exp(-1j * 2.0 * np.pi * frequency_hz * time_vector)

    window_norm = np.sum(window)
    spectral_u1 = 2.0 * np.sum(w_ch1 * carrier_reference) / window_norm
    spectral_u2 = 2.0 * np.sum(w_ch2 * carrier_reference) / window_norm

    amp1 = float(np.abs(spectral_u1))
    amp2 = float(np.abs(spectral_u2))

    gain_db = float(20.0 * np.log10(max(amp2, 1e-6) / max(amp1, 1e-6)))
    phase_deg = float(np.rad2deg(np.angle(spectral_u2) - np.angle(spectral_u1)))
    phase_deg = (phase_deg + 180.0) % 360.0 - 180.0

    return gain_db, phase_deg, amp2


# ==============================================================================
# Sweep Configuration Dialog
# ==============================================================================

user_settings = ask_sweep_settings()
if user_settings is None:
    print("[ABORT] Sweep cancelled by user.")
    sys.exit(0)

SWEEP_FREQ_START_HZ, SWEEP_FREQ_STOP_HZ, SWEEP_NUM_POINTS, DDS_STIMULUS_VPP, AUTO_VDIV = user_settings
print(
    f"[CONFIG] {SWEEP_NUM_POINTS} points, {format_frequency(SWEEP_FREQ_START_HZ)} – "
    f"{format_frequency(SWEEP_FREQ_STOP_HZ)}, stimulus {DDS_STIMULUS_VPP:.2f} Vpp, "
    f"Auto-V/div {'ON' if AUTO_VDIV else 'OFF'}"
)


# ==============================================================================
# Hardware Initialization
# ==============================================================================

resource_manager = pyvisa.ResourceManager("@py")
scope: pyvisa.resources.MessageBasedResource = resource_manager.open_resource(VISA_RESOURCE_ID)
scope.timeout = 2500
scope.query_delay = 0.05
scope.write_termination = "\n"
scope.read_termination = "\n"

device_identity = scope.query("*IDN?").strip()
print(f"[HW-INIT] Connected to: {device_identity}")

scope.write(f":DDS:TYPE SINE;:DDS:AMP {DDS_STIMULUS_VPP:.3f};:DDS:SWIT ON")
time.sleep(0.1)

# Vertical scale as set by the user. The returned V/div already includes the
# probe attenuation (1X/10X). These values are restored after an auto-ranged sweep.
original_vdiv_ch1 = float(scope.query(":CHANnel1:SCALe?").strip())
original_vdiv_ch2 = float(scope.query(":CHANnel2:SCALe?").strip())
vdiv_ch1 = original_vdiv_ch1
vdiv_ch2 = original_vdiv_ch2
print(f"[HW-INIT] Vertical scale read from scope: CH1 = {format_vdiv(vdiv_ch1)}/div, CH2 = {format_vdiv(vdiv_ch2)}/div")

vdiv_steps_ch1: List[float] = []
vdiv_steps_ch2: List[float] = []
if AUTO_VDIV:
    probe_ch1 = float(scope.query(":CHANnel1:PROBe?").strip())
    probe_ch2 = float(scope.query(":CHANnel2:PROBe?").strip())
    vdiv_steps_ch1 = build_vdiv_steps(probe_ch1)
    vdiv_steps_ch2 = build_vdiv_steps(probe_ch2)
    print(
        f"[AUTO-VDIV] Probe CH1 {probe_ch1:g}X, CH2 {probe_ch2:g}X | CH2 range "
        f"{format_vdiv(vdiv_steps_ch2[0])} … {format_vdiv(vdiv_steps_ch2[-1])}/div"
    )
else:
    stimulus_screen_div = DDS_STIMULUS_VPP / vdiv_ch1
    if stimulus_screen_div > 8.0:
        print(
            f"[WARN] Stimulus {DDS_STIMULUS_VPP:.2f} Vpp spans {stimulus_screen_div:.1f} div on CH1 "
            f"(screen = 8 div). Increase CH1 V/div to avoid clipping."
        )


def setup_ch1_autorange() -> None:
    """Pick the finest CH1 step where the stimulus spans at most CH1_TARGET_MAX_DIV,
    then verify with a test acquisition (coarser if it still clips / exceeds 7 div)."""
    global vdiv_ch1
    target_vdiv = DDS_STIMULUS_VPP / CH1_TARGET_MAX_DIV
    candidates = [step for step in vdiv_steps_ch1 if step >= target_vdiv * 0.999]
    vdiv_ch1 = set_vdiv(scope, 1, candidates[0] if candidates else vdiv_steps_ch1[-1])

    test_freq = SWEEP_FREQ_START_HZ
    test_timebase = get_optimal_timebase(test_freq)
    scope.write(f":TIMebase:SCALe {test_timebase:.6e};:DDS:FREQ {test_freq:.2f}")
    time.sleep(get_settle_time(test_freq))
    wait_after_vdiv_switch(test_timebase)

    for _ in range(AUTORANGE_MAX_SWITCHES + 1):
        acquisition = acquire_channels(scope, discard_first=True)
        if acquisition is None:
            print("[AUTO-VDIV] CH1 test acquisition failed – keeping the calculated V/div.")
            return
        ch1_samples, _ = acquisition
        ptp_counts = int(np.ptp(ch1_samples))
        print(
            f"[AUTO-VDIV] CH1 = {format_vdiv(vdiv_ch1)}/div: test measurement "
            f"{ptp_counts * vdiv_ch1 / ADC_COUNTS_PER_DIV:.2f} Vpp ({ptp_counts / ADC_COUNTS_PER_DIV:.1f} div)"
        )
        if autorange_direction(ch1_samples) <= 0:
            return
        next_idx = nearest_step_index(vdiv_steps_ch1, vdiv_ch1) + 1
        if next_idx >= len(vdiv_steps_ch1):
            print("[AUTO-VDIV] CH1 largest V/div reached – signal may clip.")
            return
        vdiv_ch1 = set_vdiv(scope, 1, vdiv_steps_ch1[next_idx])
        wait_after_vdiv_switch(test_timebase)


def autorange_ch2(
    acquisition: Tuple[np.ndarray, np.ndarray], timebase_s: float
) -> Tuple[Optional[Tuple[np.ndarray, np.ndarray]], int]:
    """Adjust CH2 V/div step by step until the signal spans ~2.4 … 7 div.
    Returns the acquisition taken with the final V/div and the number of switches."""
    global vdiv_ch2
    switches = 0
    last_direction = 0
    while switches < AUTORANGE_MAX_SWITCHES:
        direction = autorange_direction(acquisition[1])
        # Stop if in range, or if the direction reverses (hysteresis / no ping-pong)
        if direction == 0 or (last_direction != 0 and direction != last_direction):
            break
        next_idx = nearest_step_index(vdiv_steps_ch2, vdiv_ch2) + direction
        if not 0 <= next_idx < len(vdiv_steps_ch2):
            break
        applied = set_vdiv(scope, 2, vdiv_steps_ch2[next_idx])
        if abs(applied - vdiv_ch2) <= 1e-9 * max(applied, vdiv_ch2):
            break  # scope refused the step (range limit)
        vdiv_ch2 = applied
        switches += 1
        last_direction = direction
        wait_after_vdiv_switch(timebase_s)
        new_acquisition = acquire_channels(scope, discard_first=True)
        if new_acquisition is None:
            return None, switches
        acquisition = new_acquisition
    return acquisition, switches


sweep_frequencies = np.geomspace(SWEEP_FREQ_START_HZ, SWEEP_FREQ_STOP_HZ, SWEEP_NUM_POINTS)
sweep_records: List[Dict[str, Any]] = []

raw_freqs: List[float] = []
raw_mag_db: List[float] = []
raw_phase_deg: List[float] = []
raw_amps_ch2: List[float] = []
raw_vdiv_ch2: List[float] = []
clipping_events: List[Tuple[float, str]] = []
total_vdiv_switches = 0


# ==============================================================================
# Frequency Sweep Acquisition Loop
# ==============================================================================

sweep_start_timestamp = time.time()

try:
    if AUTO_VDIV:
        setup_ch1_autorange()
        # CH2 starts at the CH1 range (passband gain ~0 dB is the common case);
        # afterwards every point starts with the V/div of the previous point.
        vdiv_ch2 = set_vdiv(scope, 2, vdiv_steps_ch2[nearest_step_index(vdiv_steps_ch2, vdiv_ch1)])

    print(f"\n[SWEEP] Commencing automated sweep ({SWEEP_NUM_POINTS} points)...")

    for point_idx, current_freq in enumerate(sweep_frequencies):
        timebase_val = get_optimal_timebase(current_freq)
        scope.write(f":TIMebase:SCALe {timebase_val:.6e};:DDS:FREQ {current_freq:.2f}")

        time.sleep(get_settle_time(current_freq))

        # The initial CH2 switch above needs the same stale-data handling as any other switch
        acquisition = acquire_channels(scope, discard_first=(AUTO_VDIV and point_idx == 0))
        if acquisition is None:
            continue

        if AUTO_VDIV:
            acquisition, point_switches = autorange_ch2(acquisition, timebase_val)
            total_vdiv_switches += point_switches
            if acquisition is None:
                continue

        ch1_offset_binary, ch2_offset_binary = acquisition
        scale_factor_ch1 = vdiv_ch1 / ADC_COUNTS_PER_DIV
        scale_factor_ch2 = vdiv_ch2 / ADC_COUNTS_PER_DIV

        clipped_ch1 = detect_clipping(ch1_offset_binary)
        clipped_ch2 = detect_clipping(ch2_offset_binary)
        if clipped_ch1 or clipped_ch2:
            clipped_channels = " + ".join(
                name for name, flag in (("CH1", clipped_ch1), ("CH2", clipped_ch2)) if flag
            )
            clipping_events.append((current_freq, clipped_channels))

        ch1_centered = ch1_offset_binary.astype(float)
        ch2_centered = ch2_offset_binary.astype(float)

        ch1_v = (ch1_centered - np.mean(ch1_centered)) * scale_factor_ch1
        ch2_v = (ch2_centered - np.mean(ch2_centered)) * scale_factor_ch2

        if timebase_val <= 500e-9:
            sample_interval = 2.0e-9
            onscreen_samples = int(round((14.0 * timebase_val) / sample_interval))
            zoom_limit = int(np.clip(onscreen_samples, 25, len(ch1_v)))
        else:
            sample_interval = (14.0 * timebase_val) / len(ch1_v)
            zoom_limit = len(ch1_v)

        gain, phase, amp_ch2 = compute_spectral_point(ch1_v, ch2_v, sample_interval, current_freq)

        raw_freqs.append(current_freq)
        raw_mag_db.append(gain)
        raw_phase_deg.append(phase)
        raw_amps_ch2.append(amp_ch2)
        raw_vdiv_ch2.append(vdiv_ch2)

        ch1_onscreen = ch1_v[:zoom_limit]
        ch2_onscreen = ch2_v[:zoom_limit]
        time_vector = np.arange(zoom_limit) * sample_interval
        total_time_s = time_vector[-1]

        if total_time_s < 1e-6:
            unit_multiplier, unit_label = 1e9, "ns"
        elif total_time_s < 1e-3:
            unit_multiplier, unit_label = 1e6, "µs"
        elif total_time_s < 1.0:
            unit_multiplier, unit_label = 1e3, "ms"
        else:
            unit_multiplier, unit_label = 1.0, "s"

        sweep_records.append({
            "freq_hz": current_freq,
            "timebase_s": timebase_val,
            "time_axis_scaled": time_vector * unit_multiplier,
            "total_time_scaled": total_time_s * unit_multiplier,
            "time_unit": unit_label,
            "unit_multiplier": unit_multiplier,
            "ch1_waveform": ch1_onscreen,
            "ch2_waveform": ch2_onscreen,
            "vpp_ch1": float(np.ptp(ch1_v)),
            "vpp_ch2": float(np.ptp(ch2_v)),
            "gain_db": gain,
            "phase_deg": phase,
            "clipped_ch1": clipped_ch1,
            "clipped_ch2": clipped_ch2,
            "vdiv_ch1": vdiv_ch1,
            "vdiv_ch2": vdiv_ch2,
        })

        freq_str = (
            f"{current_freq:6.1f} Hz" if current_freq < 1e3 else (
                f"{current_freq/1e3:6.2f} kHz" if current_freq < 1e6 else f"{current_freq/1e6:6.2f} MHz"
            )
        )
        print(
            f"[{point_idx+1:3d}/{len(sweep_frequencies)}] f = {freq_str:>10s} | Gain = {gain:6.1f} dB"
            f" | CH2 {format_vdiv(vdiv_ch2):>7s}/div   ",
            end="\r",
        )

except Exception as err:
    print(f"\n[ERROR] Acquisition fault: {err}")
finally:
    if AUTO_VDIV:
        try:
            scope.write(f":CHANnel1:SCALe {original_vdiv_ch1:g}")
            scope.write(f":CHANnel2:SCALe {original_vdiv_ch2:g}")
            print(
                f"\n[AUTO-VDIV] {total_vdiv_switches} CH2 range switches. Restored V/div: "
                f"CH1 {format_vdiv(original_vdiv_ch1)}/div, CH2 {format_vdiv(original_vdiv_ch2)}/div"
            )
        except Exception as restore_err:
            print(f"\n[WARN] Could not restore original V/div settings: {restore_err}")
    try:
        scope.write(":KEY:LOCK DISABLE")
    except Exception:
        pass
    scope.close()

elapsed_duration = time.time() - sweep_start_timestamp
print(f"\n[SWEEP] Finished. {len(sweep_records)} points acquired in {elapsed_duration:.1f} s.")

if not sweep_records:
    print("[FATAL] No measurement records captured. Halting.")
    sys.exit(1)

if clipping_events:
    clip_hint = (
        "auto-ranging could not avoid it (range limit or too many switches)"
        if AUTO_VDIV else "increase V/div on the scope"
    )
    print(f"[WARN] Signal clipped (ADC limit) at {len(clipping_events)} point(s) – {clip_hint}:")
    for clip_freq, clip_channels in clipping_events:
        print(f"         {format_frequency(clip_freq):>11s}: {clip_channels}")


# ==============================================================================
# Analytical Post-Processing: Topology & Boundary Classification
# ==============================================================================

freqs_array = np.array(raw_freqs)
mag_db_array = np.array(raw_mag_db)
phase_deg_array = np.array(raw_phase_deg)
amps2_array = np.array(raw_amps_ch2)
unwrapped_phase_array = np.rad2deg(np.unwrap(np.deg2rad(phase_deg_array)))

# --- Bidirectional Parasitic Limit Detection (trend-based) ---
# 1. Passband anchor: the highest point that is also FLAT. A parasitic rise can reach
#    passband level at the band edge (e.g. +0.5 dB at 25 MHz behind a 1st-order RC),
#    but it is steep there, while a real passband (or resonance peak) has slope ~0.
# 2. From the anchor the curve is followed outwards to BOTH sides (low-pass: right,
#    high-pass: left, band-pass: both). Within the stopband, parasitic coupling
#    (ground inductance, crosstalk, ESL) is detected where
#      - the magnitude rises again by more than PARASITIC_MAG_RISE_DB above its minimum,
#      - the phase reverses its trend by more than PARASITIC_PHASE_REVERSAL_DEG, or
#      - CH2 drops below the ADC resolution (per-point V/div).
# 3. The dip (partial cancellation of filter and parasitic path) marks the takeover.
#    The parasitic level is extrapolated from the trend of the curve beyond the dip back
#    towards the passband. Points are valid if they lie at least PARASITIC_MARGIN_DB
#    above it (10 dB -> max. error approx. +-3 dB / +-18 deg).
PASSBAND_CANDIDATE_DB: float = 6.0             # anchor candidates: within 6 dB of the maximum
PARASITIC_STOPBAND_GATE_DB: float = 10.0       # checks only active this far below the passband
PARASITIC_MAG_RISE_DB: float = 3.0
PARASITIC_PHASE_REVERSAL_DEG: float = 10.0
PARASITIC_PHASE_TREND_MIN_DEG: float = 10.0    # minimum phase travel to establish a trend
PARASITIC_MIN_AMPLITUDE_LSB: float = 0.9       # 0.9 LSB = 18 mV at 0.5 V/div
PARASITIC_LEVEL_WINDOW_DEC: float = 0.5        # parasitic trend is fitted over 0.5 decade …
PARASITIC_FIT_START_DEC: float = 0.2           # … starting 0.2 decade beyond the dip
PARASITIC_MAX_SLOPE_DB_DEC: float = 40.0       # clamp for the extrapolated parasitic slope
PARASITIC_MARGIN_DB: float = 10.0              # required distance of valid points above the parasitic level


def find_passband_anchor(log_f: np.ndarray, mag_db: np.ndarray) -> int:
    slope_db_per_dec = np.abs(np.gradient(mag_db, log_f))
    candidates = np.where(mag_db >= np.max(mag_db) - PASSBAND_CANDIDATE_DB)[0]
    return int(candidates[np.argmin(slope_db_per_dec[candidates])])


def estimate_parasitic_level(
    log_f: np.ndarray, mag_db: np.ndarray, dip_idx: int, step: int,
) -> Callable[[int], float]:
    """Fit the parasitic trend beyond the dip and return level(index) in dB.
    Points right after the dip are still lowered by the cancellation, so the fit uses
    the window PARASITIC_FIT_START_DEC … +PARASITIC_LEVEL_WINDOW_DEC beyond the dip,
    where the parasitic path clearly dominates. Falls back to closer points if needed."""
    all_idx = np.arange(len(mag_db))
    outward_dist = (log_f[all_idx] - log_f[dip_idx]) * step
    beyond = all_idx[(outward_dist >= PARASITIC_FIT_START_DEC)
                     & (outward_dist <= PARASITIC_FIT_START_DEC + PARASITIC_LEVEL_WINDOW_DEC)]
    if len(beyond) < 2:
        beyond = all_idx[((all_idx - dip_idx) * step >= 2) & (outward_dist <= PARASITIC_LEVEL_WINDOW_DEC)]
    if len(beyond) == 0:
        beyond = np.array([i for i in (dip_idx + step,) if 0 <= i < len(mag_db)], dtype=int)
    if len(beyond) == 0:
        flat_level = float(mag_db[dip_idx])
        return lambda _i: flat_level

    center_log_f = float(np.mean(log_f[beyond]))
    center_level = float(np.mean(mag_db[beyond]))
    outward_slope = 0.0
    if len(beyond) >= 2:
        outward_slope = float(np.polyfit(log_f[beyond], mag_db[beyond], 1)[0]) * step
        # The parasitic path must not grow towards the passband; limit the extrapolation
        outward_slope = float(np.clip(outward_slope, 0.0, PARASITIC_MAX_SLOPE_DB_DEC))
    return lambda i: center_level + outward_slope * (log_f[i] - center_log_f) * step


def find_reliable_boundary(
    log_f: np.ndarray, mag_db: np.ndarray, phase_deg: np.ndarray, below_floor: np.ndarray,
    anchor_idx: int, step: int, stopband_level_db: float,
) -> int:
    """Walk from anchor_idx in direction step (+1 = higher f, -1 = lower f) and
    return the last index that still belongs to the real filter response."""
    n_points = len(mag_db)
    entry_idx: Optional[int] = None
    min_idx = extreme_idx = anchor_idx
    phase_trend = 0
    trigger_idx: Optional[int] = None
    floor_hit = False

    idx = anchor_idx + step
    while 0 <= idx < n_points:
        if below_floor[idx]:
            trigger_idx, floor_hit = idx, True
            break
        if entry_idx is None:
            if mag_db[idx] < stopband_level_db:
                entry_idx = min_idx = extreme_idx = idx
                travel = phase_deg[idx] - phase_deg[anchor_idx]
                phase_trend = int(np.sign(travel)) if abs(travel) >= PARASITIC_PHASE_TREND_MIN_DEG else 0
        else:
            if mag_db[idx] < mag_db[min_idx]:
                min_idx = idx
            if phase_trend and (phase_deg[idx] - phase_deg[extreme_idx]) * phase_trend > 0:
                extreme_idx = idx
            magnitude_rises = mag_db[idx] > mag_db[min_idx] + PARASITIC_MAG_RISE_DB
            phase_reverses = bool(phase_trend) and (
                (phase_deg[extreme_idx] - phase_deg[idx]) * phase_trend > PARASITIC_PHASE_REVERSAL_DEG
            )
            if magnitude_rises or phase_reverses:
                trigger_idx = idx
                break
        idx += step

    if trigger_idx is None:
        return idx - step                      # reached the end of the sweep undisturbed
    if entry_idx is None:
        return max(trigger_idx - step, anchor_idx) if step > 0 else min(trigger_idx - step, anchor_idx)

    # Dip = lowest point from the stopband entry up to half a decade beyond the trigger
    region = np.arange(len(mag_db))
    region = region[((region - entry_idx) * step >= 0)
                    & ((log_f[region] - log_f[trigger_idx]) * step <= PARASITIC_LEVEL_WINDOW_DEC)]
    dip_idx = int(region[np.argmin(mag_db[region])])
    start_idx = dip_idx - step
    if floor_hit and (trigger_idx - dip_idx) * step <= 0:
        start_idx = trigger_idx - step         # never beyond the ADC-floor point

    parasitic_level = estimate_parasitic_level(log_f, mag_db, dip_idx, step)
    idx = start_idx
    while idx != anchor_idx and mag_db[idx] < parasitic_level(idx) + PARASITIC_MARGIN_DB:
        idx -= step
    return idx


log_freqs_array = np.log10(freqs_array)
noise_floor_amplitude_v = PARASITIC_MIN_AMPLITUDE_LSB * np.array(raw_vdiv_ch2) / ADC_COUNTS_PER_DIV
below_noise_floor = amps2_array < noise_floor_amplitude_v

passband_anchor_idx = find_passband_anchor(log_freqs_array, mag_db_array)
stopband_level_db = float(mag_db_array[passband_anchor_idx]) - PARASITIC_STOPBAND_GATE_DB
left_boundary_idx = find_reliable_boundary(
    log_freqs_array, mag_db_array, unwrapped_phase_array, below_noise_floor,
    passband_anchor_idx, -1, stopband_level_db,
)
right_boundary_idx = find_reliable_boundary(
    log_freqs_array, mag_db_array, unwrapped_phase_array, below_noise_floor,
    passband_anchor_idx, +1, stopband_level_db,
)

valid_slice = slice(left_boundary_idx, right_boundary_idx + 1)
f_valid = freqs_array[valid_slice]
mag_valid = mag_db_array[valid_slice]
phase_valid = unwrapped_phase_array[valid_slice]

# Parasitic segments share their boundary point with the valid curve so the lines connect
parasitic_segments: List[slice] = []
if left_boundary_idx > 0:
    parasitic_segments.append(slice(0, left_boundary_idx + 1))
if right_boundary_idx < len(freqs_array) - 1:
    parasitic_segments.append(slice(right_boundary_idx, len(freqs_array)))

print(
    f"[ANALYSIS] Reliable range: {format_frequency(freqs_array[left_boundary_idx])} – "
    f"{format_frequency(freqs_array[right_boundary_idx])} "
    f"(down to {float(np.min(mag_valid)):.1f} dB)"
)

# --- Filter character & -3 dB cutoff (evaluated on the reliable range only) ---
# Low-pass has high gain at the start and low at the end; high-pass is vice versa.
start_gain_median = float(np.median(mag_valid[:5]))
end_gain_median = float(np.median(mag_valid[-5:]))
is_highpass = end_gain_median > (start_gain_median + 6.0)

passband_reference_db = float(np.percentile(mag_valid, 90))
cutoff_target_db = passband_reference_db - 3.01
detected_cutoff_hz: Optional[float] = None

crossing_indices = np.where(np.diff(mag_valid >= cutoff_target_db))[0]
if len(crossing_indices) > 0:
    idx_cutoff_band = int(crossing_indices[0] if is_highpass else crossing_indices[-1])
    f_lower, f_upper = f_valid[idx_cutoff_band], f_valid[idx_cutoff_band + 1]
    g_lower, g_upper = mag_valid[idx_cutoff_band], mag_valid[idx_cutoff_band + 1]
    if g_upper != g_lower:
        log_fc = np.interp(
            cutoff_target_db,
            [min(g_lower, g_upper), max(g_lower, g_upper)],
            [np.log10(f_lower), np.log10(f_upper)] if g_lower < g_upper else [np.log10(f_upper), np.log10(f_lower)],
        )
        detected_cutoff_hz = 10.0**log_fc


# ==============================================================================
# Enterprise GUI Layout & Graphics Initialization
# ==============================================================================

figure = plt.figure(figsize=(15.0, 8.5))

grid_spec = figure.add_gridspec(
    2, 2,
    left=0.08, right=0.95, top=0.92, bottom=0.22,
    width_ratios=[1.1, 1.0], height_ratios=[1, 1],
    hspace=0.35, wspace=0.25,
)

ax_oscilloscope = figure.add_subplot(grid_spec[:, 0])
ax_magnitude    = figure.add_subplot(grid_spec[0, 1])
ax_phase        = figure.add_subplot(grid_spec[1, 1])

for axis in (ax_oscilloscope, ax_magnitude, ax_phase):
    axis.set_navigate(False)

initial_record = sweep_records[0]

# Voltage axis scales with the largest measured peak (CH1 or CH2) across the sweep
peak_voltage = max(
    max(float(np.max(np.abs(rec["ch1_waveform"]))), float(np.max(np.abs(rec["ch2_waveform"]))))
    for rec in sweep_records
)
osc_y_limit = PLOT_Y_HEADROOM * max(peak_voltage, 1e-3)

# --- 1. Oscilloscope Time-Domain Axis ---
ax_oscilloscope.axhline(0.0, color="black", linestyle="-", lw=1.0, alpha=0.7)

trace_ch1, = ax_oscilloscope.plot(
    initial_record["time_axis_scaled"], initial_record["ch1_waveform"],
    label="CH1 (Input)", color="tab:blue", lw=2.0
)
trace_ch2, = ax_oscilloscope.plot(
    initial_record["time_axis_scaled"], initial_record["ch2_waveform"],
    label="CH2 (Output)", color="tab:orange", lw=2.0
)

line_marker_ch1 = ax_oscilloscope.axvline(0, color="tab:blue", linestyle="--", lw=1.3, alpha=0.85)
line_marker_ch2 = ax_oscilloscope.axvline(0, color="tab:orange", linestyle="--", lw=1.3, alpha=0.85)

DIMENSION_ARROW_Y: float = -0.78 * osc_y_limit
ARROW_LABEL_OFFSET: float = 0.03 * osc_y_limit
arrow_delta_t = ax_oscilloscope.annotate(
    "",
    xy=(0, DIMENSION_ARROW_Y),
    xytext=(0, DIMENSION_ARROW_Y),
    arrowprops=dict(arrowstyle="<->", color="black", lw=1.4, shrinkA=0, shrinkB=0),
)
label_arrow_delta_t = ax_oscilloscope.text(
    0, DIMENSION_ARROW_Y + ARROW_LABEL_OFFSET, r"$\Delta t$",
    color="black", fontsize=10, ha="center", va="bottom", fontweight="bold"
)

text_phase_telemetry = ax_oscilloscope.text(
    0.04, 0.05, "", transform=ax_oscilloscope.transAxes, fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85),
)

ax_oscilloscope.grid(True, linestyle="--", alpha=0.6)
ax_oscilloscope.set_ylabel("Voltage [V]")
ax_oscilloscope.set_ylim(-osc_y_limit, osc_y_limit)
ax_oscilloscope.set_xlim(0, initial_record["total_time_scaled"])
ax_oscilloscope.set_xlabel(f"Time [{initial_record['time_unit']}]")
ax_oscilloscope.legend(loc="upper right")

# --- 2. Bode Magnitude Axis ---
ax_magnitude.semilogx(f_valid, mag_valid, color="tab:blue", lw=2.0, label="Magnitude Response")
for segment_no, segment in enumerate(parasitic_segments):
    ax_magnitude.semilogx(
        freqs_array[segment], mag_db_array[segment], color="crimson", lw=2.0, linestyle="--",
        label="Parasitic / Noise" if segment_no == 0 else None,
    )

cursor_mag_marker, = ax_magnitude.plot([initial_record["freq_hz"]], [mag_db_array[0]], "ro", ms=7)
cursor_mag_vline = ax_magnitude.axvline(initial_record["freq_hz"], color="red", linestyle=":", alpha=0.8)

if detected_cutoff_hz:
    ax_magnitude.axvline(
        detected_cutoff_hz, color="forestgreen", linestyle="--", lw=1.5,
        label=f"$f_c$ = {detected_cutoff_hz/1e3:.2f} kHz (-3 dB)"
    )
    ax_magnitude.axhline(cutoff_target_db, color="forestgreen", linestyle=":", alpha=0.6)

ax_magnitude.set_ylabel("Magnitude [dB]")
ax_magnitude.set_title("Bode Diagram: Magnitude Response", fontsize=11)
ax_magnitude.grid(True, which="both", linestyle="--", alpha=0.5)
ax_magnitude.legend(loc="lower left" if is_highpass else "lower left", fontsize=8)

# --- 3. Bode Phase Axis ---
ax_phase.semilogx(f_valid, phase_valid, color="tab:orange", lw=2.0, label="Phase Response")
for segment in parasitic_segments:
    ax_phase.semilogx(
        freqs_array[segment], unwrapped_phase_array[segment], color="crimson", lw=2.0, linestyle="--"
    )

cursor_phase_marker, = ax_phase.plot([initial_record["freq_hz"]], [unwrapped_phase_array[0]], "ro", ms=7)
cursor_phase_vline = ax_phase.axvline(initial_record["freq_hz"], color="red", linestyle=":", alpha=0.8)

if detected_cutoff_hz:
    phase_at_cutoff = float(np.interp(detected_cutoff_hz, freqs_array, unwrapped_phase_array))
    ax_phase.axvline(
        detected_cutoff_hz, color="forestgreen", linestyle="--", lw=1.5,
        label=rf"$\varphi(f_c)$ = {phase_at_cutoff:.1f}°"
    )
    ax_phase.axhline(phase_at_cutoff, color="forestgreen", linestyle=":", alpha=0.6)

ax_phase.set_ylabel("Phase [°]")
ax_phase.set_xlabel("Frequency [Hz]")
ax_phase.set_title("Bode Diagram: Phase Response", fontsize=11)
ax_phase.grid(True, which="both", linestyle="--", alpha=0.5)
ax_phase.legend(loc="lower left" if not is_highpass else "lower left", fontsize=8)

# Auto-range phase response cleanly for either HP (0..+90°) or LP (0..-90°)
min_p = float(np.min(phase_valid))
max_p = float(np.max(phase_valid))
ax_phase.set_ylim(min(-15.0, min_p - 15.0), max(15.0, max_p + 15.0))


# ==============================================================================
# Helper Formatting & Interaction Handlers
# ==============================================================================

def get_scaled_timebase(record: Dict[str, Any]) -> float:
    return record["timebase_s"] * record["unit_multiplier"]


def build_oscilloscope_title(record: Dict[str, Any], gain_db: float) -> str:
    title = (
        f"Frequency: {format_frequency(record['freq_hz'])} | "
        f"Timebase: {get_scaled_timebase(record):.1f} {record['time_unit']}/Div | "
        f"V/Div: {format_vdiv(record['vdiv_ch1'])} / {format_vdiv(record['vdiv_ch2'])}\n"
        f"CH1: {record['vpp_ch1']:.2f} Vpp | CH2: {record['vpp_ch2']*1000:.1f} mVpp | "
        f"Gain: {gain_db:.1f} dB"
    )
    clipped = [name for name, key in (("CH1", "clipped_ch1"), ("CH2", "clipped_ch2")) if record[key]]
    if clipped:
        title += f" | ⚠ {' + '.join(clipped)} clipped"
    return title


def find_interpolated_zero_crossing(
    t_arr: np.ndarray, y_arr: np.ndarray, target_t: float
) -> Optional[float]:
    """Finds the rising zero-crossing closest to target_t using sub-sample linear interpolation."""
    cross_idx = np.where((y_arr[:-1] <= 0.0) & (y_arr[1:] > 0.0))[0]
    if len(cross_idx) == 0:
        return None

    # Linear interpolation: t_zero = t0 + (0 - y0) * (t1 - t0) / (y1 - y0)
    exact_times = []
    for idx in cross_idx:
        t0, t1 = t_arr[idx], t_arr[idx + 1]
        y0, y1 = y_arr[idx], y_arr[idx + 1]
        if y1 != y0:
            t_interp = t0 - y0 * (t1 - t0) / (y1 - y0)
            exact_times.append(t_interp)

    if not exact_times:
        return None

    exact_times_arr = np.array(exact_times)
    best_idx = np.argmin(np.abs(exact_times_arr - target_t))
    return float(exact_times_arr[best_idx])


title_handle = ax_oscilloscope.set_title(
    build_oscilloscope_title(initial_record, mag_db_array[0]),
    fontsize=10,
)

ax_slider_handle = figure.add_axes([0.15, 0.08, 0.70, 0.04])
slider_control = Slider(
    ax=ax_slider_handle,
    label="Frequency Step ",
    valmin=0,
    valmax=len(sweep_records) - 1,
    valinit=0,
    valstep=1,
    color="tab:blue",
)


def on_slider_update(val: float) -> None:
    selected_idx = int(slider_control.val)
    record = sweep_records[selected_idx]
    current_freq = record["freq_hz"]
    current_phase = unwrapped_phase_array[selected_idx]

    # 1. Update oscilloscope waveform traces
    t_axis = record["time_axis_scaled"]
    trace_ch1.set_data(t_axis, record["ch1_waveform"])
    trace_ch2.set_data(t_axis, record["ch2_waveform"])
    ax_oscilloscope.set_xlim(0, record["total_time_scaled"])
    ax_oscilloscope.set_xlabel(f"Time [{record['time_unit']}]")

    # 2. True zero-crossing measurement between CH1 and CH2
    screen_center_t = record["total_time_scaled"] / 2.0
    t_zero_ch1 = find_interpolated_zero_crossing(t_axis, record["ch1_waveform"], screen_center_t)

    # Calculate physical expected zero crossing of CH2 based on phase
    if t_zero_ch1 is not None:
        period_scaled = (1.0 / current_freq) * record["unit_multiplier"]
        # delta_t is positive if CH1 leads (lagging CH2 / lowpass), negative if CH2 leads (highpass)
        ideal_shift_scaled = -(current_phase / 360.0) * period_scaled
        target_ch2_search = t_zero_ch1 + ideal_shift_scaled

        t_zero_ch2 = find_interpolated_zero_crossing(t_axis, record["ch2_waveform"], target_ch2_search)

        if t_zero_ch2 is not None:
            t_left = min(t_zero_ch1, t_zero_ch2)
            t_right = max(t_zero_ch1, t_zero_ch2)
            measured_delta_t = abs(t_zero_ch2 - t_zero_ch1)

            line_marker_ch1.set_xdata([t_zero_ch1, t_zero_ch1])
            line_marker_ch1.set_visible(True)

            line_marker_ch2.set_xdata([t_zero_ch2, t_zero_ch2])
            line_marker_ch2.set_visible(True)

            arrow_delta_t.xy = (t_right, DIMENSION_ARROW_Y)
            arrow_delta_t.set_position((t_left, DIMENSION_ARROW_Y))
            arrow_delta_t.set_visible(True)

            label_arrow_delta_t.set_position(((t_left + t_right) / 2.0, DIMENSION_ARROW_Y + ARROW_LABEL_OFFSET))
            label_arrow_delta_t.set_visible(True)

            text_phase_telemetry.set_text(
                f"$\\Delta t$ = {measured_delta_t:.2f} {record['time_unit']}  ($\\Delta \\varphi$ = {current_phase:+.1f}°)"
            )
        else:
            line_marker_ch1.set_visible(False)
            line_marker_ch2.set_visible(False)
            arrow_delta_t.set_visible(False)
            label_arrow_delta_t.set_visible(False)
            text_phase_telemetry.set_text(f"$\\Delta \\varphi$ = {current_phase:+.1f}°")
    else:
        line_marker_ch1.set_visible(False)
        line_marker_ch2.set_visible(False)
        arrow_delta_t.set_visible(False)
        label_arrow_delta_t.set_visible(False)
        text_phase_telemetry.set_text(f"$\\Delta \\varphi$ = {current_phase:+.1f}°")

    title_handle.set_text(build_oscilloscope_title(record, mag_db_array[selected_idx]))
    title_handle.set_color("crimson" if (record["clipped_ch1"] or record["clipped_ch2"]) else "black")
    trace_ch1.set_label(f"CH1: {record['vpp_ch1']:.2f} Vpp")
    trace_ch2.set_label(f"CH2: {record['vpp_ch2']*1000:.1f} mVpp")
    ax_oscilloscope.legend(loc="upper right")

    # 3. Synchronize cursor telemetry markers in frequency-domain plots
    cursor_mag_marker.set_data([current_freq], [mag_db_array[selected_idx]])
    cursor_mag_vline.set_xdata([current_freq, current_freq])

    cursor_phase_marker.set_data([current_freq], [current_phase])
    cursor_phase_vline.set_xdata([current_freq, current_freq])

    figure.canvas.draw_idle()


on_slider_update(0)
slider_control.on_changed(on_slider_update)


# ==============================================================================
# Bode Plot Export ("Save Bode Plot" button)
# ==============================================================================

EXPORT_DPI: int = 200
EXPORT_PADDING_IN: float = 0.15
EXPORT_FILETYPES = (
    ("PNG image", "*.png"),
    ("PDF document", "*.pdf"),
    ("SVG vector graphic", "*.svg"),
)
slider_cursor_artists = (cursor_mag_marker, cursor_mag_vline, cursor_phase_marker, cursor_phase_vline)


def export_bode_plot(path: str) -> None:
    """Save only the two Bode diagrams (magnitude + phase, incl. titles, labels and
    legends). The slider cursor (red dot / dotted line) is hidden for the export."""
    previous_visibility = [artist.get_visible() for artist in slider_cursor_artists]
    for artist in slider_cursor_artists:
        artist.set_visible(False)
    try:
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        bode_bbox = Bbox.union([axis.get_tightbbox(renderer) for axis in (ax_magnitude, ax_phase)])
        bode_bbox_inches = bode_bbox.transformed(figure.dpi_scale_trans.inverted()).padded(EXPORT_PADDING_IN)
        figure.savefig(path, dpi=EXPORT_DPI, bbox_inches=bode_bbox_inches, facecolor="white")
    finally:
        for artist, visible in zip(slider_cursor_artists, previous_visibility):
            artist.set_visible(visible)
        figure.canvas.draw_idle()


def ask_export_path() -> Optional[str]:
    """File dialog (Tk) with a timestamped default name; falls back to the working
    directory if no Tk dialog is available (e.g. other Matplotlib backends)."""
    default_name = time.strftime("bode_plot_%Y%m%d_%H%M%S.png")
    try:
        from tkinter import filedialog
        parent = figure.canvas.get_tk_widget() if hasattr(figure.canvas, "get_tk_widget") else None
        chosen = filedialog.asksaveasfilename(
            parent=parent,
            title="Save Bode Plot",
            initialfile=default_name,
            defaultextension=".png",
            filetypes=EXPORT_FILETYPES,
        )
        return chosen or None
    except Exception:
        return default_name


def on_save_clicked(_event: Any) -> None:
    path = ask_export_path()
    if not path:
        return
    try:
        export_bode_plot(path)
        print(f"[EXPORT] Bode plot saved to {path}")
        save_button.label.set_text("Saved ✓")
    except Exception as export_err:
        print(f"[EXPORT] Saving failed: {export_err}")
        save_button.label.set_text("Save failed")
    figure.canvas.draw_idle()


ax_save_button = figure.add_axes([0.80, 0.015, 0.15, 0.045])
save_button = Button(ax_save_button, "Save Bode Plot", color="0.92", hovercolor="0.80")
save_button.on_clicked(on_save_clicked)

plt.show()