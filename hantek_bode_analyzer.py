"""Automated Frequency Response Analyzer (Bode Plotter) for Hantek DSO2D10.

========================================================================
Architecture Overview:
----------------------
This production-grade measurement tool controls a Hantek DSO2D10 digital
storage oscilloscope and its integrated DDS waveform generator to execute
a fully automated frequency sweep (10 Hz to 25 MHz).

Key Architectural Components:
1. Low-Level Instrument I/O (SCPI / PyVISA):
   - High-throughput raw binary waveform acquisition via proprietary SCPI commands.
   - Non-blocking error recovery and automatic communication timeout tuning.
2. Digital Signal Processing (DSP):
   - Coherent Single-Bin Discrete Fourier Transform (DFT) with Hann windowing
     to eliminate spectral leakage without assumptions about filter type or order.
   - Automated -3 dB cut-off frequency detection via sub-bin log-interpolation.
   - Dynamic noise floor and parasitic coupling boundary classification.
   - Intelligent stopband phase handling: boundary detection triggers where parasitic
     crosstalk begins, and phase axis scaling is strictly bounded to the physical
     filter passband to prevent visual distortion.
3. Interactive Visualization Engine (Matplotlib):
   - Real-time time-domain signal reconstruction synced with Bode magnitude and phase.
   - Center-seeking zero-crossing phase-lag dimensioning (Delta-t).
   - High-performance idle redraw triggers without memory leaks.

Author: Mechatronics Engineering Pipeline
License: MIT License
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import pyvisa
from scipy.signal import savgol_filter

# Global UI configuration: eliminate standard toolbar for an embedded feel
plt.rcParams["toolbar"] = "None"

# Standard 1-2-5 oscilloscope timebase increments (seconds per division)[cite: 5]
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
SWEEP_FREQ_START_HZ: float = 10.0
SWEEP_FREQ_STOP_HZ: float = 25_000_000.0
SWEEP_NUM_POINTS: int = 100
DDS_STIMULUS_VPP: float = 2.0
FIXED_VDIV_VOLTS: float = 0.5


# ==============================================================================
# Signal Processing and Instrument Control Routines
# ==============================================================================

def get_optimal_timebase(freq_hz: float) -> float:
    """Calculates optimal timebase scale to display ~3.5 cycles across 14 divisions.

    Args:
        freq_hz: Excitation frequency in Hertz.

    Returns:
        Closest discrete timebase setting supported by the DSO hardware in s/div.
    """
    target_timebase = 3.5 / (14.0 * freq_hz)
    return min(TIMEBASE_STEPS, key=lambda tb: abs(tb - target_timebase))


def fetch_waveform_raw(scope_inst: pyvisa.resources.MessageBasedResource) -> bytes:
    """Acquires raw 8000-byte acquisition record directly from hardware buffer[cite: 5].

    Args:
        scope_inst: Active PyVISA instrument connection.

    Returns:
        Combined binary waveform payload, or empty bytes on communication failure.
    """
    samples_total = -1
    samples_data: Optional[bytearray] = None
    samples_received = 0

    # Retry loop with micro-sleeps to absorb internal ARM bus contention
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


def compute_spectral_point(
    ch1_voltage: np.ndarray,
    ch2_voltage: np.ndarray,
    sample_interval_s: float,
    frequency_hz: float,
) -> Tuple[float, float, float]:
    """Computes transfer function magnitude and phase using windowed Fourier projection.

    A symmetric Hann window is applied prior to projection to suppress spectral leakage
    and prevent noise contamination from incomplete terminal periods.

    Args:
        ch1_voltage: Voltage array of reference excitation (CH1).
        ch2_voltage: Voltage array of device under test output (CH2).
        sample_interval_s: Delta time between consecutive digitizer samples.
        frequency_hz: Excitation frequency in Hertz.

    Returns:
        A tuple of (gain_db, phase_degrees, output_amplitude_volts).
    """
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

# Configure integrated DDS signal generator[cite: 5]
scope.write(f":DDS:TYPE SINE;:DDS:AMP {DDS_STIMULUS_VPP:.1f};:DDS:SWIT ON")

# Apply deterministic analog input scaling
scope.write(f":CHANnel1:SCALe {FIXED_VDIV_VOLTS:.2f}")
scope.write(f":CHANnel2:SCALe {FIXED_VDIV_VOLTS:.2f}")
time.sleep(0.1)

vdiv_ch1 = float(scope.query(":CHANnel1:SCALe?").strip())
vdiv_ch2 = float(scope.query(":CHANnel2:SCALe?").strip())
scale_factor_ch1 = vdiv_ch1 / 25.0
scale_factor_ch2 = vdiv_ch2 / 25.0

sweep_frequencies = np.geomspace(SWEEP_FREQ_START_HZ, SWEEP_FREQ_STOP_HZ, SWEEP_NUM_POINTS)
sweep_records: List[Dict[str, Any]] = []

raw_freqs: List[float] = []
raw_mag_db: List[float] = []
raw_phase_deg: List[float] = []
raw_amps_ch2: List[float] = []


# ==============================================================================
# Frequency Sweep Acquisition Loop
# ==============================================================================

print(f"\n[SWEEP] Commencing automated sweep ({SWEEP_NUM_POINTS} points)...")
sweep_start_timestamp = time.time()

try:
    for point_idx, current_freq in enumerate(sweep_frequencies):
        timebase_val = get_optimal_timebase(current_freq)
        scope.write(f":TIMebase:SCALe {timebase_val:.6e};:DDS:FREQ {current_freq:.2f}")

        settle_time_s = max(0.015, min(0.12, 2.5 / current_freq))
        time.sleep(settle_time_s)

        raw_bytes = fetch_waveform_raw(scope)
        if len(raw_bytes) < 8000:
            time.sleep(0.02)
            raw_bytes = fetch_waveform_raw(scope)
            if len(raw_bytes) < 8000:
                continue

        raw_buffer = np.frombuffer(raw_bytes, dtype=np.uint8)
        block_len = 2000

        ch1_interleaved = np.concatenate([raw_buffer[0:block_len], raw_buffer[2*block_len:3*block_len]])
        ch2_interleaved = np.concatenate([raw_buffer[block_len:2*block_len], raw_buffer[3*block_len:4*block_len]])

        ch1_centered = (ch1_interleaved ^ 0x80).astype(float)
        ch2_centered = (ch2_interleaved ^ 0x80).astype(float)

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
        })

        freq_str = (
            f"{current_freq:6.1f} Hz" if current_freq < 1e3 else (
                f"{current_freq/1e3:6.2f} kHz" if current_freq < 1e6 else f"{current_freq/1e6:6.2f} MHz"
            )
        )
        print(f"[{point_idx+1:3d}/{len(sweep_frequencies)}] f = {freq_str:>10s} | Gain = {gain:6.1f} dB", end="\r")

except Exception as err:
    print(f"\n[ERROR] Acquisition fault: {err}")
finally:
    try:
        scope.write(":KEY:LOCK DISABLE")
    except Exception:
        pass
    scope.close()

elapsed_duration = time.time() - sweep_start_timestamp
print(f"\n[SWEEP] Finished. {len(sweep_records)} points acquired in {elapsed_duration:.1f} s.")

if not sweep_records:
    print("[FATAL] No measurement records captured. Halting.")
    exit(1)


# ==============================================================================
# Analytical Post-Processing & Limit Classification
# ==============================================================================

freqs_array = np.array(raw_freqs)
mag_db_array = np.array(raw_mag_db)
phase_deg_array = np.array(raw_phase_deg)
amps2_array = np.array(raw_amps_ch2)

# --- 3.1 Cut-Off Frequency Detection (-3 dB relative to passband level) ---
passband_reference_db = float(np.percentile(mag_db_array, 90))
cutoff_target_db = passband_reference_db - 3.01
detected_cutoff_hz: Optional[float] = None
idx_cutoff_band: int = 0

crossing_indices = np.where(np.diff(mag_db_array < cutoff_target_db))[0]
if len(crossing_indices) > 0:
    idx_cutoff_band = crossing_indices[0]
    f_lower, f_upper = freqs_array[idx_cutoff_band], freqs_array[idx_cutoff_band + 1]
    g_lower, g_upper = mag_db_array[idx_cutoff_band], mag_db_array[idx_cutoff_band + 1]
    if g_upper != g_lower:
        log_fc = np.interp(
            cutoff_target_db,
            [min(g_lower, g_upper), max(g_lower, g_upper)],
            [np.log10(f_upper), np.log10(f_lower)] if g_lower > g_upper else [np.log10(f_lower), np.log10(f_upper)],
        )
        detected_cutoff_hz = 10.0**log_fc

# --- 3.2 Dynamic Noise Floor & Parasitic Limit Detection ---
unwrapped_phase_array = np.rad2deg(np.unwrap(np.deg2rad(phase_deg_array)))

# Parasitic limit criteria:
# 1. Magnitude drops below -38 dB or output voltage drops below 18 mV
past_cutoff_selector = np.arange(len(mag_db_array)) >= idx_cutoff_band
hit_noise_indices = np.where(past_cutoff_selector & ((mag_db_array <= -38.0) | (amps2_array < 0.018)))[0]
mag_limit_idx = hit_noise_indices[0] if len(hit_noise_indices) > 0 else len(mag_db_array) - 1

# 2. Phase reaches physical minimum and turns back up (+4° reversal from minimum)
phase_limit_idx = len(unwrapped_phase_array) - 1
if np.any(past_cutoff_selector):
    past_indices = np.where(past_cutoff_selector)[0]
    min_phase_idx_local = int(np.argmin(unwrapped_phase_array[past_indices]))
    global_min_phase_idx = past_indices[min_phase_idx_local]
    phase_reversal = np.where(
        (np.arange(len(unwrapped_phase_array)) > global_min_phase_idx) &
        (unwrapped_phase_array > unwrapped_phase_array[global_min_phase_idx] + 4.0)
    )[0]
    if len(phase_reversal) > 0:
        phase_limit_idx = global_min_phase_idx

# Boundary transitions exactly where crosstalk begins to dominate
boundary_idx = min(mag_limit_idx, phase_limit_idx)
parasitic_limit_hz = freqs_array[boundary_idx] if boundary_idx < len(freqs_array) - 2 else None

# Partition datasets strictly at boundary_idx
if parasitic_limit_hz:
    f_valid = freqs_array[:boundary_idx+1]
    mag_valid = mag_db_array[:boundary_idx+1]
    phase_valid = unwrapped_phase_array[:boundary_idx+1]

    f_parasitic = freqs_array[boundary_idx:]
    mag_parasitic = mag_db_array[boundary_idx:]
    phase_parasitic = unwrapped_phase_array[boundary_idx:]
else:
    f_valid, mag_valid, phase_valid = freqs_array, mag_db_array, unwrapped_phase_array
    f_parasitic, mag_parasitic, phase_parasitic = [], [], []


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

DIMENSION_ARROW_Y: float = -1.25
arrow_delta_t = ax_oscilloscope.annotate(
    "",
    xy=(0, DIMENSION_ARROW_Y),
    xytext=(0, DIMENSION_ARROW_Y),
    arrowprops=dict(arrowstyle="<->", color="black", lw=1.4, shrinkA=0, shrinkB=0),
)
label_arrow_delta_t = ax_oscilloscope.text(
    0, DIMENSION_ARROW_Y + 0.05, r"$\Delta t$",
    color="black", fontsize=10, ha="center", va="bottom", fontweight="bold"
)

text_phase_telemetry = ax_oscilloscope.text(
    0.04, 0.05, "", transform=ax_oscilloscope.transAxes, fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.85),
)

ax_oscilloscope.grid(True, linestyle="--", alpha=0.6)
ax_oscilloscope.set_ylabel("Voltage [V]")
ax_oscilloscope.set_ylim(-1.6, 1.6)
ax_oscilloscope.set_xlim(0, initial_record["total_time_scaled"])
ax_oscilloscope.set_xlabel(f"Time [{initial_record['time_unit']}]")
ax_oscilloscope.legend(loc="upper right")

# --- 2. Bode Magnitude Axis ---
ax_magnitude.semilogx(f_valid, mag_valid, color="tab:blue", lw=2.0, label="Magnitude Response")
if len(f_parasitic) > 0:
    ax_magnitude.semilogx(
        f_parasitic, mag_parasitic, color="crimson", lw=2.0, linestyle="--", label="Parasitic Coupling"
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
ax_magnitude.legend(loc="lower left", fontsize=8)

# --- 3. Bode Phase Axis ---
ax_phase.semilogx(f_valid, phase_valid, color="tab:orange", lw=2.0, label="Phase Response")
if len(f_parasitic) > 0:
    ax_phase.semilogx(f_parasitic, phase_parasitic, color="crimson", lw=2.0, linestyle="--")

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
ax_phase.legend(loc="lower left", fontsize=8)

# Bound phase scale strictly to valid physical filter range to prevent noise blowout
min_p = float(np.min(phase_valid))
max_p = float(np.max(phase_valid))
y_phase_lower = min(-95.0, min_p - 15.0)
y_phase_upper = max(10.0, max_p + 10.0)
ax_phase.set_ylim(y_phase_lower, y_phase_upper)


# ==============================================================================
# Helper Formatting & Interaction Handlers
# ==============================================================================

def format_frequency(freq_val: float) -> str:
    """Formats frequency with adaptive SI unit prefixes."""
    if freq_val < 1e3:
        return f"{freq_val:.1f} Hz"
    if freq_val < 1e6:
        return f"{freq_val/1e3:.2f} kHz"
    return f"{freq_val/1e6:.2f} MHz"


def get_scaled_timebase(record: Dict[str, Any]) -> float:
    """Scales numeric timebase according to the active unit multiplier."""
    return record["timebase_s"] * record["unit_multiplier"]


title_handle = ax_oscilloscope.set_title(
    f"Frequency: {format_frequency(initial_record['freq_hz'])} | "
    f"Timebase: {get_scaled_timebase(initial_record):.1f} {initial_record['time_unit']}/Div\n"
    f"CH1: {initial_record['vpp_ch1']:.2f} Vpp | CH2: {initial_record['vpp_ch2']*1000:.1f} mVpp | "
    f"Gain: {mag_db_array[0]:.1f} dB",
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
    """Synchronizes UI panes when inspecting individual sweep frequency points."""
    selected_idx = int(slider_control.val)
    record = sweep_records[selected_idx]
    current_freq = record["freq_hz"]
    current_phase = unwrapped_phase_array[selected_idx]

    # 1. Update oscilloscope waveform plots
    trace_ch1.set_data(record["time_axis_scaled"], record["ch1_waveform"])
    trace_ch2.set_data(record["time_axis_scaled"], record["ch2_waveform"])
    ax_oscilloscope.set_xlim(0, record["total_time_scaled"])
    ax_oscilloscope.set_xlabel(f"Time [{record['time_unit']}]")

    # 2. Precise phase-lag dimensioning (Anchored near screen center)
    delta_t_seconds = abs(current_phase) / (360.0 * current_freq)
    delta_t_scaled = delta_t_seconds * record["unit_multiplier"]

    ch1_samples = record["ch1_waveform"]
    zero_crossings = np.where((ch1_samples[:-1] <= 0) & (ch1_samples[1:] > 0))[0]

    if len(zero_crossings) > 0:
        screen_center_t = record["total_time_scaled"] / 2.0
        crossing_times = record["time_axis_scaled"][zero_crossings]
        best_crossing_idx = zero_crossings[np.argmin(np.abs(crossing_times - screen_center_t))]

        t_start = record["time_axis_scaled"][best_crossing_idx]
        t_target = t_start + delta_t_scaled

        line_marker_ch1.set_xdata([t_start, t_start])
        line_marker_ch1.set_visible(True)

        line_marker_ch2.set_xdata([t_target, t_target])
        line_marker_ch2.set_visible(True)

        arrow_delta_t.xy = (t_target, DIMENSION_ARROW_Y)
        arrow_delta_t.set_position((t_start, DIMENSION_ARROW_Y))
        arrow_delta_t.set_visible(True)

        label_arrow_delta_t.set_position(((t_start + t_target) / 2.0, DIMENSION_ARROW_Y + 0.05))
        label_arrow_delta_t.set_visible(True)

        text_phase_telemetry.set_text(
            f"$\\Delta t$ = {delta_t_scaled:.2f} {record['time_unit']}  ($\\Delta \\varphi$ = {current_phase:.1f}°)"
        )
    else:
        line_marker_ch1.set_visible(False)
        line_marker_ch2.set_visible(False)
        arrow_delta_t.set_visible(False)
        label_arrow_delta_t.set_visible(False)
        text_phase_telemetry.set_text(f"$\\Delta \\varphi$ = {current_phase:.1f}°")

    timebase_display = get_scaled_timebase(record)
    title_handle.set_text(
        f"Frequency: {format_frequency(current_freq)} | Timebase: {timebase_display:.1f} {record['time_unit']}/Div\n"
        f"CH1: {record['vpp_ch1']:.2f} Vpp | CH2: {record['vpp_ch2']*1000:.1f} mVpp | Gain: {mag_db_array[selected_idx]:.1f} dB"
    )
    trace_ch1.set_label(f"CH1: {record['vpp_ch1']:.2f} Vpp")
    trace_ch2.set_label(f"CH2: {record['vpp_ch2']*1000:.1f} mVpp")
    ax_oscilloscope.legend(loc="upper right")

    # 3. Synchronize cursor telemetry markers in frequency-domain plots
    cursor_mag_marker.set_data([current_freq], [mag_db_array[selected_idx]])
    cursor_mag_vline.set_xdata([current_freq, current_freq])

    cursor_phase_marker.set_data([current_freq], [current_phase])
    cursor_phase_vline.set_xdata([current_freq, current_freq])

    figure.canvas.draw_idle()


# Initialize state for index 0
on_slider_update(0)

slider_control.on_changed(on_slider_update)

plt.show()
