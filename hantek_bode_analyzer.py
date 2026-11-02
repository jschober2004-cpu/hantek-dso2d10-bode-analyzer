"""Automated Frequency Response Analyzer (Bode Plotter) for Hantek DSO2D10.

========================================================================
Key Enhancements:
1. Universal Filter Topology Detection:
   - Automatically classifies Low-Pass vs. High-Pass vs. Flat passbands.
   - Accurately classifies parasitic noise boundaries on the appropriate side
     (left for high-pass, right for low-pass).
2. Sub-Sample Interpolated Zero-Crossing Dimensioning (Delta-t):
   - Computes real rising zero-crossings for both CH1 and CH2 via linear interpolation.
   - Correctly anchors phase lead (Delta-t < 0 / CH2 ahead) and phase lag (Delta-t > 0).
3. Dynamic Passband Axis Scaling:
   - Symmetrically bounds the phase axis to prevent visual blowout from stopband noise.
"""

from __future__ import annotations

import time
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
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
SWEEP_FREQ_START_HZ: float = 10.0
SWEEP_FREQ_STOP_HZ: float = 25_000_000.0
SWEEP_NUM_POINTS: int = 100
DDS_STIMULUS_VPP: float = 2.0
FIXED_VDIV_VOLTS: float = 0.5


def get_optimal_timebase(freq_hz: float) -> float:
    target_timebase = 3.5 / (14.0 * freq_hz)
    return min(TIMEBASE_STEPS, key=lambda tb: abs(tb - target_timebase))


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

scope.write(f":DDS:TYPE SINE;:DDS:AMP {DDS_STIMULUS_VPP:.1f};:DDS:SWIT ON")
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
# Analytical Post-Processing: Topology & Boundary Classification
# ==============================================================================

freqs_array = np.array(raw_freqs)
mag_db_array = np.array(raw_mag_db)
phase_deg_array = np.array(raw_phase_deg)
amps2_array = np.array(raw_amps_ch2)
unwrapped_phase_array = np.rad2deg(np.unwrap(np.deg2rad(phase_deg_array)))

# Determine filter character (Low-Pass vs. High-Pass)
# Low-pass has high gain at start and low at end; High-pass is vice-versa
start_gain_median = float(np.median(mag_db_array[:5]))
end_gain_median = float(np.median(mag_db_array[-5:]))
is_highpass = end_gain_median > (start_gain_median + 6.0)

passband_reference_db = float(np.percentile(mag_db_array, 90))
cutoff_target_db = passband_reference_db - 3.01
detected_cutoff_hz: Optional[float] = None
idx_cutoff_band: int = 0

crossing_indices = np.where(np.diff(mag_db_array >= cutoff_target_db))[0]
if len(crossing_indices) > 0:
    # Pick the appropriate transition point
    idx_cutoff_band = crossing_indices[0] if is_highpass else crossing_indices[-1]
    f_lower, f_upper = freqs_array[idx_cutoff_band], freqs_array[idx_cutoff_band + 1]
    g_lower, g_upper = mag_db_array[idx_cutoff_band], mag_db_array[idx_cutoff_band + 1]
    if g_upper != g_lower:
        log_fc = np.interp(
            cutoff_target_db,
            [min(g_lower, g_upper), max(g_lower, g_upper)],
            [np.log10(f_lower), np.log10(f_upper)] if g_lower < g_upper else [np.log10(f_upper), np.log10(f_lower)],
        )
        detected_cutoff_hz = 10.0**log_fc

# --- Bidirectional Noise Floor / Parasitic Limit Detection ---
# Identifies points dominated by ADC quantization floor or crosstalk (< -36 dB or < 18 mV)
noise_floor_mask = (mag_db_array <= -36.0) | (amps2_array < 0.018)

if is_highpass:
    # Stopband is at LOW frequencies (left side)
    # Search backwards from cutoff to find where reliable passband/transition ends
    stopband_indices = np.where(noise_floor_mask & (np.arange(len(mag_db_array)) <= idx_cutoff_band))[0]
    left_boundary_idx = int(np.max(stopband_indices)) if len(stopband_indices) > 0 else 0

    f_valid = freqs_array[left_boundary_idx:]
    mag_valid = mag_db_array[left_boundary_idx:]
    phase_valid = unwrapped_phase_array[left_boundary_idx:]

    f_parasitic = freqs_array[:left_boundary_idx+1]
    mag_parasitic = mag_db_array[:left_boundary_idx+1]
    phase_parasitic = unwrapped_phase_array[:left_boundary_idx+1]
else:
    # Low-pass: Stopband is at HIGH frequencies (right side)
    stopband_indices = np.where(noise_floor_mask & (np.arange(len(mag_db_array)) >= idx_cutoff_band))[0]
    right_boundary_idx = int(np.min(stopband_indices)) if len(stopband_indices) > 0 else len(mag_db_array) - 1

    f_valid = freqs_array[:right_boundary_idx+1]
    mag_valid = mag_db_array[:right_boundary_idx+1]
    phase_valid = unwrapped_phase_array[:right_boundary_idx+1]

    f_parasitic = freqs_array[right_boundary_idx:]
    mag_parasitic = mag_db_array[right_boundary_idx:]
    phase_parasitic = unwrapped_phase_array[right_boundary_idx:]


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
        f_parasitic, mag_parasitic, color="crimson", lw=2.0, linestyle="--", label="Parasitic / Noise"
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
ax_phase.legend(loc="lower left" if not is_highpass else "lower left", fontsize=8)

# Auto-range phase response cleanly for either HP (0..+90°) or LP (0..-90°)
min_p = float(np.min(phase_valid))
max_p = float(np.max(phase_valid))
ax_phase.set_ylim(min(-15.0, min_p - 15.0), max(15.0, max_p + 15.0))


# ==============================================================================
# Helper Formatting & Interaction Handlers
# ==============================================================================

def format_frequency(freq_val: float) -> str:
    if freq_val < 1e3:
        return f"{freq_val:.1f} Hz"
    if freq_val < 1e6:
        return f"{freq_val/1e3:.2f} kHz"
    return f"{freq_val/1e6:.2f} MHz"


def get_scaled_timebase(record: Dict[str, Any]) -> float:
    return record["timebase_s"] * record["unit_multiplier"]


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

            label_arrow_delta_t.set_position(((t_left + t_right) / 2.0, DIMENSION_ARROW_Y + 0.05))
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


on_slider_update(0)
slider_control.on_changed(on_slider_update)

plt.show()