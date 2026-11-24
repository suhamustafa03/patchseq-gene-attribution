"""Minimal single-compartment conductance-based model, used ONLY as a
biophysical plausibility screen for candidate genes -- never as
validation of a statistical attribution result.

Requires Brian2 (`pip install brian2`); importing this module without it
installed raises a clear ImportError rather than silently no-op'ing, so a
caller can't mistake a skipped simulation for a null plausibility result.
"""

from __future__ import annotations

from dataclasses import dataclass

try:
    import brian2 as b2
    from brian2 import NeuronGroup, StateMonitor, SpikeMonitor, ms, mV, nS, pF, nA, run, start_scope

    _BRIAN2_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only when brian2 missing
    _BRIAN2_AVAILABLE = False


def _require_brian2() -> None:
    if not _BRIAN2_AVAILABLE:
        raise ImportError(
            "simulation.py requires Brian2 for the biophysical plausibility "
            "screen. Install it with `pip install brian2` and re-run. This "
            "is an explicit requirement, not a silently skipped feature, "
            "because a missing simulation must never be mistaken for a "
            "null plausibility result."
        )


# Illustrative, simplified gene -> ion-channel-conductance mapping.
# THIS IS NOT A VALIDATED MAPPING -- each entry must be literature-checked
# per real candidate gene before any scientific use. It exists only so
# this plausibility screen has a concrete conductance to perturb.
GENE_CHANNEL_MAP: dict[str, str] = {
    "Kcnc1": "Kv3",  # Kv3.1 channel subunit, fast-spiking-associated
    "Scn1a": "Nav",  # Nav1.1 channel subunit
}


@dataclass
class PlausibilityScreenResult:
    gene: str
    channel: str
    conductance_multiplier: float
    baseline_firing_rate_hz: float
    perturbed_firing_rate_hz: float
    baseline_ap_width_ms: float
    perturbed_ap_width_ms: float
    firing_rate_shift_direction: str
    ap_width_shift_direction: str
    note: str = (
        "Biophysical plausibility screen only -- not validation of any "
        "statistical attribution result."
    )


def _run_single_compartment(
    g_na_multiplier: float = 1.0, g_kv3_multiplier: float = 1.0, input_current_nA: float = 1.2, duration_ms: float = 500.0
):
    _require_brian2()
    start_scope()

    # Minimal Hodgkin-Huxley-style single compartment, deliberately
    # simplified -- this is a plausibility screen, not a biophysically
    # exhaustive model of any specific interneuron subclass.
    eqs = """
    dv/dt = (gL*(EL-v) + gNa*m**3*h*(ENa-v) + gKv3*n**4*(EK-v) + I)/C : volt
    dm/dt = alpha_m*(1-m) - beta_m*m : 1
    dh/dt = alpha_h*(1-h) - beta_h*h : 1
    dn/dt = alpha_n*(1-n) - beta_n*n : 1
    alpha_m = 0.1/mV*(v+40*mV)/(1-exp(-(v+40*mV)/(10*mV)))/ms : Hz
    beta_m = 4*exp(-(v+65*mV)/(18*mV))/ms : Hz
    alpha_h = 0.07*exp(-(v+65*mV)/(20*mV))/ms : Hz
    beta_h = 1/(1+exp(-(v+35*mV)/(10*mV)))/ms : Hz
    alpha_n = 0.01/mV*(v+55*mV)/(1-exp(-(v+55*mV)/(10*mV)))/ms : Hz
    beta_n = 0.125*exp(-(v+65*mV)/(80*mV))/ms : Hz
    gNa : siemens
    gKv3 : siemens
    I : amp
    """
    gL = 3 * nS
    EL = -65 * mV
    ENa = 50 * mV
    EK = -90 * mV
    C = 100 * pF

    neuron = NeuronGroup(
        1,
        eqs,
        threshold="v > -10*mV",
        refractory="v > -10*mV",
        method="exponential_euler",
        namespace={"gL": gL, "EL": EL, "ENa": ENa, "EK": EK, "C": C},
    )
    neuron.v = EL
    neuron.m = 0.05
    neuron.h = 0.6
    neuron.n = 0.32
    neuron.gNa = 120 * nS * g_na_multiplier
    neuron.gKv3 = 36 * nS * g_kv3_multiplier
    neuron.I = input_current_nA * nA

    spikes = SpikeMonitor(neuron)
    state = StateMonitor(neuron, "v", record=True)
    run(duration_ms * ms)

    firing_rate_hz = spikes.count[0] / (duration_ms / 1000.0)
    ap_width_ms = _estimate_ap_width_ms(state.v[0], float(b2.defaultclock.dt / ms))
    return firing_rate_hz, ap_width_ms


def _estimate_ap_width_ms(v_trace, dt_ms: float) -> float:
    import numpy as np

    v = v_trace / mV if hasattr(v_trace, "__truediv__") else v_trace
    v = np.asarray(v)
    threshold = -20.0
    above = v > threshold
    if not above.any():
        return 0.0
    # width of the longest contiguous above-threshold run, as a crude AP-width proxy
    runs, current = [], 0
    for flag in above:
        current = current + 1 if flag else 0
        runs.append(current)
    return max(runs) * dt_ms


def run_plausibility_screen(gene: str, conductance_multiplier: float) -> PlausibilityScreenResult:
    """Perturb the conductance mapped to `gene` and report the direction
    of the resulting firing-rate / AP-width shift. Plausibility screen
    only; see GENE_CHANNEL_MAP docstring and PlausibilityScreenResult.note.
    """
    if gene not in GENE_CHANNEL_MAP:
        raise ValueError(
            f"No illustrative channel mapping for gene={gene!r}. "
            f"Known genes: {sorted(GENE_CHANNEL_MAP)}. Add a literature-"
            "checked mapping before using a new gene here."
        )
    channel = GENE_CHANNEL_MAP[gene]

    baseline_rate, baseline_width = _run_single_compartment()
    if channel == "Nav":
        perturbed_rate, perturbed_width = _run_single_compartment(g_na_multiplier=conductance_multiplier)
    elif channel == "Kv3":
        perturbed_rate, perturbed_width = _run_single_compartment(g_kv3_multiplier=conductance_multiplier)
    else:  # pragma: no cover - guarded by GENE_CHANNEL_MAP contents
        raise ValueError(f"Unhandled channel {channel!r}")

    def direction(before: float, after: float) -> str:
        if after > before:
            return "increase"
        if after < before:
            return "decrease"
        return "no_change"

    return PlausibilityScreenResult(
        gene=gene,
        channel=channel,
        conductance_multiplier=conductance_multiplier,
        baseline_firing_rate_hz=float(baseline_rate),
        perturbed_firing_rate_hz=float(perturbed_rate),
        baseline_ap_width_ms=float(baseline_width),
        perturbed_ap_width_ms=float(perturbed_width),
        firing_rate_shift_direction=direction(baseline_rate, perturbed_rate),
        ap_width_shift_direction=direction(baseline_width, perturbed_width),
    )
