#!/usr/bin/env python3
"""
How to use it
-------------
1) Validate all systems
   python3 tool.py data.json validate

   Validate a single system (by name or BH label), strict mode:
   python3 tool.py data.json validate --system Biosphere_NPP_diurnal_eta0p01 --strict

2) Dump a CSV for a time-series system (e.g., the 12h diurnal biosphere)
   python3 tool.py data.json dump-csv --system Biosphere_NPP_diurnal_eta0p01 --out diurnal.csv

3) Build the analytic notch delta summary
   python3 tool.py data.json notch-summary --out notch_delta_summary.csv \
       --Tsink 255 \
       --Pon_baseline_12h 2.6618540345721604e14 \
       --Pon_winter_9h   3.5491387127628805e14 \
       --Pon_summer_15h  2.1294832276577283e14
"""

import json
import math
import argparse
import sys
import csv
from typing import Dict, Any, List, Tuple, Optional

# ===== Constants (override from JSON if present) =====
K_B = 1.380649e-23
LN2 = 0.6931471805599453
KBLN2 = K_B * LN2                  # J/K per bit^-1
BITS_PER_J_PER_K = 1.0 / KBLN2
SECONDS_PER_YEAR = 31557600.0

# ===== Carrier aliases → canonical names =====
CARRIER_ALIASES = {
    "electricity": "work",
    "mechanical": "work",
    "electric": "work",
}
CANONICAL_CARRIERS = {"work", "chemical", "heat", "radiation_bb", "radiation_spectral"}

# ===== Utilities =====
def finite(x) -> bool:
    try:
        return (x is not None) and math.isfinite(float(x))
    except Exception:
        return False

def close(a: float, b: float, rel: float = 1e-3, abs_tol: float = 1e-12) -> bool:
    return abs(a - b) <= max(abs_tol, rel * max(1.0, abs(a), abs(b)))

def carrier_canon(l: Dict[str, Any]) -> str:
    raw = str(l.get("carrier", "")).lower()
    return CARRIER_ALIASES.get(raw, raw)

def get_first_key(d: Dict[str, Any], keys: List[str]) -> Optional[str]:
    for k in keys:
        if k in d:
            return k
    return None

def get_float(d: Dict[str, Any], k: str, default: float = 0.0) -> float:
    try:
        return float(d.get(k, default))
    except Exception:
        return float(default)

# ===== JSON I/O =====
def load_data(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        data = json.load(f)
    # Override constants
    global K_B, LN2, KBLN2, BITS_PER_J_PER_K, SECONDS_PER_YEAR
    const = data.get("meta", {}).get("constants", {})
    if "k_B_J_per_K" in const:   K_B = float(const["k_B_J_per_K"])
    if "ln2" in const:           LN2 = float(const["ln2"])
    if "kBln2_J_per_K" in const: KBLN2 = float(const["kBln2_J_per_K"])
    else:                        KBLN2 = K_B * LN2
    BITS_PER_J_PER_K = 1.0 / KBLN2 if KBLN2 else 1.0442737824274138e23
    SECONDS_PER_YEAR = float(data.get("meta", {}).get("t_s", 31557600))
    return data

# ===== Core math =====
def leg_alpha(l: Dict[str, Any]) -> float:
    # alpha or alpha_eff for spectral
    return float(l.get("alpha", l.get("alpha_eff", 0.0)))

def leg_ndot_from_internals(l: Dict[str, Any]) -> float:
    """
    Compute instantaneous Ndot for a leg using current P_W and temperatures.
    Ndot = P * (1/Tsink - alpha/Tsrc), with alpha=0 handling (work/chem ideal).
    """
    if "P_W" not in l or "T_sink_K" not in l:
        raise ValueError(f'{l.get("leg_id","?")}: missing P_W or T_sink_K')
    P = float(l["P_W"])
    Tsink = float(l["T_sink_K"])
    if Tsink <= 0:
        raise ValueError(f'{l.get("leg_id","?")}: T_sink_K must be > 0')
    a = leg_alpha(l)
    Tsrc = l.get("T_source_K", None)

    if a == 0.0:  # work/electric/chemical ideal bound
        return P / Tsink

    if Tsrc is None:
        raise ValueError(f'{l.get("leg_id","?")}: T_source_K null but alpha!=0')

    Tsrc = float(Tsrc)
    if Tsrc <= 0:
        raise ValueError(f'{l.get("leg_id","?")}: T_source_K must be > 0')

    return P * (1.0 / Tsink - a / Tsrc)

# ===== Validation =====
_SOFT_NOTE_KEYS = {
    "Landauer erase sum differs (filtered by erasure_candidate)",
    "Landauer sink sum differs (informative)",
    "sink_label missing",
    "totals missing",
    "non-canonical carrier",
}

def validate_system(sysobj: Dict[str, Any], rel: float = 1e-3, abs_tol: float = 1e-12, strict: bool = False) -> Tuple[bool, List[str]]:
    ok = True
    notes: List[str] = []
    legs = sysobj.get("carrier_legs", []) or []
    totals = sysobj.get("totals", {}) or {}

    # Finite checks + canonical carriers
    for l in legs:
        legid = l.get("leg_id","?")
        for k in ("P_W", "T_sink_K"):
            if not finite(l.get(k, None)):
                ok = False; notes.append(f'{legid}: non-finite {k}')
        if "alpha" not in l and "alpha_eff" not in l:
            ok = False; notes.append(f'{legid}: alpha missing')

        # Canonicalize carrier names (soft warn if non-canonical)
        cc = carrier_canon(l)
        if cc and cc not in CANONICAL_CARRIERS:
            notes.append(f'{legid}: non-canonical carrier "{l.get("carrier")}" treated as "{cc}"')

    # Sums vs totals (robust, presence-guarded)
    if totals:
        Nd_sum = sum(float(l["Ndot_leg_J_per_K_s"]) for l in legs if "Ndot_leg_J_per_K_s" in l)
        Wx_sum = sum(float(l["W_ex_leg_W"])         for l in legs if "W_ex_leg_W" in l)

        if "Ndot_J_per_K_s" in totals and not close(Nd_sum, float(totals["Ndot_J_per_K_s"]), rel, abs_tol):
            ok = False; notes.append("Ndot sum mismatch")
        if "W_exergy_max_W" in totals and not close(Wx_sum, float(totals["W_exergy_max_W"]), rel, abs_tol):
            ok = False; notes.append("W_ex sum mismatch")

        # Landauer sums: allow key variants, filter by erasure_candidate for erase totals
        erase_leg_key = get_first_key(legs[0] if legs else {}, [
            "bits_per_s_landauer_at_erase_leg", "bits_per_s_landauer_erase_leg"
        ]) or "bits_per_s_landauer_erase_leg"
        sink_leg_key = get_first_key(legs[0] if legs else {}, [
            "bits_per_s_landauer_at_sink_leg", "bits_per_s_landauer_sink_leg"
        ]) or "bits_per_s_landauer_sink_leg"

        erase_total_key = get_first_key(totals, ["bits_per_s_landauer_at_erase", "bits_per_s_landauer_erase"])
        sink_total_key  = get_first_key(totals, ["bits_per_s_landauer_at_sink",  "bits_per_s_landauer_sink"])

        if erase_total_key:
            bitsL_erase_sum = 0.0
            for l in legs:
                if l.get("erasure_candidate", False) and (erase_leg_key in l):
                    bitsL_erase_sum += float(l[erase_leg_key])
            if bitsL_erase_sum and not close(bitsL_erase_sum, float(totals[erase_total_key]), rel, abs_tol):
                msg = "Landauer erase sum differs (filtered by erasure_candidate)"
                notes.append(msg)
                if strict:
                    ok = False

        if sink_total_key:
            bitsL_sink_sum = sum(float(l[sink_leg_key]) for l in legs if sink_leg_key in l)
            if bitsL_sink_sum and not close(bitsL_sink_sum, float(totals[sink_total_key]), rel, abs_tol):
                notes.append("Landauer sink sum differs (informative)")
    else:
        notes.append("totals missing")

    # Conservation bookkeeping
    Pin = sysobj.get("P_input_W", None)
    if sysobj.get("validation", {}).get("power_conservation_scope") == "ok" or Pin is not None:
        if Pin is None:
            ok = False; notes.append("Conservation scope=ok but P_input_W missing")
        else:
            Plegs = sum(float(l["P_W"]) for l in legs if "P_W" in l)
            if not close(Plegs, float(Pin), rel=1e-6, abs_tol=1e-9):
                ok = False; notes.append(f'Power conservation mismatch: legs={Plegs} vs input={Pin}')

    # Per-leg checks
    for l in legs:
        legid = l.get("leg_id","?")
        a = leg_alpha(l)
        Ts = l.get("T_source_K", None)
        if a != 0.0 and Ts is None:
            notes.append(f'{legid}: T_source_K null but alpha!=0')
        if "alpha_eff" in l:
            aeff = float(l["alpha_eff"])
            if not (0.0 - 1e-12 <= aeff <= 4.0/3.0 + 1e-12):
                ok = False; notes.append(f'{legid}: alpha_eff out of [0,4/3]')
        if not l.get("sink_label"):
            notes.append(f'{legid}: sink_label missing')
        # Non-negativity of reported exergy (if present)
        if "W_ex_leg_W" in l and float(l["W_ex_leg_W"]) < -abs_tol:
            ok = False; notes.append(f'{legid}: W_ex_leg_W negative')

        # Identities for alpha=0 and eta_use=1 with ideal bound (work/chemical)
        if a == 0.0 and float(l.get("eta_use", 1.0)) == 1.0 and carrier_canon(l) in ("work","chemical"):
            Tk = float(l["T_sink_K"])
            Nd = l.get("Ndot_leg_J_per_K_s")
            Wx = l.get("W_ex_leg_W")
            if Nd is not None and Wx is not None:
                if not close(Tk * float(Nd), float(Wx), rel, abs_tol):
                    ok = False; notes.append(f'{legid}: T_sink*Ndot != W_ex (alpha=0 ideal)')

        # Landauer erase vs sink ratio check if both available
        sink_key_leg  = get_first_key(l, ["bits_per_s_landauer_at_sink_leg",  "bits_per_s_landauer_sink_leg"])
        erase_key_leg = get_first_key(l, ["bits_per_s_landauer_at_erase_leg", "bits_per_s_landauer_erase_leg"])
        if sink_key_leg and erase_key_leg:
            Tsink  = float(l["T_sink_K"])
            Terase = float(l.get("T_erase_K", Tsink))
            sink_bits  = float(l[sink_key_leg])
            erase_bits = float(l[erase_key_leg])
            if sink_bits > 0 and Terase > 0:
                expect_ratio = Tsink / Terase
                got_ratio = erase_bits / sink_bits
                if not close(got_ratio, expect_ratio, rel=1e-2, abs_tol=1e-3):
                    notes.append(f'{legid}: Landauer erase/sink ratio deviates (got {got_ratio:.3g}, expect ~{expect_ratio:.3g})')

    # Exergy ledger sanity
    ledger = sysobj.get("exergy_ledger", {}) or {}
    if ledger:
        Win  = get_float(ledger, "W_ex_in_total_W")
        Wout = get_float(ledger, "W_ex_out_total_W")
        Wuse = get_float(ledger, "W_useful_W")
        D    = get_float(ledger, "D_W")
        if not close(D, Win - Wout - Wuse, rel=1e-6, abs_tol=1e-6):
            ok = False; notes.append(f'Exergy ledger mismatch: D != Win - Wout - Wuse ({D} vs {Win - Wout - Wuse})')
        if D < -1e-6:
            ok = False; notes.append("Exergy destruction D < 0")

    # Strict mode: treat soft notes as failures too
    if strict:
        for n in notes:
            if n in _SOFT_NOTE_KEYS:
                ok = False
                break

    return ok, notes

def validate_all(data: Dict[str, Any], system_filter: Optional[str] = None, strict: bool = False) -> Dict[str, Any]:
    results: Dict[str, Any] = {}
    schema = data.get("meta", {}).get("schema_version", "unknown")

    def do_validate(name: str, sysobj: Dict[str, Any]):
        ok, notes = validate_system(sysobj, strict=strict)
        results[name] = {"ok": ok, "notes": notes, "schema": schema}

    # Filtered validation
    if system_filter:
        sysobj = data.get("real_systems", {}).get(system_filter)
        if sysobj is not None:
            do_validate(system_filter, sysobj)
            return results
        # Try black holes by label
        for bh in data.get("black_holes", []) or []:
            if bh.get("label") == system_filter:
                do_validate(system_filter, bh)
                return results
        # Not found
        results[system_filter] = {"ok": False, "notes": [f"System '{system_filter}' not found"], "schema": schema}
        return results

    # All systems
    for name, sysobj in (data.get("real_systems", {}) or {}).items():
        do_validate(name, sysobj)
    for bh in (data.get("black_holes", []) or []):
        do_validate(bh.get("label", "BH"), bh)
    return results

# ===== Time-series CSV dumper =====
def _pad_repeat_last(arr: List[float], steps: int) -> List[float]:
    if len(arr) == steps: return arr
    if len(arr) == 0:     return [0.0] * steps
    if len(arr) > steps:  return arr[:steps]
    return arr + [arr[-1]] * (steps - len(arr))

def compute_timeseries_csv(sysobj: Dict[str, Any], outfile: str) -> None:
    """
    Robust series handling:
      - derive steps from longest available series; fallback to declared steps or 24
      - pad shorter series by repeating last value
      - clip/pad P_in_series to steps
      - emit raw_negative column (1 if Ndot_raw < 0)
      - guard leg_ndot exceptions with a single warning per leg
    """
    st = sysobj.get("time_series", {}) or {}
    default_dt = float(st.get("dt_s", 3600.0))

    legs = sysobj.get("carrier_legs", []) or []

    # Collect series and infer steps
    leg_series: List[List[float]] = []
    leg_meta: List[Dict[str, Any]] = []
    inferred_steps = 0
    for l in legs:
        series = l.get("time_series", {}).get("P_W_series", None)
        if isinstance(series, list):
            arr = [float(x) for x in series]
            inferred_steps = max(inferred_steps, len(arr))
        else:
            arr = []  # will be filled with constant later
        leg_series.append(arr)
        leg_meta.append(l)

    steps = inferred_steps if inferred_steps > 0 else int(st.get("steps", 24))
    dt = default_dt  # single dt for the CSV; mixed dts are not supported

    # Expand constants or pad/repeat-last
    for i, (l, s) in enumerate(zip(leg_meta, leg_series)):
        if len(s) == 0:
            leg_series[i] = [get_float(l, "P_W")] * steps
        else:
            leg_series[i] = _pad_repeat_last(s, steps)

    # Optional input series (clip/pad to steps)
    P_in_series = sysobj.get("inputs", {}) \
                        .get("radiation_capture", {}) \
                        .get("time_series", {}) \
                        .get("P_W_series", None)
    if isinstance(P_in_series, list):
        P_in_series = [float(x) for x in P_in_series]
        P_in_series = _pad_repeat_last(P_in_series, steps)
    else:
        P_in_series = None

    # Track which legs had ndot calculation errors so we only warn once
    warned_bad_leg = set()

    # Compute per-step
    rows = []
    J_cum = 0.0
    bits_cum = 0.0

    for i in range(steps):
        P_in = P_in_series[i] if isinstance(P_in_series, list) else ""
        Ndot_raw = 0.0
        P_chem = 0.0
        P_heat = 0.0

        for l, Parr in zip(leg_meta, leg_series):
            P_i = float(Parr[i])
            temp = dict(l)
            temp["P_W"] = P_i
            try:
                nd = leg_ndot_from_internals(temp)
            except Exception as e:
                leg_id = l.get("leg_id", f"leg{len(warned_bad_leg)}")
                if leg_id not in warned_bad_leg:
                    print(f"[warn] {leg_id}: ndot calc failed at step {i} ({e}); ndot=0 for this leg", file=sys.stderr)
                    warned_bad_leg.add(leg_id)
                nd = 0.0
            Ndot_raw += nd

            cc = carrier_canon(l)
            if cc == "chemical": P_chem += P_i
            elif cc == "heat":   P_heat += P_i

        Ndot_clamped = max(0.0, Ndot_raw)
        J_cum += Ndot_clamped * dt
        bits_cum += Ndot_clamped * dt * (1.0 / KBLN2)

        rows.append({
            "hour": i,
            "P_in_W": P_in,
            "P_chem_W": P_chem if P_chem > 0 else "",
            "P_heat_W": P_heat if P_heat > 0 else "",
            "Ndot_raw_J_per_K_s": Ndot_raw,
            "Ndot_clamped_J_per_K_s": Ndot_clamped,
            "J_cum_J_per_K": J_cum,
            "bits_entropy_cum": bits_cum,
            "raw_negative": 1 if Ndot_raw < 0.0 else 0
        })

    # Write CSV
    with open(outfile, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

# ===== Notch delta summary =====
def notch_delta_summary(outfile: str, T_sink: float, P_on_map: Dict[str, float], hours_removed: float = 2.0) -> None:
    seconds_removed = hours_removed * 3600.0
    rows = []
    for scen, Pon in P_on_map.items():
        Pon = float(Pon)
        dJ = seconds_removed * (Pon / float(T_sink))
        dBits = dJ * (1.0 / KBLN2)
        rows.append({
            "scenario": scen,
            "hours_removed": hours_removed,
            "T_sink_K": T_sink,
            "P_chem_on_W": Pon,
            "Delta_J_expected_J_per_K": dJ,
            "Delta_bits_expected": dBits
        })
    with open(outfile, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

# ===== CLI =====
def main():
    ap = argparse.ArgumentParser(description="Validator and tools for thermodynamic negentropy JSON")
    ap.add_argument("json", help="Path to data.json")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_val = sub.add_parser("validate", help="Validate systems")
    p_val.add_argument("--system", help="Validate only this system (key in real_systems) or BH label")
    p_val.add_argument("--strict", action="store_true", help="Fail on soft Landauer aggregation notes, etc.")

    p_dump = sub.add_parser("dump-csv", help="Dump a system's time series to CSV")
    p_dump.add_argument("--system", required=True, help="System name (key in real_systems) or BH label")
    p_dump.add_argument("--out", required=True, help="Output CSV path")

    p_notch = sub.add_parser("notch-summary", help="Write notch delta summary CSV")
    p_notch.add_argument("--out", required=True, help="Output CSV path")
    p_notch.add_argument("--Tsink", type=float, default=255.0, help="Sink temperature K")
    p_notch.add_argument("--Pon_baseline_12h", type=float, default=2.6618540345721604e14)
    p_notch.add_argument("--Pon_winter_9h",   type=float, default=3.5491387127628805e14)
    p_notch.add_argument("--Pon_summer_15h",  type=float, default=2.1294832276577283e14)
    p_notch.add_argument("--hours_removed",   type=float, default=2.0)

    args = ap.parse_args()
    data = load_data(args.json)
    schema = data.get("meta", {}).get("schema_version", "unknown")

    if args.cmd == "validate":
        print(f"[schema] {schema}")
        results = validate_all(data, system_filter=args.system, strict=args.strict)
        ok_all = True
        for name, res in results.items():
            status = "OK" if res["ok"] else "FAIL"
            note_str = ", ".join(res["notes"]) if res["notes"] else "no issues"
            print(f"{status:4s} | {name} | {note_str}")
            ok_all &= res["ok"]
        sys.exit(0 if ok_all else 2)

    elif args.cmd == "dump-csv":
        sysname = args.system
        sysobj = data.get("real_systems", {}).get(sysname)
        if not sysobj:
            # Try BH by label
            for bh in data.get("black_holes", []) or []:
                if bh.get("label") == sysname:
                    sysobj = bh
                    break
        if not sysobj:
            print(f"System '{sysname}' not found (real_systems or black_holes labels)", file=sys.stderr)
            sys.exit(1)
        if "time_series" not in sysobj:
            print(f"[warn] System '{sysname}' has no top-level time_series; using leg series or constants.", file=sys.stderr)
        compute_timeseries_csv(sysobj, args.out)
        print(f"Wrote {args.out}")

    elif args.cmd == "notch-summary":
        Pmap = {
            "baseline_12h": args.Pon_baseline_12h,
            "winter_9h":    args.Pon_winter_9h,
            "summer_15h":   args.Pon_summer_15h
        }
        notch_delta_summary(args.out, args.Tsink, Pmap, args.hours_removed)
        print(
            f"Wrote notch delta CSV {args.out} for {len(Pmap)} scenarios: "
            f"{', '.join(sorted(Pmap.keys()))}"
        )

if __name__ == "__main__":
    main()
