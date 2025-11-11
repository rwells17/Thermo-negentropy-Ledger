# Thermodynamic Negentropy Budget Toolkit

Author: **Ryan Wells**  
Status: **v1.0.0 (working prototype)**  
Language: Python 3

This repo contains:

- A **JSON schema** for describing thermodynamic systems in terms of power flows, entropy production, exergy, and Landauer-limited information rates.
- A **validator / analysis tool** (`tool.py`) that:
  - Verifies internal thermodynamic accounting.
  - Computes **negentropy budgets** (J/K), **Landauer bits**, and **exergy ledgers**.
  - Handles both **steady-state** and **time-series** (diurnal, seasonal, notch) systems.
- Example systems spanning:
  - Stars (Sun, neutron star)
  - Black holes (Hawking-limited, clamp behavior)
  - Nuclear power plant
  - AI / data center loads
  - **Earth’s biosphere NPP** (annual and diurnal models, winter/summer variants, and outage “notch” scenarios)
  - A toy negative-bracket demo for clamp sanity

This is a **carrier-aware entropy/exergy budget**, wired for Landauer reasoning and clean unit consistency.

---

## 1. Conceptual Overview

### 1.1 Core quantity: negentropy flux \( \dot{N} \)

For each carrier **leg** \( i \):

\[
\dot{N}_i = P_i \left( \frac{1}{T_{\text{sink},i}} - \alpha_i \cdot \frac{1}{T_{\text{source},i}} \right)
\]

The **system-level** negentropy flux at a given moment is:

\[
\dot{N} = \left[ \sum_i \dot{N}_i \right]_+
\]

Where \([x]_+ = \max(x, 0)\) is the **clamp**.

- **Steady systems**: the bracket is typically ≥ 0, so the clamp rarely “bites”.
- **Time-series systems**: clamp matters; negative raw brackets are explicitly recorded.

#### Sum-then-clamp

Time integration uses **sum-then-clamp**:

1. For each timestep, compute each leg’s \( \dot{N}_i(t) \).
2. Sum them: \( \dot{N}_{\text{raw}}(t) = \sum_i \dot{N}_i(t) \).
3. Clamp: \( \dot{N}(t) = \max(0, \dot{N}_{\text{raw}}(t)) \).
4. Accumulate:

\[
J_{\text{new}} = \int \dot{N}(t)\, dt
\]

The **Toy_NegBracket_demo** system exists specifically to show this behavior: the raw sum goes negative for a couple of hours, the clamp zeroes it, and the time-series reports those indices.

---

### 1.2 Exergy and carriers

Per leg \( i \), we track an **exergy rate** \( W_{\text{ex},i} \) (in watts):

- **Work / electricity / mechanical / chemical** (ideal bound, \( \alpha = 0\)):
  - \( W_{\text{ex},i} = P_i \).
- **Heat**:
  - \( W_{\text{ex},i} = P_i (1 - T_{\text{sink}}/T_{\text{source}}) \).
- **Blackbody radiation** (`radiation_bb`), Petela:
  - \( r = T_{\text{sink}} / T_{\text{source}} \)
  - \( W_{\text{ex}} = P \left[ 1 - \frac{4}{3}r + \frac{1}{3}r^4 \right] \).
- **Spectral radiation** (`radiation_spectral`):
  - General form: \( W_{\text{ex}} = P - T_{\text{sink}} \dot{S} \).
  - In this JSON, we use `alpha_eff` as a compact proxy with guardrails.

The **exergy ledger** enforces:

\[
D = W_{\text{ex,in}} - W_{\text{ex,out}} - W_{\text{useful}} \ge 0
\]

The validator checks that `D_W` matches `W_ex_in_total_W - W_ex_out_total_W - W_useful_W` and that `D_W ≥ 0`.

---

### 1.3 Landauer, bits, and temperatures

Given a leg with exergy \( W_{\text{ex},i} \) and erasure temperature \( T_{\text{erase},i} \):

\[
\text{bits/s}_{\text{ideal},i} \le \frac{W_{\text{ex},i}}{k_B T_{\text{erase},i} \ln 2}
\]

- The JSON stores **per-leg** Landauer bits as fields like:
  - `bits_per_s_landauer_sink_leg` (using `T_sink_K`)
  - `bits_per_s_landauer_erase_leg` (using `T_erase_K`)
- The validator sums over **`erasure_candidate: true`** legs for totals.

Usable bits can be further reduced by `eta_use`:

\[
\text{bits/s}_{\text{usable},i} = \eta_{\text{use},i} \cdot \frac{W_{\text{ex},i}}{k_B T_{\text{erase},i} \ln 2}
\]

---

## 2. JSON Structure

Top-level keys:

- `meta`: schema, units, constants.
- `real_systems`: named systems (Sun, biosphere, data center, etc.).
- `black_holes`: Hawking-limited examples, demonstrating clamp behavior.

### 2.1 `meta`

Key pieces:

- `schema_version`: e.g. `"1.5-hardened"`.
- `constants`:
  - `k_B_J_per_K`
  - `ln2`
  - `kBln2_J_per_K`
  - `bits_per_J_per_K`
- `time_model`:
  - `steady`: values are constant over `t_s`.
  - `time_series`: integration rules (left Riemann, sum-then-clamp).
- `units`: explicit units for field names.
- `carriers`: documentation-only descriptions for each carrier type.

---

### 2.2 Systems and legs

Each entry in `real_systems` (or `black_holes`) looks roughly like:

```json
"My_System": {
  "P_input_W": 1.0e6,         // optional for conservation checks
  "inputs": { ... },          // optional structured source description
  "carrier_legs": [ ... ],    // required
  "totals": { ... },          // recommended
  "exergy_ledger": { ... },   // recommended when in/out exergy is split
  "validation": { ... },      // bookkeeping hints
  "time_series": { ... },     // when using diurnal / notched profiles
  "time_model": "steady"      // or "time_series"
}

