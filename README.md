# Modular Multiplication for Shor's Algorithm in $2n+3$ qubits

A Qiskit implementation of the modular exponentiation circuit for Shor's algorithm following the construction described in Beauregard (2003). The notebook also includes a full demo of Shor's algorithm via Quantum Phase Estimation (QPE).

## Overview

Shor's algorithm factors an integer $N$ in polynomial time on a quantum computer by reducing the factoring problem to quantum period finding. The most expensive subroutine is **modular exponentiation** — computing $a^x \bmod N$ in superposition — which requires a carefully constructed quantum oracle.

This project implements that oracle using **Beauregard's construction**, which achieves an efficient qubit count by building the circuit from four nested layers:

| Layer | Gate | What it does |
|-------|------|--------------|
| 1 | $\varphi\text{ADD}(a)$ | Adds a classical constant $a$ in the frequency domain. |
| 2 | $\varphi\text{ADD}(a)\text{MOD}(N)$ | Modular addition using a single ancilla to detect/correct overflow. |
| 3 | $\text{CMULT}(a)\text{MOD}(N)$ | Controlled multiplication, mapping $\|x, b\rangle \mapsto \|x, (b + ax) \bmod N\rangle$. |
| 4 | $C\text{-}U_a$ | Full oracle: controlled $\|x\rangle \mapsto \|ax \bmod N\rangle$ via multiply–swap–unmultiply. |

## Background

The implementation follows:

> S. Beauregard, *Circuit for Shor's algorithm using 2n+3 qubits*, Quantum Information & Computation, 3(2), 2003. [[arXiv:quant-ph/0205095]](https://arxiv.org/abs/quant-ph/0205095)

Key design choices from Beauregard's paper:
- **QFT adder**: Addition is performed in the frequency domain, avoiding carry qubits and reducing ancilla overhead.
- **Controlled modular addition**: Uses a reversible overflow trick with a single ancilla qubit (see below).
- **2n+3 qubits total**: Register layout is `c(1) | x(n) | b(n+1) | anc(1)`, where $n = \lceil \log_2 N \rceil$.

## Contents

| File | Description |
|------|-------------|
| `shor_oracle.ipynb` | Oracle construction, step-by-step derivation, correctness tests, circuit visualisations, and QPE demo for $N=15$. |

## Implementation Notes

### Circuit Structure

The oracle is implemented as a `ShorOracle(a, N)` class. Its circuit is assembled from private builder methods that correspond directly to Beauregard's figures:

1. **`_qft()`** — QFT subcircuit (used internally for arithmetic in the frequency domain)
2. **`_adder()`** — $\varphi\text{ADD}(a)$: adds a classical constant via phase rotations; supports 0, 1, or 2 control qubits
3. **`_mod_adder()`** — $\varphi\text{ADD}(a)\text{MOD}(N)$: modular addition using a single ancilla (see below)
4. **`_multiplier()`** — $\text{CMULT}(a)\text{MOD}(N)$: controlled modular multiplication via $n$ calls to `_mod_adder`
5. **`_build()`** — $C\text{-}U_a$: the full oracle, composed as `CMULT` → controlled-SWAPs → `CMULT†`

### Modular Adder: Reversible Overflow Detection

The modular adder is the most subtle piece. It must compute $b \mapsto (b + a) \bmod N$ **reversibly** using a single ancilla qubit — no measurement allowed. The procedure (all operations in the frequency domain unless otherwise noted):

1. **Add $a$**: Doubly-controlled $\varphi\text{ADD}(a)$ maps $b \mapsto b + a$.
2. **Subtract $N$**: Uncontrolled $\varphi\text{ADD}(N)^\dagger$ yields $b + a - N$.
3. **Set ancilla**: Apply $\text{QFT}^\dagger$; if $b + a - N < 0$ the MSB is $|1\rangle$. Copy it into the ancilla via CNOT, then reapply QFT. The ancilla is now $|1\rangle$ iff $b + a < N$.
4. **Conditionally restore**: Ancilla-controlled $\varphi\text{ADD}(N)$ adds $N$ back when needed, yielding $(b + a) \bmod N$.
5. **Uncompute ancilla**: Apply $\varphi\text{ADD}(a)^\dagger$ and $\text{QFT}^\dagger$. By the identity $(a + b) \bmod N \geq a \Leftrightarrow a + b < N$, the MSB is $|0\rangle$ precisely when the ancilla is $|1\rangle$. A NOT–CNOT–NOT sequence resets the ancilla, and a final QFT and $\varphi\text{ADD}(a)$ restore the register.

### Correctness Testing

`test_oracle(a, N)` checks the oracle against classical modular exponentiation using Qiskit's statevector simulator. For each input $|x\rangle$ and both `ctrl=0` and `ctrl=1`, it verifies that the output state is the correct computational basis state with probability 1 and that all scratch registers ($b$, overflow qubit, ancilla) are reset to $|0\rangle$. All values of $a$ coprime to $N = 15$ are tested.

```
a= 2, N=15 | 11 qubits (2n+3=11) | All passed ☑️
a= 4, N=15 | 11 qubits (2n+3=11) | All passed ☑️
a= 7, N=15 | 11 qubits (2n+3=11) | All passed ☑️
...
```

### Example Circuit Summary $(a = 7,\ N = 15)$

| Property | Value |
|----------|-------|
| Qubits | 11 |
| Depth | 403 |
| Total gates | 750 |

| Gate | Count |
|------|-------|
| `cp` | 400 |
| `h` | 154 |
| `mcphase` | 120 |
| `p` | 32 |
| `u2` | 24 |
| `cx` | 16 |
| `cswap` | 4 |

### Full Shor's Algorithm Demo

Beyond the oracle, the notebook implements a complete QPE-based factoring of $N = 15$:

- **Circuit**: $4n + 3$ qubits total — $2n$ counting qubits and the $2n+3$ total qubit oracle registers.
- **Setup**: The $x$ register is initialised to $|1\rangle$; each counting qubit $j$ controls $U_{a^{2^j}}$, implemented by constructing a fresh `ShorOracle` for the pre-computed classical value $a^{2^j} \bmod N$.
- **Measurement**: After inverse QFT on the counting register, outcomes cluster at multiples of $2^{2n}/r$, where $r$ is the multiplicative order of $a$ mod $N$.
- **Classical post-processing**: The continued fraction algorithm (`Fraction(φ).limit_denominator(N)`) recovers $r$ from the measured phase $\phi = c / 2^{2n}$. Factors are then extracted via $\gcd(a^{r/2} \pm 1,\ N)$.

## Requirements

```
qiskit
qiskit-aer
matplotlib
numpy
pandas
jupyter
```

The oracle can be instantiated directly:

```python
oracle = ShorOracle(a=7, N=15)
print(oracle)          # ShorOracle(a=7, N=15; n=4, qubits=11)
oracle.qc.draw("mpl")  # Draw the full C-U_a circuit
```

## References

- Beauregard, S. (2003). Circuit for Shor's algorithm using 2n+3 qubits. *QIC*, 3(2), 175–185.
- [Qiskit Documentation](https://docs.quantum.ibm.com/)

## License

MIT
