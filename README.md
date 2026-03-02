# Modular Multiplication for Shor's Algorithm in $2n+3$ Qubits

A Qiskit implementation of the modular exponentiation circuit for Shor's algorithm following the method described in Beauregard (2003). The notebook also includes a short factorization demo utilizing this circuit.

## Overview

Shor's algorithm factors an integer $N$ in polynomial time on a quantum computer by reducing the factoring problem to quantum period finding. The most expensive subroutine is **modular exponentiation** — computing $a^x \bmod N$ in superposition — which requires a carefully constructed quantum oracle.

This project implements that oracle using Beauregard's construction, which achieves an efficient qubit count by building the circuit from four nested layers:

| Layer | Gate | What it does |
|-------|------|--------------|
| 1 | $\varphi\text{ADD}(a)$ | Adds a classical constant $a$ in the frequency domain. |
| 2 | $\varphi\text{ADD}(a)\text{MOD}(N)$ | Modular addition using a single ancilla to detect/correct overflow. |
| 3 | $\text{CMULT}(a)\text{MOD}(N)$ | Controlled multiplication, mapping $\|x, b\rangle \mapsto \|x, (b + ax) \bmod N\rangle$. |
| 4 | $C\text{-}U_a$ | Full oracle: controlled $\|x\rangle \mapsto \|ax \bmod N\rangle$ via multiply–swap–unmultiply. |

Register layout: `c(1) | x(n) | b(n+1) | anc(1)` $\Rightarrow 2n+3$ qubits, where $n = \lceil \log_2 N \rceil$.

## Background

> Beauregard, S. "Circuit for Shor's algorithm using 2n+3 qubits", *Quantum Information & Computation*, 3(2), 175-185, 2003. [[arXiv:quant-ph/0205095]](https://arxiv.org/abs/quant-ph/0205095)

## Implementation

The oracle is implemented as a `ShorOracle(a, N)` class with builder methods corresponding to figures in Beauregard's paper:

1. **`_qft()`** — QFT subcircuit via Qiskit's built-in QFT (used for arithmetic in the frequency domain).
2. **`_adder()`** — $\varphi\text{ADD}(a)$: phase-rotation adder with 0, 1, or 2 control qubits.
3. **`_mod_adder()`** — $\varphi\text{ADD}(a)\text{MOD}(N)$: modular addition with a single ancilla.
4. **`_multiplier()`** — $\text{CMULT}(a)\text{MOD}(N)$: controlled modular multiplication via $n$ calls to `_mod_adder`.
5. **`_build()`** — $C\text{-}U_a$: the full oracle, composed as `CMULT` → controlled-SWAPs → `CMULT†`

### Correctness Testing

`test_oracle(a, N)` checks the oracle against classical modular exponentiation using Qiskit's statevector simulator. For each input $|x\rangle$ and both `ctrl=0` and `ctrl=1`, it verifies the output (with high probability) and that all scratch registers ($b$, overflow qubit, ancilla) are reset to $|0\rangle$. All values of $a$ coprime to $N = 15$ are tested.

```
a= 2, N=15 | 11 qubits (2n+3=11) | All passed ☑️
a= 4, N=15 | 11 qubits (2n+3=11) | All passed ☑️
a= 7, N=15 | 11 qubits (2n+3=11) | All passed ☑️
...
```

### Circuit Summary $(a = 7,\ N = 15)$

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
oracle.qc.draw("mpl")  # Draw the base C-U_a circuit. Deeper levels can be drawn as shown in the notebook
```

## References
- Beauregard, S. "Circuit for Shor's algorithm using 2n+3 qubits", *Quantum Information & Computation*, 3(2), 175-185, 2003. [[arXiv:quant-ph/0205095]](https://arxiv.org/abs/quant-ph/0205095)
- [Qiskit Documentation](https://docs.quantum.ibm.com/)

## License

MIT
