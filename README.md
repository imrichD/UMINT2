# UMINT2 – Umelá inteligencia 2 (FEI STU, 2026)

Lab assignments from the course *Umelá inteligencia 2* (I-UI2), master's program Robotics and Cybernetics, FEI STU Bratislava.
Everything is written in MATLAB and focuses on evolutionary and swarm optimization algorithms, tested on standard multimodal benchmark functions (Schwefel, Eggholder) in 10 dimensions.

## Contents

| Folder / file | Topic | Test function |
|---|---|---|
| `zadanie1UMINT2.m` | Genetic algorithm – parametrization (selection, crossover, mutation settings), 5 runs per variant | Schwefel (`testfn3`) |
| `zadnie2/` | Parallel GA (island model, star topology, 5 islands) vs. standard GA | Eggholder |
| `zadanie3/` | Parallel GA with migration vs. standard GA, plus CPGA variant | Eggholder |
| `zadanie4/` | Particle Swarm Optimization (PSO) vs. GA, incl. variant with "scout" particles and swarm-movement visualization | Schwefel |
| `zadanie5/` | Differential Evolution (DE/rand/1/bin) – F × CR parameter grid search, comparison with GA | Schwefel |
| `Dokumentacie/` | Lab reports (in Slovak) | |
| `genetic/` | Genetic toolbox by Ivan Sekaj (FEI STU), used by all assignments | |

Each assignment folder contains the script and a `vystupy*/` folder with the generated graphs (`.png`) and result tables (`.csv`, `.mat`).

## Highlights

- **DE (zadanie 5):** low crossover rate (CR ≤ 0.3) is the key on the separable Schwefel function – all such settings reach the global minimum −4189.83. Best setting F = 0.3, CR = 0.1 reaches it in ~180 generations, GA in ~200–320. GA is faster early on, DE fine-tunes faster thanks to its self-scaling difference-vector step.
- **PSO (zadanie 4):** linearly decreasing inertia weight; scout particles help escape local minima of the Schwefel function.
- **Parallel GA (zadania 2–3):** island model with explorer islands and a central exploiting island, compared against a single-population GA of equal size.

## How to run

1. Open MATLAB (tested on R2024b).
2. Run the script of the assignment, e.g. `zadanie5/zadanie5UMINT2.m`. Paths to `genetic/` and `testfn3.m` are added automatically.
3. Outputs are written into the assignment's `vystupy/` folder.

## Author

Imrich Dzvoňár – FEI STU Bratislava
