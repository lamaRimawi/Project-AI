# Job Shop Scheduling - Genetic Algorithm

Solves Job Shop Scheduling Problem using Genetic Algorithm to minimize makespan (total completion time).

## Quick Start

1. **Install dependencies:**
   ```bash
   pip install pandas numpy matplotlib plotly openpyxl
   ```

2. **Prepare data:** Create `datasheet.xlsx` with two sheets:
   - **"Processing Time"**: Job processing times for each operation
   - **"Machines Sequence"**: Machine order for each job

3. **Run:**
   ```bash
   python AI.py
   ```

## Input Parameters
- Population size (default: 20)
- Crossover rate (default: 0.9) 
- Mutation rate (default: 0.1)
- Mutation selection rate (default: 0.1)
- Number of iterations (default: 1000)

## Output
- Optimal job sequence and makespan value
- Makespan convergence plot
- Interactive Gantt chart (`gantt_chart.html`)

## Data Format Example

**Processing Time sheet:**
```
Job | Op1 | Op2 | Op3
J1  | 10  | 20  | 15
J2  | 25  | 10  | 30
```

**Machines Sequence sheet:**
```
Job | Op1 | Op2 | Op3  
J1  | 1   | 3   | 2
J2  | 2   | 1   | 3
```

## Algorithm
Uses genetic algorithm with 2-point crossover, repair mechanism for infeasible solutions, circular shift mutation, and roulette wheel selection.
