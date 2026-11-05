import numpy as np


def wright_fisher(N, p0, generations, rng):
    """
    Simulate neutral drift in one population.

    N           : number of diploid individuals (so 2N gene copies)
    p0          : starting frequency of allele A
    generations : how many generations to run
    rng         : a numpy random generator

    Returns an array of allele frequencies, length generations + 1.
    """
    p = p0
    trajectory = [p]

    for _ in range(generations):
        p = rng.binomial(2 * N, p) / (2 * N)
        trajectory.append(p)

    return np.array(trajectory)


def wright_fisher_selection(N, p0, s, generations, rng):
    """
    Wright-Fisher with selection.

    s : selection coefficient. s > 0 favours allele A.

    Selection shifts the expected frequency each generation; drift then
    samples around that shifted expectation. Once the allele fixes or is
    lost the state is absorbing, so the trajectory is padded and the loop
    exits early.
    """
    p = p0
    trajectory = [p]

    for _ in range(generations):
        p_sel = (p * (1 + s)) / (1 + p * s)
        p = rng.binomial(2 * N, p_sel) / (2 * N)
        trajectory.append(p)

        if p == 0.0 or p == 1.0:
            trajectory.extend([p] * (generations - len(trajectory) + 1))
            break

    return np.array(trajectory)


def run_replicates(N, p0, generations, n_reps, rng):
    """
    Run many independent populations under identical conditions.

    A single trajectory tells you nothing in a stochastic model — every
    result here comes from the distribution across replicates.

    Returns a 2D array, shape (n_reps, generations + 1).
    """
    return np.array([
        wright_fisher(N, p0, generations, rng)
        for _ in range(n_reps)
    ])


if __name__ == "__main__":
    rng = np.random.default_rng(seed=42)

    N = 1000
    s = 0.01
    p0 = 1 / (2 * N)          # a single new mutant copy

    n_reps = 5000
    fixed = 0

    for _ in range(n_reps):
        traj = wright_fisher_selection(N, p0, s, generations=2000, rng=rng)
        if traj[-1] == 1.0:
            fixed += 1

    print(f"Fixed: {fixed} / {n_reps}")
    print(f"Observed fixation probability: {fixed / n_reps:.4f}")
    print(f"Haldane's prediction (2s):     {2 * s:.4f}")