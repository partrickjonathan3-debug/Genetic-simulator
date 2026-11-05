import numpy as np
import matplotlib.pyplot as plt

from drift import run_replicates, wright_fisher_selection


def plot_trajectories(N, p0, generations, n_reps, rng,
                      outfile="figures/trajectories.png"):
    """Plot a sample of neutral drift trajectories."""
    trajs = run_replicates(N, p0, generations, n_reps, rng)

    fig, ax = plt.subplots(figsize=(8, 5))

    for traj in trajs:
        ax.plot(traj, alpha=0.5, linewidth=1)

    ax.axhline(p0, color="black", linestyle="--", alpha=0.6,
               label=f"Starting frequency ({p0})")

    ax.set_xlabel("Generation")
    ax.set_ylabel("Allele frequency")
    ax.set_ylim(0, 1)
    ax.set_title(f"Genetic drift (N = {N}, {n_reps} populations)")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_population_sizes(p0, generations, n_reps, sizes, rng,
                          outfile="figures/population_sizes.png"):
    """
    Compare drift in populations of different sizes.

    sharey=True locks all panels to the same axis — without it matplotlib
    autoscales each one and the comparison becomes meaningless.
    """
    fig, axes = plt.subplots(1, len(sizes), figsize=(4 * len(sizes), 4),
                             sharey=True)

    for ax, N in zip(axes, sizes):
        trajs = run_replicates(N, p0, generations, n_reps, rng)

        for traj in trajs:
            ax.plot(traj, alpha=0.5, linewidth=1)

        ax.axhline(p0, color="black", linestyle="--", alpha=0.6)
        ax.set_title(f"N = {N}")
        ax.set_xlabel("Generation")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Allele frequency")
    fig.suptitle("Drift is stronger in small populations")

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_selection_vs_drift(N, p0, generations, n_reps, s_values, rng,
                            outfile="figures/selection.png"):
    """
    Compare trajectories under different selection coefficients.

    Population size and starting frequency are held constant — only the
    strength of selection changes between panels.
    """
    fig, axes = plt.subplots(1, len(s_values),
                             figsize=(4 * len(s_values), 4), sharey=True)

    for ax, s in zip(axes, s_values):
        for _ in range(n_reps):
            traj = wright_fisher_selection(N, p0, s, generations, rng)
            ax.plot(traj, alpha=0.4, linewidth=1)

        ax.axhline(p0, color="black", linestyle="--", alpha=0.6)
        ax.set_title(f"s = {s}")
        ax.set_xlabel("Generation")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)

    axes[0].set_ylabel("Allele frequency")
    fig.suptitle(f"Selection against drift (N = {N}, starting at {p0})")

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved {outfile}")


if __name__ == "__main__":
    rng = np.random.default_rng(seed=1)

    plot_trajectories(N=50, p0=0.5, generations=200, n_reps=30, rng=rng)

    plot_population_sizes(p0=0.5, generations=100, n_reps=30,
                          sizes=[10, 100, 1000], rng=rng)

    plot_selection_vs_drift(N=100, p0=0.1, generations=200, n_reps=40,
                            s_values=[0.0, 0.02, 0.1], rng=rng)