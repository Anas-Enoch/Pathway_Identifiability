import time
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# TODO: replace this with your real alignment call
# -------------------------
def run_alignment_once(N: int, seed: int = 0) -> None:
    """
    Placeholder for JL-FGW alignment runtime.
    Replace the body with your actual call, e.g.:
        jl_fgw_align(G1, G2, ...)
    """
    rng = np.random.default_rng(seed)
    # proxy: O(N^2) matrix ops to mimic transport-plan-like work
    A = rng.normal(size=(N, N))
    _ = A @ A.T  # O(N^3) actually; keep it O(N^2) instead:
    # Better O(N^2) proxy:
    B = rng.normal(size=(N, N))
    _ = A * B  # O(N^2)


def measure_time(N: int, repeats: int = 5) -> float:
    times = []
    for r in range(repeats):
        t0 = time.perf_counter()
        run_alignment_once(N, seed=r)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return float(np.median(times))


def main():
    # pathway sizes (typical curated pathways)
    Ns = [20, 30, 40, 60, 80, 100, 120, 150]
    repeats = 7

    runtimes = [measure_time(N, repeats=repeats) for N in Ns]

    # Plot
    plt.figure(figsize=(6, 4))
    plt.plot(Ns, runtimes, marker="o")
    plt.xlabel("Number of nodes (N)")
    plt.ylabel("Median runtime per alignment (seconds)")
    plt.title("Runtime scaling of JL-FGW alignment")

    plt.tight_layout()
    plt.savefig("FigS1_RuntimeScaling.pdf")
    plt.close()

    print("Saved FigS1_RuntimeScaling.pdf")
    print("N:", Ns)
    print("runtime(s):", runtimes)


if __name__ == "__main__":
    main()
