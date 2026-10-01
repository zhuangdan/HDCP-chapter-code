"""Reproduce the illustrative Gaussian experiment in the associated book chapter.

These are simple, known-covariance Max, Quadratic and Bonferroni-combined
CUSUM benchmarks, NOT reimplementations of DMS, INSPECT, or other named methods.

Run: python reproduce_illustrations.py
Existing completed cells are reused. Use a clean results directory to change settings.
Dependencies: numpy, pandas, matplotlib. All paths are relative to this script.
"""
from __future__ import annotations

from pathlib import Path
import json
import platform
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FIG = ROOT / 'figures'
RES = ROOT / 'results'
N, P = 200, 200
TRIM, ALPHA = 20, 0.05
B_NULL, R_NULL, R_ALT = 4999, 2000, 1000
BASE_SEED = 20260923
SPARSITIES = (1, 5, 10, 20, 50, 100, 200)
DELTAS = (0.6, 1.0, 1.4)
K = np.arange(TRIM, N - TRIM + 1)
METHODS = ('Max', 'Quadratic', 'Adaptive-Bonferroni')


def seed_rng(*parts: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([BASE_SEED, *parts]))


def contrasts(x: np.ndarray) -> np.ndarray:
    """Rows are replicate x time x coordinate; CUSUM is after minus before."""
    if x.ndim != 3 or x.shape[1:] != (N, P):
        raise ValueError(f'Expected shape (replicate, {N}, {P}).')
    z = np.cumsum(x, axis=1)
    left = z[:, K - 1, :]
    total = z[:, -1:, :]
    kk = K[None, :, None]
    return (kk * total - N * left) / np.sqrt(N * kk * (N - kk))


def scores(c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    max_process = np.abs(c).max(axis=2)
    quad_process = np.square(c).sum(axis=2)
    process = np.stack([max_process, quad_process], axis=2)
    loc_idx = process.argmax(axis=1)
    values = process.max(axis=1)
    locations = K[loc_idx]
    return values, locations


def score_samples(rng: np.random.Generator, count: int,
                  tau: int = 100, sparsity: int = 0,
                  delta: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    if not 1 <= tau < N:
        raise ValueError('tau must lie between 1 and N-1.')
    if not 0 <= sparsity <= P or delta < 0:
        raise ValueError('Invalid signal specification.')
    vals, locs = [], []
    for start in range(0, count, 50):
        size = min(50, count - start)
        x = rng.standard_normal((size, N, P))
        if sparsity:
            x[:, tau:, :sparsity] += delta / np.sqrt(sparsity)
        v, l = scores(contrasts(x))
        vals.append(v)
        locs.append(l)
    return np.concatenate(vals), np.concatenate(locs)


def evaluate(values: np.ndarray, locations: np.ndarray,
             null_sorted: np.ndarray, tau: int | None) -> list[dict]:
    p = np.column_stack([
        (1 + len(null_sorted) - np.searchsorted(null_sorted[:, j], values[:, j],
                                               side='left')) / (len(null_sorted)+1)
        for j in range(2)
    ])
    p_all = np.column_stack([p, np.minimum(1, 2*p.min(axis=1))])
    chosen = p.argmin(axis=1)
    adaptive_loc = locations[np.arange(len(locations)), chosen]
    loc_all = np.column_stack([locations, adaptive_loc])
    result = []
    for j, method in enumerate(METHODS):
        reject = p_all[:, j] <= ALPHA
        rate = float(reject.mean())
        row = {'method': method, 'repetitions': len(values), 'rejections': int(reject.sum()),
               'rejection_rate': rate,
               'mcse_conditional_on_reference': float(np.sqrt(rate*(1-rate)/len(values)))}
        if tau is not None:
            err = np.abs(loc_all[:, j]-tau)
            row.update({'unconditional_mae': float(err.mean()),
                        'conditional_mae': float(err[reject].mean()) if reject.any() else np.nan,
                        'joint_detect_within_5': float((reject & (err <= 5)).mean())})
        result.append(row)
    return result


def savefig(fig, basename: str) -> None:
    fig.tight_layout()
    fig.savefig(FIG / f'{basename}.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'{basename}.png', dpi=180, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    FIG.mkdir(exist_ok=True); RES.mkdir(exist_ok=True)
    settings={'n':N,'p':P,'trim':TRIM,'alpha':ALPHA,'null_reference':B_NULL,
              'null_repetitions':R_NULL,'alternative_repetitions':R_ALT,
              'seed':BASE_SEED,'sparsities':list(SPARSITIES),'deltas':list(DELTAS),
              'implementation_version':1}
    config=RES/'experiment_config.json'
    if config.exists() and json.loads(config.read_text()) != settings:
        raise RuntimeError('Settings differ from saved results. Use a clean results directory before rerunning.')
    config.write_text(json.dumps(settings,indent=2),encoding='utf-8')
    start = time.perf_counter()
    if (RES/'null_reference.npz').exists() and (RES/'null_size.csv').exists():
        ref=np.load(RES/'null_reference.npz')['sorted_statistics']
        null_df=pd.read_csv(RES/'null_size.csv')
    else:
        print('Simulating independent null reference...', flush=True)
        ref, _ = score_samples(seed_rng(0), B_NULL)
        ref = np.sort(ref, axis=0)
        np.savez_compressed(RES/'null_reference.npz', sorted_statistics=ref)
        v, l = score_samples(seed_rng(1), R_NULL)
        null_df = pd.DataFrame(evaluate(v, l, ref, tau=None))
        null_df.to_csv(RES/'null_size.csv', index=False)
    print(null_df.to_string(index=False), flush=True)
    power_path=RES/'power_and_localization.csv'
    rows=pd.read_csv(power_path).to_dict('records') if power_path.exists() else []
    for di, delta in enumerate(DELTAS):
        for si, s in enumerate(SPARSITIES):
            if any(row['Delta']==delta and row['s']==s for row in rows):
                continue
            v, l = score_samples(seed_rng(2, di, si), R_ALT, tau=100, sparsity=s, delta=delta)
            for row in evaluate(v, l, ref, tau=100):
                rows.append({'n':N, 'p':P, 'tau':100, 's':s, 'Delta':delta, **row})
            pd.DataFrame(rows).to_csv(power_path,index=False)
            print(f'Finished Delta={delta}, s={s}.', flush=True)
    power = pd.DataFrame(rows)
    boundary_path=RES/'boundary_results.csv'
    if boundary_path.exists():
        rows=pd.read_csv(boundary_path).to_dict('records')
    else:
        rows=power[(power.Delta==1.0)&power.s.isin([5,100])].to_dict('records')
    for ti, tau in enumerate((40, 160)):
        for si, s in enumerate((5, 100)):
            if any(row['tau']==tau and row['s']==s for row in rows):
                continue
            v, l = score_samples(seed_rng(3, ti, si), R_ALT, tau=tau, sparsity=s, delta=1.0)
            for row in evaluate(v,l,ref,tau=tau):
                rows.append({'n':N,'p':P,'tau':tau,'s':s,'Delta':1.0,**row})
            pd.DataFrame(rows).to_csv(boundary_path,index=False)
            print(f'Finished boundary tau={tau}, s={s}.',flush=True)
    boundary=pd.DataFrame(rows).sort_values(['s','tau','method'])
    boundary.to_csv(boundary_path,index=False)

    # Figure 1: theoretical drift per unit nonzero mean change, not observations.
    fig, ax = plt.subplots(figsize=(6.2,3.5))
    kk = np.arange(1,N)
    for tau, ls in ((40,'--'),(100,'-'),(160,':')):
        h = np.where(kk <= tau, (N-tau)*np.sqrt(kk/(N*(N-kk))),
                     tau*np.sqrt((N-kk)/(N*kk)))
        ax.plot(kk, h, linestyle=ls, label=fr'$\tau={tau}$')
    ax.set(xlabel='Candidate change-point location k',
           ylabel=r'CUSUM drift factor $d_\tau(k)$', xlim=(1,N-1))
    ax.legend(frameon=False)
    savefig(fig,'cusum_drift')

    # Figure 2: simulation with Delta=1; all methods share data per replicate.
    fig, ax = plt.subplots(figsize=(6.2,3.6))
    for method, marker, ls in zip(METHODS, ('o','s','^'), ('-','--',':')):
        d = power[(power.method==method)&(power.Delta==1.0)].sort_values('s')
        ax.errorbar(d.s,d.rejection_rate,yerr=d.mcse_conditional_on_reference,
                    marker=marker,linestyle=ls,capsize=3,label=method)
    ax.set(xscale='log',xlabel='Number of affected coordinates s (log scale)',
           ylabel='Empirical rejection probability',ylim=(-0.02,1.04))
    ax.set_xticks(SPARSITIES,labels=[str(s) for s in SPARSITIES])
    ax.legend(frameon=False,loc='best')
    savefig(fig,'power_by_sparsity')

    # Figure 3: joint detection/localization, no conditioning on rejection.
    fig, ax = plt.subplots(figsize=(6.2,3.6))
    for method, marker, ls in zip(METHODS, ('o','s','^'), ('-','--',':')):
        d=power[(power.method==method)&(power.Delta==1.0)].sort_values('s')
        ax.plot(d.s,d.joint_detect_within_5,marker=marker,linestyle=ls,label=method)
    ax.set(xscale='log',xlabel='Number of affected coordinates s (log scale)',
           ylabel='P(reject and localization error <= 5)',ylim=(-0.02,1.04))
    ax.set_xticks(SPARSITIES,labels=[str(s) for s in SPARSITIES]);ax.legend(frameon=False)
    savefig(fig,'joint_detection_localization')

    info = {'n':N,'p':P,'trim':TRIM,'alpha':ALPHA,'null_reference_draws':B_NULL,
            'independent_null_repetitions':R_NULL,'alternative_repetitions_per_cell':R_ALT,
            'seed':BASE_SEED,'deltas':DELTAS,'sparsities':SPARSITIES,
            'calibration':'Known N(0,I_p) null, independent Monte Carlo rank p-values; adaptive Bonferroni',
            'reference_reused':True,
            'uncertainty_note':'Displayed MCSE conditions on the simulated null reference; it excludes reference Monte Carlo error.',
            'python':platform.python_version(),'numpy':np.__version__,
            'pandas':pd.__version__,'matplotlib':matplotlib.__version__,
            'wall_seconds_this_invocation':time.perf_counter()-start}
    (RES/'simulation_metadata.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    print('All illustrative outputs saved.',flush=True)

if __name__ == '__main__':
    main()
