"""Throughput scaling test: independent forward-model loops on N GPUs.

Every MPI rank creates its own BeAR instance on its own GPU and evaluates the
*same* sequence of random prior draws (same seed on all ranks), so every GPU
gets an identical workload. The ranks start and stop at MPI barriers; the
aggregate throughput is N * n_calls / (slowest rank's wall time). No data is
exchanged during the timed loop, so this measures the hardware throughput
ceiling, not the scaling of a sampler-driven retrieval.

Usage (one task per GPU, no master rank):
    srun python tools/scaling_test_mpi.py <retrieval_folder> [options]

Options:
    --calls N        timed likelihood calls per rank      (default 100000)
    --warmup N       untimed calls per rank beforehand    (default 2000);
                     -1 = one full untimed pass over the timed set (cleanest)
    --seed S         RNG seed for the prior draws         (default 1)
    --results FILE   append one summary line per run to FILE (CSV)
    --label TEXT     free-text tag written into the results line

Warm-up matters: the first evaluations in a new temperature range load native
opacity files (~100 ms each), and that must not bleed into the timed loop. The
script reports the rate of the second half of the loop separately; if it is
noticeably higher than the full-loop rate, increase --warmup.
"""
import argparse
import os
import socket
import subprocess
import sys
import time

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except ImportError:
    # single-process fallback (workstation without MPI): same code path, N = 1
    class _SerialComm:
        def Get_rank(self): return 0
        def Get_size(self): return 1
        def Barrier(self): pass
        def gather(self, obj, root=0): return [obj]

    class MPI:  # noqa: N801
        COMM_WORLD = _SerialComm()
        Wtime = staticmethod(time.perf_counter)

    comm = MPI.COMM_WORLD
    print("mpi4py not found: running as a single process")

# ----------------------------------------------------------------------------
# GPU assignment: must happen before the first CUDA call, i.e. before the bear
# module is imported. If the scheduler already restricted CUDA_VISIBLE_DEVICES
# to one device per task (e.g. --gpus-per-task=1), that is respected.
# ----------------------------------------------------------------------------
rank = comm.Get_rank()
size = comm.Get_size()

local_rank = int(
    os.environ.get('OMPI_COMM_WORLD_LOCAL_RANK')
    or os.environ.get('MV2_COMM_WORLD_LOCAL_RANK')
    or os.environ.get('SLURM_LOCALID')
    or rank
)

visible = os.environ.get('CUDA_VISIBLE_DEVICES')
if visible:
    gpu_ids = visible.split(',')
else:
    gpu_list = subprocess.run(
        ['nvidia-smi', '-L'], capture_output=True, text=True, check=True
    ).stdout
    gpu_ids = [str(i) for i in range(len(gpu_list.strip().splitlines()))]

gpu_id = gpu_ids[local_rank % len(gpu_ids)]
os.environ['CUDA_VISIBLE_DEVICES'] = gpu_id

try:
    gpu_name = subprocess.run(
        ['nvidia-smi', '-i', gpu_id, '--query-gpu=name', '--format=csv,noheader'],
        capture_output=True, text=True, check=True
    ).stdout.strip()
except Exception:
    gpu_name = 'unknown'

try:
    nb_cpus = len(os.sched_getaffinity(0))
except AttributeError:
    nb_cpus = os.cpu_count()

# ----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('retrieval_folder')
parser.add_argument('--calls', type=int, default=100000)
parser.add_argument('--warmup', type=int, default=2000)
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--results', default=None)
parser.add_argument('--label', default='')
args = parser.parse_args()

bear_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(bear_root, 'python'))

from lib import bear  # noqa: E402
import numpy as np   # noqa: E402

retrieval_folder = os.path.abspath(args.retrieval_folder) + '/'

# ----------------------------------------------------------------------------
# Setup (untimed): model, identical parameter sequence on every rank
# ----------------------------------------------------------------------------
t_setup0 = time.perf_counter()

model_config = bear.Config(retrieval_folder)
model = bear.Retrieval(model_config)
nb_param = model.nbParameters()

rng = np.random.default_rng(args.seed)
n_warmup = max(args.warmup, 0)
cubes = rng.uniform(0.05, 0.95, size=(n_warmup + args.calls, nb_param))
phys = [np.array(model.convertCubeParameters(c)[1], dtype=float) for c in cubes]

timed = phys[n_warmup:]
# --warmup -1: warm up on the timed set itself, so every opacity file the
# timed loop will touch is already loaded (doubles the runtime, cleanest)
warmup_set = phys[:n_warmup] if args.warmup >= 0 else timed

t_setup = time.perf_counter() - t_setup0

# warm-up: loads opacity files, fills caches, JIT etc.
t_warm0 = time.perf_counter()
for p in warmup_set:
    model.computeLikelihood(p)
t_warm = time.perf_counter() - t_warm0

# ----------------------------------------------------------------------------
# Timed loop between two barriers
# ----------------------------------------------------------------------------
half = len(timed) // 2

comm.Barrier()
t_global0 = MPI.Wtime()

t0 = time.perf_counter()
logl = np.empty(len(timed))
for i, p in enumerate(timed):
    logl[i] = model.computeLikelihood(p)
    if i == half:
        t_half = time.perf_counter()
t_local = time.perf_counter() - t0

comm.Barrier()
t_global = MPI.Wtime() - t_global0

t_second_half = (t0 + t_local) - t_half
finite = int(np.isfinite(logl).sum())

record = dict(
    rank=rank,
    host=socket.gethostname(),
    gpu=gpu_id,
    gpu_name=gpu_name,
    cpus=nb_cpus,
    t_setup=t_setup,
    t_warm=t_warm,
    t_local=t_local,
    rate_local=args.calls / t_local,
    rate_second_half=(len(timed) - half) / t_second_half,
    finite=finite,
    logl_median=float(np.nanmedian(logl)),
)

records = comm.gather(record, root=0)

if rank == 0:
    rates = np.array([r['rate_local'] for r in records])
    rates_2nd = np.array([r['rate_second_half'] for r in records])
    t_locals = np.array([r['t_local'] for r in records])
    aggregate = size * args.calls / t_global

    print(f"\nBeAR forward-model scaling test: {retrieval_folder}")
    print(f"ranks (GPUs): {size}   parameters: {nb_param}   "
          f"calls/rank: {args.calls}   warm-up/rank: {args.warmup}   seed: {args.seed}\n")

    print(f"{'rank':>4} {'host':<16} {'gpu':>3} {'name':<24} {'cpus':>4} "
          f"{'setup[s]':>9} {'warm[s]':>8} {'loop[s]':>9} "
          f"{'calls/s':>9} {'2nd half':>9} {'ms/call':>8} {'finite':>7} {'med logL':>10}")
    for r in records:
        print(f"{r['rank']:>4} {r['host'][:16]:<16} {r['gpu']:>3} {r['gpu_name'][:24]:<24} "
              f"{r['cpus']:>4} {r['t_setup']:>9.1f} {r['t_warm']:>8.1f} {r['t_local']:>9.1f} "
              f"{r['rate_local']:>9.1f} {r['rate_second_half']:>9.1f} "
              f"{1e3 / r['rate_local']:>8.3f} {r['finite']:>7} {r['logl_median']:>10.1f}")

    print()
    print(f"per-rank rate   : mean {rates.mean():.1f}  min {rates.min():.1f}  "
          f"max {rates.max():.1f} calls/s   (spread {100 * (rates.max() / rates.min() - 1):.1f} %)")
    print(f"second-half rate: mean {rates_2nd.mean():.1f} calls/s   "
          f"({100 * (rates_2nd.mean() / rates.mean() - 1):+.1f} % vs full loop; "
          f"large positive => increase --warmup)")
    print(f"slowest rank    : {t_locals.max():.1f} s   barrier-to-barrier: {t_global:.1f} s")
    print(f"AGGREGATE       : {aggregate:.1f} calls/s on {size} GPU(s)   "
          f"= {aggregate / size:.1f} calls/s per GPU")
    print("Ideal scaling = N x (per-GPU rate of the N=1 run); compare across jobs.")

    if args.results:
        new_file = not os.path.exists(args.results)
        with open(args.results, 'a') as f:
            if new_file:
                f.write('# n_gpus, calls_per_rank, aggregate_calls_per_s, '
                        'per_gpu_calls_per_s, mean_rank_rate, min_rank_rate, '
                        'max_rank_rate, second_half_rate, wall_s, gpu_name, '
                        'cpus_per_rank, retrieval_folder, label\n')
            f.write(f"{size}, {args.calls}, {aggregate:.2f}, {aggregate / size:.2f}, "
                    f"{rates.mean():.2f}, {rates.min():.2f}, {rates.max():.2f}, "
                    f"{rates_2nd.mean():.2f}, {t_global:.2f}, "
                    f"{records[0]['gpu_name']}, {records[0]['cpus']}, "
                    f"{retrieval_folder}, {args.label}\n")
        print(f"appended summary to {args.results}")
