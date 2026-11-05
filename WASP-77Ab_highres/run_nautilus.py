import os
import sys

bear_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(bear_root, 'python'))

from lib import bear
import numpy as np
from nautilus import Sampler


retrieval_folder = os.path.dirname(os.path.realpath(__file__)) + '/'

print(f"Loading config from: {retrieval_folder}")
model_config = bear.Config(retrieval_folder)
model = bear.Retrieval(model_config)

nb_param = model.nbParameters()
print(f"Number of parameters: {nb_param}")


def prior(cube):
    converted = model.convertCubeParameters(np.array(cube, dtype=float))
    return np.array(converted[1], dtype=float)  # physical_parameters


def loglike(physical_params):
    return float(model.computeLikelihood(np.array(physical_params, dtype=float)))


print("Starting retrieval\n")

sampler = Sampler(
    prior,
    loglike,
    n_dim=nb_param,
    n_live=3000,
    filepath=retrieval_folder + 'nautilus.hdf5',
)

sampler.run(verbose=True)

points, log_w, log_l = sampler.posterior(equal_weight=True)

posteriors = [model.convertToParameterUnits(list(s)) for s in points]
log_likelihoods = list(log_l)

with open(retrieval_folder + 'post_equal_weights_nautilus.dat', 'w') as f:
    for params, logl in zip(posteriors, log_likelihoods):
        f.write('  '.join(f'{v:.6e}' for v in params) + f'  {logl:.6e}\n')

post_process = bear.PostProcess(model_config)
print("Starting post process\n")
post_process.run(posteriors, log_likelihoods)
