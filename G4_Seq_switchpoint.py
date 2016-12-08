import pymc as pm
import numpy as np
import sys
from tqdm import tqdm

# def get_switchpoint(data):

#     def make_model(data, lower, upper):
#         switchpoint = pm.DiscreteUniform('switchpoint', lower=lower,
#                                          upper=upper)
#         early_rate = pm.Exponential('early_rate', 1)
#         late_rate = pm.Exponential('late_rate', 1)

#         @pm.deterministic(plot=False)
#         def rate(s=switchpoint, early=early_rate, late=late_rate):
#             out = np.empty(len(data))
#             out[:s] = early
#             out[s:] = late
#             return out

#         phredscore = pm.Poisson('phredscore', mu=rate, value=data,
#                                 observed=True)
#         return locals()

#     bases = np.array(list(range(len(data))))

#     M = pm.MCMC(make_model(data, bases.min(), bases.max()))
#     M.sample(iter=6000, burn=1000, progress_bar=False)
#     return np.array([
#         M.trace('switchpoint').gettrace()[-1000:],
#         M.trace('early_rate').gettrace()[-1000:],
#         M.trace('late_rate').gettrace()[-1000:],
#     ])


def get_switchpoint(data):

    def make_model(data):
        switchpoint = pm.DiscreteUniform('switchpoint', 0, len(data))
        early_rate = pm.Beta('early_rate', 0.5, 0.5)
        late_rate = pm.Beta('late_rate', 0.5, 0.5)

        @pm.deterministic(plot=False)
        def rate(s=switchpoint, early=early_rate, late=late_rate):
            out = np.empty(len(data))
            out[:s] = early
            out[s:] = late
            return out

        phredscore = pm.Bernoulli('phredscore', p=rate, value=data,
                                  observed=True)
        return locals()

    M = pm.MCMC(make_model(data))
    M.sample(iter=10000, burn=1000, progress_bar=False)
    return np.array([
        M.trace('switchpoint').gettrace()[-2000:],
        M.trace('early_rate').gettrace()[-2000:],
        M.trace('late_rate').gettrace()[-2000:],
    ])
    return M


if __name__ == "__main__":
    data_file = sys.argv[1]
    out_file = sys.argv[2]
    traces = []
    for data in tqdm(np.load(data_file)):
        traces.append(get_switchpoint(data))

    np.save(out_file, np.array(traces))
