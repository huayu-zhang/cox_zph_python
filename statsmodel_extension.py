import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats.distributions import chi2
from statsmodels.duration.hazard_regression import PHRegResults


def cox_ph_test(self):
    """
    The function is the equivalent of R function survival::cox.zph
    :param self: PHRegResults
    :return: dict of {'rho':[one per exog], 'z': [one per exog], 'p_value': [one per exog]}
    """

    with_events = ~np.isnan(self.schoenfeld_residuals).any(axis=1)
    n_events = sum(with_events)

    sresid = self.schoenfeld_residuals[with_events, :]
    ttimes = self.model.endog[with_events].argsort().argsort()

    xx = ttimes - np.mean(ttimes)
    r2 = sresid @ self.cov_params() * n_events
    test = xx @ r2

    corel = [rho for rho, _ in
             [pearsonr(col, xx) for col in r2.T]]

    z = test ** 2 / (np.diag(self.cov_params()) * n_events * sum(xx ** 2))
    pz_chi2 = chi2.sf(z, 1)

    self.cox_ph_test_results = {'rho': corel, 'z': z, 'p_value': pz_chi2}

    return self.cox_ph_test_results


PHRegResults.cox_ph_test = cox_ph_test
