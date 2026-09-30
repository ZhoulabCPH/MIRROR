"""Optional Cox analysis and public survival metrics.

Median-based groups describe the supplied risk scores. The Cox wrapper returns
coefficient estimates and corresponding statistical summaries.
"""

import numpy as np
import pandas as pd

from utils.survival import CIndex, accuracy, accuracy_cox, cox_log_rank

__all__ = ["accuracy", "accuracy_cox", "cox_log_rank", "CIndex", "coxph_log_rank"]


def coxph_log_rank(survtime_all, labels, covariates):
    """Fit a univariate Cox model and return C-index, Wald p-value, HR, and CI.

    ``covariates`` supplies one continuous or binary covariate per patient.
    The reported confidence limits apply to its hazard ratio. A convergence
    failure yields NaN values, since no fitted estimate is available. The
    p-value corresponds to the Cox coefficient Wald statistic.
    """
    from lifelines import CoxPHFitter
    from lifelines.exceptions import ConvergenceError

    data = pd.DataFrame({
        "covariates": np.asarray(covariates, dtype=float).reshape(-1),
        "survtime": np.asarray(survtime_all, dtype=float).reshape(-1),
        "event": np.asarray(labels).reshape(-1),
    })
    if not data["event"].isin([0, 1]).all():
        raise ValueError("Event labels must contain only zero and one.")
    if not np.isfinite(data[["covariates", "survtime"]].to_numpy()).all():
        raise ValueError("Covariates and survival times must be finite.")
    try:
        model = CoxPHFitter()
        model.fit(data, duration_col="survtime", event_col="event")
    except ConvergenceError:
        return (float("nan"),) * 5
    intervals = model.confidence_intervals_.loc["covariates"]
    return (
        float(model.concordance_index_),
        float(model.summary.loc["covariates", "p"]),
        float(model.hazard_ratios_.loc["covariates"]),
        float(np.exp(intervals.iloc[0])),
        float(np.exp(intervals.iloc[1])),
    )
