"""Compute classification and survival metrics from supplied predictions."""

import numpy as np


def accuracy(output, labels):
    """Return top-class accuracy for a tensor of class scores and integer labels."""
    if len(labels) == 0:
        return float("nan")
    predictions = output.argmax(dim=1).type_as(labels)
    return predictions.eq(labels).double().mean()


def _survival_arrays(hazards, labels, survival_times=None):
    """Convert paired survival inputs to one-dimensional numeric arrays."""
    risks = np.asarray(hazards, dtype=float).reshape(-1)
    events = np.asarray(labels).reshape(-1)
    if risks.size != events.size:
        raise ValueError("Risk scores and event indicators must have equal lengths.")
    if not np.isfinite(risks).all() or not np.isin(events, [0, 1]).all():
        raise ValueError("Risk scores must be finite and event indicators must be binary.")
    if survival_times is None:
        return risks, events
    times = np.asarray(survival_times, dtype=float).reshape(-1)
    if times.size != risks.size or not np.isfinite(times).all():
        raise ValueError("Survival times must be finite and match the risk-score length.")
    return risks, events, times


def accuracy_cox(hazardsdata, labels):
    """Compare median-dichotomized risks with observed event indicators.

    This descriptive statistic does not account for censoring or follow-up time.
    """
    risks, events = _survival_arrays(hazardsdata, labels)
    if risks.size == 0:
        return float("nan")
    predicted_events = risks > np.median(risks)
    return float(np.mean(predicted_events == events))


def cox_log_rank(hazardsdata, labels, survtime_all):
    """Return the log-rank p-value between median-defined risk groups."""
    from lifelines.statistics import logrank_test

    risks, events, times = _survival_arrays(hazardsdata, labels, survtime_all)
    if risks.size == 0:
        return float("nan")
    low_risk = risks <= np.median(risks)
    if low_risk.all() or not low_risk.any():
        return float("nan")
    result = logrank_test(
        times[low_risk], times[~low_risk],
        event_observed_A=events[low_risk], event_observed_B=events[~low_risk],
    )
    return float(result.p_value)


def concordance_index_manual(hazards, labels, survtime_all):
    """Score strictly ordered comparable pairs, with half credit for tied risks.

    Higher risk predicts an earlier event. A pair is comparable when its earlier
    observation records an event. Equal survival times are excluded. A dataset
    with no comparable pairs returns NaN.
    """
    risks, events, times = _survival_arrays(hazards, labels, survtime_all)
    concordant = 0.0
    comparable = 0
    for index in np.flatnonzero(events == 1):
        later = times > times[index]
        comparable += int(later.sum())
        concordant += float(np.sum(risks[later] < risks[index]))
        concordant += 0.5 * float(np.sum(risks[later] == risks[index]))
    return concordant / comparable if comparable else float("nan")


def concordance_index_lifelines(hazards, labels, survtime_all):
    """Return lifelines concordance with higher scores indicating greater risk."""
    from lifelines.utils import concordance_index

    risks, events, times = _survival_arrays(hazards, labels, survtime_all)
    if risks.size == 0:
        return float("nan")
    try:
        return float(concordance_index(times, -risks, events))
    except ZeroDivisionError:
        return float("nan")


# These public aliases expose the metric names used by existing model callers.
CIndex = concordance_index_manual
CIndex_lifeline = concordance_index_lifelines
