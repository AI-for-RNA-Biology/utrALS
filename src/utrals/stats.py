"""Small statistics helpers used by the notebooks."""
import numpy as np


def find_elbow(values):
    values = np.asarray(values, dtype=float)
    x = np.arange(1, len(values) + 1, dtype=float)

    first_point = np.array([x[0], values[0]])
    last_point = np.array([x[-1], values[-1]])

    line_vec = last_point - first_point
    line_vec = line_vec / np.linalg.norm(line_vec)

    points = np.column_stack((x, values))
    vec_from_first = points - first_point
    projection = np.outer(vec_from_first @ line_vec, line_vec) + first_point
    distances = np.linalg.norm(points - projection, axis=1)

    elbow_idx = int(np.argmax(distances))
    return elbow_idx, distances


## extract the outliers
def get_outliers(serie):
    """Get outliers of a ps.Series.
    Args:
        serie (pd.Series): pd.Serie with numerical values

    Returns:
        ps.Series: Two series with the lower outliers and upper outliers.
    """
    Q1, Q3 = serie.quantile([0.25, 0.75]).values
    lower_bound = Q1 - 1.5*(Q3-Q1)
    upper_bound = Q3 + 1.5*(Q3-Q1)
    return serie[serie < lower_bound], serie[serie > upper_bound]
