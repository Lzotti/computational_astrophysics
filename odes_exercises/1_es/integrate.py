import numpy as np

def integrate(f, t_range, IC, h, args=None, method=None, extra_par=None):

    t0, t_end = float(t_range[0]), float(t_range[1])
    t = t0
    y = np.array(IC)

    if h == 0:
        raise ValueError("h must be non-zero")

    direction = np.sign(t_end - t0)
    if direction == 0:
        return np.array([t]), np.array([y])
    if np.sign(h) != direction:
        raise ValueError("h has wrong sign for the given t_range")

    t_array = [t]
    y_array = [y]

    # Tolleranza, rivedo l'espressione se è sensata o se è meglio il paragone con h
    tol = 10.0 * np.finfo(float).eps * max(1.0, abs(t0), abs(t_end))

    while (t_end - t) * direction > tol:

        # Use a reduced last step so the final time is exactly t_end.
        h_step = min(abs(h), abs(t_end - t)) * direction

        y_next = method(t, y, h_step, f, args, extra_par=extra_par)

        t += h_step
        if abs(t_end - t) <= tol:
            t = t_end
        y = np.array(y_next)

        t_array.append(t)
        y_array.append(y)
    
    t_array = np.array(t_array)
    y_array = np.array(y_array)

    return t_array, y_array