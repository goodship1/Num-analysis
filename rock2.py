#Implementation of the ROCK2 paper on the 1D brusselator with convergence studies 
import numpy as np
import time
import tracemalloc
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Define the 2D Brusselator system
def brusselator_2d(t, y):
    x, y_ = y
    a, b = 1.0, 1.0
    dx_dt = a - (b + 1) * x + x**2 * y_
    dy_dt = b * x - x**2 * y_
    return np.array([dx_dt, dy_dt], dtype=np.float64)

# ROCK2 Method with adjustable M (number of stages)
def rock2(f, t_span, y0, h, m=2):
    """
    ROCK2 Method with adjustable M (number of stages).
    
    Parameters:
    f : callable
        The ODE system to solve.
    t_span : tuple
        The time span (t0, tf) for integration.
    y0 : ndarray
        The initial conditions.
    h : float
        The time step size.
    m : int
        Number of stages (default is 2).
    """
    t0, tf = t_span
    t = t0
    y = y0.copy()
    t_values = [t]
    y_values = [y.copy()]
    
    # Stability-enhancing coefficients
    a1, a2 = 1 - 1 / m, 1 / m
    b1, b2 = 0.5 * (1 - a1), 0.5 * (1 - a2)
    
    num_steps = int(np.ceil((tf - t0) / h))
    
    for _ in range(num_steps):
        k1 = f(t, y)
        y_star = y + h * a1 * k1
        k2 = f(t + h * a2, y_star)
        y += h * (b1 * k1 + b2 * k2)
        t += h
        t_values.append(t)
        y_values.append(y.copy())
    
    return np.array(t_values), np.array(y_values)

# L2 norm error calculation
def l2_norm_error(y_numerical, y_reference):
    error = y_numerical - y_reference
    l2_norm = np.linalg.norm(error, axis=1)
    return np.sqrt(np.mean(l2_norm**2))

# Function to calculate order of convergence
def calculate_order_of_convergence(errors, hs):
    orders = []
    for i in range(1, len(errors)):
        if errors[i-1] == 0:
            orders.append(np.nan)
        else:
            order = np.log(errors[i] / errors[i-1]) / np.log(hs[i] / hs[i-1])
            orders.append(order)
    return orders

# Interpolate solution for error calculation
def interpolate_solution(t_values, y_values, t_values_ref):
    y_interpolated = []
    for i in range(y_values.shape[1]):
        f_interp = interp1d(t_values, y_values[:, i], kind='cubic', fill_value="extrapolate")
        y_interpolated.append(f_interp(t_values_ref))
    return np.array(y_interpolated).T

# Convergence study for ROCK2
def rock2_convergence_study(m_values, t_span, y0, h_values):
    for m in m_values:
        errors = []
        iterations_list = []
        computation_times = []
        memory_usage = []

        print(f"\nGenerating reference solution for ROCK2 with M={m} and h=0.001...")
        t_ref, y_ref = rock2(brusselator_2d, t_span, y0, 0.0001, m=m)
        print("Reference solution generated.")

        for h in h_values:
            print(f"Running ROCK2 with M={m} and h={h}...")
            tracemalloc.start()
            start_time = time.time()

            # Run ROCK2
            t, y = rock2(brusselator_2d, t_span, y0, h, m=m)
            
            computation_times.append(time.time() - start_time)
            current, peak = tracemalloc.get_traced_memory()
            memory_usage.append(peak)
            tracemalloc.stop()
            
            # Interpolate solution for error calculation
            y_interp = interpolate_solution(t, y, t_ref)
            error = l2_norm_error(y_interp, y_ref)
            errors.append(error)
            iterations_list.append(len(t) - 1)

        # Calculate order of convergence
        orders = calculate_order_of_convergence(errors, h_values)

        # Print results
        print(f"\nROCK2 Method with M={m} on 2D Brusselator")
        print(f"Step sizes: {h_values}")
        print(f"Errors: {errors}")
        print(f"Orders of convergence: {orders}")
        print(f"Iterations per step size: {iterations_list}")
        print(f"Computation times (s): {computation_times}")
        print(f"Memory usage (bytes): {memory_usage}")

        # Plot convergence
        plt.figure(figsize=(10, 6))
        plt.loglog(h_values, errors, 'o-', label=f"Numerical Error reference")
        plt.loglog(h_values, [errors[0] * (h / h_values[0])**2 for h in h_values], 'k--', label="Second Order")
        plt.xlabel("Step Size (h)")
        plt.ylabel("Error Norm")
        plt.title(f"Convergence Study: ROCK2 on 1D Brusselator")
        plt.legend()
        plt.grid(True, which="both", ls="--")
        plt.show()

# Main function
if __name__ == "__main__":
    t_span = (0, 100000)
    y0 = np.array([1.2, 2.5])
    h_values = np.linspace(0.1, 0.001, 100)  # Step sizes for convergence study
    m_values = [2]  # Different numbers of stages to test

    # Perform convergence study for ROCK2 with different M values
    rock2_convergence_study(m_values, t_span, y0, h_values)

