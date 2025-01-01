import numpy as np
import time
import tracemalloc
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Configure NumPy to raise exceptions on floating-point errors
np.seterr(over='raise', invalid='raise', divide='raise', under='warn')

# Kahan summation for improved floating-point precision
def kahan_sum(values):
    sum_ = np.float128(0.0)
    compensation = np.float128(0.0)
    for value in values:
        y = np.float128(value) - compensation
        temp = sum_ + y
        compensation = (temp - sum_) - y
        sum_ = temp
    return sum_

# 1D Brusselator PDE function using method of lines
def brusselator_1d(t, y, N, dx, D_u, D_v, a_param, b_param):
    # Split the combined y vector into u and v components
    u = y[:N]
    v = y[N:]

    # Initialize derivatives
    du_dt = np.zeros(N, dtype=np.float128)
    dv_dt = np.zeros(N, dtype=np.float128)

    # Second spatial derivatives (Laplacian) with Neumann boundary conditions
    # Improved approximation for second derivatives at boundaries
    u_xx = np.zeros(N, dtype=np.float128)
    v_xx = np.zeros(N, dtype=np.float128)

    # Internal points
    u_xx[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dx**2
    v_xx[1:-1] = (v[2:] - 2 * v[1:-1] + v[:-2]) / dx**2

    # Neumann boundary conditions (zero flux)
    u_xx[0] = (u[2] - 2 * u[1] + u[0]) / dx**2
    u_xx[-1] = (u[-1] - 2 * u[-2] + u[-3]) / dx**2
    v_xx[0] = (v[2] - 2 * v[1] + v[0]) / dx**2
    v_xx[-1] = (v[-1] - 2 * v[-2] + v[-3]) / dx**2

    # Reaction terms
    reaction_u = a_param - (b_param + 1) * u + u**2 * v
    reaction_v = b_param * u - u**2 * v

    # Combine diffusion and reaction
    du_dt = D_u * u_xx + reaction_u
    dv_dt = D_v * v_xx + reaction_v

    # Combine derivatives into a single vector
    dydt = np.concatenate([du_dt, dv_dt])

    return dydt

# General Runge-Kutta integration function
def runge_kutta_general(f, t_span, y0, h, a_values, b_values, args=()):
    t0, tf = t_span
    t = np.float128(t0)
    y = y0.copy()
    t_values = [t0]
    y_values = [y.copy()]

    s = len(b_values)  # Number of stages
    num_steps = int(np.ceil((tf - t0) / h))
    c_values = np.array([kahan_sum(a_values[i][:i]) for i in range(s)], dtype=np.float128)

    for step in range(num_steps):
        # Adjust the last step to end exactly at tf
        if t + h > tf + 1e-12:
            h = tf - t

        # If h is zero or negative due to floating-point errors, break the loop
        if h <= 0:
            break

        k = [np.zeros_like(y, dtype=np.float128) for _ in range(s)]
        for i in range(s):
            y_stage = y.copy()
            for j in range(i):
                y_stage += h * a_values[i][j] * k[j]
            t_stage = t + c_values[i] * h
            try:
                k[i] = f(t_stage, y_stage, *args)
            except FloatingPointError as e:
                print(f"Floating point error at time {t_stage}: {e}")
                return np.array(t_values, dtype=np.float128), np.array(y_values, dtype=np.float128)

        for i in range(s):
            y += h * b_values[i] * k[i]

        t += h

        # Only append if t is greater than the last t_value to avoid duplicates
        if t > t_values[-1] + 1e-12:
            t_values.append(t)
            y_values.append(y.copy())
        else:
            # Skip appending duplicate t_values
            pass

    return np.array(t_values, dtype=np.float128), np.array(y_values, dtype=np.float128)

# L2 norm error calculation
def l2_norm_error(y_numerical, y_reference):
    error = y_numerical - y_reference
    l2_norm = np.linalg.norm(error, axis=1)
    return np.sqrt(np.mean(l2_norm**2, dtype=np.float128))

# Function to calculate order of convergence
def calculate_order_of_convergence(errors, hs):
    orders = []
    for i in range(1, len(errors)):
        if errors[i-1] == 0 or hs[i-1] == hs[i]:
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

if __name__ == "__main__":
    # Parameters for convergence study
    t_span = (0, 10)  # Simulation time
    L = 10.0          # Length of the spatial domain
    N = 100           # Number of spatial points
    x = np.linspace(0, L, N)
    dx = x[1] - x[0]

    # Diffusion coefficients
    D_u = 0.1
    D_v = 0.05

    # Reaction parameters
    a_param = 1.0
    b_param = 3.0

    # Initial conditions with small perturbations
    np.random.seed(0)  # For reproducibility
    u0 = a_param + 0.1 * np.random.rand(N)
    v0 = b_param / a_param + 0.1 * np.random.rand(N)

    y0 = np.concatenate([u0, v0]).astype(np.float128)

    # Calculate maximum allowable time step for stability
    D_max = max(D_u, D_v)  # D_u = 0.1, D_v = 0.05
    dt_max = dx**2 / (2 * D_max)
    print(f"Maximum allowable Δt for stability: {dt_max}")

    # Adjust time step sizes accordingly
    hs = np.linspace(dt_max * 0.9, 0.001, 120)  # Start from 90% of dt_max

    errors = []
    iterations_list = []
    memory_usage = []
    computation_times = []

    # Define the coefficients for the 16-stage ESRK method
    s = 16  # Number of stages for your 16-stage ESRK method

    # Insert your b_values and a_values here
    b=[0.047706708943468935, 0.047706708943468935, 0.047706708943468935, 0.047706708943468935, 0.047706708943468935, 0.047706708943468935, 0.19110588425278838, -0.10217917299269486, -0.12795089861786924, 0.15817620320928832, 0.08069004693897402, -0.028653274616902827, 0.18231131798295502, 0.18109086580209066, 0.08612145658862838, 0.09304731779192842]

    b_values = np.array(b)
    
    a =[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.156231948809828, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0850384950103904, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.207712154404172, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, -1.66713860251474e-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, -0.204466488978545, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.226439386595127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, -0.228285639718228, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, 0.123390982023330, 0, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, 0.217065648281025, 0, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, -0.127950898617869, -0.207491798514979, 0, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, -0.127950898617869, 0.158176203209288, 0.000735658290826195, 0, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, -0.127950898617869, 0.158176203209288, 0.0806900469389740, 0.207272622131221, 0, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, -0.127950898617869, 0.158176203209288, 0.0806900469389740, -0.0286532746169028, 0.00150345391912373, 0, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, -0.127950898617869, 0.158176203209288, 0.0806900469389740, -0.0286532746169028, 0.182311317982955, 0.144902888444827, 0, 0], [0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.0477067089434689, 0.191105884252788, -0.102179172992695, -0.127950898617869, 0.158176203209288, 0.0806900469389740, -0.0286532746169028, 0.182311317982955, 0.181090865802091, 0.223463917740452, 0]]

    a_values = np.array(a)
   



    # Generate the reference solution with high precision using a very small step size
    print("Generating reference solution with h=0.0001...")
    t_values_ref, y_values_ref = runge_kutta_general(
        brusselator_1d, t_span, y0, 0.0001, a_values, b_values, args=(N, dx, D_u, D_v, a_param, b_param))
    print("Reference solution generated.")

    for h in hs:
        print(f"Running ESRK method with h={h}...")
        # Start memory tracking
        tracemalloc.start()

        # Start time measurement
        start_time = time.time()

        # Run numerical method
        try:
            t_values, y_values = runge_kutta_general(
                brusselator_1d, t_span, y0, h, a_values, b_values, args=(N, dx, D_u, D_v, a_param, b_param))
        except FloatingPointError as e:
            print(f"Floating point error encountered at h={h}: {e}")
            errors.append(np.nan)
            iterations_list.append(np.nan)
            computation_times.append(np.nan)
            memory_usage.append(np.nan)
            continue

        # Record computation time and memory usage
        computation_times.append(time.time() - start_time)
        current, peak = tracemalloc.get_traced_memory()
        memory_usage.append(peak)  # Peak memory usage
        tracemalloc.stop()

        # Interpolate solution
        y_values_interpolated = interpolate_solution(t_values, y_values, t_values_ref)

        # Calculate L2 norm error
        error = l2_norm_error(y_values_interpolated, y_values_ref)
        errors.append(error)
        iterations_list.append(len(t_values) - 1)  # Number of iterations

        print(f"h={h}: Error={error}, Iterations={len(t_values)-1}, Time={computation_times[-1]:.6f}s, Memory={peak} bytes")

    # Calculate order of convergence
    orders = calculate_order_of_convergence(errors, hs)

    # Print results
    print("\n2D Brusselator Simulation with 16-Stage Reduced ESRK Method with non-negative Coeffients")
    print("Step sizes:", hs)
    print("Errors:", errors)
    print("Orders of convergence:", orders)
    print("Iterations per step size:", iterations_list)
    print("Computation times (s):", computation_times)
    print("Memory usage (bytes):", memory_usage)

    # Plotting the convergence
    plt.figure(figsize=(8, 6))
    plt.loglog(hs, errors, 'o-', label='Numerical Error')
    # Reference line for expected order of convergence (update exponent as needed)
    expected_order = 4  # Update if your ESRK method has a different order
    plt.loglog(hs, [errors[0] * (h / hs[0]) ** expected_order for h in hs], 'k--', label=f'{expected_order}th Order')
    plt.xlabel('Time Step Size (Δt)')
    plt.ylabel('Error Norm')
    plt.title('Convergence Study of 2D Brusselator with 16-Stage ESRK reduced scheme postive Coeffients ')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()

