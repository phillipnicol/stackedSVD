import numpy as np
from typing import List, Tuple
from tqdm import tqdm
import pickle
from datetime import datetime
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

from typing import Tuple
import torch
import pandas as pd
import os
import argparse



def largest_sv_vector_power_method(
    X: np.ndarray,
    max_iter: int = 10_000,
    tol: float = 1e-8,
    dtype: torch.dtype = torch.bfloat16,  # Changed to bfloat16
    device=None,
) -> Tuple[np.ndarray, float]:
    """
    GPU-optimized power iteration that returns the dominant
    right singular vector v and corresponding singular value s.
    
    This version uses bfloat16 and a pre-transposed matrix for efficiency.
    """
    
    
    # --- device / dtype housekeeping ----------------------------------------
    if device is None:
        device = "cuda:0" if torch.cuda.is_available() else "cpu"
    
    ######## if min dim < 1000, use torch svd
    if X.shape[0] <= 5000 or X.shape[1] <= 5000:
        X_torch = torch.from_numpy(X).to(device=device, dtype=torch.float32).contiguous()
        # Use torch's SVD for small matrices
        U, S, Vt = torch.linalg.svd(X_torch, full_matrices=False)
        return Vt[0].cpu().numpy(), S[0].item()

    
    
    # Move data to GPU once and create a transposed copy
    X_torch = torch.from_numpy(X).to(device=device, dtype=dtype).contiguous()
    # Pre-transpose X for efficient X.T @ u operation
    X_t_torch = X_torch.T.contiguous()
    
    m, n = X_torch.shape
    
    # --- pre-allocate tensors -----------------------------------------------
    v = torch.randn(n, device=device, dtype=dtype)
    v.div_(torch.norm(v))
    
    # --- power iteration (JIT-compiled) -------------------------------------
    @torch.jit.script
    def power_iter_step(X: torch.Tensor, Xt: torch.Tensor, v: torch.Tensor):
        # 1. u = X @ v
        u = torch.matmul(X, v)
        # 2. v_new = X.T @ u
        v_new = torch.matmul(Xt, u)
        
        norm_v_new = torch.norm(v_new)
        if norm_v_new > 1e-14:
            v_new.div_(norm_v_new)
        return v_new, norm_v_new

    # --- main iteration loop ------------------------------------------------
    with torch.no_grad():
        for i in range(max_iter):
            v_old = v.clone()
            v, norm_v_new = power_iter_step(X_torch, X_t_torch, v)

            if norm_v_new < 1e-14:
                # print("Numerical underflow, stopping iteration.")
                break
            
            # Check for convergence
            if torch.norm(v - v_old) < tol:
                # print(f"Converged after {i+1} iterations.")
                break
    
    # --- Final singular value calculation -----------------------------------
    with torch.no_grad():
        # u = X @ v is already the first step of the iteration
        u_final = torch.matmul(X_torch, v)
        s = torch.norm(u_final).item()
        
    return v.to(dtype=torch.float32).cpu().numpy(), s
    
    

def generate_data(d: int, c_arr: List[float], theta_arr: List[float]) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Generate synthetic data according to the model X_i = theta_i * u_i * v^T + E_i
    
    Args:
        d: Column dimension (shared across matrices)
        c_arr: List of ratios n_i / d for each matrix
        theta_arr: List of signal strengths
        
    Returns:
        X_list: List of generated matrices
        v: True shared component
    """
    M = len(c_arr)
    assert len(theta_arr) == M, "theta_arr must have same length as c_arr"
    
    assert all(c > 0 for c in c_arr), "All elements of c_arr must be positive"
    
    # Generate true v
    v = np.random.normal(size=d)
    v /= np.linalg.norm(v)
    
    X_list = []
    for i in range(M):
        n_i = int(c_arr[i] * d)
        u_i = np.random.normal(size=n_i)
        u_i /= np.linalg.norm(u_i)
        
        # Generate noise matrix
        E_i = np.random.normal(0, 1/np.sqrt(d), size=(n_i, d))
        
        # Generate data matrix
        X_i = theta_arr[i] * np.outer(u_i, v) + E_i
        X_list.append(X_i)
    
    return X_list, v


def generate_data_exponential(d: int, c_arr: List[float], theta_arr: List[float]) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Generate synthetic data according to the model X_i = theta_i * u_i * v^T + E_i
    
    Args:
        d: Column dimension (shared across matrices)
        c_arr: List of ratios n_i / d for each matrix
        theta_arr: List of signal strengths
        
    Returns:
        X_list: List of generated matrices
        v: True shared component
    """
    M = len(c_arr)
    assert len(theta_arr) == M, "theta_arr must have same length as c_arr"
    
    assert all(c > 0 for c in c_arr), "All elements of c_arr must be positive"
    
    # Generate true v
    v = np.random.normal(size=d)
    v /= np.linalg.norm(v)
    
    X_list = []
    for i in range(M):
        n_i = int(c_arr[i] * d)
        u_i = np.random.normal(size=n_i)
        u_i /= np.linalg.norm(u_i)
        
        # Generate noise matrix with centered exponential distribution
        E_i = np.random.exponential(scale=1/np.sqrt(d), size=(n_i, d)) - 1/np.sqrt(d)
        
        # Generate data matrix
        X_i = theta_arr[i] * np.outer(u_i, v) + E_i
        X_list.append(X_i)
    
    return X_list, v


def compute_optimal_weights_stacksvd(theta_arr: List[float], c_arr: List[float]) -> np.ndarray:
    """
    Compute optimal weights for the Stack-SVD method, proportional to:
        w_i ∝ theta_i / sqrt(theta_i^2 + c_i)

    Args:
        theta_arr: List of signal strengths for each matrix.
        c_arr: List of ratios n_i / d for each matrix.

    Returns:
        Normalized weight vector as a NumPy array.
    """
    weights = np.array([theta / np.sqrt(theta**2 + c) for theta, c in zip(theta_arr, c_arr)])
    return weights / np.linalg.norm(weights)


def compute_optimal_weights_svdstack(theta_arr: List[float], c_arr: List[float]) -> np.ndarray:
    """
    Compute optimal weights for the SVD-Stack method, proportional to:
        w_i ∝ theta_i * sqrt((theta_i^2 + 1) / (theta_i^2 + c_i))

    Args:
        theta_arr: List of signal strengths for each matrix.
        c_arr: List of ratios n_i / d for each matrix.

    Returns:
        Normalized weight vector as a NumPy array.
    """
    weights = np.array([theta * np.sqrt((theta**2 + 1)/(theta**2 + c)) 
                        for theta, c in zip(theta_arr, c_arr)])
    return weights / np.linalg.norm(weights)


def estimate_stacksvd(X_list: List[np.ndarray], weights: np.ndarray = None) -> np.ndarray:
    """
    Estimate the principal right singular vector using the Stack-SVD method.
    Optionally applies weights to each matrix before stacking.

    Args:
        X_list: List of data matrices.
        weights: Optional array of weights for each matrix. If None, uses equal weights.

    Returns:
        Estimated principal right singular vector as a NumPy array.
    """
    if weights is None:
        weights = np.ones(len(X_list)) / np.sqrt(len(X_list))
        
    # Stack matrices with weights
    X_stacked = np.vstack([w * X for w, X in zip(weights, X_list)])
    
    # Compute principal right singular vector using power method
    v_est, _ = largest_sv_vector_power_method(X_stacked)
    return v_est

def estimate_svdstack(X_list: List[np.ndarray], weights: np.ndarray = None) -> np.ndarray:
    """
    Estimate the principal right singular vector using the SVD-Stack method.
    Optionally applies weights to each singular vector before stacking.

    Args:
        X_list: List of data matrices.
        weights: Optional array of weights for each matrix. If None, uses equal weights.

    Returns:
        Estimated principal right singular vector as a NumPy array.
    """
    if weights is None:
        weights = np.ones(len(X_list)) / np.sqrt(len(X_list))
    
    # Get principal right singular vectors for each matrix
    V_list = []
    for X in X_list:
        v_i, _ = largest_sv_vector_power_method(X)
        V_list.append(v_i)
    
    # Stack vectors with weights
    V_stacked = np.vstack([w * v.reshape(1, -1) for w, v in zip(weights, V_list)])
    
    # Compute principal right singular vector of stacked vectors
    v_est, _ = largest_sv_vector_power_method(V_stacked)
    return v_est

def run_simulation(c_arr: List[float], d: int, theta_arr: List[float], num_trials: int = 100) -> dict:
    """
    Run simulation comparing the four methods.
    
    Args:
        c_arr: List of ratios n_i / d for each matrix
        d: Column dimension
        theta_arr: List of signal strengths
        num_trials: Number of simulation trials
        
    Returns:
        Dictionary containing simulation results
    """

    
    # Initialize arrays to store results
    results = {
        'unweighted_stacksvd': np.zeros(num_trials),
        'weighted_stacksvd': np.zeros(num_trials),
        'unweighted_svdstack': np.zeros(num_trials),
        'weighted_svdstack': np.zeros(num_trials)
    }
    
    # Compute optimal weights
    w_stacksvd = compute_optimal_weights_stacksvd(theta_arr, c_arr)
    w_svdstack = compute_optimal_weights_svdstack(theta_arr, c_arr)
    
    for trial in tqdm(range(num_trials), desc="Running trials"):
        # Generate data
        X_list, v_true = generate_data(d, c_arr, theta_arr)
        
        # Run all four methods
        v_stacksvd = estimate_stacksvd(X_list)
        v_stacksvd_w = estimate_stacksvd(X_list, w_stacksvd)
        v_svdstack = estimate_svdstack(X_list)
        v_svdstack_w = estimate_svdstack(X_list, w_svdstack)
        
        # Store results (squared correlation)
        results['unweighted_stacksvd'][trial] = (v_stacksvd @ v_true)**2
        results['weighted_stacksvd'][trial] = (v_stacksvd_w @ v_true)**2
        results['unweighted_svdstack'][trial] = (v_svdstack @ v_true)**2
        results['weighted_svdstack'][trial] = (v_svdstack_w @ v_true)**2
    
    return results

def run_simulation_exponential(c_arr: List[float], d: int, theta_arr: List[float], num_trials: int = 100) -> dict:
    """
    Run simulation comparing the four methods with exponential noise.
    
    Args:
        c_arr: List of ratios n_i / d for each matrix
        d: Column dimension
        theta_arr: List of signal strengths
        num_trials: Number of simulation trials
        
    Returns:
        Dictionary containing simulation results
    """
    M = len(c_arr)
    n_arr = np.array(c_arr) * d
    
    # Initialize arrays to store results
    results = {
        'unweighted_stacksvd': np.zeros(num_trials),
        'weighted_stacksvd': np.zeros(num_trials),
        'unweighted_svdstack': np.zeros(num_trials),
        'weighted_svdstack': np.zeros(num_trials)
    }
    
    # Compute optimal weights
    w_stacksvd = compute_optimal_weights_stacksvd(theta_arr, c_arr)
    w_svdstack = compute_optimal_weights_svdstack(theta_arr, c_arr)
    
    for trial in tqdm(range(num_trials), desc="Running trials"):
        # Generate data
        X_list, v_true = generate_data_exponential(d, c_arr, theta_arr)
        
        # Run all four methods
        v_stacksvd = estimate_stacksvd(X_list)
        v_stacksvd_w = estimate_stacksvd(X_list, w_stacksvd)
        v_svdstack = estimate_svdstack(X_list)
        v_svdstack_w = estimate_svdstack(X_list, w_svdstack)
        
        # Store results (squared correlation)
        results['unweighted_stacksvd'][trial] = (v_stacksvd @ v_true)**2
        results['weighted_stacksvd'][trial] = (v_stacksvd_w @ v_true)**2
        results['unweighted_svdstack'][trial] = (v_svdstack @ v_true)**2
        results['weighted_svdstack'][trial] = (v_svdstack_w @ v_true)**2
    
    return results


def compute_asymptotic_power(method: int, theta_arr: List[float], c_arr: List[float]) -> float:
    """For one of the four methods, compute the asymptotic power.

    Args:
        method (int): {0: 'unweighted_stacksvd', 1: 'weighted_stacksvd', 2: 'unweighted_svdstack', 3: 'weighted_svdstack'}
        theta_arr (List[float]): Length M array of signal strengths
        c_arr (List[float]): Length M array of ratios n_i / d for each matrix

    Returns:
        float: Asymptotic power of the method
    """
    M = len(c_arr)
    c_arr = np.array(c_arr).astype(float)
    theta_arr = np.array(theta_arr).astype(float)
    theta_norm_sq = np.linalg.norm(theta_arr)**2
    
    if method == 0:  # unweighted stack-svd
        return max((theta_norm_sq**2 - c_arr.sum()) / (theta_norm_sq * (theta_norm_sq + 1)),0)
    
    if method == 1:  # weighted stack-svd
        theta_red = theta_arr[theta_arr > 0]
        c_red = c_arr[theta_arr > 0]
        def f(x):
            return sum(theta**2 * (1 - x) / (c*theta**-2 + x) for theta, c in zip(theta_red, c_red)) - 1
        
        # Binary search for the unique solution x in (0, 1)
        low, high = 0, 1
        while high - low > 1e-7:
            mid = (low + high) / 2
            if f(mid) > 0:
                low = mid
            else:
                high = mid
        x = (low + high) / 2
        return x
    
    ####### for svd-stack methods
    denom = theta_arr**4 + theta_arr**2
    inv_denom = np.ones(M)
    inv_denom[denom > 0] = 1 / denom[denom > 0]
    beta_vec = np.sqrt(np.maximum((theta_arr**4 - c_arr) * inv_denom, 0))

    if method == 2: # unweighted svd-stack    
        A = np.outer(beta_vec,beta_vec) + np.diag(1 - beta_vec**2)
        eigenvalues, eigenvectors = np.linalg.eigh(A)
        idx = np.argmax(np.abs(eigenvalues))
        
        lam = eigenvalues[idx]
        principal_eigenvector = eigenvectors[:, idx]
        return (principal_eigenvector @ beta_vec) ** 2 / lam
    
    if method == 3: # weighted svd-stack
        S = np.sum(beta_vec**2 / (1 - beta_vec**2))
        return S / (S+1)
    
    raise ValueError("Invalid method number")


def save_simulation_results_pickle(results: dict, c_arr: List[float], d: int, theta_arr: List[float], folder_path: str, filestr: str = ''):
    """
    Save the results of a run simulation to a pickle file.
    
    Args:
        results: Dictionary containing simulation results
        c_arr: List of ratios n_i / d for each matrix
        d: Column dimension
        theta_arr: List of signal strengths
        folder_path: Path to the folder where the results should be saved
        filestr: Optional string to include in the filename. If not provided, a default name will be used.
    """
    data = {
        'results': results,
        'c_arr': c_arr,
        'd': d,
        'theta_arr': theta_arr,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }
    
    # filename = f"{folder_path}/sim_res_{filestr}_{datetime.now().strftime('%H%M%S')}.pkl"
    filename = f"{folder_path}/sim_res_{filestr}.pkl"
    
    with open(filename, 'wb') as f:
        pickle.dump(data, f)
        
        
def generate_simulations(M_arr: np.ndarray, c_arr_arr: List[List[float]], theta_arr_arr: List[List[float]], d_arr: np.ndarray, num_trials: int, folder_path: str, filestr_prefix: str):
    for i, (M, d) in enumerate(zip(M_arr, d_arr)):
        print(f"Generating simulations for M={M}, d={d}")
        c_arr = c_arr_arr[i]
        theta_arr = theta_arr_arr[i]
        
        results = run_simulation(c_arr, d, theta_arr, num_trials)
        # Save results
        filestr = f'{filestr_prefix}_M_{M}_d_{d}'
        print(f"Saving results to {filestr}")
        save_simulation_results_pickle(results, c_arr, d, theta_arr, folder_path, filestr)
        
            
def generate_simulations_exponential(M_arr: np.ndarray, c_arr_arr: List[List[float]], theta_arr_arr: List[List[float]], d_arr: np.ndarray, num_trials: int, folder_path: str, filestr_prefix: str):
    for i, (M, d) in enumerate(zip(M_arr, d_arr)):
        print(f"Generating simulations exponential for M={M}, d={d}")
        c_arr = c_arr_arr[i]
        theta_arr = theta_arr_arr[i]
        
        results = run_simulation_exponential(c_arr, d, theta_arr, num_trials)
        # Save results
        filestr = f'{filestr_prefix}_M_{M}_d_{d}'
        print(f"Saving results to {filestr}")
        save_simulation_results_pickle(results, c_arr, d, theta_arr, folder_path, filestr)


def plot_results(x_arr: np.ndarray, mean_arr: np.ndarray, std_err_arr: np.ndarray, asymp_opt_power_arr: np.ndarray, folder_path: str, filestr_prefix: str, plot_var: str, title_str: str = ""):
    fig, ax = plt.subplots()
    for j in range(mean_arr.shape[1]):
        ax.errorbar(x_arr, mean_arr[:, j], yerr=std_err_arr[:, j], label=method_names[j])
        ax.plot(x_arr, asymp_opt_power_arr[:, j], linestyle=(0, dash_pattern[j]), color=ax.get_lines()[-1].get_color())
    
    if plot_var == 'd':
        ax.set_xlabel('d')
        plt.xscale('log')
    else:
        ax.set_xlabel('M')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    ax.set_ylabel('Squared inner product')
    ax.legend()
    ax.set_title(title_str)
    os.makedirs(folder_path, exist_ok=True)
    fig.savefig(f"{folder_path}/{filestr_prefix}_plot.pdf")
    
    

def load_simulation_results(filestr: str, folder_path: str, M_arr: np.ndarray, d_arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load simulation results from pickle files and compute mean, standard error, and asymptotic optimal power arrays.
    
    Args:
        filestr: File string prefix used in the filenames.
        folder_path: Path to the folder where the results are saved.
        M_arr: Array of M values.
        d_arr: Array of d values.
        
    Returns:
        mean_arr: Array of mean values for each method.
        std_err_arr: Array of standard error values for each method.
        asymp_opt_power_arr: Array of asymptotic optimal power values for each method.
    """
    mean_arr = np.zeros((len(M_arr), 4))
    std_err_arr = np.zeros((len(M_arr), 4))
    asymp_opt_power_arr = np.zeros((len(M_arr), 4))
    
    
    ######### construct a df that has columns for M, d, method, asymp_opt_power, trial number, and value
    df = pd.DataFrame(columns=['M', 'd', 'Method', 'asymp_opt_power', 'Trial', 'Value'])
    df_arr = []
    
    for i, (M, d) in enumerate(zip(M_arr, d_arr)):
        with open(f"{folder_path}/sim_res_{filestr}_M_{M}_d_{d}.pkl", 'rb') as f:
            data = pickle.load(f)
        
        results = data['results']
        num_trials = len(next(iter(results.values())))
        
        for j, (method, values) in enumerate(results.items()):
            mean = values.mean()
            std_err = values.std() / np.sqrt(num_trials)
            mean_arr[i, j] = mean
            std_err_arr[i, j] = std_err
            
            # Compute asymptotic optimal power
            asymp_opt_power = compute_asymptotic_power(j, data['theta_arr'], data['c_arr'])
            asymp_opt_power_arr[i, j] = asymp_opt_power
            
            ###### append to df
            for trial in range(num_trials):
                df_arr.append({
                    'M': M,
                    'd': d,
                    'Method': method_names[j],
                    'asymp_opt_power': asymp_opt_power,
                    'Trial': trial,
                    'Value': values[trial]
                })
    df = pd.DataFrame(df_arr)
    df.to_csv(f"{folder_path}/{filestr}_overall_results.tsv", sep='\t', index=False)
    
    
    return mean_arr, std_err_arr, asymp_opt_power_arr

def run_d_trial_exponential(rerun_trial=True, num_trials=10, folder_path='sim_results', d_max=10000):
    """Run simulations for M=3, increasing d, with c_i=1 and theta_i=[1.7, 1.6, 1.5] under exponential noise.

    Args:
        rerun_trial (bool): Whether to rerun the simulations.
        num_trials (int): Number of trials.
        folder_path (str): Folder to save results.
        d_max (int): Maximum value of d (default 10000).
    """
    M = 3
    d_arr = np.logspace(np.log10(100), np.log10(d_max), num=5, dtype=int)
    c_arr = [1.0] * M
    theta_arr = [1.7, 1.6, 1.5]
    c_arr_arr = [c_arr for _ in d_arr]
    theta_arr_arr = [theta_arr for _ in d_arr]
    M_arr = np.array([M for _ in d_arr])
    if rerun_trial:
        generate_simulations_exponential(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'd_trial_exponential')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('d_trial_exponential', folder_path, M_arr, d_arr)
    fig_title = "All c_i equal to 1, theta_i = [1.7, 1.6, 1.5], M=3, increasing d, exponential noise"
    plot_results(d_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'd_trial_exponential', 'd', fig_title)
    
    
def run_M_trial_exponential(rerun_trial=True, num_trials=10, folder_path='sim_results', d=2000):
    """Runs simulation trials for increasing M with fixed d and exponential noise, c_i = 1, theta_i = 1+2/(i+1)."""

    M_arr = np.array([2, 3, 5, 10])
    d_arr = np.array([d] * len(M_arr))
    
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr_arr = [np.array([1 + 2/(j+1) for j in range(M)]) for M in M_arr]
    
    if rerun_trial:
        generate_simulations_exponential(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'trial2_exp')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('trial2_exp', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 2: theta_i = 1 + 2/i \n c all equal to 1, d = {d}, exponential noise"
    plot_results(M_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'trial2_exp', 'M', title_str)
    
    
    
def run_fig2a_rate(rerun_trial=True, num_trials=10, folder_path='sim_results', dmax=10000, dmin=100, num_points=5):
    """
    Runs simulation for Figure 2a, with theta values set to [0.97, 0.8] and c values set to [1, 1], varying the dimension parameter `d` while keeping other parameters fixed.
    """
    M = 2
    d_arr = np.logspace(np.log10(dmin), np.log10(dmax), num=num_points, dtype=int)
    M_arr = np.array([M for _ in d_arr])
    
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr = np.array([.97, 0.8])
    theta_arr_arr = [theta_arr for _ in d_arr]
    
    print(f"Running trial 2a: thetas = [{theta_arr[0]},{theta_arr[1]}], c = [1,1], M=2, increasing d")

    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'fig2a_rate')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('fig2a_rate', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 2a: theta_i = [{theta_arr[0]}, {theta_arr[1]}] \n c_i = [1, 1], M=2, increasing d"
    plot_results(d_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'fig2a_rate', 'd', title_str)


def run_fig2b_rate(rerun_trial=True, num_trials=10, folder_path='sim_results', dmax=10000, num_points=5):
    """
    Runs simulation for Figure 2b with c = [1, 1] and theta = [1.2, 1.05], increasing d.
    """
    print("Running trial 2b: thetas = [1.2, 1.05], c = [1,1], M=2, increasing d")
    M = 2
    d_arr = np.logspace(np.log10(100), np.log10(dmax), num=num_points, dtype=int)
    M_arr = np.array([M for _ in d_arr])
    
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr = np.array([1.2, 1.05])
    theta_arr_arr = [theta_arr for _ in d_arr]
    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'fig2b_rate')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('fig2b_rate', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 2b: theta_i = [1.2, 1.05] \n c_i = [1, 1], M=2, increasing d"
    plot_results(d_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'fig2b_rate', 'd', title_str)
    
    
def run_fig2c_rate(rerun_trial=True, num_trials=10, folder_path='sim_results', dmax=10000, num_points=5):
    """Runs and plots simulation results for Figure 2c with theta = [2, 1.3], c = [1, 1], and varying dimension d."""
    
    print("Running trial 2c: thetas = [2, 1.3], c = [1,1], M=2, increasing d")
    M = 2
    d_arr = np.logspace(np.log10(100), np.log10(dmax), num=num_points, dtype=int)
    M_arr = np.array([M for _ in d_arr])
    
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr = np.array([2, 1.3])
    theta_arr_arr = [theta_arr for _ in d_arr]
    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'fig2c_rate')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('fig2c_rate', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 2c: theta_i = [2, 1.3] \n c_i = [1, 1], M=2, increasing d"
    plot_results(d_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'fig2c_rate', 'd', title_str)


######## fig 3a: \theta_i =.9, c_i = 1, increasing M from 2 to M_max
def run_fig3a_selectedEx(rerun_trial=True, num_trials=10, folder_path='sim_results', d=5000, M_max=10):
    """Runs simulation for Figure 3a with theta_i = 0.7, c_i = 1, and increasing M from 2 to M_max."""
    
    print("Running trial 3a: theta_i=.7, c_i=1, increasing M")
    M_arr = np.arange(2, M_max+1, 2)
    d_arr = np.array([d] * len(M_arr))
    
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr_arr = [np.array([0.7] * M) for M in M_arr]
    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'fig3a_selectedEx')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('fig3a_selectedEx', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 3a: theta_i = [0.9] \n c_i = [1], increasing M"
    plot_results(M_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'fig3a_selectedEx', 'M', title_str)

def run_fig3b_selectedEx(rerun_trial=True, num_trials=10, folder_path='sim_results', d=5000, M_max=10):
    """Runs simulation for Figure 3b with theta_i = [1.2, 1.2, 0, ...], c_i = 1, and increasing M from 2 to M_max."""
    print("Running trial 3b: thetas = [1.2, 1.2,0...], c_i=1, increasing M")
    M_arr = np.arange(2, M_max+1, 2)
    d_arr = np.array([d] * len(M_arr))
    
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr_arr = [np.array([1.2] * 2 + [0] * (M - 2)) for M in M_arr]
    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'fig3b_selectedEx')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('fig3b_selectedEx', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 3b: theta_i = [1.2] \n c_i = [1], increasing M"
    plot_results(M_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'fig3b_selectedEx', 'M', title_str)

def run_fig3c_selectedEx(rerun_trial=True, num_trials=10, folder_path='sim_results', dmax=10000, num_points=5):
    """Runs simulation for Figure 3c with theta_i = [.95, .95, 0], c_i = [1, 1, 2], and increasing d."""
    print("Running trial 3c: thetas = [.95,.95,0], c = [1,1,2], increasing d")
    M = 3
    d_arr = np.logspace(np.log10(100), np.log10(dmax), num=num_points, dtype=int)
    M_arr = np.array([M for _ in d_arr])
    
    c_arr = np.array([1, 1, 2]).astype(float)
    theta_arr = np.array([.95, .95, 0]).astype(float)
    
    c_arr_arr = [c_arr for _ in d_arr]
    theta_arr_arr = [theta_arr for _ in d_arr]
    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'fig3c_selectedEx')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('fig3c_selectedEx', folder_path, M_arr, d_arr)
    
    title_str = f"Trial 3c: theta_i = [0.95, 0.95, 0] \n c_i = [1, 1, 2], increasing d"
    plot_results(d_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'fig3c_selectedEx', 'd', title_str)



def run_svdstack_edgecase(rerun_trial=True, num_trials=10, folder_path='sim_results', d=5000,M_max=10,num_points=5):
    """Runs simulation for an edge case of SVD-Stack where beta_1>0 and beta_2=0."""
    print("Running SVD-Stack edge case: d=5000, M=3")
    M_arr = np.linspace(2, M_max, num_points, dtype=int)
    d_arr = np.array([d] * len(M_arr))
        
    c_arr_arr = [np.ones(M) for M in M_arr]
    theta_arr_arr = [np.array([2.0] + [0.0] * (M - 1)) for M in M_arr]
    
    if rerun_trial:
        generate_simulations(M_arr, c_arr_arr, theta_arr_arr, d_arr, num_trials, folder_path, 'svdstack_edgecase')
    mean_arr, std_err_arr, asymp_opt_power = load_simulation_results('svdstack_edgecase', folder_path, M_arr, d_arr)

    title_str = f"Trial svdstack_edgecase: theta_i = [2.0, 0.0, ...] \n c_i = [1, 1, ...], increasing M"
    plot_results(M_arr, mean_arr, std_err_arr, asymp_opt_power, folder_path, 'svdstack_edgecase', 'M', title_str)




method_names = ['Unweighted Stack-SVD', 'Weighted Stack-SVD', 'Unweighted SVD-Stack', 'Weighted SVD-Stack']
dash_pattern = [[4,4], [5,4],[6,4],[7,4]] 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script for figure generation for \"Stacked SVD or SVD stacked? Optimal data integration via Random Matrix Theory\" by Baharav*, Nicol* et. al. "
                    "Call with --test to use small parameter choices for efficient testing, or larger custom values with command line specification. "
                    "Script by default uses GPUs if available, otherwise falls back to CPU."
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help='Use test settings: num_trials=10, dmax=1000, num_points=5, M_max=10 (overrides other parameters; default: False)'
    )
    parser.add_argument('--num_trials', type=int, default=100, help='Number of trials (default: 100)')
    parser.add_argument('--dmax', type=int, default=20000, help='Maximum d value (default: 20000)')
    parser.add_argument('--num_points', type=int, default=10, help='Number of points (default: 10)')
    parser.add_argument('--M_max', type=int, default=16, help='Maximum M value (default: 16)')
    parser.add_argument('--rerun_trials', type=lambda x: (str(x).lower() == 'true'), default=True,
                        help='Whether to rerun the simulations (default: True)')
    parser.add_argument('--folder_path', type=str, default='sim_results_test',
                        help='Folder to save results (default: sim_results_test)')

    
    args = parser.parse_args()

    if args.test:
        num_trials = 5
        dmax = 1000
        num_points = 3
        M_max = 6
        rerun_trials = True
        folder_path = 'sim_results_test'
    else:
        num_trials = args.num_trials
        dmax = args.dmax
        num_points = args.num_points
        M_max = args.M_max
        rerun_trials = args.rerun_trials
        folder_path = args.folder_path
        
    print(f"Running simulations with num_trials={num_trials}, dmax={dmax}, num_points={num_points}, M_max={M_max}, rerun_trials={rerun_trials}, folder_path='{folder_path}'")
    
    os.makedirs(folder_path, exist_ok=True)


    ###### run fig2a, fig2b, fig2c
    run_fig2a_rate(rerun_trials, num_trials, folder_path, dmax=dmax, num_points=num_points)
    run_fig2b_rate(rerun_trials, num_trials, folder_path, dmax=dmax, num_points=num_points)
    run_fig2c_rate(rerun_trials, num_trials, folder_path, dmax=dmax, num_points=num_points)
    
    ######## run fig 3 sims
    run_fig3a_selectedEx(rerun_trials, num_trials, folder_path, d=dmax, M_max=M_max)
    run_fig3b_selectedEx(rerun_trials, num_trials, folder_path, d=dmax, M_max=M_max)
    run_fig3c_selectedEx(rerun_trials, num_trials, folder_path, dmax=dmax, num_points=num_points)
    
    #### run exponential noise trials
    run_d_trial_exponential(rerun_trials, num_trials, folder_path, d_max=dmax)
    run_M_trial_exponential(rerun_trials, num_trials, folder_path, d=dmax)
    
    ##### run svdstack edge case
    run_svdstack_edgecase(True, num_trials, folder_path, d=dmax,M_max=M_max,num_points=num_points)


