import numpy as np
import scipy.io
import matplotlib.pyplot as plt


def compute_3d_energy_spectrum(u, v, w, L):
    """
    Computes the 1D energy spectrum E(k) from a 3D velocity field.

    Parameters:
    u, v, w : 3D numpy arrays (N, N, N) representing velocity components
    L       : Physical domain length (float)

    Returns:
    k_shells : Wavenumber magnitudes
    E_k      : Energy spectrum values
    """
    N = u.shape[0]
    dx = L / N

    # 1. Perform 3D Real FFT and normalize by the number of grid points
    # We use rfftn for efficiency if the input is real
    u_hat = np.fft.fftn(u)/N**3
    v_hat = np.fft.fftn(v)/N**3
    w_hat = np.fft.fftn(w)/N**3

    # 2. Compute Energy Density in Fourier space: 1/2 * (|u^|^2 + |v^|^2 + |w^|^2)
    # Note: For HIT, E(k) is the trace of the spectral tensor
    energy_density = 0.5 * (np.abs(u_hat)**2 + np.abs(v_hat)**2 + np.abs(w_hat)**2)

    # 3. Create wavenumber coordinates
    # k_n = 2*pi*n / L
    k = np.fft.fftfreq(N, d=dx) * 2 * np.pi

    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # 4. Binning: Group energy by wavenumber magnitude k
    k_min = 2 * np.pi / L
    k_max = np.pi * N / L  # Nyquist limit

    # Define bins (shells)
    num_bins = N // 2
    bins = np.linspace(k_min, k_max, num_bins)
    k_shells = 0.5 * (bins[1:] + bins[:-1])

    # Sum energy density within each shell
    E_k, _ = np.histogram(K_mag.flatten(), bins=bins, weights=energy_density.flatten())
    dk = bins[1] - bins[0]
    E_k = E_k / dk
    # Adjust for spherical shell volume density (if normalization is required for 1D spectrum)
    # The histogram weights sum up the discrete modes.

    return k_shells, E_k


def read_data(N, sys_size, filename):
    # Read binary data as double-precision floats
    # np.fromfile is the equivalent of fread
    with open(filename, 'rb') as fileID:
        A = np.fromfile(fileID, dtype=np.float64)

    # Reassign density & velocity components
    # MATLAB's reshape/permute logic maps to NumPy's reshape with 'F' order
    # MATLAB: reshape(A, mp, np, pp, sys_size) -> (mp, np, pp, sys_size)
    # MATLAB: permute(..., [4 1 2 3]) -> (sys_size, mp, np, pp)

    qc = A.reshape((N+1, N+1, N+1, sys_size), order='F')
    qc = np.transpose(qc, (3, 0, 1, 2))

    rho = qc[0, :, :, :]
    u = qc[1, :, :, :] / rho
    v = qc[2, :, :, :] / rho
    w = qc[3, :, :, :] / rho

    return u, v, w


# Parameters
L = 0.09*2*np.pi
sys_size = 1 + 3 + 1 + 1

# Post-processing
# n31
N = 31
klim31 = np.pi/(L/N)
u, v, w = read_data(N, sys_size, "n31/restart_data/lustre_0.dat")
k31, Ek31 = compute_3d_energy_spectrum(u, v, w, L)
# n63
N = 63
klim63 = np.pi/(L/N)
u, v, w = read_data(N, sys_size, "n63/restart_data/lustre_0.dat")
k63, Ek63 = compute_3d_energy_spectrum(u, v, w, L)
# n127
N = 127
klim127 = np.pi/(L/N)
u, v, w = read_data(N, sys_size, "n127/restart_data/lustre_0.dat")
k127, Ek127 = compute_3d_energy_spectrum(u, v, w, L)
# n255
N = 255
klim255 = np.pi/(L/N)
u, v, w = read_data(N, sys_size, "n255/restart_data/lustre_0.dat")
k255, Ek255 = compute_3d_energy_spectrum(u, v, w, L)

# von Karman-Pao spectrum from Saad et al. (2017, AIAA)
rho0 = 1000.0  # [kg/m3]
mu = 1.002e-03  # [kg/m/s]
hit_alpha = 1.453  # [1]
hit_kappa_e = 40*(5/12)**0.5  # [1/m]
hit_uprime = 0.25  # [m/s]
hit_lscale = 0.746834/hit_kappa_e  # [m]
hit_epsilon = hit_uprime**3/hit_lscale  # [m2/s3]
hit_nu = mu / rho0  # [m2/s]
hit_kappa_eta = hit_epsilon**(1/4)*hit_nu**(-3/4)  # [1/m]
hit_domain = 0.09*2*np.pi  # [m]
k = np.linspace(10, 2e3, 10000)
E_true = hit_alpha*(hit_uprime**2/hit_kappa_e)*(k/hit_kappa_e)**4/(1+(k/hit_kappa_e)**2)**(17/6)*np.exp(-2*(k/hit_kappa_eta)**2)

# Plot
plt.figure(figsize=(5, 5))
plt.rcParams['text.usetex'] = True
plt.loglog(k, E_true, 'r-', linewidth=2, label=r'$\mbox{Analytic}$')
plt.loglog(k31, Ek31, 'ko-', linewidth=1, label=r'$N=32^3$')
plt.loglog([klim31, klim31], [1e-6, 1e-3], 'k--', linewidth=1)
plt.loglog(k63, Ek63, 'b^-', linewidth=1, label=r'$N=64^3$')
plt.loglog([klim63, klim63], [1e-6, 1e-3], 'b--', linewidth=1)
plt.loglog(k127, Ek127, 'g>-', linewidth=1, label=r'$N=128^3$')
plt.loglog([klim127, klim127], [1e-6, 1e-3], 'g--', linewidth=1)
plt.loglog(k255, Ek255, 'm<-', linewidth=1, label=r'$N=256^3$')
plt.loglog([klim255, klim255], [1e-6, 1e-3], 'm--', linewidth=1)
plt.legend()
plt.xlabel(r'$\kappa (m^{-1})$')
plt.ylabel(r'$E(\kappa) (m^3 / s^2)$')
plt.xlim(left=10, right=2e3)
plt.ylim(top=1e-3, bottom=1e-6)
plt.grid(color='k', linestyle='-', linewidth=0.3, which='major')
plt.grid(color='k', linestyle=':', linewidth=0.2, which='minor')
plt.savefig('energy_spectrum.png')
plt.close()
