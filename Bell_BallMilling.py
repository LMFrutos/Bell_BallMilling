#!/usr/bin/env python3
"""
Ball-Milling-Bell-model
-----------------------
-----------------------

This program evaluates the work term W(φ,θ) applied to a molecular geometry 
under external pressure.

General
----------------------------------------------
1.  The two geometries and the input parameters are *not* hard-coded: they
    are loaded from external text files:
        • R_min.txt
        • R_TS.txt
        • input.txt
   (formats are described below).
2.  Output files with the results:
        • output.txt  – main numbers and diagnostics
        • W.mtx       – W(φ,θ) with coordinates (φ,θ,W(φ,θ))

File formats
------------
R_min.txt and R_TS.txt
    Line 1         : integer N (number of atoms)
    Lines 2 … N+1  : Z  x  y  z   (atomic number and Cartesian coordinates, Å)

input.txt
    Line 1         : float     external pressure, P_ext, in GPa
    Line 2         : integer   ntheta (grid resolution for θ coordinate)
    Line 3         : string    'YES' or 'NO' to include van der Waals radii

The script performs the following major steps
---------------------------------------------
1.  Load input files and verify consistency.
2.  Remove overall translation by shifting both geometries to their
    respective centres of mass (mass-weighted).
3.  Remove overall rotation with a *mass-weighted* Kabsch alignment
    so that R_TS is optimally superimposed onto R_min.
4.  Compute delta_R = R_TS_aligned – R_min_c (internal displacements).
5.  Determine the radius of the smallest sphere that contains either
    geometry; use it to compute the external force for the requested pressure.
6.  Loop over a (φ,θ) grid (nφ × nθ) and evaluate W(φ,θ).
7.  Write summary numbers to *output.txt* and the full W function to *W.mtx*.

"""
from pathlib import Path
import numpy as np
from scipy.io import mmwrite
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# --------------------------------------------------------------------------
# ------------------------- 1.  Utility functions --------------------------
# --------------------------------------------------------------------------

def get_vdw_radius(z: int) -> float:
    """
    Get van der Waals radius for atomic number Z.
    Returns radius in Ångstroms.
    
    Parameters
    ----------
    z : int
        Atomic number
        
    Returns
    -------
    float
        van der Waals radius in Å
    """
    # Van der Waals radii database (Å)
    vdw_radii = {
        1: 1.20,   # H
        6: 1.70,   # C
        7: 1.55,   # N
        8: 1.52,   # O
        9: 1.47,   # F
        15: 1.80,  # P
        16: 1.80,  # S
        17: 1.75,  # Cl
        # Add more atomic numbers and their vdW radii here as needed        
    }
    
    if z not in vdw_radii:
        raise ValueError(f"van der Waals radius not defined for atomic number {z}")
    
    return vdw_radii[z]


def get_atomic_mass(z: int) -> float:
    """
    Get atomic mass for atomic number Z.
    Returns mass in atomic units (a.u.).
    
    Parameters
    ----------
    z : int
        Atomic number
        
    Returns
    -------
    float
        atomic mass in a.u.
    """
    # Atomic masses database (a.u.)
    atomic_masses = {
        1: 1.008,    # H
        6: 12.011,   # C
        7: 14.007,   # N
        8: 15.999,   # O
        9: 18.998,   # F
        15: 30.974,  # P
        16: 32.065,  # S
        17: 35.453,  # Cl
        # Add more atomic numbers and their masses here as needed        
    }
    
    if z not in atomic_masses:
        raise ValueError(f"Atomic mass not defined for atomic number {z}")
    
    return atomic_masses[z]


def read_geometry(filename: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Read a geometry file with format:
        N
        Z1 x1 y1 z1
        ...
        ZN xN yN zN
    Returns
    -------
    coords : (N, 3) ndarray of floats (Å)
    zatom : (N,) ndarray of ints (atomic numbers)
    """
    with open(filename, "r", encoding="utf-8") as f:
        try:
            N = int(f.readline().strip())
        except ValueError as exc:
            raise ValueError(f"First line of {filename} must be an integer.") from exc
        data = np.loadtxt(f)
    if data.shape != (N, 4):
        raise ValueError(
            f"{filename}: expected {N} lines with 4 columns (Z x y z), found shape {data.shape}"
        )
    
    zatom = data[:, 0].astype(int)
    coords = data[:, 1:4]
    return coords, zatom


def read_input_parameters(filename: str) -> tuple[float, int, str]: # <--- CAMBIO
    """
    Read the input file with format:
        P_ext_GPa
        ntheta
        YES/NO (for vdW radii)
    Returns
    -------
    P_ext_GPa : float
    ntheta : int
    add_vdw_str : str
    """
    with open(filename, "r", encoding="utf-8") as f:
        try:
            P_ext_GPa = float(f.readline().strip())
        except ValueError as exc:
            raise ValueError(f"First line of {filename} must be a float.") from exc
        try:
            ntheta = int(f.readline().strip())
        except ValueError as exc:
            raise ValueError(f"Second line of {filename} must be an integer.") from exc
        # <--- CAMBIO: Bloque añadido para leer la tercera línea
        try:
            add_vdw_str = f.readline().strip().upper()
            if add_vdw_str not in ["YES", "NO"]:
                raise ValueError("Third line must be 'YES' or 'NO'.")
        except IndexError as exc:
            raise ValueError("Third line ('YES' or 'NO') is missing in input file.") from exc

    return P_ext_GPa, ntheta, add_vdw_str


def kabsch_weighted(A: np.ndarray, B: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    Mass-weighted Kabsch algorithm.
    Given two coordinate sets A and B (N×3) already translated to their
    respective centres of mass, and a weight vector w (masses, length N),
    return the optimal rotation matrix R such that
        B_rot = B @ R
    minimises RMSD between A and B_rot.
    """
    # Correlation matrix (3 × 3)
    C = A.T @ np.diag(w) @ B
    # Singular-value decomposition
    U, _, Vt = np.linalg.svd(C)
    V = Vt.T
    # Correct for improper rotation if necessary
    d = np.sign(np.linalg.det(V @ U.T))
    D = np.diag([1.0, 1.0, d])
    # Optimal rotation
    R = V @ D @ U.T
    return R


def generate_spherical_plot(W: np.ndarray, nphi: int, ntheta: int, min_W: float, max_W: float) -> None:
    """
    Generate a spherical representation of W(φ,θ) using a rainbow colormap.
    
    Parameters
    ----------
    W : np.ndarray
        W(φ,θ) matrix (nphi × ntheta)
    nphi : int
        Number of phi points
    ntheta : int
        Number of theta points
    min_W : float
        Minimum W value
    max_W : float
        Maximum W value
    """
    # Create coordinate grids
    phi = np.linspace(0, 2*np.pi, nphi)
    theta = np.linspace(0, np.pi, ntheta)
    PHI, THETA = np.meshgrid(phi, theta, indexing='ij')
    
    # Convert spherical to Cartesian coordinates for 3D plotting
    X = np.sin(THETA) * np.cos(PHI)
    Y = np.sin(THETA) * np.sin(PHI)
    Z = np.cos(THETA)
    
    # Create rainbow colormap (red to blue)
    colors = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue']
    n_bins = 256
    rainbow_cmap = LinearSegmentedColormap.from_list('rainbow', colors, N=n_bins)
    
    # Create the plot
    fig = plt.figure(figsize=(12, 10))
    
    # 3D spherical plot
    ax1 = fig.add_subplot(221, projection='3d')
    surf = ax1.plot_surface(X, Y, Z, facecolors=rainbow_cmap((W - min_W) / (max_W - min_W)),
                           alpha=0.9, linewidth=0, antialiased=True, rcount=ntheta, ccount=nphi)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title('W(φ,θ) on Sphere - 3D View')
    ax1.set_box_aspect([1,1,1])

#    # 2D rectangular plot
    ax3 = fig.add_subplot(223)
    im3 = ax3.imshow(W.T, extent=[0, 2*np.pi, 0, np.pi], aspect='auto', 
                     cmap=rainbow_cmap, origin='lower')
    ax3.set_xlabel('φ (radians)')
    ax3.set_ylabel('θ (radians)')
    ax3.set_title('W(φ,θ) - Rectangular Plot')
    
    # Colorbar
    ax4 = fig.add_subplot(224)
    cbar = plt.colorbar(im3, ax=ax4, orientation='vertical')
    cbar.set_label('W (kcal/mol)')
    ax4.axis('off')
    ax4.text(0.1, 0.9, f'Min W: {min_W:.3f} kcal/mol', transform=ax4.transAxes, fontsize=10)
    ax4.text(0.1, 0.8, f'Max W: {max_W:.3f} kcal/mol', transform=ax4.transAxes, fontsize=10)
    ax4.text(0.1, 0.7, f'Range: {max_W-min_W:.3f} kcal/mol', transform=ax4.transAxes, fontsize=10)
    ax4.text(0.1, 0.5, 'Red: Lower energy', transform=ax4.transAxes, fontsize=10, color='red')
    ax4.text(0.1, 0.4, 'Blue: Higher energy', transform=ax4.transAxes, fontsize=10, color='blue')
    
    plt.tight_layout()
    plt.savefig('W_spherical.png', dpi=300, bbox_inches='tight')
    plt.close()


# --------------------------------------------------------------------------
# ------------------------- 2.  Main computation ---------------------------
# --------------------------------------------------------------------------

def main() -> None:
    # ---- 2.1  Read inputs -------------------------------------------------
    geom_min, zatom_min = read_geometry("R_min.txt")
    geom_TS, zatom_TS = read_geometry("R_TS.txt")
    P_ext_GPa, ntheta, add_vdw_str = read_input_parameters("input.txt") 

    if geom_min.shape != geom_TS.shape:
        raise RuntimeError("R_min and R_TS contain different numbers of atoms.")
    if not np.array_equal(zatom_min, zatom_TS):
        raise RuntimeError("R_min and R_TS have different atomic numbers.")

    N = geom_min.shape[0]                   # number of atoms
    zatom = zatom_min                       # atomic numbers array
    
    # Generate masses from atomic numbers
    masses = np.array([get_atomic_mass(z) for z in zatom])
    masses = masses.astype(float)
    mass_total = masses.sum()

    # Get van der Waals radii for all atoms
    vdw_radii = np.array([get_vdw_radius(z) for z in zatom])
    
    # Logic to use or not the van der Waals radii
    use_vdw = (add_vdw_str == "YES")
    if not use_vdw:
        vdw_radii.fill(0.0) # Set all vdW radii to zero if input is 'NO'

    # ---- 2.2  Remove overall translation (shift to centre of mass) -------
    COM_min = (masses @ geom_min) / mass_total     # 1 × 3
    COM_TS  = (masses @ geom_TS)  / mass_total
    R_min_c = geom_min - COM_min                   # centred geometries
    R_TS_c  = geom_TS  - COM_TS

    # ---- 2.3  Remove overall rotation (weighted Kabsch) ------------------
    R_opt = kabsch_weighted(R_min_c, R_TS_c, masses)
    R_TS_aligned = R_TS_c @ R_opt                # now aligned with R_min_c

    # ---- 2.4  Internal displacement --------------------------------------
    delta_R = R_TS_aligned - R_min_c             # N × 3

    # ---- 3.   Sphere radius enclosing both geometries --------------------
    # Find farthest atom distance from origin and add its vdW radius
    dists_Rmin = np.linalg.norm(geom_min, axis=1)
    dists_RTS = np.linalg.norm(geom_TS, axis=1)
    
    idx_max_Rmin = np.argmax(dists_Rmin)
    idx_max_RTS = np.argmax(dists_RTS)
    
    # The vdw_radii array will be zero if 'NO' was selected
    max_dist_Rmin = dists_Rmin[idx_max_Rmin] + vdw_radii[idx_max_Rmin]
    max_dist_RTS = dists_RTS[idx_max_RTS] + vdw_radii[idx_max_RTS]
    
    R_sphere = max(max_dist_Rmin, max_dist_RTS)  # Å

    # ---- 4.   External force corresponding to the requested pressure -----
    # Force (nN) according to the original formula:
    #   F_ext = P_ext_GPa * 0.01 * p * R_sphere²
    F_ext = P_ext_GPa * 0.01 * np.pi * R_sphere**2    # nN
    Fext_i = F_ext / N                                # per-atom force

    # ---- 5.   Loop over (φ,θ) grid and evaluate W ------------------------
    nphi   = 2 * ntheta  # Keep 2:1 ratio
    dphi   = 2.0 * np.pi / nphi
    dtheta = np.pi / ntheta
    factor_conv = 1.88973 * 627.5 / 82.387            # conversion to a.u.

    W      = np.zeros((nphi, ntheta))
    R_mat  = np.zeros_like(W)
    sum_Ea = 0.0

    # Pre-compute norms and unit vectors of R_min 
    norms_Rmin = np.linalg.norm(geom_min, axis=1)
    # Avoid division by zero for any possible origin-coincident atom
    with np.errstate(divide="ignore", invalid="ignore"):
        dR_min_unit = np.where(
            norms_Rmin[:, None] > 0,
            geom_min / norms_Rmin[:, None],
            0.0
        )

    for iphi in range(nphi):
        phi = iphi * dphi
        cos_phi, sin_phi = np.cos(phi), np.sin(phi)

        for itheta in range(ntheta):
            theta = itheta * dtheta
            sin_theta, cos_theta = np.sin(theta), np.cos(theta)

            # Outward normal vector n(φ,θ)
            n_vec = np.array([sin_theta * cos_phi,
                              sin_theta * sin_phi,
                              cos_theta])

            # ---- 5.1  Largest distance of projection of R_min_c onto
            #           plane perpendicular to n_vec  -----------------------
            proj_along_n = R_min_c @ n_vec
            Rproj        = R_min_c - proj_along_n[:, None] * n_vec
            r_i          = np.linalg.norm(Rproj, axis=1)
            
            # Find the index of the atom with the largest projection
            idx_max_proj = np.argmax(r_i)
            # The base distance is the largest of the projections
            r_max_base = r_i[idx_max_proj]
            # The effective radius is the base distance plus the van der Waals radius of the atom
            # The vdw_radii array already contains zeros if the user selected 'NO'.
            r_max = r_max_base + vdw_radii[idx_max_proj]


            R_mat[iphi, itheta] = (r_max / R_sphere)**2

            # ---- 5.2  Work term W(φ,θ) -----------------------------------
            dot_dR_n = dR_min_unit @ n_vec
            Feff = (-dot_dR_n[:, None] * n_vec) * Fext_i * (r_max / R_sphere)**2
            dotFeff_deltaR = np.einsum("ij,ij->i", Feff, delta_R)

            suma = -factor_conv * dotFeff_deltaR.sum()
            W[iphi, itheta] = suma

            sum_Ea += (-factor_conv
                       * (dphi * dtheta * sin_theta)
                       * dotFeff_deltaR.sum())

    sum_Ea /= (4.0 * np.pi)
    min_W, max_W = W.min(), W.max()

    # ---- 6.   Write outputs ---------------------------------------------
    with open("output.txt", "w", encoding="utf-8") as f:
        f.write("# =====  Results generated by compute_W.py  =====\n")
        f.write(f"Number of atoms, N                 : {N:d}\n")
        f.write(f"External pressure, P_ext (GPa)     : {P_ext_GPa:.6f}\n")
        f.write(f"Include vdW radii?               : {add_vdw_str}\n")  
        f.write(f"Sphere radius, R_sphere (Å)        : {R_sphere:.6f}\n")
        f.write(f"Total external force, F_ext (nN)   : {F_ext:.6f}\n")
        f.write(f"Per-atom external force, F_i (nN)  : {Fext_i:.6f}\n")
        f.write(f"Average ?Ea over sphere (kcal/mol) : {sum_Ea:.6f}\n")
        f.write(f"Minimum W value (kcal/mol)         : {min_W:.6f}\n")
        f.write(f"Maximum W value (kcal/mol)         : {max_W:.6f}\n")
        f.write("# W(phi,theta) is stored separately in W.mtx (phi,theta,W(phi,theta))\n")

    # Write W.mtx 
    with open("W.mtx", "w", encoding="utf-8") as f:
        f.write(f"% W(phi,theta) generated by compute_W.py. Format: (phi,theta,W(phi,theta))\n")
        f.write(f"{nphi} {ntheta} {nphi * ntheta}\n")
        
        # Write data: phi theta W(phi,theta)
        for iphi in range(nphi):
            phi = iphi * dphi
            for itheta in range(ntheta):
                theta = itheta * dtheta
                f.write(f"{phi:.6f} {theta:.6f} {W[iphi, itheta]:.6f}\n")

    # ---- 7.   Generate spherical visualization of W(φ,θ) ----------------
    generate_spherical_plot(W, nphi, ntheta, min_W, max_W)

    print("Computation Ball-Milling-Bell-model finished.")
    print("  · Numerical summary : output.txt")
    print("  · Coordinate file   : W.mtx")
    print("  · Spherical plot    : W_spherical.png")


# --------------------------------------------------------------------------
# --------------------------- 3.  Script entry -----------------------------
# --------------------------------------------------------------------------

if __name__ == "__main__":
    main()