<script>
MathJax = {
  tex: {
    inlineMath: [['$', '$'], ['\\(', '\\)']]
  }
};
</script>
<script id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js">
</script>

# Specification for Controllable Output Density Synthetic Axis-Aligned Hyper-rectangle Generator

## 1. Background and Objectives

### 1.1 Problem Definition

In the performance evaluation of multidimensional spatial database systems, the spatial join is one of the most computationally challenging operators. Let $R$ and $S$ be two sets containing $d$-dimensional axis-aligned hyper-rectangles (hereinafter referred to as "boxes"). Each object $b$ is defined as the Cartesian product of $d$ half-open intervals:

$$
b = \prod_{k=1}^{d}[\ell_k,\ u_k) \quad \subseteq \mathbb{R}^d
$$

The result set $J$ of a spatial join under the "intersects" predicate is defined as:

$$
J(R, S) = \{(r, s) \in R \times S \mid r \cap s \neq \varnothing\}
$$

To ensure the credibility and fairness of benchmark results, the synthetic data generation mechanism must possess **statistical controllability**. That is, beyond the appearance of the generated geometric objects, the computational load of the join operation must be strictly controlled.

### 1.2 Output Density Metric $\alpha_{\mathrm{out}}$

The cardinality $\vert J \vert$ of the join result is typically non-linear and highly dependent on input scale and data distribution. To establish a unified load metric across experiments of different scales, the normalized output density $\alpha_{\mathrm{out}}$ is defined as follows:

$$
\alpha_{\mathrm{out}} = \frac{|J(R, S)|}{|R| + |S|}
$$

This metric quantifies the average number of result edges contributed by each node (input object) in the bipartite join graph.

### 1.3 Generation Goal

The core objective of this generator is to solve the **inverse problem**: Given input cardinalities $(n_R, n_S)$, dimension $d$, and a target output density $\alpha_{\mathrm{out}}^{\star}$, automatically derive and generate datasets $R, S$ that satisfy the following statistical expectation:

$$
\mathbb{E}\left[ \frac{|J(R, S)|}{n_R + n_S} \right] \approx \alpha_{\mathrm{out}}^{\star}
$$

This process must simultaneously output complete metadata, including random seeds and parameter trajectories used in the generation process, to ensure strict reproducibility of the experiments.

## 2. Generation Principles and Methodology

### 2.1 Space Definition and Intersection Test

Define the spatial universe $\mathcal{U}$ as a $d$-dimensional hyper-rectangle:

$$
\mathcal{U} = \prod_{k=1}^{d}[u_k^{\min},\ u_k^{\max})
$$

Let $W_k = u_k^{\max} - u_k^{\min}$ be the span of the $k$-th dimension. The total volume of the universe is $V_{\mathcal{U}} = \prod_{k=1}^d W_k$.

For any two boxes $r, s$, the necessary and sufficient condition for their intersection is that there is an overlap in the intervals across all dimensions. Using the half-open interval representation, this condition can be formalized as:

$$
\forall k \in \{1, \dots, d\}: \max(\ell_k^{(r)}, \ell_k^{(s)}) < \min(u_k^{(r)}, u_k^{(s)})
$$

### 2.2 Coverage-based Scale Control

To control the average scale of geometric objects, we introduce a **coverage** parameter $C \in (0, \infty)$. For any set of objects $T \in \{R, S\}$, its coverage is defined as the ratio of the sum of the volumes of all objects $\sum_{i=1}^{N} \mathrm{vol}(r_i)$ to the volume of the universe $V_{\mathcal{U}}$.

Based on a given $C$ and set cardinality $n_T$, the **expected volume** $\bar{v}_T$ of objects in that set is derived as:

$$
\bar{v}_T = \frac{C \cdot V_{\mathcal{U}}}{n_T}
$$

This formula ensures **scale invariance** of the generation model: as the sample size $n_T$ increases, the average volume of individual objects automatically scales by a factor of $1/n_T$, thereby maintaining relatively stable spatial filling characteristics for a fixed $C$.

### 2.3 Stochastic Geometric Generation Model

For each object in set $T$, its geometric attributes are generated through three orthogonal independent sampling steps:

#### 2.3.1 Volume Sampling

First, sample the scalar volume $V$ of the object. The system supports various probability distributions $\mathcal{D}_V(\bar{v}_T, \theta)$, where $\mathbb{E}[V] \approx \bar{v}_T$:

- **Fixed**: $V = \bar{v}_T$ (Constant value).
- **Exponential**: $V \sim \mathrm{Exp}(\lambda=1/\bar{v}_T)$.
- **Normal / LogNormal**: Parameterized by mean $\bar{v}_T$ and Coefficient of Variation ($\mathrm{CV}$).

#### 2.3.2 Shape Factor and Side Length Derivation

To decouple volume from aspect ratio, we introduce a $d$-dimensional shape factor vector $\mathbf{g} = (g_1, \dots, g_d)$, subject to the constraint $\prod_{k=1}^d g_k = 1$.

The shape factors follow a log-normal distribution. First, latent variables $z_k \sim \mathcal{N}(0, \sigma_{\mathrm{shape}}^2)$ are sampled, followed by geometric mean normalization:

$$
g_k = \frac{e^{z_k}}{\left(\prod_{j=1}^d e^{z_j}\right)^{1/d}}
$$

Finally, the side length $\lambda_k$ for the $k$-th dimension of the object is derived:

$$
\lambda_k = V^{1/d} \cdot g_k
$$

To ensure physical feasibility, side lengths are truncated to ensure $\lambda_k < W_k$.

#### 2.3.3 Spatial Location Distribution

Objects are placed using a **Uniform Distribution** within the universe $\mathcal{U}$. For the $k$-th dimension, the lower bound $\ell_k$ is sampled from:

$$
\ell_k \sim \mathrm{Uniform}(u_k^{\min}, u_k^{\max} - \lambda_k)
$$

The upper bound is determined by $u_k = \ell_k + \lambda_k$. This location distribution forms the basis for the subsequent probabilistic derivation.

### 2.4 Inverse Parameter Solving

The core logic of the generator is to solve for $C$ such that the expected output density under this coverage equals the target value $\alpha_{\mathrm{out}}^{\star}$.

#### 2.4.1 Intersection Probability Model

Define the function $p(C)$ as the probability that a randomly selected pair of objects $(r, s)$ from $R$ and $S$ intersect, given coverage parameter $C$.

Based on the assumption of uniform spatial distribution, the probability $P_{\mathrm{1D}}$ that two intervals of lengths $a$ and $b$ intersect within a domain of length $W$ has an analytical solution:

$$
P_{\mathrm{1D}}(a, b; W) = \begin{cases} 1 & \text{if } a+b \ge W 
\\
1 - \frac{(W - a - b)^2}{(W - a)(W - b)} & \text{if } a+b < W \end{cases}
$$

In $d$-dimensional space, given side length vectors $\boldsymbol{\lambda}^{(r)}, \boldsymbol{\lambda}^{(s)}$, the intersection probability is the product of the probabilities in each dimension. Therefore, the overall intersection probability $p(C)$ is the expectation over the side length distributions:

$$
p(C) = \mathbb{E}_{\boldsymbol{\lambda}^{(r)}, \boldsymbol{\lambda}^{(s)}} \left[ \prod_{k=1}^d P_{\mathrm{1D}}(\lambda_k^{(r)}, \lambda_k^{(s)}; W_k) \right]
$$

#### 2.4.2 Expected Density Equation

The expected output density $\alpha_{\mathrm{exp}}(C)$ is linearly related to $p(C)$:

$$
\alpha_{\mathrm{exp}}(C) = p(C) \cdot \frac{n_R n_S}{n_R + n_S}
$$

Since $p(C)$ is monotonically increasing with respect to $C$ but lacks a closed-form analytical expression (due to the integration over volume and shape distributions), the system employs a **numerical root-finding method**.

#### 2.4.3 Solving Algorithm

The system executes the following steps to solve the equation $\alpha_{\mathrm{exp}}(C) = \alpha_{\mathrm{out}}^{\star}$:

1. **Monte Carlo Estimation**: Quickly and unbiasedly estimate $p(C)$ by sampling a batch of side length vectors and using the analytical $P_{\mathrm{1D}}$ formula.
2. **Interval Bracketing**: Search for the lower and upper bounds $[C_{\mathrm{lo}}, C_{\mathrm{hi}}]$ for $C$.
3. **Log-space Binary Search**: Perform an iterative search in the $\log C$ space until the relative error $\vert \alpha_{\mathrm{exp}}(C) - \alpha_{\mathrm{out}}^{\star}\vert / \alpha_{\mathrm{out}}^{\star}$ is less than a preset threshold (e.g., $2\%$).

## 3. Interface Specification and Usage

### 3.1 Data Structure

The generated data is encapsulated in a `BoxSet` structure, which contains the following fields:

- **lower**: An $\mathbb{R}^{n \times d}$ matrix, representing the lower bound coordinates of $n$ objects in each dimension.
- **upper**: An $\mathbb{R}^{n \times d}$ matrix, representing the upper bound coordinates of $n$ objects in each dimension.
- **universe**: An $\mathbb{R}^{d \times 2}$ matrix, defining the extent of the universe.

Data precision can be configured as `float32` or `float64`. The system performs numerical stability corrections at the output layer to ensure that: 
$$
\forall i, k: \mathrm{upper}_{i,k} > \mathrm{lower}_{i,k}
$$
holds strictly under the specified precision.

### 3.2 Core Function Interface

```python
def make_rectangles_R_S(
    nR: int,
    nS: int,
    alpha_out: float,
    d: int = 2,
    universe: Optional[np.ndarray] = None,
    volume_dist: Literal["fixed", "exponential", "normal", "lognormal"] = "fixed",
    volume_cv: float = 0.25,
    shape_sigma: float = 0.0,
    tune_samples: int = 200_000,
    tune_tol_rel: float = 0.02,
    seed: int = 0,
    dtype: np.dtype = np.float32
) -> Tuple[BoxSet, BoxSet, Dict]
```

**Parameter Description**:

- `nR`, `nS`: Cardinalities of sets $R$ and $S$.
- `alpha_out`: Target output density $\alpha_{\mathrm{out}}^{\star}$.
- `d`: Spatial dimension.
- `volume_dist`, `volume_cv`: Type of volume distribution and its coefficient of variation.
- `shape_sigma`: Log-standard deviation $\sigma_{\mathrm{shape}}$ for the shape factor distribution. $\sigma=0$ indicates hypercubes.
- `tune_samples`: Number of Monte Carlo samples used for $p(C)$ estimation.

### 3.3 Output Metadata

The third element returned by the function, the `info` dictionary, contains critical data required for experiment reproducibility:

| **Field**                    | **Symbol**                      | **Description**                                              |
| ---------------------------- | ------------------------------- | ------------------------------------------------------------ |
| `coverage`                   | $C$                             | The final coverage parameter obtained via numerical solving. |
| `alpha_target`               | $\alpha_{\mathrm{out}}^{\star}$ | The target density specified by the user.                    |
| `alpha_expected_est`         | $\alpha_{\mathrm{exp}}(C)$      | Estimated expected density based on the final $C$.           |
| `pair_intersection_prob_est` | $p(C)$                          | Estimated probability of intersection for a single pair corresponding to $C$. |
| `tune_history`               | -                               | Trajectory of parameter iteration during the solving process (used for auditing convergence). |
| `params`                     | -                               | A complete snapshot of input parameters (including the random seed). |

### 3.4 Usage Example

```bash
pip install alacarte-rectgen
```

The following code demonstrates how to generate two sets of 2D rectangle data with 500,000 objects each, a target output density of 10, and how to verify the generated metadata.

```python
import numpy as np
import alacarte_rectgen as ar

# Parameter settings
N_R, N_S = 500_000, 500_000
TARGET_ALPHA = 10.0

# 1. Generate data and solve parameters
R, S, info = ar.make_rectangles_R_S(
    nR=N_R, 
    nS=N_S, 
    alpha_out=TARGET_ALPHA,
    d=2,
    universe=None,          # Default is [0, 1)^2
    volume_dist="normal",   # Use normal distribution for volume
    volume_cv=0.25,         # Coefficient of variation for volume
    shape_sigma=0.5,        # Enable aspect ratio variation
    seed=42,
    tune_tol_rel=0.01       # Solving tolerance 1%
)

# 2. Access generation results
print(f"Generated R size: {R.n}, S size: {S.n}")
print(f"Coordinates shape: {R.lower.shape}")

# 3. Audit generation parameters
print("\n--- Generation Audit ---")
print(f"Solved Coverage (C): {info['coverage']:.6e}")
print(f"Target Alpha:        {info['alpha_target']:.4f}")
print(f"Expected Alpha:      {info['alpha_expected_est']:.4f}")
print(f"Intersection Prob:   {info['pair_intersection_prob_est']:.6e}")

```


