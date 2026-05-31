# Code and Erratum for *Growth in affine Hecke categories*

This repository contains code accompanying the paper

> *Growth in affine Hecke categories*.

The repository is the one cited in the paper as

> *Code and more for “Growth in affine Hecke categories”*, 2026.

The files collect the computer calculations used for the examples and numerical checks discussed in the paper, in particular:

- tensor-power calculations for Kazhdan--Lusztig basis elements, including the affine $\widetilde A_2$ wall element $b_{123}$;
- polynomial-ray calculations for $\nu_\delta(b_w)$ in affine types $\widetilde A_1$ and $\widetilde A_2$;
- a Python/Colab notebook reproducing the long affine $\widetilde A_2$ tensor-power sequence and the numerical plots supporting the predicted $n^{-3/2}8^n$ growth.

## Contact

If you find any errors in the paper **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.

## Files

### `hecke_tensor_powers_clean.m`

This Magma script computes tensor powers

$$
    b^n = \sum_x a_x(n) b_x
$$

in the Kazhdan--Lusztig basis and prints:

- $\nu(b^n)=\sum_x a_x(n)$, evaluated at $v=1$;
- $\nu_\delta(b^n)$, the standard-basis coefficient sum at $v=1$;
- the number of KL basis terms with nonzero coefficient at $v=1$;
- optionally, normalized values such as
  $$
      \frac{\nu(b^n)n^{3/2}}{\beta^n},
      \qquad \beta=\nu_\delta(b),
  $$
  for the affine $\widetilde A_2$ test case.

The default parameters are

```magma
AFFINE_TYPE := "A~2";
WORD := [1,2,3];
ELEMENT_MODE := "KL";
```

so that the script studies

$$
    b=b_{123}.
$$

With these defaults, the sequence $\nu(b^n)$ begins

$$
    1,\ 3,\ 17,\ 83,\ 472,\ldots
$$

as recorded in the paper.

The file also contains an optional pretty-printer for small affine $\widetilde A_2$ expansions, using the notation $x_j,y_j,z_j,\theta(m,n)$ from the paper.

### `hecke_polynomial_rays_clean.m`

This Magma script computes the polynomial-growth quantity

$$
    \nu_\delta(b_w),
$$

i.e. the sum of standard-basis coefficients of $b_w$ at $v=1$, along the ray families considered in the paper.

It contains two blocks.

#### Affine $\widetilde A_1$

The two alternating rays

$$
    1,12,121,1212,\ldots
$$

and

$$
    2,21,212,2121,\ldots
$$

are tested against the formula

$$
    \nu_\delta(b_w)=2\ell(w)+1.
$$

#### Affine $\widetilde A_2$

The script computes $\nu_\delta$ for the wall families

$$
    x_j,
    \qquad
    y_j,
    \qquad
    z_j,
$$

using the notation of the paper. These are the singular wall rays whose polynomial growth is analyzed in the affine $\widetilde A_2$ section.

The output includes the raw values and, optionally, the diagnostic ratios $p_j/j^2$, useful for seeing the expected quadratic growth numerically.

### `KL_affineA2.ipynb`

This Python/Colab notebook gives a separate affine $\widetilde A_2$ computation for the guiding example $b_{123}$. It:

- implements affine $\widetilde A_2$ using affine permutations in window notation;
- uses the Libedinsky--Patimo families and triangular reduction to compute the tensor-power sequence recursively;
- reproduces $b_n$
  up to the longer range used for numerical experimentation;
- plots the normalized sequence
  $b_n/(n^{-3/2}8^n)$,
  together with nearby comparison normalizations, to visualize the predicted exponent.

The notebook is mainly intended as a transparent numerical companion to the affine $\widetilde A_2$ discussion in the paper. It can be opened directly in Jupyter or Google Colab.

## Requirements

### Magma files

The two `.m` scripts use:

1. **Magma**;
2. Joel Gibson's **IHecke** package.

Both files begin with

```magma
AttachSpec("IHecke/IHecke.spec");
```

Hence they should be run from a directory in which Magma can see

```text
IHecke/IHecke.spec
```

For example, if `IHecke/` lives inside `~/Documents/ASLoc`, then run:

```bash
cd ~/Documents/ASLoc
magma /path/to/growthaffinehecke/hecke_tensor_powers_clean.m
magma /path/to/growthaffinehecke/hecke_polynomial_rays_clean.m
```

### Notebook

The notebook uses standard Python together with `matplotlib` for plotting. It runs naturally in Google Colab or in a local Jupyter environment.

## Parameters

Both Magma files have a short user-parameter block near the top. The most useful options are:

### Tensor powers

```magma
NMAX := 12;
PRINT_KL_EXPANSIONS := false;
PRINT_NORMALIZED_RATIO := true;
NORMALIZING_EXPONENT := 3/2;
```

### Polynomial rays

```magma
RUN_A1 := true;
RUN_A2_WALL := true;
A1_MAX_LENGTH := 20;
A2_MAX_J := 40;
PRINT_A2_QUADRATIC_RATIOS := true;
```

## Conventions

- The Magma scripts use the generator numbering provided by Magma/IHecke.
- In affine $\widetilde A_2$, the word `[1,2,3]` denotes $s_1s_2s_3$.
- `ELEMENT_MODE := "KL"` means the single KL basis element $b_w$, not the product of simple KL generators $b_{s_{i_1}}\cdots b_{s_{i_k}}$. The tensor-power example in the paper uses the former.

## Scope

These files are intended to make the computations used in the paper transparent and reproducible. They are deliberately focused on the main computational tasks needed there, rather than on the larger collection of exploratory files used during development.


## Erratum

Empty so far.
