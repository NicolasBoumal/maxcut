# Max-Cut
Optimization algorithms to compute approximate solutions to the Max-Cut problem

These are codes from 2016 used for the experiments in the following paper:

> N. Boumal, V. Voroninski and A. S. Bandeira,
> 
> The non-convex Burer-Monteiro approach works on smooth semidefinite programs,
> 
> in the [proceedings of NIPS 2016]([url](https://papers.nips.cc/paper_files/paper/2016/hash/3de2334a314a7a72721f1f74a6cb4cee-Abstract.html)).

We followed up with more theory in another paper:

> N. Boumal, V. Voroninski and A. S. Bandeira,
> 
> Deterministic guarantees for Burer-Monteiro factorizations of smooth semidefinite programs,
> 
> in [Communications on Pure and Applied Mathematics, 2019]([url](https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.21830)https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.21830).

Usage:

1. Install [Manopt](https://www.manopt.org), [SDPLR](https://yalmip.github.io/solver/sdplr/) and [CVX](https://cvxr.com/cvx/) (or a subset, depending on which methods you want to test)
2. Make sure they are all on your Matlab path
3. Run rudytest.m for the actual experiments (edit it first to choose which graphs to run on etc.)
4. Run rudytest_to_table.m (this produces rudytable.tex with latex code for a table of results; without header though: see the paper)

Uploaded to github on March 25, 2024.
