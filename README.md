# Hyperplanes_and_Matroids

This is code doing computations with hyperplane arrangements and oriented matroids.

See examples of usage in ``examples.jl``.

In particular, it can generate random central and affine hyperplane arrangements, and produce the central arrangement corresponding to an affine one.


This happens in ``generate_arrangement.jl``.

It can compute the oriented matroid from a hyperplane arrangement, this happens in ``hyperplane_to_matroid.jl``.

It can also tranform between different presentations of an oriented matroid (topes, circuits, cocircuits, covectors, chirotope), this is done in ``conversions.jl``.

It can check whether a matroid satisfies the circuit/cocircuit axioms, this is done in ``check_axioms.jl``.
