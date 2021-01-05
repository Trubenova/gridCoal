# Simulating species range expansion

In this worked out example, we create simple input data representing an expansion of a fictitious species over 10 time points spaced 20 years apart, with generation time of 2 years, on a rectangular grid of 5 x 4. We will sample two middle lines a calculate coalescence times for pairs of samples from each pair of sampled cells.


Then we run 100 coalescence simulations using the main gridCoal simulator (Main.py), generating 100 output files with coalescence times, one with a summary of inputs, as well as very detailed demography debugger produced by msprime.

Then we analyse the outputs, calculating global Fst anf F^*.
