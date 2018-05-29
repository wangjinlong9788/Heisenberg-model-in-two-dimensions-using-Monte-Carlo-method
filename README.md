# Heisenberg model in two dimensions using Monte Carlo method of Metropolis–Hastings algorithm

For the Monte Carlo step, generate trial configurations by the Barker-Watts spin rotation [Chem. Phys. Lett. 3 , 144
(1969)] of randomly selected spins with a magnitude of the maximum spin rotation adjusted to ensure half of the trial
configurations are rejected in the equilibrium state.

 Hamiltonian of the  Heisenberg model
 ![image](https://github.com/wangjinlong9788/Heisenberg-model-in-two-dimensions-using-Monte-Carlo-method/blob/master/Model.PNG)

See also: Physical Review B, 75, 014425 (2007)

i and j are lattice indexes in 2D, <ij> stands for nearest neighbors, r_{ij} = r_j − r_i are lattice vectors with lattice spacing a,

and <i,j>_{r_c} stands for all spins with distance r_{ij} ≤ r_c.
 
See also: Physical Review B, 75, 014425 (2007). (*)

The parameters of the model are J-exchange interaction, A - anisotropy, and D - dipolar coupling. The boundary conditions shall be periodic.
