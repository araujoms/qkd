### Supporting files for [Quantum key distribution rates from semidefinite programming](https://arxiv.org/abs/2211.05725)

#### Mateus Araújo, Marcus Huber, Miguel Navascués, Matej Pivoluska, Armin Tavakoli

Code works in MATLAB and Octave. Requires the development version of [YALMIP](https://yalmip.github.io/) to run, and a SDP solver. Good choices are [SeDuMi](https://github.com/sqlp/sedumi) (open source, reliable, slow, runs on Octave) or [MOSEK](https://www.mosek.com/) (proprietary, less reliable, fast, requires MATLAB).

The main files are

 - [mub_qkd.m](https://github.com/araujoms/qkd/blob/main/mub_qkd.m): Computes the key rate for the MUB protocol from section 4.1 assuming the probabilities of the isotropic state.
 - [subspace_qkd.m](https://github.com/araujoms/qkd/blob/main/subspace_qkd.m): Computes the key rate for the subspace protocol from section 4.2 assuming the probabilities of the isotropic state.
 - [overlap_qkd.m](https://github.com/araujoms/qkd/blob/main/overlap_qkd.m): Computes the key rate for the overlap protocol from section 4.3 assuming the probabilities of the isotropic state.
 - [matej_experiment_overlap.m](https://github.com/araujoms/qkd/blob/main/matej_experiment_overlap.m): Computes the key rate for the time-bin entanglement experiment described in section 5.1 with the overlap protocol.
 - [matej_experiment_subspace.m](https://github.com/araujoms/qkd/blob/main/matej_experiment_overlap.m): Computes the key rate for the time-bin entanglement experiment described in section 5.1 with the subspace protocol.
 - [mehul_experiment_full.m](https://github.com/araujoms/qkd/blob/main/mehul_experiment_full.m): Computes the key rate for the pixel entanglement experiment described in section 5.1 with the full MUB protocol.
 - [mehul_experiment_restricted.m](https://github.com/araujoms/qkd/blob/main/mehul_experiment_full.m): Computes the key rate for the pixel entanglement experiment described in section 5.1 with the restricted MUB protocol.
