# TCCM Homework 1: Computation of the MP2 Energy

## 📌 Project Description
This project implements a C program to compute the **Hartree-Fock (HF) energy** and **MP2 correlation energy** for a closed-shell system. It reads molecular integrals and orbital energies from **TREXIO** files (HDF5 format) and calculates the total energy.

The current implementation focuses on the Hartree-Fock self-consistent field (SCF) energy calculation, as described in the `mp2.pdf` assignment.

## 📂 Directory Structure
The repository is organized as follows:

* **`Heredia_Cazzanti_HF_MP2/`**: Main source code directory.
  * **`HF/`**: Contains the source code for the Hartree-Fock calculation (e.g., `HF.c`).
  * **`MP2/`**: Contains the source code for the MP2 energy correction.
* **`data/`**: Contains the input data files (e.g., `c2h4.h5`) required for the calculation.
* **`mp2.pdf`**: The official assignment documentation and theoretical background.
* **`INSTALL.md`**: Instructions for compiling and running the program.
