# Installation and Usage Instructions

## 1. Prerequisites

To compile and run this program, you must have the **TREXIO** library installed. TREXIO depends on the **HDF5** library.

### Step 1: Install HDF5
* **Ubuntu/Debian:**
  ```bash
  sudo apt install libhdf5-dev
  ```

* **macOS:**
  ```bash
  brew install hdf5
  ```

 ### Step 2: Install TREXIO
  ```bash
  wget [https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz](https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz)
  tar -zxvf trexio-2.5.0.tar.gz
  cd trexio-2.5.0
  ./configure
  make
  sudo make install
 ```
This installs the library to `/usr/local/lib` and headers to `/usr/local/include`.

## 2. Compilation

First navigate to the Hartree-Fock source directory and compile the code using `gcc`. You must link the TREXIO library.

```bash
cd Heredia_Cazzanti_Sujal_HF_MP2/HF
gcc -I/usr/local/include -L/usr/local/lib -ltrexio HF.c -o hf_calc
```
After the complilation of HF is done, navigate to MP2 source directory and complie the code using `gcc` and do not forget to link TREXIO library.

```bash
cd Heredia_Cazzanti_Sujal_HF_MP2/MP2
gcc -I/usr/local/include -L/usr/local/lib -ltrexio MP2.c -o mp2_calc
```

**Note:** If you encounter an error about loading shared libraries, add the library path to your environment: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib`


## 3. How to run

The program currently expects the input file (e.g., `c2h4.h5`) to be present in the same directory where you run the executable.


* Copy the data file from the data/ folder to your executable folder:
  ```bash
  cp ../../data/c2h4.h5 .
  ```
* Run the program
  ```bash
  ./hf_calc
  ```

For the `c2h4.h5` (Ethylene) molecule, the HF code will output:

* Nuclear repulsion energy

* One-electron energy

* Coulomb and Exchange energy contributions

* Final Total Energy

and MP2 code will output correlation energy (MP2 correction) that could be added to the HF energy.

