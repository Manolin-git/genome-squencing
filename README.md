# genome-squencing

This repository contains the code and supplementary materials for our study on applying variational quantum algorithms (VQAs) to the genome assembly problem.

## Before Using it

If you clone the repo, please install the dependencies in order to use them, and simply do below if you are using `pip`.

To construct the cost function or Hamiltonian for VQA in HOBO encoding, we use [PyHOBO](https://github.com/Manolin-git/pyhobo) which in turn rely on [OpenFermion](https://quantumai.google/openfermion). These can be installed using the following command

```bash
pip install openfermion
pip install pyhobo
```

For running VQA optimization loop, we use various Qiskit libraries -

```bash
pip install qiskit[visualization]==1.0.2
pip install qiskit_aer
pip install qiskit_ibm_runtime
pip install matplotlib
pip install pylatexenc
pip install qiskit-algorithms
```

## Content

### DBG Approach

These results are presented in `vqa-DBG-7-nodes` for 7 nodes system (as discussed in report). The folder contains `data` file which contains the graph details such as weight matrix and edge frequencies.

There are three anstaz in consideration : product-state, block-linear-entanglement and full-linear-entanglement. These contain `.py` and `.ipynb` files which contain simulation for various optimizers. The interactive python book contains further description about the code.

<!-- LICENSE -->

## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

## Contact

Created by [Himanshu](mailto:himanshu.sahu@partner.ibm.com) - feel free to contact me!
