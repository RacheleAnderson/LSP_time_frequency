# Time-frequency feature extraction for classification of episodic memory

This repository contains the code to reproduce the simulation study presented in the paper

Anderson, R., Sandsten, M., "Time-frequency feature extraction for classification of episodic memory". EURASIP J. Adv. Signal Process. 2020, 19 (2020).
https://doi.org/10.1186/s13634-020-00681-8

### Keywords
Time-frequency features, Classification, Non-stationary signals, Neural networks, EEG signals, Locally stationary processes (LSP), Optimal spectral estimation

### How to
Run the scripts main_MSE_methods_comparison.m to reproduce the results presented in the mean square error (MSE) evaluation simulation study. In the folder "functions" you find the implementation of the LSP inference method (HATS, see folder LSP_inference for more details on the inference method), as well as the implementation of the MSE optimal LSP kernel and computation of the corresponding multitapers. 

### Abstract of the paper
This paper investigates the extraction of time-frequency (TF) features for classification of electroencephalography (EEG) signals and episodic memory. We propose a model based on the definition of locally stationary processes (LSPs), estimate the model parameters, and derive a mean square error (MSE) optimal Wigner-Ville spectrum (WVS) estimator for the signals. The estimator is compared with state-of-the-art TF representations: the spectrogram, the Welch method, the classically estimated WVS, and the Morlet wavelet scalogram. First, we evaluate the MSE of each spectrum estimate with respect to the true WVS for simulated data, where it is shown that the LSP-inference MSE optimal estimator clearly outperforms other methods. Then, we use the different TF representations to extract the features which feed a neural network classifier and compare the classification accuracies for simulated datasets. Finally, we provide an example of real data application on EEG signals measured during a visual memory encoding task, where the classification accuracy is evaluated as in the simulation study. The results show consistent improvement in classification accuracy by using the features extracted from the proposed LSP-inference MSE optimal estimator, compared to the use of state-of-the-art methods, both for simulated datasets and for the real data example.

#### Note: 
The episodic memory EEG data used in the example is not included in this repository.
