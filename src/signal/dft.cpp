// Include required libraries.
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>


// Define a Discrete Fourier Transform function.
std::vector<std::complex<double>> dft(const std::vector<std::complex<double>> &signal) {
    // Store the length of the signal.
    const int num_samples = signal.size();

    // Create a vector to store the Fourier coefficients.
    std::vector<std::complex<double>> fourier_coefficients;
    fourier_coefficients.reserve(num_samples);

    // Calculate each Fourier coefficient (k).
    for (int k = 0; k < num_samples; k++) {
        // Create a variable to calculate the Fourier coefficient.
        std::complex<double> fourier_coefficient = (0.0, 0.0);

        // Calculate each sample's (n) contribution.
        for (int n = 0; n < num_samples; n++) {
            // Calculate the angular term.
            const double theta = (2.0 * M_PI * k * n) / num_samples;

            // Calculate the cosine and sine terms.
            const double real = std::cos(theta);
            const double imag = std::sin(theta);

            // Update the Fourier coefficient.
            const std::complex<double> coeff(real, -imag);
            fourier_coefficient += signal[n] * coeff;
        }

        // Store the Fourier coefficient.
        fourier_coefficients.push_back(fourier_coefficient);
    }

    return fourier_coefficients;
}


// Define a function that shows a given number of Fourier coefficients.
void show_fourier_coefficients(
    const std::vector<std::complex<double>> &fourier_coefficients, 
    const int &num_coefficients
) {
    // Store the number of samples.
    const int num_samples = fourier_coefficients.size();

    // Output general information.
    std::cout << "----" << std::endl;
    std::cout << "First " << num_coefficients << " Fourier Coefficients" << std::endl;
    std::cout << std::endl;

    // Output headers.
    std::cout << "k\t" << std::setw(12) 
            << "Real\t" << std::setw(12) 
            << "Imag\t" << std::setw(12) 
            << "Magnitude" << std::endl;
    
    // Output each Fourier coefficient (k).
    for (int k = 0; k < num_coefficients; k++) {
        std::cout << k << "\t" << std::setw(12) 
                << fourier_coefficients[k].real() / num_samples << "\t" << std::setw(12) 
                << fourier_coefficients[k].imag() / num_samples << "\t" << std::setw(12)
                << std::abs(fourier_coefficients[k]) / num_samples << std::endl;
    }

    // Output end.
    std::cout << "----" << std::endl;
}


int main(void) {
    // Create a signal.
    // Define the number of samples.
    const int num_samples = 1000;

    // Define the frequency and phase.
    const double freq = 3.0;
    const double phase = 0.0;

    // Create a vector to store the signal.
    std::vector<std::complex<double>> signal;
    signal.reserve(num_samples);

    // Calculate each sample (n).
    for (int n = 0; n < num_samples; n++) {
        // Calculate the cosine and sine terms.
        const double real = std::cos(((2.0 * M_PI * freq * n) / num_samples) + phase);
        const double imag = 0.0;

        // Store the sample.
        const std::complex<double> sample(real, -imag);
        signal.push_back(sample);
    }
    
    // Calculate the Discrete Fourier Transform of the signal.
    const std::vector<std::complex<double>> fourier_coefficients = dft(signal);
    
    // Output the Fourier coefficients.
    show_fourier_coefficients(fourier_coefficients, 5);

    return EXIT_SUCCESS;
}
