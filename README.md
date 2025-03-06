Scalar Curvature Calculator
This repository contains a pure Python implementation designed for calculating scalar curvature based on user-defined metrics. The project was developed as part of an engineering thesis and focuses on symbolic tensor computations using the Sympy library.
![Uploading screencapture-itensor-online-2025-03-06-16_39_32.png…]()

Overview
The Scalar Curvature Calculator is a command-line tool that accepts a metric in a simple text format and computes various tensor quantities, including:

Metric tensor components,
Christoffel symbols,
Riemann tensor,
Ricci tensor,
Einstein tensor,
And, finally, the scalar curvature.
The main focus of the project is on the precise symbolic computation of these quantities, making it suitable for research and academic work in the field of differential geometry and general relativity.

Features
Pure Python Implementation:
Developed entirely in Python without reliance on external backend frameworks.

Symbolic Computation:
Uses Sympy to manipulate symbolic expressions, ensuring high precision in mathematical calculations.

Modular Design:
The code is organized into distinct modules:

metric_models.py: Contains data models and parsing logic for the input metric.
tensor_calculations.py: Implements functions for computing tensor components and scalar curvature.
LaTeX Conversion: Additional functions generate LaTeX-formatted representations of the computed tensors for easy inclusion in academic documents.
Ease of Use:
Input metrics are provided in a simple text format, and the tool outputs both plain text and LaTeX-formatted results.

Technologies
Python 3.10+
Sympy: For symbolic mathematics and tensor computations.
Standard Python libraries: (e.g., math, re) for input processing and formatting.

Security & Reliability
Since the calculator is intended for academic use as a standalone Python tool:

Input validation is performed using controlled parsing with sympy.sympify to prevent execution of arbitrary code.
The project is designed with simplicity in mind, reducing the risk of security vulnerabilities inherent in web backends.
Future Work
Numerical Simulation Integration:
Possibility to integrate with numerical simulation data from tools like MATLAB or Mathematica.

Enhanced LaTeX Export:
Further refinements in the LaTeX output to support copy-paste functionality for academic papers.

Interactive Visualization:
Adding interactive, animated visualizations of tensor components using libraries such as Three.js.
Feel free to modify this README.md to better match your project's specifics and any additional functionalities you implement.
