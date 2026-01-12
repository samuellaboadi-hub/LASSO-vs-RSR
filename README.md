# LASSO-vs-RSR
Implementation and comparison of LASSO and Ranked Sparsity Regularization (RSR) using coordinate descent, evaluating structured sparsity and prediction accuracy on real data.
This project implements and compares LASSO and Ranked Sparsity Regularization (RSR) for feature selection and predictive modeling using a coordinate descent optimization framework. While LASSO enforces unstructured sparsity by penalizing all predictors equally, RSR introduces structured sparsity by assigning different penalty weights to features based on their rank (e.g., main effects vs. interactions), allowing for more flexible and interpretable models.
The goal of this project is to evaluate how structured regularization impacts model sparsity, interpretability, and prediction accuracy when applied to real-world data.
Methods
Implemented coordinate descent algorithms for both LASSO and RSR.
Constructed rank-based feature groups (e.g., main effects, interactions, polynomial terms).
Applied rank-weighted soft-thresholding for RSR.
Tuned regularization parameters using cross-validation.
Compared models using:
Prediction error (MSE)
Sparsity (number of nonzero coefficients)
Model interpretability
Data
The methods were evaluated on real-world datasets, including the mtcars dataset, which contains automobile characteristics such as weight, horsepower, and number of cylinders, with fuel efficiency (mpg) as the response variable.
Key Findings
LASSO produces highly sparse and interpretable models dominated by main effects.
RSR selects structured interactions and nonlinear terms while controlling complexity.
RSR captures richer relationships and can improve prediction accuracy in some settings.
Structured sparsity provides a principled way to incorporate domain knowledge into regularization.

