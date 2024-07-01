This R package was developed using the theoretical foundations presented in the article "Latent map Gaussian processes for mixed variable metamodeling" by Nicholas Oune and Ramin Bostanabad.
We provide two main functions described below, aiming to fit Latent Map Gaussian Process (LMGP) models to mixed data (categorical and quantitative).

Project Participants:
########################

- Julien Crimotel
- Gaulthier Voland
- Lucas LALOUE

Dependencies:
#############

library(Rcpp)
library(dplyr)

Function Documentation:
##############################

Function LMGP.fit
#################

The LMGP.fit function is used to fit an LMGP model to hybrid data (both character
and numeric). This function takes the data and several optional parameters as input
and returns a fitted model.

Function Parameters

• data_frame_type (data.frame, size n × p): The raw dataframe containing
explanatory variables, with columns of type "numeric" or "character".
• Y_vec (numeric, size n): The response vector.
• N_try (integer): The number of initializations to perform during optimization.
• kernel_str (character): The kernel to use ("exp" for exponential or "mat" for
Matern).
• F_matrix (matrix, size n × h): The base matrix F described in the paper.
Default is the n × 1 matrix (1)i∈[[1,n]].
44
• ncol_A (integer): The number of columns in matrix A describing the dimension
of the projection from the qualitative matrix to the latent space.
• reg (numeric, size 1): The constant to add to the diagonal of matrix R
(smoothing parameter). Default is NA. This means that the nugget is optimized
by maximum likelihood.
• silent (logical): Flag to enable/disable messages during execution.
• scale (logical): Flag to enable/disable data scaling.
• amplitude_A (numeric, size 2): Initialization amplitude for elements of matrix A.
• amplitude_w (numeric, size 2): Initialization amplitude for elements of vector ω.
• amplitude_reg (numeric, size 2): Initialization amplitude for the nugget.
• optimiseur (character): The optimization method ("L-BFGS-B","BFGS","bobyka","cmaes").
• max_it_optimization (integer): Maximum number of iterations in one optimization try.
• relative_tol_optimization (integer): Relative tolerance for the cost function optimization.
• pow_w_10 (logical): Should optimization be performed on a logarithmic scale?
• type_init (character): The initialization method ("lhs","sobol","unif").

Function Result

• model (list): The fitted model, encoded in a list, ready for direct use in the
LMGP.predict function.
Potential Errors
• No Convergence: If none of the N_try initializations converges, an error is returned.
• Ineffective Model: If the model predicts less accurately than the simple mean
of the Y_vec vector in terms of MSE, an error is returned.

Function LMGP.predict
#####################

The LMGP.predict function is used to make predictions from new data using a pretrained model. It predicts the mean and variance of the response conditionally on the
predictors.

Function Parameters

• fitted_LMGP (list): The model returned by the LMGP.fit function.
• new_data_frame_type (data.frame, size k×p): The new raw dataframe containing explanatory variables, with columns of type "numeric" or "character".
• new_F_matrix (matrix, size k × h): The new base matrix F described in the
paper. Note that the functional basis used to construct F should be the same as
the one used in the LMGP.fit function.
Function Result
• predictions (list): The first element of the list contains the predicted mean,
and the second element contains the predicted variances.

