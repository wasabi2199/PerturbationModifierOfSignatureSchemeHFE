Master thesis 2023. Implementation of the HFE signature scheme without/with perturbation modifier.

All possible command line options are:
-h By using this switch, we can display all defined switches along with their descriptions.
-k is used to initiate the generation of a private and public key.
-s is used to initiate the generation of a signature using the corresponding private key.
-v is used to initiate the verification of a signature using the corresponding public key.
-n sets the parameter n, which is the size of the polynomial system.
-d sets the parameter d, which is the degree of the HFE polynomial.
-t sets the parameter t, which is the perturbation dimension.
-a sets the parameter a of the "minus" modifier.
-i sets the path to the directory with keys (where they are located during signing/verification or where the generated keys should be stored).
-f sets the path to the file we are duplicating (or for which we are verifying the signature).
-g sets the path to the directory where the signature will be saved (or where we can find the signature we are verifying).
-p - this switch activates the perturbed version of the HFE scheme.
