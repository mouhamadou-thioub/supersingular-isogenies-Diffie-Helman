# supersingular-isogenies-Diffie-Helman
it is a post-quantum key exchange protocol that is based on supersingular elliptic curves
sike1.c is the protocol without optimization.
sike2.c uses optimal strategies
exec. is a small program that allows you to choose the protocol to execute among sike1 and sike2
In this work, the body we are working on is an extension of Fp. The arithmetic of this extension is not available on the GMP library, we have used some functions of the GMP library to be able to implement the arithmetic of Fp / (x ^ 2 + 1).
To have the correct settings during compilation,  consult the original document submitted to NIST:http://sike.org
