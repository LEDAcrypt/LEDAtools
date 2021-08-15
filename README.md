# LEDAtools
This repository contains a set of tools to automatically generate
LEDA cryptosystem parameters for a given security level, or
evaluate the security of a given parameter set.

---
Available tools

 * `parameter_generator`: generates a set of code parameters for the LEDA cryptosystem
 * `work_factor_computation`: calculates the required computational effort to break
   perform either a message recovery or a key recovery attack against a LEDA[kem|pkc]
   instance
 * `enumeration_complexity`: calculates the computational effort required to enumerate all
   the possible H and Q matrices (parts of LEDA[kem|pkc] private key)
 * `constant_weight_encodable_bits`: computes the maximum length of the bitstring encodable
   by the constant weight encoding in LEDApkc

---

Installation requirements
-------------------------

The tools are written in C++ and employ Victor Shoup's NTL library.
The building automation is performed with a plain makefile, hence GNU Make,
or any compatible tool should be available on the target system.
To compile the tools it is sufficient to type

`make`

in the directory where the sources are present.

---

LEDAkem and LEDApkc code parameter generator
--------------------------------------------
The `parameter_generator` tool computes the set of code parameters for a LEDA
cryptosystems, given the desired security level.
The tool takes as input:
  * the binary logarithm of the classical and quantum
computational complexity to be taken as target
  * the number of circulant blocks n_0 and the safety factor epsilon to obtain parameters with
a satisfactory DFR (for all the proposed parameters epsilon=0.3)
  * a lower bound on the value of the size of the circulant block from which
the parameter search should be started. Such a value allows to speed
up the search in case a rough approximation of the final value is known.

The tool can be run as

`./parameter_generator security_level_classic security_level_pq n_0 epsilon starting_prime_lb`

and will print out the security level of the considere ISDs when exploring the
parameter space, ending with a printout of the required parameters, namely
the circulant block size p, the number of errors t, the number of ones in a circulant
block of H, d_v, the number of ones per row of circulant blocks in Q.


Calculator for the complexity of an ISD attack for message/key recovery
-----------------------------------------------------------------------
The tool computes the binary logarithm of the computational time complexities
of performing an information set decoding attack against the quasi cyclic
code of which the parameters are taken as input.

The tool takes as inputs:
  * the size of the code and its dimension
  * the number of errors to be decoded
  * the size of the circulant block for quasi cyclic codes (can be set to 1 to
  indicate a non QC code)
  * a boolean value which indicates if the attack of which the complexity is computed
  is a decoding attack (`is_kra=0`) or a key recovery attack (`is_kra=1`). This takes
  into account the speedup provided by the Decode One Out of Many technique for
  decoding attacks, and the QC speedup factor for key recovery attacks.

The tool can be run as
`./work_factor_computation parameter_file.csv`

and will output the Quantum complexities for the quantum Lee-Brickell and quantum Stern ISDs, followed
by the complexities for the ISD algorithms by Prange, Lee and Brickell, Leon, Stern, May Meurer and Thomae,
and Becker Joux May and Meyer.
The tool accepts as input a comma-separated-value file, where each line counts
as a parameter set. The parameter_sets.csv file is provided as an example for
the format. Lines starting with a `#` will be ignored.

Constant weight encoding capacity computer
------------------------------------------

This tool computes the length of the longest bistring which can be encoded in a
constant weight error vector with length n and number of errors t.

The tool can thus be run from the commandline as:

`./constant_weight_encodable_bits <codeword_size> <number_of_errors>`

The tool outputs a C++ `define` statement ready for use in the `qc_ldpc_parameters.h` file of LEDApkc looking like:

`#define MAX_ENCODABLE_BIT_SIZE_CW_ENCODING (700)`

---

H and Q matrix enumeration complexity computer
----------------------------------------------
This tool calculates the computational effort required to enumerate
all the possible instances of the H or the Q matrix of a LEDA cryptosystem
private key.
The tool requires the sizes and number of asserted terms in a row of a circulant
block of H (i.e.,d_v) and in the rows of the circulant blocks of Q (i.e. m_0, m_1,
m_2, m_3), together with the circulant block size p and the number of blocks n_0.
In case the values m_2 and m_3 are not available for the chosen number of circulant
blocks, they should be passed as zeroes.

The tool can be run as follows:

`./enumeration_complexity <p> <d_v> <n_0> <m_0> <m_1> <m_2> <m_3>`

---

