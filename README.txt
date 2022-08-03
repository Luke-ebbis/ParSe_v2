To use Parse_v2:
1. Compile Parse_v2.f using a fortran compiler. The program was tested using GNU Fortran (Homebrew GCC 11.3.0_2) 11.3.0.
2. Run the program by using the executable followed by the protein primary sequence, for example "./a.out SEQUENCESEQUENCESEQUENCESEQUENCE" without the quotes.

The primary sequence length must be at least 25 residues, and no longer than 10000. 

Program output:

Converted sequence
Summed classifier distance of P-labeled windows (∑ P)
Summed classifier distance of P-labeled windows + Uπ + Uq correction (∑ P+Uπ+Uq)
CSV file of Residue number, Amino Acid type, Residue label, ∑ P, Residue label w/ Uπ+Uq corrections, ∑ P+Uπ+Uq