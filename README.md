# Parse

Predicts protein regions that are disordered, and which subset of those can undergo phase separation.

# Installation

```
# install pixi
curl -fsSL https://pixi.sh/install.sh | bash
# from within the root of the repo
pixi install
```


# Usage


```
pixi run predict <sequence>

```

The Parse.f program limits primary sequence lengths to at least 25 residues, which is required for the sliding-window scheme, and no longer than 10000. 

# Output of the programme


Converted sequence
Summed classifier distance of P-labeled windows
Summed classifier distance of P-labeled windows with Uπ and Uq extension
CSV file of Residue number, Amino Acid type, Residue label, Classifier Distance, Residue label (w/ Uπ Uq extension), Classifier Distance (w/ Uπ Uq extension)

# See also

[webserver](stevewhitten.github.io/Parse_v2_web/)

Ibrahim, A. Y., Khaodeuanepheng, N. P., Amarasekara, D. L., Correia, J. J., Lewis, K. A., Fitzkee, N. C., Hough, L. E., & Whitten, S. T. (2023). Intrinsically disordered regions that drive phase separation form a robustly distinct protein class. Journal of Biological Chemistry, 299(1), 102801. https://doi.org/10.1016/j.jbc.2022.102801



