# ArduinoLMI
Linear Matrix Inequality (LMI) solver for Arduino. Allows the resolution of a set of simple standard LMI problems such as robust pole placement, robust stabilisation, gain-scheduled control.
## References
If you use this code for research purposes, please cite the article  [*Hiba Houmsi, Paolo Massioni, Federico Bribiesca Argomedo, Romain Delpoux. "Embedded controller optimization for efficient electric motor drive," 2023 IEEE Vehicle Power and Propulsion (VPPC 2023), Oct 2023, Milan, Italy.`*](https://ec-lyon.hal.science/hal-04176290/). The article contains a description of the embedded implementation of the interior-point method solver used in this code.

For a general introduction on Linaear Matrix inequalities (LMI) problems, refer to  [*M. Chilali, P. Gahinet and P. Apkarian, "Robust pole 
 placement in LMI regions," IEEE Transactions on Automatic Control, vol. 44, no. 12,
 pp. 2257-2270, Dec. 1999*](https:/doi.org/10.1109/9.81120), which also provides an introduction to the interior point algorithm.

For the regional pole placement method, refer to  [Stephen Boyd, Laurent El Ghaoui, E. Feron, 
 and V. Balakrishnan, "Linear Matrix Inequalities in System and Control Theory",
 Volume 15 of Studies in Applied Mathematics, Society for Industrial and Applied 
 Mathematics (SIAM), 1994.*](https://web.stanford.edu/~boyd/lmibook/lmibook.pdf).
