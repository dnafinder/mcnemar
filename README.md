# mcnemar
Permorm McNemar's chi square on a 2x2 matrix<br/>
In statistics, McNemar's test is a non-parametric method used on nominal
data to determine whether the row and column marginal frequencies are
equal. It is named after Q. McNemar, who introduced it in 1947. It is
applied to 2x2 contingency tables with a dichotomous trait with matched
pairs of subjects.

Syntax: 	mcnemar(x,alpha)
     
    Inputs:
          X - 2x2 data matrix 
          ALPHA (default 0.05) 
    Outputs:
          - Chi Square critical value
          - Chi square value
          - p-value
          - Power
  Example:
In the following example, a researcher attempts to determine if a drug
has an effect on a particular disease. 

                     Drug
                 +         -
            --------------------
        +   |   101   |   59   |
            |-------------------  Placebo
        -   |   121   |   33   |
            --------------------
                                      

  x=[101 59; 121 33];

  Calling on Matlab the function: 
            mcnemar(x)

  Answer is:

Critical value at 95% fiducial level = 3.8415<br/>
McNemar chi-square (with Yates'es correction) = 20.672222    p = 0.000005<br/>
alpha = 0.0500  Zb = 2.7566  Power (2-tails) = 0.0058<br/>

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2007) McNemar test: perform the McNemar test on a 2x2
matrix. 
http://www.mathworks.com/matlabcentral/fileexchange/15472
