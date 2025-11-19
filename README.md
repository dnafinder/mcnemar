[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/mcnemar)

📘 Overview
mcnemar is a MATLAB implementation of McNemar's chi-square test for matched pairs, designed for 2×2 contingency tables with paired binary outcomes. The test evaluates whether the marginal proportions of two related conditions are equal (for example, response before vs after a treatment, or under two different diagnostic methods).

The input is a 2×2 table:

             Condition 2
               +      -
Condition 1  +  a      b
            -  c      d

The test is based on the discordant pairs b and c, and uses McNemar's chi-square statistic with Yates' continuity correction and 1 degree of freedom.

✨ Features
- Performs McNemar's chi-square test on a 2×2 matched-pairs table
- Computes:
  • Critical chi-square value at a chosen significance level  
  • Chi-square test statistic with Yates' correction (df = 1)  
  • Two-sided p-value  
  • Approximate two-sided power  
- Accepts a simple 2×2 numeric input
- Provides both printed output and a structured result (STATS)

📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/mcnemar

2. Add the folder to your MATLAB path:
      addpath('path_to_mcnemar')

3. Verify that MATLAB can find the function:
      which mcnemar

⚙️ Requirements
- MATLAB (any recent version)
- No additional toolboxes are required

📈 Usage
Drug vs placebo example:

    % 2x2 matched-pairs table
    %               Drug
    %             +     -
    % Placebo +  101   59
    %         -  121   33
    x = [101 59; 121 33];

    % Run McNemar's test with default alpha = 0.05
    STATS = mcnemar(x);

Without output argument:

    mcnemar(x);

the function prints something like:

    Critical value at 95% significance level = 3.8415
    McNemar chi-square (with Yates' correction) = 20.672222    p = 0.000005
    alpha = 0.0500  Zb = 2.7566  Power (2-tails) = 0.0058

Custom significance level:

    alpha = 0.01;
    STATS = mcnemar(x, alpha);

🔢 Inputs
mcnemar(x)
mcnemar(x, alpha)

- x     : 2×2 numeric matrix of nonnegative integers representing paired outcomes.
          x(1,1) = a (both positive),
          x(1,2) = b (positive only under condition 1),
          x(2,1) = c (positive only under condition 2),
          x(2,2) = d (both negative).

- alpha : Optional significance level in (0,1). Default is 0.05.

📤 Outputs
If called with an output argument, mcnemar returns a structure STATS with fields:

- STATS.chisq      : McNemar chi-square statistic with Yates' correction
- STATS.df         : degrees of freedom (1)
- STATS.crit       : chi-square critical value at the specified alpha
- STATS.pvalue     : p-value of the test (two-sided)
- STATS.alpha      : significance level
- STATS.Zb         : Z_beta, used in approximate power calculation
- STATS.power      : approximate two-sided power
- STATS.N          : total number of pairs (sum of all cells in x)
- STATS.table      : original 2×2 table x
- STATS.discordant : [b c], counts of discordant pairs

🧠 Interpretation
- A large chi-square value with a small p-value suggests a significant difference between the two matched conditions.
- If p < alpha, you reject the null hypothesis that the marginal proportions are equal.
- The power value is an approximation and should be treated as indicative rather than as a full design-stage power analysis.

📌 Notes
- McNemar's test is recommended for matched-pair designs with binary outcomes.
- When there are no discordant pairs (b = c = 0), the test cannot be computed; the function will issue a warning.
- For small discordant counts, McNemar's test can be compared with exact alternatives such as Liddell's test (see separate implementation).

🧾 Citation
If you use this function in publications or analyses, please cite:

Cardillo G. (2025). mcnemar: McNemar's chi-square test for matched pairs in MATLAB.  
Available at: https://github.com/dnafinder/mcnemar

👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

📄 License
mcnemar is distributed under the terms specified in the LICENSE file available at:
https://github.com/dnafinder/mcnemar
