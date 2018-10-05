## Code and data for On the predictability of infectious disease outbreaks

### Citation
Scarpino, S. V., & Petri, G. (2017). On the predictability of infectious disease outbreaks. arXiv preprint arXiv:1703.07317. https://arxiv.org/abs/1703.07317

### Notes on the code
The .R files contain the code to recreate the figures and main statistical analyses. 

### Data
Empirical data for all diseases–aside from Zika and dengue–were obtained from the US. National
Notifiable Diseases Surveillance System as digitized by Project Tycho [1]. Zika data were obtained from public health reports from Colombia and Mexico as digitized by [2].  Dengue data were obtained from the Pandemic Prediction and Forecasting Science and Technology Interagency Working Group under the National Science and Technology Council [3].

[1] van Panhuis, W. G. et al. Contagious diseases in the United States from 1888 to the
present. The New England journal of medicine 369, 2152–2152 (2013).

[2] Rodriguez, D. M. et al. 10.5281zenodo.344913 (2017).

[3]  Dengue Forecasting project website. Retrieved from http://dengueforecasting.noaa.gov/  (last accessed: Oct 5th  2018, Oct 16).

### Supplement
The supplement is available at http://scarpino.github.io/files/supplementary-information-predictability.pdf.

### Acknowledgements
We thank Joshua Garland and Pej Rohani for productive conversations on permutation entropy and helpful
comments on an earlier version of the manuscript. S.V.S. received funding support from the
University of Vermont and Northeastern University. G.P. received funding support from Fondazione Compagnia San Paolo. S.V.S. and G.P. conducted performed the study as fellows at IMeRA and drafted the manuscript at [Four Corners of the Earth](https://www.instagram.com/fourcornersoftheearthdeli/) in Burlington Vermont.

### Study Abstract
Infectious disease outbreaks recapitulate biology: they emerge from the multi-level interaction of hosts, pathogens, and their shared environment.  As a result, predicting when, where, and how far diseases will spread requires a complex systems approach to modeling.  Recent studies have demonstrated that predicting different components of outbreaks--e.g., the expected number of cases, pace and tempo of cases needing treatment, demand for prophylactic equipment, importation probability etc.--is feasible.  Therefore, advancing both the science and practice of disease forecasting now requires testing for the presence of fundamental limits to outbreak prediction.  To investigate the question of outbreak prediction, we study the information theoretic limits to forecasting across a broad set of infectious diseases using permutation entropy as a model independent measure of predictability.  Studying the predictability of a diverse collection of historical outbreaks--including, chlamydia, dengue, gonorrhea, hepatitis A, influenza, measles, mumps, polio, and whooping cough--we identify a fundamental entropy barrier for infectious disease time series forecasting.  However, we find that for most diseases this barrier to prediction is often well beyond the time scale of single outbreaks, implying prediction is likely to succeed.  We also find that the forecast horizon varies by disease and demonstrate that both shifting model structures and social network heterogeneity are the most likely mechanisms for the observed differences in predictability across contagions.  Our results highlight the importance of moving beyond time series forecasting, by embracing dynamic modeling approaches to prediction, and suggest challenges for performing model selection across long disease time series.  We further anticipate that our findings will contribute to the rapidly growing field of epidemiological forecasting and may relate more broadly to the predictability of complex adaptive systems.