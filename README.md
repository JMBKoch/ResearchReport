# One- vs. two-step Approach in Bayesian Regularized Structual Equation Modeling (SEM)

This repository contains the code of my Masters Thesis in Applied Statistics at Utrecht University (September '21 - May '22). 

While the concepts of *penalization* or *regularization* have been common place in regression and machine learning for a long time, their usage has only more recently been proposed in the context of structural equation modeling (SEM, Jacobucci, Grimm, McArdle, 2016). For instance, cross-loadings in a factor model can be penalized, such that only large cross-loadings are included in the model, whereas small cross loadings are shrunken to zero. 

In a Bayesian context, regularization can be achieved through so-called shrinkage priors (see van Erp, Oberski, & Mulder, 2019 for an overview). One such shrinkage prior is the *small-variance normal prior*, which corresponds to the classical ridge penalty (Muthén & Asparouhov, 2012). The issue with this shrinkage prior, however, is that it not only shrinks small coefficients (to zero), but also shrinks large coefficients (to a smaller estimate). Consequently, the estimates of the larger coefficients are biased substantially. This approach therefore requires a second step, where the model is re-estimated without the shrinkage prior, while fixing those coefficients that were shrunk to zero in the first step to zero. 

The aim of this project is to assess through simulation whether other, more sophisticated shrinkage priors (van Erp, Oberski, & Mulder, 2019) are able to outperform the small-variance normal prior, within a simpler, one-step approach.

### References

Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized structural equation modeling. Structural Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. doi: 10.1080/10705511.2016.1154793 

van Erp, S., Oberski, D. L., & Mulder, J. (2019). Shrinkage priors for bayesian penalized regression. Journal of Mathematical Psychology, 89, 31–50. doi: 10.1016/j.jmp.2018.12.004 

Muthén, B. O., & Asparouhov, T. (2012). Bayesian structural equation modeling: A more flexible representation of substantive theory. Psychological Methods, 17(3), 313–335. doi: 10.1037/a0026802
