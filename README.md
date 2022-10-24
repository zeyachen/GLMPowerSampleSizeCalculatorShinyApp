# GLMPowerSampleSizeCalculatorShinyApp
Visualization of Score/LRT based Power Calculator for Generalized Linear Models via Shiny App.

h2('$$g(E[Y]) = g(\\mu) = \\eta = Z^{T}\\psi + X^{T}\\lambda$$'),
                           h4(HTML("<ul><li>\\(\\mu\\) is the mean of Y, \\(\\eta\\) is the linear predictor,
                                    </li><li>where Z is a vector of p covariates, X is a vector of q covariates, </li>
                                    <li>and \\(\\psi\\)  and \\(\\lambda\\) denote their respective regression coefficients.</li></ul>")),
                           h4('Our calculater returns power or sample size based on score test with \\(H_0:\\psi = \\psi_0\\). \\(\\lambda\\) is/are treated as nuisance parameter'),
                           h4('Here are the steps to use the Calculator.'),
                           h4(HTML("<ul>
                                    <li> 1. Start by selecting power or sample size you would like to compute based on Score or Likelihood Ratio test. </li>
                                    <li> 2. Z is/are the variable/s of interest. Minimum 1 z and maximum 2 each. Input numbers of Z and X, ie \\(z_{1}\\), \\(z_{2}\\)... </li>
                                    <li> 3. Set the weights of design matrix. If choose unequal, input weights by comma separated form and make sure they adds up to 1.</li>
                                    <li> 4. Input the coefficient from your glm with \\(\\psi_i\\) correponds to \\(z_{i}\\), \\(\\lambda_i\\) correponds to \\(x_i\\).</li>
                                    <li> 5. Input \\(\\hat{\\psi_i}\\) from your glm, desired \\(\\alpha_{0}\\), and sample size or nominal power.</li>
                                    <li> 6. Select corresponding family.</li>
                                    </ul>")),
                           h4('Here are instructions for plotting.'),
                           h4(HTML("<ul>
                                    <li> Our plot is done by linear interpolation of several estimates. Each polyline corresponds to a \\(\\psi_1\\) value.  </li>
                                    <li> You could use default setup to decide the number of estimates, default has 5 etimates with 3rd one being the input power or sample size and  1st, 2nd, 4th, 5th being the 0.8, 0.9, 1.1, 1.2 times of the input.  </li>
                                    <li> Or you can set it up mannually by specifying the upper bound, lower bound, and number of estimates. Calculator will generate a list of inputs evenly spread out in the range. Number of estimates is capped at 15. </li>
                                    <li> \\(\\delta\\):Difference of \\(\\psi_1\\) value is used to produce multiple polyline, \\(\\psi_1\\) values are input\\(\\psi_1\\) - 2\\(\\delta\\), input\\(\\psi_1\\) - \\(\\delta\\), input\\(\\psi_1\\), input\\(\\psi_1\\) + \\(\\delta\\), input\\(\\psi_1\\) + 2\\(\\delta\\).</li>
                                    </ul>"))
