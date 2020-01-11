# R programming exercise (Group 5): Frequency dependent fitness

- __Project__: OncoSimulR, frequency-dependent fitness

- __Subject__: Programming and Statistics with R (Master's Degree in Bioinformatics and Computational Biology, Universidad Autónoma de Madrid).

- __Authors__: Sara Dorado Alfaro, Miguel Hernández del Valle, Álvaro Huertas García, Diego Mañanes Cayero, Alejandro Martín Muñoz. 

According to the specifications given, our group has worked on the second asked exercise, frequency-dependent fitness, whose results are available in the [OncoSimulR_Ejercicio2 GitHub repository](https://github.com/SaraLite/OncoSimulR_Ejercicio2) . 
 Modified the original OncoSimulR.Rmd vignette
The vignette_group5.Rmd is a markdown file, which contains the code and comments required to generate the corresponding vignette. As it can be seen, once we understood and got used to the OncoSimulR package and its functionality, we proceeded to develop some examples: Rock-paper-scissors model in bacterial community, a simple example of relationships between different bacterial species; the classic evolutionary game Hawk and Dove, a didactical example of Lotka-Volterra competition models to understand better the concepts and situations treated by the package; and three examples of cancer, situation intended to be treated by the package particularly. In the told vignette, these examples and their corresponding codes can be analysed in higher detail. After running this markdown file, a html file is created, which is also present in the repository. 

Regarding the package functionality, we have modified some functions in order to display some plots in a more comfortable way. Moreover, we have created some functions for showing a summary of the information given by multiple simulations of each situation, thus reaching a global vision of each scenario after simulating it several times. Those functions can be analysed in the file Funciones_plot_markdown.R. 

Derived from the creation of these new functions, we faced and solved different problems, which are displayed and showed in more detail in the file Problems_faced.R. 

We have also written some tests for the generated code, which can be run and analysed in the file test_functions.R. 

Finally, the repository contains three different folders: css, required for the vignette file style; Images_markdown, which contains the images used in the told vignette; and presentation, which contains the different files related to the pdf file for the exposition of this work.

### Modified the original OncoSimulR.Rmd vignette
Furthermore, the integration of the new code created in the original OncoSimulR package vignette has also been carried out. The additional examples added now correspond to section 11 of the vignette. Both codes have been successfully integrated. The functions used by our examples load the necessary functions at the beginning of the section, but are removed later at the end of the section. Thus, the results of the original vignette are not changed and they have been merged nicely. 

In the original OncoSimulR repository the new markdown has been added, its corresponding html and the Images_markdown folder where the images needed to generate the new vignette are located. The test document test_functions_group5.R has also been included in the testthat test folder, the scripts functions_plot_markdown.R and legend_outside.R in the R folder with the rest of the scripts.  
