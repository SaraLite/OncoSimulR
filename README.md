# R programming exercise (Group 5): Frequency dependent fitness

- __Project__: OncoSimulR, frequency-dependent fitness

- __Subject__: Programming and Statistics with R (Master's Degree in Bioinformatics and Computational Biology, Universidad Autónoma de Madrid).

- __Authors__: Sara Dorado Alfaro, Miguel Hernández del Valle, Álvaro Huertas García, Diego Mañanes Cayero, Alejandro Martín Muñoz. 

According to the specifications given, our group has worked on the second asked exercise, frequency-dependent fitness, whose results are available in the following github repositories:

* __[OncoSimulR_Ejercicio2 GitHub repository](https://github.com/SaraLite/OncoSimulR_Ejercicio2)__
* __[OncoSimul Github repository](https://github.com/SaraLite/OncoSimul)__

Both repositories contain the same modifications and examples added with respect to the original code. Nevertheless, the difference between them is that [OncoSimulR_Ejercicio2 GitHub repository](https://github.com/SaraLite/OncoSimulR_Ejercicio2) exclusively contains the new modifications and examples added, while the second repository ([OncoSimul Github repository](https://github.com/SaraLite/OncoSimul)) contains the succesful integration of the original OncoSimulR package code with the modifications and examples proposed. 

Before explaining in more detail the content of the repositories, __as a summary, the achievements of this work are__:
      
1. __Understand how to use and familiarize with the frequency dependent fitness functionality in OncoSimulR pacakge__
2. __Add examples from the literature and map them to code:__
    + Rock-paper-scissors model in bacterial community
    + Classic evolutionary game: Hawk and Dove 
    + Didactical example of Lotka-Volterra competition models
    + Game Theory with social dilemmas of tumour acidity and vasculature
    + Prostate cancer tumour-stroma interactions
    + Evolutionary Dynamics of Tumor-Stroma Interactions in Multiple Myeloma
3. __Run simulations of different scenarios and cases where frequency-dependent fitness can make a difference__
4. __Modify the original documentation:__
    + Addition of boxplot and strichart plot functions to show the distribution of results among several simulation
    + New functionality to place the legend outside the plot in the desired position
    + Implementation of the new legend parameters in the package help documentation (man/plot.oncosimul)
    + Merging the developed code into the original vignette
    + LaTeX exposition document
    + Addition of the references used
5. __Test the correct functionality of the new functions added__

### [OncoSimulR_Ejercicio2 GitHub repository](https://github.com/SaraLite/OncoSimulR_Ejercicio2)

The vignette_group5.Rmd is a markdown file, which contains the code and comments required to generate the corresponding vignette. As it can be seen, once we understood and got used to the OncoSimulR package and its functionality, we proceeded to develop some examples: Rock-paper-scissors model in bacterial community, a simple example of relationships between different bacterial species; the classic evolutionary game Hawk and Dove, a didactical example of Lotka-Volterra competition models to understand better the concepts and situations treated by the package; and three examples of cancer, situation intended to be treated by the package particularly. In the told vignette, these examples and their corresponding codes can be analysed in higher detail. After running this markdown file, a html file is created, which is also present in the repository. 

Regarding the package functionality, we have modified some functions in order to display some plots in a more comfortable way. Moreover, we have created some functions for showing a summary of the information given by multiple simulations of each situation, thus reaching a global vision of each scenario after simulating it several times. Those functions can be analysed in the boxplot_stripchart_pop.R. The new functionality of plot.oncosimul to place the legend outside the plot is in the legend_outside.R script. Both scripts (boxplot_stripchart_pop.R and legend_outside.R) are loaded into vignette_group5.Rmd. This allows you to see the new features added without changing the initial features of the package (discussed below). 

Derived from the creation of these new functions, we faced and solved different problems, which are displayed and showed in more detail in the file Problems_faced.R. 

We have also written some tests for the generated code, which can be run and analysed in the file test_functions_group5.R. This file shows different situations in which the functions showed a wrong behavior, which led to the use of vapply to ensure the correct functionality and the considerations to be taken. 

Finally, the repository contains three different folders: css, required for the vignette file style; Images_markdown, which contains the images used in the told vignette; and presentation, which contains the different files related to the pdf file for the exposition of this work.

### [OncoSimul Github repository](https://github.com/SaraLite/OncoSimul)

### Modified the original OncoSimulR.Rmd vignette
Furthermore, the integration of the new code created in the original OncoSimulR package vignette has also been carried out. The additional examples added now correspond to section 11 of the vignette. Both codes have been successfully integrated. The functions used by our examples load the necessary functions at the beginning of the section, but are removed later at the end of the section. Thus, the results of the original vignette are not changed and they have been merged nicely. 

In the original OncoSimulR repository the new markdown has been added, its corresponding html and the Images_markdown folder where the images needed to generate the new vignette are located. The test document test_functions_group5.R has also been included in the testthat test folder, the scripts Functions_plot_markdown.R and legend_outside.R in the R folder with the rest of the scripts.  
