# R programming exercise (Group 5)

- Project: OncoSimulR, frequency-dependent fitness

- Subject: Programming and Statistics with R (Master's Degree in Bioinformatics and Computational Biology, Universidad Autónoma de Madrid).

- Authors: Sara Dorado Alfaro, Miguel Hernández del Valle, Álvaro Huertas García, Diego Mañanes Cayero, Alejandro Martín Muñoz. 

According to the specifications given, our group has worked on the second asked exercise, frequency-dependent fitness, whose results are available in the OncoSimulR_Ejercicio2 GitHub repository. The prostate-cancer_example.rmd is a markdown file, which contains the code and comments required to generate the corresponding vignette. As it can be seen, once we understood and got used to the OncoSimulR package and its functionality, we proceeded to develop some examples: Rock-paper-scissors model in bacterial community, a simple example of relationships between different bacterial species; the classic evolutionary game Hawk and Dove, a didactical example of Lotka-Volterra competition models to understand better the concepts and situations treated by the package; and three examples of cancer, situation intended to be treated by the package particularly. In the told vignette, these examples and their corresponding codes can be analysed in higher detail. After running this markdown file, a html file is created, which is also present in the repository. 

Regarding the package functionality, we have modified some functions in order to display some plots in a more comfortable way. Moreover, we have created some functions for showing a summary of the information given by multiple simulations of each situation, thus reaching a global vision of each scenario after simulating it several times. Those functions can be analysed in the file Funciones_plot_markdown.R. 

Derived from the creation of these new functions, we faced and solved different problems, which are displayed and showed in more detail in the file Problems_faced.R. 

We have also written some tests for the generated code, which can be run and analysed in the file test_functions.R. 

Finally, the repository contains three different folders: css, required for the vignette file style; Images_markdown, which contains the images used in the told vignette; and presentation, which contains the different files related to the pdf file for the exposition of this work.
