TB NanoString Shiny App

Project Description

This project contains a Shiny application for analyzing and visualizing the TB NanoString gene expression dataset. It uses dimension reduction and clustering methods to explore differences between active TB and latent TB infection (LTBI) samples.

Files Included

app.R: R code for the Shiny application.

TBnanostring.rds: Dataset required to run the application.

Required R Packages

Install the required packages in R:

install.packages(c("shiny", "tidyverse", "ggplot2", "umap", "NMF", "mclust"))

How to Run the Application

Download or clone this repository.

Open the project folder in RStudio.

Install the required R packages.

Open app.R.

Run the following command in the R console:

shiny::runApp()

The application will open in your web browser.

Dataset

The application uses the TB NanoString gene expression dataset containing samples from individuals with active TB and LTBI.
