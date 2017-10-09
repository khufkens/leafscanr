# leafscanr

The **leafscanr** R package provides functions to (batch) process leaf scans made on a flatbed scanner and extract individual leaf traits such as leaf area, greenness (Gcc) and redness (Rcc) values.

## Installation

clone the project to your home computer using the following command (with git installed) or use the following set of command in your R console.

```R
if(!require(devtools)){install.package(devtools)}
devtools::install_github("khufkens/leafscanr")
library(leafscanr)
```

## Use

The algorithm uses a kmeans clustering and a priori knowledge of the scanner's background matte as a way to segment the image with limited user interaction.

To run the algorithm on a batch of images put all images in a single directory with no other files and make sure to either manually exclude any reference strips or if set to the left remove it by deleting x number of pixels (see reference argument). Make sure **no leaves touch each other** or they will be grouped together. Remove any other objects which are not leaves from the scan. Results will be returned enumerated from right to left, top to bottom.

The code below looks for images in the home directory for any images. These images are then processed with the assumption of a white scanner background, and a dpi (dots per inch) scanner setting of 300. No columns of the image are removed as the reference is set to NULL and the argument is ignored. When the out_path is specified resulting data is saved in a file called **leaf\_trait\_image\_analysis.rds**. If the plot argument is set to TRUE a plot of every processed image will be shown for quality control purposes (see figure below). A test image is provided in the **analysis** folder of this github repository.

```R
leaf_traits = calculate_leaf_traits(path = "~",
                                    dpi = 300,
                                    background = "white",
                                    reference = NULL,
                                    out_path = NULL,
                                    plot = TRUE)
```

![](https://raw.githubusercontent.com/khufkens/leafscanr/master/analysis/example_analysis.png)

[an example image with enumerated leaves and summary statistics]

## Acknowledgements

This software was written for the MicroMic project.
