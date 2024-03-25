# HaploCestry

## Overview

HaploCestry is a web application designed to facilitate the exploration of paternal lineages using haplogroups. It utilizes data from the [Allen Ancient DNA Resource (AADR)](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data) database and the International Society of Genetic Genealogy (ISOGG) [Y-DNA Haplogroup Tree](https://isogg.org/tree/) to dynamically generate ancestry reports and visualizations.

## Features

- **Ancestry Exploration:** Users can explore their paternal lineage history by inputting their haplogroup ID (in ISOGG format) into the application.
- **Ancestry Reports:** HaploCestry generates detailed reports, including information on when the haplogroup formed, its parent haplogroup, immediate descendant haplogroups, and the number of DNA-tested descendants associated with each haplogroup.
- **Country of Origin:** The application provides information on the countries of DNA-tested descendants associated with each haplogroup.
- **Visualizations:** Users can visualize their ancestry through simple ancestry plots generated by the application.

## Getting Started

To use HaploCestry, follow these steps:

1. Visit the [HaploCestry web app](https://elijahugoh.shinyapps.io/HaploCestry/).
2. Enter your haplogroup ID (in ISOGG format) into the text input field.
3. Click the "Submit" button to generate your ancestry report and visualization.

## Running HaploCestry Locally
You can also clone the repository to a local difrectory and run HaploCestry from [RStudio](https://posit.co/download/rstudio-desktop/). 

```bash
# Run the command from a terminal in a chosen directory
git clone https://github.com/Elijah-Ugoh/HaploCestry.git
```
Make sure all the app files, including tree data are saved in the same directory as the app.  

## Requirements

HaploCestry requires a modern web browser with JavaScript enabled to function properly.

## Limitations

- The accuracy of ancestry reports and visualizations depends on the completeness and accuracy of the underlying data sources.
- The application currently only analyzes paternal lineage ancestry and may not provide insights into maternal lineage or admixture ancestry.

## Contributing

Contributions to HaploCestry are welcome! If you encounter any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request on GitHub.

## License

HaploCestry is licensed under under the [GPL-3.0 license](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Acknowledgements

HaploCestry utilizes data from the Allen Ancient DNA Resource (AADR) database and the International Society of Genetic Genealogy (ISOGG) Y-DNA Haplogroup Tree. Special thanks to the contributors of these datasets.
