# boom-bust-genetics

Boom-bust population dynamics drive rapid genetic change

**There are five analysis script files associated with this analysis**

*There are also four function R scripts and one script containing the libraries used for the analysis. These file are called from the main analysis scripts.*

Analysis scripts are labeled 01\_\* to 05\_\* and should be run sequentially.

For neEstimator and stairway analysis it will be necessary to download the associated software to run the analysis.

**There are four data files associated with this analysis.**

-   pherm_genotypes.rds
-   syoung_genotypes.rds
-   ecological_data.csv
-   em_gridcoortable.csv

### Details for: pherm_genotypes.rds

-   Description: a rdata file containing the single-nucleotide polymorphism (SNP) genotypes for 1399 *Peusdomys hermannsburgensis* samples and their associated metadata.

-   Format(s): .rds

-   Size(s): 88.6 MB

-   Missing data codes: missing genotypes are coded as NA

-   Other encoding details: 0 = homozygote, 1 = heterozygote, 2 = homozygote (for alternate allele)

-   Metadata variables:

    -   id: individual identification number
    -   pop: site location
    -   lat: latitude
    -   lon: longitude
    -   gridId: grid location
    -   trip: month captured
    -   species: *Peusdomys hermannsburgensis*

### Details for: syoung_genotypes.rds

-   Description: a rdata file containing the single-nucleotide polymorphism (SNP) genotypes for 677 *Sminthopsis youngsoni* samples and their associated metadata.

-   Format(s): .rds

-   Size(s): 19.3 MB

-   Missing data codes: missing genotypes are coded as NA

-   Other encoding details: 0 = homozygote, 1 = heterozygote, 2 = homozygote (for alternate allele)

-   Metadata variables:

    -   id: individual identification number
    -   pop: site location
    -   lat: latitude
    -   lon: longitude
    -   gridId: grid location
    -   trip: month captured
    -   species: *Sminthopsis youngsoni*

### Details for: ecological_data.csv

-   Description: a comma-delimited file containing rainfall and capture data per month. Also included are years since our defined population booms.

-   Format(s): .csv

-   Size(s): 15.7 KB

-   Dimensions: 228 rows x 8 columns

-   Variables:

    -   trip: month
    -   rain: average rainfall over our study grids per trip
    -   captures: *Pseudomys hermannsburgensis* mean captures per trip
    -   capturesSy: *Sminthopsis youngsoni* mean captures per trip
    -   year: year of trip
    -   sinceEvent: Period (how we subset our data based on rainfall and captures)
    -   monthsSince: months since the beginning of the boom
    -   yearsSince: years since the beginning of the boom

-   Missing data codes: NA

### Details for: em_gridcoortable.csv

-   Description: a comma-delimited file containing the coordinates of our sampling grids

-   Format(s): .csv

-   Size(s): 3.2 KB

-   Dimensions: 56 rows x 6 columns

-   Variables:

    -   gridId: grid id
    -   site: site code
    -   lat: latitude
    -   lon: longitude
    -   coor: True coordinates (grid) or an average site coordinate (site mean)
    -   siteName: Full site name

-   Missing data codes: blank cell
