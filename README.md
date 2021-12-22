# cdbn-home-away

This Github repository has data and scripts needed to replicate the analyses in the paper:
"Local to continental-scale variation in fitness and heritability in common bean (*Phaseolus vulgaris*)"
DOI: 10.1002/csc2.20694

## Authors
Patrick Ewing, Alice H. MacQueen, Jeffrey W. White, Phillip E. McClean, Bryan Runck, Michael Kantar, Colin Khoury, Phil Miklas

## Questions
What is the role of plant breeding on local adaptation in a selfing species?

What is the role of plant breeding on local adaptation in different domestication clades?
How does the amount of local adaptation over time?

## Data in data/ folder

phenotypes.rda is an R object that contains biomass/yield data (column: SY) for 35 years of the CDBN (column:Year) over multiple locations (column: Location_code). Bean gene pool, Race, and a CDBN ID and Seq_ID to tie to previously published GBS data is also available.

## Results in results/ folder
Contains the figures and tables generated for the published paper.

## Data in data-raw/
Contains code to prepare the datasets used in this publication.

## Code in R/
Contains code to run the analyses used in this publication.
