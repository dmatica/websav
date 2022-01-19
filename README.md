# webSAV: a browser-based approach to Sequencing Analysis Viewer

Illumina’s Sequencing Analysis Viewer is an application that allows for reviewing the quality metrics of a Illumina sequencing run. This is a vital step in the Next Generation Sequencing workflow, to make sure the sequencing run itself performed as expected before proceeding with the downstream data analysis. Each sequencing run generates an InterOp folder, which contains a set of a binary files that provide run performance metrics, from broad metrics, such as the overall amount of data generated from the run and the percentage of reads with a Q-score greater than 30, along with more granular information pertaining to each specific lane and cycle. 

A key issue with Sequencing Analysis Viewer (SAV) is that the compatibility is limited to computers running Windows. With a growing number of scientists using either Linux or Mac OS as their primary computer, this creates a notable gap in terms of access to the software.

One tool available is the InterOp package released by Illumina, which is a command-line parser for Illumina sequencing runs, but does not provide the same user interface/experience as SAV. Here, I would like to introduce webSAV as an OS-agnostic web browser-based approach to review Illumina sequencing performance data.
