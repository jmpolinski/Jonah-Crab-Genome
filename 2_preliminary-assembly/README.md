
# Preliminary assembly of the Jonah Crab Genome

## Testing assemblers

As there are many genome assembly algorithms available and we were unsure which would perform the best for our data, multiple were tested. Assembly statistics for all tested programs are shown below:  

| Assembler                   | WTDBG2 (DTG)   | Shasta (0.7.0) | Miniasm (0.3)  | Flye (2.9)     |
| --------------------------- |:--------------:|:--------------:|:--------------:|:--------------:|
| # Scaffolds                 | 16001          | 15822          | 7407           | 5456           |
| Scaffold total (Mbp)        | 720.682        | 554.622        | 825.136        | 701.064        |
| N50                         | 168.886 Mbp    | 152.089 Mbp    | 143.008 Mbp    | 240.256 Kbp    |
| L50                         | 969            | 979            | 1477           | 719            |
| N90                         | 15.934 Kbp     | 34.359 Kbp     | 52.11 Kbp      | 62.122 Kbp     |
| L90                         | 6991           | 3895           | 5367           | 3028           |
| Max scaffold length         | 3.553 Mbp      | 3.448 Mbp      | 2.323 Mbp      | 3.637 Mbp      |
| BUSCO dataset               | arthopod_odb10 | arthopod_odb10 | arthopod_odb10 | arthopod_odb10 |
| BUSCO (% complete/% single) | 76.2% / 75.0%  | 85.2% / 83.7%  | 2.5% / 2.4%    | 90.1% / 85.5%  |
  
Note that the CANU assembler was also tested but the assembly was aborted after running over 30 days.  
  
The Flye assembly was moved forward to Omni-C scaffolding. Code for how this assembly was generated can be found in 1_JC-flye.sh.

    
*All assemblies and QC done by Jennifer Polinski.* 
