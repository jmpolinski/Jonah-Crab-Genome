#!/data/resources/app_modules/CAFE/release/cafe
date
version

#specify data file, p-value threshold, # of threads to use, and log file
load -i counts4cafe_2024-04.txt -p 0.01 -t 10 -l cafe_2024-04.log

#the phylogenetic tree structure with branch lengths
tree (((((((Cborealis:86.9310953,Ptrituberculatus:86.9310953):43.4769649,Esinensis:130.4080602):136.7180419,((Pclarkii:244.4875339,Hamericanus:244.4875339):12.4965591,(Pjaponicus:71.1576657,(Pchinensis:59.7075194,Lvannamei:59.7075194):11.4501463):185.8264273):10.1420092):92.6981857,Hazteca:359.8242878):120.0354805,(Lsalmonis:279.2744982,Tcalifornicus:279.2744982):200.5852701):51.3916201,(Dmagna:57.4685029,Dpulex:57.4685029):473.7828855):12.73,Dmelanogaster:543.98);

#search for 2 parameter model (insect, crustacean: brachiopoda, copepoda, malacostraca)
lambda -s -t (((((((4,4)4,4)4,((4,4)4,(4,(4,4)4)4)4)4,4)3,(3,3)3)2,(2,2)2)1,1)

# generate a report
report cafe_timetree_report.txt

