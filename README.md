# clas12-background
Tools to study detector rates and occupancies in GEMC beam-background simulations.

### Prerequisites
  * A Linux or Mac computer
  * Java Development Kit 11 or newer
  * maven

### Build and run
Clone this repository:
```  
  git clone https://raffaelladevita/clas12-background
```
Go to the folder clas12-background and compile with maven:
```
  cd clas12-background
  mvn install
```
Run the code with:
```
  ./bin/background 

     Usage : background [options] file1 file2 ... fileN 

   Options :
    -histo : read histogram file (0/1) (default = 0)
  -modules : colon-separated list of modules to be activated (default = )
        -n : maximum number of events to process (default = -1)
        -o : histogram file name prefix (default = )
     -plot : display histograms (0/1) (default = 1)
    -print : print histograms (0/1) (default = 0)
    -stats : histogram stat option (e.g. "10" will display entries) (default = )
     -time : simulated time window per event in ns (default = 250)
```
