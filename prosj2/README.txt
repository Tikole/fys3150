Kommando for å bygge programmet: ($ make build). Dette produserer en fil kalt 'main.exe'.
($ make all) vil i tillegg til å bygge programmet kjøre alle mulige tester. 
Testene kan også kjøres etter bygging med ($ make test) eller ($ ./main.exe test).
For å produsere data om Jakobimetodens skalering med matrisestørrelse brukes kommandoen
($ ./main.exe scaling <N>) hvor argumentet angir hvor store matriser den skal teste, den tester
med N=3,4,...100, resultatene lagres i to filer 'scaling_dense.txt' og 'scaling_tridiag.txt'.
Jeg har ikke prøvd med noe mer enn 100. Løsning av ligningen produseres med kommandoen
($ ./main.exe solve <n>) hvor n er en integer som angir antall steps, har også her bare prøvd opp
til n=100, løsningen tilhørende de tre laveste egenverdiene lagres, sammen med analytiske
løsninger i filen 'solution_<n>.txt. Alle mulige plots produseres ved kjøring av 'plot.py'
uten argumenter, etter at dataene er produsert med main.exe.