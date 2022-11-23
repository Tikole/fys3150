Man trenger 2 make kommandoer
make all lager 3 kjørbare filer:
    - isingdb er kompilert med -g og -O0 og kan brukes
      med gdb.
    - ising er optimalisert me -O3 og kompileres med parallellisering
    - isingST er er også optimalisert -O3, men er uten parallellisering

make scrub: sletter alle filene output_*.txt i 'data' mappen, dette er den mest
      måten å bruke 'plot.py', siden denne antar alt i denne mappen skal med i plottet.

Alle ising-filene har de samme kommandoene:
    -defaults som viser hvilke verdier som er default
    -h som viser en ikke implementert hjelpemelding
     kjørekommandoen begynner med en integer som angir L, og man må angi T enten som en 
     enkelt verdi av T med -T <float>, eller ved linspace argumenter: -Tspace T0 T1 steps
    I tillegg finnes følgende:
    -random og -ordered, random er default bestemmer initialkonfigurasjonen
    -runs er antall uavhengige kjøringer som skal gjøres for hver temperatur per default
     er den 1. Alle runs skrives hver for seg til data mappen
    -sample hvor mange spins som skal flippes for hver sample, reduserer minne og diskbruk
     default er 1.
    -chunk en størrelse som brukes internt for automatisk å bedømme ekvilibrasjon.
     Bør typisk være 10 000 eller noe slikt.
    -max antall chunks å beregne. Antall samples blir altså chunk * max, gitt at man
     ikke bruker automatisk ekvilibrasjonsfunksjonen. Bruker man denne fungerer max 
     som en limiter som hindrer programmet fra evig eksekusjon.
    -eps angir "epsilonpølsa" som kvalifiserer ekvilibrium, aner ikke hva denne burde
     være, har ikke rukket å finne ut av dette. På null per default, hvilket skrur av
     autogreia
    -target minimum antall chunks systemet skal beregne, ignoreres hvis ikke -eps er
     > 0.

Eksempler:
    >>>./ising 100 -Tspace 1.8 2.4 20 -runs 8 -chunk 10000 -max 10 -sample 100
    >>>./ising 20 -T 1.0 -chunk 10000 -max 20 -target 10 -epsilon 0.1 -ordered

All plottingen gjøres av 'plot.py', men dette skriptet er bare rot og vanskelig å gi noen
instruksjoner i. Dog instrueres den 