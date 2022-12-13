Prosjekt 5:
---
Den kjørbare filen 'slitex' lages med kommandoen 'make exe'. Den krever hdf5 er tilgjengelig på systemet.
På mitt system, vanlig Ubuntu i WSL, finnes headeren hdf5.h i /usr/include/hdf5/serial/ og arkivet hdf5.a
i /usr/lib/x86_64-linux-gnu/hdf5/serial/ linker og inkluderer disse mappene.

slitex kjørbaren produserer en fil 'slitex.hdf5' med resultatene. En typisk kjøring bruker snaut to minutter.
Mulige argumenter til slitex finner man med kommandoen "./slitex -help". Kjøringene som er gjort ila. prosjektet
er: 
    - ./slitex (default uten spalter, uten vegg)
    - ./slitex -slits 2 -sigma 0.05 0.10 (2 spalter, økt bredde langs y-aksen for initialtilstanden)
    - ./slitex -slits 2 -sigma 0.05 0.20 (2 spalter, enda større sigma_y)
    - ./slitex -slits 1 -sigma 0.05 0.20 (1 spalte)
    - ./slitex -slits 3 -sigma 0.05 0.20 (3 spalter)

Plottene produseres med skriptet 'plot.py'. Tar ingen argumenter, antar at det finnes en fil 'slitex.hdf5' tilgjengelig.
Gjør alle plott, men jeg har kommentert ut animasjonsdelen, da denne tar en del tid. 2 animasjonseksempler ligger i repo'et.
'plot.py' er avhengig av pythonpakken 'h5py' for å les HDF5 filer. Finnes på pip som 'h5py'.