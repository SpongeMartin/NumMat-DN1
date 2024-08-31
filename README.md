# Martin Starič

## QR razcep zgornje hessenbergove matrike


Projekt vsebuje funkcije za QR razcep hessenbergove matrike z Givens rotacijami in QR iteracijo in primerjavo časovne zahtevnosti napisanih funkcij z vgrajenimi.
### Zaganjanje kode

Kodo zaženemo tako da:

1. **Aktiviramo okolje:**
   - Odpri Julia REPL in pojdi v način pkg tako da napišeš `]`.
   - Aktiviraj okolje z ukazom:
     ```julia
     activate DN1
     ```

2. **Uporaba kode:**
   - Sledi primeru v `\test\demo.jl`.

3. **Testi:**
    - Testi so napisani v `\test\runtests.jl`.

### Generiranje .tex datoteke

`.tex` datoteka je generirana s skripto `\build\makedoc.jl`, ki uporablja paket `Weave.jl`.



