# Física analítica de la cadena 1D — referencia para validación y consolidación

> Documento de trabajo para el refactor de LinQT. Todas las expresiones de abajo
> están **verificadas simbólicamente con sympy** (`derive_core.py`, `derive_spin.py`)
> y reproducidas en el módulo de referencia `chain_reference.py`, que se usa
> directamente como golden master. Donde hay una elección de convención de
> normalización se indica explícitamente, porque es justo donde un test puede
> "fallar" sin que haya un error físico.

---

## 0. Qué hay aquí y por qué

El objetivo es tener **todos los resultados analíticos cerrados** de la cadena lineal
monoatómica que sirven para (a) validar la implementación y (b) decidir la
consolidación del refactor. El punto central que pediste —que la distinción entre los
dos *flavors* (`state-projection` de Aron / `density-projection` de José) se vea
reflejada en los cálculos— se materializa en la **Sección 3**: en la cadena de canje
ambos flavors dan **exactamente lo mismo** (test de equivalencia → consolidación), y
en cuanto el campo efectivo de espín depende de `k` **divergen de forma analítica**
(test de discriminación). Eso es exactamente lo que afirmé en la discusión previa,
ahora con forma cerrada.

Hay dos ejes ortogonales que conviene no mezclar:

- **Eje "representación de la conductividad"** (Kubo-Greenwood / VAC / MSD): cómo se
  expande la misma σ. Sector de carga, Sección 2.
- **Eje "flavor de proyección"** (state vs density): cómo se prepara el conjunto
  resuelto en energía para la polarización de espín. Sector de espín, Sección 3.

---

## 1. Convenciones y modelo

Cadena monoatómica, un orbital por sitio, salto a primeros vecinos, condiciones
periódicas, energía in-situ nula.

| Símbolo | Significado |
|---|---|
| `γ` | amplitud de salto (`> 0`); es la `t` de la tesis |
| `a` | distancia interatómica |
| `t` | **tiempo** (sin colisión de nombre con el salto) |
| `ε(k)` | `-2γ cos(ka)`, con `k ∈ (-π/a, π/a]` |
| `ΔE` | semiancho de banda `= 2γ`; centro `Ē = 0` |
| `ε̃(k)` | `ε/ΔE = -cos(ka)` (espectro reescalado a `[-1,1]`) |
| `v(k)` | `(1/ħ) ∂ₖε = (2γa/ħ) sin(ka)` |

Todas las cantidades son **por sitio**.

Identidad útil (la usa todo lo de abajo): como `ε̃(k) = -cos(ka)`,

$$ T_m(\tilde\varepsilon(k)) = T_m(-\cos ka) = (-1)^m \cos(m\,ka). $$

---

## 2. Sector de carga (cadena sin espín)

### 2.1 Banda, velocidad, DOS

$$ \rho(E) = \frac{1}{\pi\sqrt{4\gamma^2 - E^2}}\quad (|E|<2\gamma), \qquad
   v^2(E) = \frac{a^2}{\hbar^2}\,(4\gamma^2 - E^2). $$

Singularidades de van Hove en los bordes `E = ±2γ` (`ρ → ∞`, `v → 0`).

### 2.2 Momentos de Chebyshev — formas cerradas exactas

**DOS.** El momento es una delta de Kronecker (verificado para `m = 0..8`):

$$ \bar C^{\mathrm{DOS}}_m \;=\; \frac1N\,\mathrm{Tr}\,T_m(\tilde H)
   \;=\; \frac{1}{2\pi}\int_{-\pi}^{\pi} T_m(-\cos\theta)\,d\theta \;=\; \delta_{m,0}. $$

Reinsertado en la reconstrucción KPM reproduce **exactamente** la DOS cerrada de
arriba (solo sobrevive `m = 0`). Esto hace de la DOS de la cadena limpia un test
analítico trivial y muy sensible: cualquier fuga de peso a `m ≥ 1` delata un bug en la
traza estocástica o en el reescalado del espectro.

**Kubo-Greenwood.** El elemento de matriz de velocidad es diagonal en `k`,
`⟨k|V̂|k'⟩ = (2γa/ħ) sin(ka) δ_{kk'}` (tesis B.9). La matriz de momentos
$C^{KG}_{mn} = \mathrm{Tr}[T_m V T_n V]/[N(1+\delta_{m0})(1+\delta_{n0})]$ tiene la
forma cerrada **dispersa** (tesis B.15), que sympy confirma término a término hasta
`m,n = 6`:

$$ C^{KG}_{m,n} \;=\; \frac{\gamma^2 a^2}{\hbar^2}
   \left[\frac{\delta_{mn}}{1+\delta_{m0}+\delta_{m1}}
   \;-\; \tfrac12\,(1-\delta_{mn})\,\delta_{|m-n|,2}\right]. $$

En unidades de `γ²a²/ħ²` los primeros bloques son:

```
        n=0   n=1   n=2   n=3   n=4
 m=0 [  1/2    0   -1/2    0     0  ]
 m=1 [   0    1/2    0   -1/2    0  ]
 m=2 [ -1/2    0    1     0   -1/2 ]
 m=3 [   0   -1/2    0    1     0  ]
 m=4 [   0     0   -1/2   0    1   ]
```

Solo diagonal y segunda subdiagonal son no nulas → test barato y muy específico.

> **Nota de normalización (importante para el test).** El factor
> `(1+δ_{m0})(1+δ_{n0})` está en el **denominador** de la definición de la tesis.
> La integral cruda `Tr[T_m V T_n V]/N` da `C^{KG}_{00} = 2`, no `1/2`; al dividir por
> `(2)(2)=4` recupera `1/2`. El golden master debe usar **la misma** convención que
> compute LinQT en la reconstrucción (ec. 4.6 de la tesis). `chain_reference.kg_moment`
> usa la convención B.15.

### 2.3 Régimen balístico dependiente del tiempo

En la cadena limpia `[H, V] = 0` (ambos diagonales en `k`), de modo que la velocidad
en Heisenberg no decae, `V(t) = V`. De ahí, en forma cerrada:

$$ C_{vv}(E,t) = v^2(E)\ \text{(constante)},\qquad
   \Delta X^2(E,t) = v^2(E)\,t^2,\qquad
   \sigma(E,t) = e^2\rho(E)\,v^2(E)\,t. $$

Sustituyendo `ρ v²`:

$$ \boxed{\ \sigma(E,t) = \frac{e^2 a^2}{\pi\hbar^2}\,\sqrt{4\gamma^2 - E^2}\;t\ } $$

**Crece linealmente en `t` y diverge**: la cadena limpia no tiene límite DC. Esto no es
un defecto numérico sino el contenido físico del régimen balístico (review eq. 53–54 y
Fig. 5). Las representaciones VAC y MSD deben dar **idéntico** `σ(E,t)` aquí
(coincidencia exacta) — ese es el test de equivalencia VAC↔MSD.

### 2.4 La cantidad DC finita: conductancia balística

Lo finito no es σ sino la conductancia, vía la longitud de propagación
`L(E,t) = 2√(ΔX²) = 2 v(E) t` (review ec. 56–58):

$$ g(E) = \frac{A}{L}\,\sigma(E,t)\ \xrightarrow{\text{1D}}\ \text{plateau plano},
   \quad g = \frac{2e^2}{h}\ \text{(cuanto de conductancia)}. $$

Como `ρ(E) v(E) = a/(πħ)` es **independiente de la energía**, el plateau es plano sobre
toda la banda (salvo sobre-disparo en `E = ±2γ`, donde mezcla de bandas y resolución
KPM rompen la aproximación — review Fig. 6). Test: planitud + valor del cuanto.

---

## 3. Sector de espín (cadena magnética) — los dos flavors

Modelo de la tesis: `H = H₀ - J_ex σ_z` (canje homogéneo a lo largo de `ẑ`), espín
inicial a lo largo de `+x̂`. En cada bloque-`k` el espacio es 2D y
`H_k = ε(k)·𝟙 - J_ex σ_z`.

### 3.1 Precesión — forma cerrada

La fase orbital `e^{-iε(k)t/ħ}` se cancela en `U† σ U` (la parte orbital conmuta con el
espín), de modo que la dinámica de espín está gobernada **solo** por el canje. sympy da:

$$ \langle \hat S_x\rangle(t) = \tfrac{\hbar}{2}\cos\!\Big(\tfrac{2J_{ex}t}{\hbar}\Big),
\quad
\langle \hat S_y\rangle(t) = -\tfrac{\hbar}{2}\sin\!\Big(\tfrac{2J_{ex}t}{\hbar}\Big),
\quad
\langle \hat S_z\rangle(t)=0. $$

Frecuencia `|f| = J_ex/(πħ)`. **Independiente de `E_F` y de la estructura de banda `γ`**
(tesis 4.12 / Fig. 4.3). Es el test cerrado más limpio que tenemos.

### 3.2 Los dos flavors coinciden aquí — y por qué

- **density-projection (José)**: `⟨S⟩ = Tr[U† S U · P† δ(H-E_F) P] / Tr[P† δ(H-E_F) P]`,
  con `P = P(γ,φ,s)` el proyector de espín. La polarización vive en el proyector.
- **state-projection (Aron)**: `⟨S⟩ = [½⟨ψ(t)|s δ(E-H)|ψ(t)⟩ + h.c.] / ⟨ψ(t)|δ(E-H)|ψ(t)⟩`,
  con `|ψ(0)⟩ = ½(𝟙 + ĵ·ŝ)|φ_r⟩` (review ec. 175–176). La polarización vive en el estado.

Para la cadena de canje el bloque de espín es **idéntico en todo `k`** (el campo
efectivo `-J_ex ẑ` no depende de `k`). sympy confirma que ambas expresiones colapsan al
mismo resultado:

$$ \langle \sigma_x\rangle_{\text{density}}(t) = \langle \sigma_x\rangle_{\text{state}}(t)
   = \cos\!\Big(\tfrac{2J_{ex}t}{\hbar}\Big). $$

→ **Test de equivalencia / consolidación.** Si los dos flavors no coinciden (a nivel de
ruido estocástico) en este modelo, el refactor introdujo una asimetría espuria.

### 3.3 Dónde divergen — campo efectivo dependiente de `k`

Modelo-juguete analítico que aísla la diferencia:
`H_k = ε(k)·𝟙 + B(k) σ_z`, con `B(k) = J_ex + δ cos(ka)`. La frecuencia de precesión
pasa a depender de `k`: `ω(k) = 2B(k)/ħ`. Entonces:

- **density (E_F afilada).** La capa de energía en 1D es `{+k₀, -k₀}` y `ε` es par, así
  que `B(+k₀) = B(-k₀)` → **una sola frecuencia, sin desfase**, pero ahora
  **dependiente de `E_F`**:
  $$ \omega_{\text{dens}}(E_F) = \frac{2}{\hbar}\Big(J_{ex} - \frac{\delta\,E_F}{2\gamma}\Big). $$
- **state (ventana de energía `η`).** Promedia `cos(2B(k)t/ħ)` sobre un rango de `k`
  (de energías) → **desfase** (la `s(t) = ℒ(E)∘cos(ω(E)t)` de Cummings, ec. 1).

Demostración numérica (`γ=1, J_ex=0.1, δ=0.05, η=0.3`, `E_F=0`):

```
 t        :   0      50      100     150     200
 density  : 1.000  -0.839   0.408   0.154  -0.667   (coherente, |·|≤1)
 state    : 1.000  -0.762   0.271   0.053  -0.030   (decae / desfasa)
```

Esto es **exactamente** la afirmación previa, ahora con forma cerrada: misma física de
fondo, dos preparaciones del conjunto resuelto en energía (capa iso-energía vs ventana
de energía) → mismos resultados cuando el campo no depende de `k`, distintos cuando sí.

> En 1D el contraste es especialmente nítido porque la capa iso-energía tiene solo
> `±k₀`. En 2D (grafeno/valle-Zeeman) la capa tiene una curva con un rango de `ω`, y
> entonces **también** la density-projection desfasa (tesis Fig. 5.6); la state-projection
> añade encima el ensanchamiento por `η` y, por la inyección direccional `|ψ₀⟩`, la
> **anisotropía** zig-zag/armchair de Cummings (Fig. 3). Esa anisotropía es un knob
> exclusivo del flavor state.

---

## 4. Cómo se reflejan los dos flavors en los momentos

Ambos flavors comparten el **mismo kernel caliente** (lo derivamos en la discusión
previa):

$$ \langle v\,|\;\hat s_H(t)\,T_m(\tilde H)\;|\,v\rangle,\qquad
   \hat s_H(t)=\hat U^\dagger(t)\,\hat s\,\hat U(t). $$

Lo único que cambia en los momentos es la frontera:

| | source `v` | filtro de energía | normalización | h.c. |
|---|---|---|---|---|
| **density** | `P\|φ_r⟩` (proyector de espín) | sandwich `T_m` | DOS proyectada en espín `⟨φ_r\|P†T_mP\|φ_r⟩` | implícito (traza real) |
| **state** | `½(𝟙+ĵ·ŝ)\|φ_r⟩` (textura de espín) | multiplica + h.c. | DOS pesada por el estado `⟨ψ₀\|T_m\|ψ₀⟩` | explícito (estado único) |

Para la cadena de canje, como el bloque de espín es `k`-independiente, ambas
normalizaciones y ambos sources dan el mismo cociente → la igualdad de 3.2. El `+h.c.`
de Aron es un artefacto de evaluar en **un** estado; en traza estocástica (LinQT) la
parte imaginaria promedia a cero, así que el flag de simetrización debe ir atado a
"single-state vs traza", **no** al flavor.

---

## 5. Especificación de tests (golden masters)

Todos contra `chain_reference.py`. Tolerancias indicativas para `N ~ 10⁴–10⁶`,
`N_r` pocos.

**Carga (sin espín):**

1. **DOS-moments**: `Cbar^DOS_m ≈ δ_{m,0}`; `|C_m| < tol` para `m ≥ 1` (tol ~ ruido
   estocástico `~1/√(N N_r)`).
2. **DOS reconstruida**: comparar con `1/(π√(4γ²-E²))` en una malla de `E` lejos de los
   bordes; error relativo `< 1%`.
3. **KG-matrix**: `C^KG_{mn}` vs forma cerrada B.15 (matriz dispersa); patrón exacto de
   ceros + valores `{1/2, 1, -1/2}·γ²a²/ħ²`. **Fija la convención de normalización aquí.**
4. **VAC↔MSD**: `σ(E,t)` por VAC y por MSD deben coincidir (equivalencia de
   representación) y crecer lineal en `t` con pendiente `(e²a²/πħ²)√(4γ²-E²)`.
5. **Conductancia balística**: plateau plano en energía; valor = cuanto de conductancia.

**Espín:**

6. **Precesión (equivalencia / consolidación)**: ambos flavors → `cos(2J_ex t/ħ) x̂ -
   sin(2J_ex t/ħ) ŷ`, `S_z = 0`, frecuencia `J_ex/(πħ)`, **independiente de `E_F`**.
   Igualdad flavor-a-flavor dentro del ruido.
7. **Discriminación de flavors** (modelo `B(k)=J_ex+δcos(ka)`): density permanece
   coherente con `ω(E_F)=2(J_ex - δE_F/2γ)/ħ`; state desfasa. El test debe **afirmar**
   que (a) coinciden en `δ=0` y (b) difieren de forma controlada en `δ≠0`. Si esto se
   cumple, los dos flavors capturan lo que dijimos → **consolidar**.

---

## 6. Verificación simbólica (reproducibilidad)

- `derive_core.py` — DOS-moments (`= δ_{m0}`), reconstrucción de DOS, matriz KG vs B.15.
- `derive_spin.py` — precesión cerrada, igualdad de los dos flavors (canje),
  modelo `k`-dependiente (density coherente / state desfasa).
- `chain_reference.py` — formas cerradas como funciones + autotest; úsese como golden
  master en la suite de tests.

Resultado de la verificación: todas las identidades de arriba dan `True` / coincidencia
exacta en sympy. Pendiente de tu lectura: confirmar que la convención de normalización
del punto 3 (lista 5) es la que computa LinQT en la ruta de reconstrucción, antes de
congelar el test.
