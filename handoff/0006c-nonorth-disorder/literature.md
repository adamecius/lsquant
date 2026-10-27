# On-site disorder and local potentials in a non-orthogonal LCAO basis

Literature study for LSQUANT, stage 6c. Studied 23 October 2026. The basis and overlap matrix are fixed throughout; only the potential matrix changes. Equations below use **potential energy**, in energy units: an electric scalar potential Φ contributes −eΦ for an electron.

## 1. Short answer and recommendation

**“On-site potential” has no unique, basis-independent meaning once the atomic orbitals overlap.** One must specify the operator, its matrix representation, and what “site” denotes. Three defensible definitions address different questions:

1. **A prescribed real-space scalar potential:** project it into the unchanged basis,
   \[
   V_{ij}=\langle\phi_i|V(\mathbf r)|\phi_j\rangle.
   \]
   This is the reference definition for electrostatics. Generally it includes off-diagonal elements. SIESTA evaluates local-potential matrix elements by grid integration [R8, §§IV–VI]. An impurity potential centred on a third atom produces three-centre integrals [R29, §2.1].
2. **Atom-labelled energy shifts with an inexpensive sparse construction:** use
   \[
   D_{ij}=w_{A(i)}\delta_{ij},\qquad
   V^{\rm M}=\tfrac12(DS+SD),\quad
   V^{\rm M}_{ij}=\tfrac12(w_{A(i)}+w_{A(j)})S_{ij}.
   \]
   This is a Mulliken population coupling, and a familiar approximation to scalar-potential matrix elements. Explicit examples are DFTB external fields [R16, §VI, Eq. (66)] and self-consistent extended-Hückel transport [R17, §II, Eq. (4)]. It is Hermitian for real shifts, gives exactly \(wS\) for a uniform shift, and preserves the numerical sparsity pattern of S.
3. **Anderson disorder defined in a particular orthonormal orbital basis:** choose that basis explicitly. In the symmetric Löwdin basis \(\chi=\phi S^{-1/2}\), impose \(V^\chi=D\). In the original LCAO representation this means
   \[
   V^{\rm L}=S^{1/2}DS^{1/2}.
   \]
   This is exact for that model, with Hermiticity, the correct constant shift and positive projectors. It generally introduces longer-range matrix elements in the original basis. Löwdin orthogonalization does **not**, by itself, imply that a physical scalar potential is diagonal in the new basis [R2, R11].

**Recommended default for LSQUANT's fixed-SIESTA-basis, atom-labelled disorder: \(V^{\rm M}=(DS+SD)/2\), documented as a Mulliken site-shift model.** Draw one real \(w_A\) per atom and apply it to all radial, angular and spin components on that atom, unless orbital-dependent disorder is explicitly intended. This choice combines a documented electrostatic convention with a sparse matrix action and the correct energy origin. It is a modelling recommendation, not a claim that the literature uniquely defines Anderson disorder this way.

Expose Löwdin-diagonal disorder as a separate model when the objective is the canonical independent-site Anderson ensemble in that orthonormal basis. For a specified spatial field or screened impurity, prefer real-space integrals. Keep raw diagonal LCAO matrix disorder available only under an explicit definition such as “random diagonal Hamiltonian-matrix perturbation.” **Do not add a one-sided \(DS\) to the covariant Hamiltonian as a general real scalar potential.**

The recommended sparse default has a substantive limitation: Mulliken site operators need not be positive. It is therefore unsuitable as an exact representation of every positive local barrier. Section 4 gives a small counterexample.

### Property comparison

Assume linearly independent orbitals, \(S=S^\dagger>0\), real D, and \(H_{ij}=\langle\phi_i|\hat H|\phi_j\rangle\). “Uniform” means setting **all** site values to the same w, including any regions to which a global shift is meant to apply.

| Added covariant matrix V | Hermitian? | Uniform w gives wS? | S → I | Meaning and locality | Compatible with a closed, Hermitian generalized eigenproblem? |
|---|---|---|---|---|---|
| D | Yes | No, generally wI | D | Diagonal matrix disorder; its reconstructed operator uses dual orbitals and need not be spatially local | Yes, as a different Hermitian model |
| (DS + SD)/2 | Yes | Yes | D | Couples to Mulliken populations; numerical support contained in S; not a positive partition | Yes |
| DS, or SD | Only if [D,S] = 0 | Yes | D | One-sided oblique assignment; representation conventions must be checked | Generally no, if used as covariant V |
| S^(1/2) D S^(1/2) | Yes | Yes | D | Diagonal disorder in Löwdin orbitals; generally longer range in LCAO indices | Yes |
| ⟨φᵢ\|V(r)\|φⱼ⟩ | Yes for real V(r) | Yes | General Hermitian matrix, **not necessarily D** | Real-space locality; nonzero only where orbital supports and the potential intersect | Yes |

The orthogonal-limit qualification in the last row matters: orthogonality removes the unweighted overlap integral; it does not make every weighted integral vanish. Recovering the textbook diagonal Anderson model requires an additional site-model approximation.

## 2. What the matrices represent

### Covariant elements, coefficient action and the identity

For \(|\psi\rangle=\sum_i|\phi_i\rangle c_i\),
\[
\langle\psi|\psi\rangle=c^\dagger Sc,\qquad
(H+V)c=ESc.
\]
The matrix of the identity operator is S, so an exact constant potential energy w gives
\[
H\longmapsto H+wS,\qquad E_n\longmapsto E_n+w.
\]
This is an algebraic requirement for a constant scalar potential, not a convention peculiar to DFTB. It also follows from invariance of \(ES-H\) when the energy origin and H shift together.

The coefficient-action matrix is \(A=S^{-1}(H+V)\). It usually is not Hermitian in the Euclidean inner product, but satisfies
\[
A^\dagger S=SA
\]
when the covariant H+V is Hermitian. Confusing these two matrices explains some apparent “non-Hermitian Hamiltonians” in non-orthogonal formulations [R5, R21, R22, R23]. A non-Hermitian coefficient-action representation can be correct; a non-Hermitian covariant potential is a different issue.

If a paper calls D the diagonal **mixed** coefficient-action representation, then its covariant matrix is SD. For real D this is self-adjoint only when DS=SD. Conversely, DS can occur in a representation acting on a different kind of coefficient. A printed one-sided product is therefore interpretable only after checking the indices and metric. It is not evidence that adding DS to SIESTA's H is generally legitimate.

There is a simple obstruction to demanding that every overlapping original orbital be an eigenstate with its own real energy shift: Hermiticity would imply
\[
(w_i-w_j)S_{ij}=0.
\]
Overlapping eigenstates of a Hermitian operator cannot carry different eigenvalues. Thus “shift each non-orthogonal orbital independently” requires a further definition.

### Operator definitions using dual orbitals

Define \(|\phi^i\rangle=\sum_j|\phi_j\rangle(S^{-1})_{ji}\), so \(\langle\phi^i|\phi_j\rangle=\delta_{ij}\). Any covariant V reconstructs the operator in the retained space as
\[
\hat V=\sum_{ij}|\phi^i\rangle V_{ij}\langle\phi^j|.
\]
Consequently:

* Covariant **V=D** means \(\hat V=\sum_iw_i|\phi^i\rangle\langle\phi^i|\). It is Hermitian, but a constant collection of wᵢ is not the identity. The duals can extend beyond a nominal atom.
* The superficially natural expression \(\sum_iw_i|\phi_i\rangle\langle\phi_i|\) instead has covariant matrix **SDS**. A uniform w gives wS², also not wS in general.
* Oblique projectors \(P_i=|\phi_i\rangle\langle\phi^i|\) sum to the identity but are not Hermitian. Symmetrizing \(\sum_iw_iP_i\) gives precisely \(V^{\rm M}\).
* Löwdin rank-one projectors \(|\chi_i\rangle\langle\chi_i|\) are Hermitian, positive and mutually orthogonal. Weighting them gives \(V^{\rm L}\).

These statements are direct derivations using the representation theory in [R5, R22, R23]; they are not four alternative names for the same operator. For multi-orbital atoms, sum over all i belonging to an atom. Global duals, subspace duals and orthogonal projectors onto an atomic subspace are also distinct constructions [R21]. Orthogonal projectors onto overlapping atomic subspaces need not sum to the identity.

### Why the endpoint average is called Mulliken-like

For a single state the real Mulliken population on orbital i is
\[
q_i^{\rm M}=\operatorname{Re}\{c_i^*(Sc)_i\},\qquad
\sum_iq_i^{\rm M}=c^\dagger Sc.
\]
Defining the external energy as \(E_{\rm ext}=\sum_iw_iq_i^{\rm M}\) gives
\[
E_{\rm ext}=c^\dagger\tfrac12(DS+SD)c.
\]
Thus the symmetrized expression is **exact for coupling to the chosen Mulliken populations**. It is approximate when those populations substitute for a spatial charge distribution or a spatial scalar-potential integral. SCC-DFTB uses related population derivatives for its self-consistent charge correction [R7, R16, R25]. The SCC correction is not itself a random Anderson potential; its algebra supplies a precedent for an externally imposed site-energy coupling.

Koskinen and Mäkinen explicitly motivate the endpoint average by dividing overlap between the two atomic regions and assuming a smoothly varying external field [R16, Eq. (66)]. Their derivation supplies no controlled error bound for independent atom-to-atom white noise.

The connection to Löwdin disorder can be made precise. Write S=I+σ with small operator norm σ. Then
\[
V^{\rm L}=D+\tfrac12(\sigma D+D\sigma)
 +\tfrac14\sigma D\sigma-\tfrac18(\sigma^2D+D\sigma^2)
 +O(\|D\|\|\sigma\|^3).
\]
The first two terms are exactly \(V^{\rm M}\). They agree through first order in overlap, but differ at second order. This is a useful connection to overlap expansions [R11], not a guarantee of accuracy for a poorly conditioned, multi-zeta basis.

### A real-space site model is also possible

Choose real weights \(f_A(\mathbf r)\) forming a partition of unity, \(\sum_Af_A=1\), and set
\[
V(\mathbf r)=\sum_Aw_Af_A(\mathbf r),\qquad
V_{ij}=\sum_Aw_A\int\phi_i^*f_A\phi_j\,d\mathbf r.
\]
This gives a genuine scalar random field, Hermiticity and the exact uniform shift without changing the basis. Nonnegative weights also preserve positivity for nonnegative w_A. However, cell indicators, smooth atom weights and screened impurity profiles define different physical models; H and S alone do not determine their integrals. SIESTA orbital functions/grid information, or separately supplied integrals, are needed.

For \(V=\sum_AV_A(\mathbf r-\mathbf R_A)\), terms with orbital centres i, j and potential centre A all different are three-centre integrals. An arbitrary grid potential need not be decomposed into atom-centred contributions at all. “Three-centre” describes the integration structure, not a new metric rule [R8, R29].

## 3. Literature table

**Evidence codes:** **F** = abstract and relevant full-text sections inspected in an accessible paper/preprint; **A** = abstract/metadata inspected, full text unavailable in this run. **B** marks a background or representation reference rather than a paper proposing an Anderson-disorder prescription. Formulae described as “implies” are deductions in the notation of this note, not a claim that the paper prints that expression. Each DOI links to the published work; openly accessible versions used for the central readings are listed after the table.

| Ref. | Reference and DOI | Year | Method / evidence | Form of V, or relevant representation | One-line assessment |
|---|---|---:|---|---|---|
| R1 | Anderson, *Absence of Diffusion in Certain Random Lattices*. [10.1103/PhysRev.109.1492](https://doi.org/10.1103/PhysRev.109.1492) | 1958 | Anderson lattice; A, B | Random local energies in the site model | Establishes the canonical disorder idea; does not select a non-orthogonal LCAO prescription. |
| R2 | Löwdin, *On the Non-Orthogonality Problem Connected with the Use of Atomic Wave Functions in the Theory of Molecules and Crystals*. [10.1063/1.1747632](https://doi.org/10.1063/1.1747632) | 1950 | Symmetric orthogonalization; A, B | χ=φS^(−1/2); implies V=S^(1/2)DS^(1/2) if Vχ=D | Defines the basis transformation; diagonality of disorder is an additional model choice. |
| R3 | Mulliken, *Electronic Population Analysis on LCAO–MO Molecular Wave Functions. I*. [10.1063/1.1740588](https://doi.org/10.1063/1.1740588) | 1955 | Population analysis; A, B | Coupling to gross populations implies (DS+SD)/2 | The source of equal overlap sharing, rather than a direct Anderson-disorder paper. |
| R4 | Hoffmann, *An Extended Hückel Theory. I. Hydrocarbons*. [10.1063/1.1734456](https://doi.org/10.1063/1.1734456) | 1963 | Extended Hückel; A, B | Non-orthogonal empirical H; no external V prescription inferred here | Its overlap-weighted hopping parametrization must not be mistaken for the electrostatic correction. |
| R5 | Artacho and Miláns del Bosch, *Nonorthogonal basis sets in quantum mechanics: Representations and second quantization*. [10.1103/PhysRevA.43.5770](https://doi.org/10.1103/PhysRevA.43.5770) | 1991 | Tensorial formalism; A, B | Different covariant/mixed operator components | Supports checking which representation a purported diagonal operator uses. |
| R6 | Nunes and Vanderbilt, *Generalization of the density-matrix method to a nonorthogonal basis*. [10.1103/PhysRevB.50.17611](https://doi.org/10.1103/PhysRevB.50.17611) | 1994 | Linear-scaling density matrix; F, B | Metric-dependent density/energy representation; no random V specified | Reinforces the distinction between density kernels and covariant Hamiltonian elements. |
| R7 | Elstner et al., *Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties*. [10.1103/PhysRevB.58.7260](https://doi.org/10.1103/PhysRevB.58.7260) | 1998 | SCC-DFTB; A; equations checked in R16/R25 | Sᵢⱼ times the mean SCC site potential | Foundational self-consistent Mulliken correction; not a unique definition of random disorder. |
| R8 | Soler et al., *The SIESTA method for ab initio order-N materials simulation*. [10.1088/0953-8984/14/11/302](https://doi.org/10.1088/0953-8984/14/11/302) | 2002 | Numerical LCAO DFT; F | Grid quadrature of φᵢ*V(r)φⱼ | Direct reference for a local scalar potential in the actual SIESTA basis. |
| R9 | Brandbyge et al., *Density-functional method for nonequilibrium electron transport*. [10.1103/PhysRevB.65.165401](https://doi.org/10.1103/PhysRevB.65.165401) | 2002 | TranSIESTA DFT–NEGF; F | Projected self-consistent potential; G uses ES−H | Solves electrostatics in a non-orthogonal basis; does not prescribe universal diagonal Anderson noise. |
| R10 | Taylor, Guo and Wang, *Ab initio modeling of quantum transport properties of molecular electronic devices*. [10.1103/PhysRevB.63.245407](https://doi.org/10.1103/PhysRevB.63.245407) | 2001 | DFT–NEGF under bias; A | Self-consistent DFT potential; exact element formula not checked here | Context for local electrostatics and open boundaries, not evidence for one chosen disorder formula. |
| R11 | Sandu, *Optical matrix elements in tight-binding models with overlap*. [10.1103/PhysRevB.72.125105](https://doi.org/10.1103/PhysRevB.72.125105) | 2005 | Overlap/operator transformations; F | Aχ=S^(−1/2)AS^(−1/2); inverse overlap expansion | Explains longer-range and intra-atomic terms generated by orthogonalization. |
| R12 | Markussen et al., *Electronic transport through Si nanowires: Role of bulk and surface disorder*. [10.1103/PhysRevB.74.245313](https://doi.org/10.1103/PhysRevB.74.245313) | 2006 | SIESTA, Kubo and recursive GF; F | Hᵢᵢ→Hᵢᵢ+δᵢ, Eq. (11); Kubo uses Löwdin H, Eq. (8) | Direct Anderson-noise precedent, but the basis in which the perturbation is applied must accompany any reuse. |
| R13 | Thygesen, *Electron transport through an interacting region: The case of a nonorthogonal basis set*. [10.1103/PhysRevB.73.035309](https://doi.org/10.1103/PhysRevB.73.035309) | 2006 | Dual-basis NEGF; A, B | Dual-basis treatment of local interactions/perturbation series | Shows why locality and lead/device partitions require metric-aware definitions. |
| R14 | Weiße et al., *The kernel polynomial method*. [10.1103/RevModPhys.78.275](https://doi.org/10.1103/RevModPhys.78.275) | 2006 | Chebyshev/KPM review; F, B | Standard Hermitian lattice operators; disorder examples | Supplies KPM methodology, not a unique extension of on-site disorder to S≠I. |
| R15 | Röder et al., *Kernel polynomial method for a nonorthogonal electronic-structure calculation of amorphous diamond*. [10.1103/PhysRevB.55.15382](https://doi.org/10.1103/PhysRevB.55.15382) | 1997 | Non-orthogonal KPM; A, B | Applies S^(−1)H using iterative solves | Direct precedent for metric-consistent KPM without a dense orthogonalization. |
| R16 | Koskinen and Mäkinen, *Density-functional tight-binding for beginners*. [10.1016/j.commatsci.2009.07.013](https://doi.org/10.1016/j.commatsci.2009.07.013) | 2009 | DFTB derivation; F | External Vᵢⱼ≈(Vᵢ+Vⱼ)Sᵢⱼ/2, Eq. (66) | The clearest explicit scalar-potential justification and its smooth-field assumption. |
| R17 | Stokbro et al., *Semiempirical model for nanoscale device simulations*. [10.1103/PhysRevB.82.075420](https://doi.org/10.1103/PhysRevB.82.075420) | 2010 | Self-consistent extended Hückel + NEGF; F | δHᵢⱼ=[δV_H(Rᵢ)+δV_H(Rⱼ)]Sᵢⱼ/2 | Explicit transport implementation of the endpoint-averaged electrostatic correction. |
| R18 | Zhou and Ceperley, *Construction of localized wave functions for a disordered optical lattice and analysis of the resulting Hubbard model parameters*. [10.1103/PhysRevA.81.013402](https://doi.org/10.1103/PhysRevA.81.013402) | 2010 | Continuum-to-lattice projection; F | Projected continuum H, then Löwdin orthogonalization | Shows distributions of both site and hopping parameters; its disorder-dependent basis differs from LSQUANT's fixed-basis assumption. |
| R19 | Gaus, Cui and Elstner, *DFTB3: Extension of the Self-Consistent-Charge Density-Functional Tight-Binding Method*. [10.1021/ct100684s](https://doi.org/10.1021/ct100684s) | 2011 | Third-order DFTB; A; context checked in R25 | Further charge-dependent Hamiltonian corrections | Higher energy expansion order does not eliminate the underlying population-model choice. |
| R20 | Giese and York, *Density-functional expansion methods: Generalization of the auxiliary basis*. [10.1063/1.3587052](https://doi.org/10.1063/1.3587052) | 2011 | Auxiliary product-density expansion; A | AO products represented by atom-centred auxiliary functions | Systematic route beyond spherical Mulliken charge fluctuations; detailed implementation not read here. |
| R21 | O’Regan, Payne and Mostofi, *Subspace representations in ab initio methods for strongly correlated systems*. [10.1103/PhysRevB.83.245124](https://doi.org/10.1103/PhysRevB.83.245124) | 2011 | Non-orthogonal Hubbard projectors; F, B | Subspace metric and Hermitian local potential construction | Highlights global-dual versus subspace-dual locality; not an endorsement of a unique Anderson model. |
| R22 | Soriano and Palacios, *Theory of projections with nonorthogonal basis sets: Partitioning techniques and effective Hamiltonians*. [10.1103/PhysRevB.90.075128](https://doi.org/10.1103/PhysRevB.90.075128) | 2014 | Projection/NEGF theory; F, B | Direct, dual and oblique projections | Establishes ambiguities in local charge/LDOS and operator partitioning. |
| R23 | Artacho and O’Regan, *Quantum mechanics in an evolving Hilbert space*. [10.1103/PhysRevB.95.115155](https://doi.org/10.1103/PhysRevB.95.115155) | 2017 | Covariant and mixed representations; F, B | Hc=ESc versus S^(−1)H action | Useful representation reference; changing-basis connection terms are irrelevant to the fixed basis here. |
| R24 | Boykin, Sarangapani and Klimeck, *Non-orthogonal tight-binding models: Problems and possible remedies for realistic nano-scale devices*. [10.1063/1.5056178](https://doi.org/10.1063/1.5056178) | 2019 | Position/charge/gauge analysis; A | Transformed position operators; first-order orthogonalization | Warns against inferring field couplings from fitted H,S alone; no detailed formula attributed from the inaccessible full text. |
| R25 | Hourahine et al., *DFTB+, a software package for efficient approximate density functional theory based atomistic simulations*. [10.1063/1.5143190](https://doi.org/10.1063/1.5143190) | 2020 | DFTB software/methods; F | Overlap-weighted charge corrections, including DFTB3, §II | Documents modern use and extensions of the Mulliken-based Hamiltonian. |
| R26 | Fan et al., *Linear scaling quantum transport methodologies*. [10.1016/j.physrep.2020.12.001](https://doi.org/10.1016/j.physrep.2020.12.001) | 2021 | Kubo/LSQT review; F | Orthogonal Anderson model with uniform site noise, Eq. (114) | Provides disorder and transport context; its usual site model does not settle the LCAO convention. |
| R27 | Martínez-Carracedo et al., *Relativistic magnetic interactions from nonorthogonal basis sets*. [10.1103/PhysRevB.108.214418](https://doi.org/10.1103/PhysRevB.108.214418) | 2023 | SIESTA local-operator partitions; F | Hermitian local block plus half intersite blocks, Eq. (13) | Related half-sharing construction satisfies sum rules; magnetic local operators are not automatically scalar disorder. |
| R28 | Vuong et al., *Multipole Expansion of Atomic Electron Density Fluctuation Interactions in the Density-Functional Tight-Binding Method*. [10.1021/acs.jctc.3c00778](https://doi.org/10.1021/acs.jctc.3c00778) | 2023 | Multipole-extended DFTB; F | Overlap-weighted monopole plus dipole/quadrupole corrections, Eq. (33) | Explicit Hamiltonian terms beyond scalar monopole site shifts; not a disorder-definition paper. |
| R29 | Van den Bossche, *Three-Center Tight-Binding Together with Multipolar Auxiliary Functions*. [10.1021/acs.jctc.4c00018](https://doi.org/10.1021/acs.jctc.4c00018) | 2024 | Three-centre TB + auxiliary products; F | ⟨φᵢ\|V_A\|φⱼ⟩ and multipolar response mappings | Explicit distinction between neglected multicentre potentials and simplified Mulliken response. |
| R30 | Santiago-García, Méndez-Sánchez and Sadurní, *Stability of transport and quantum phase transitions with nonorthogonal orbitals*. [10.1103/PhysRevB.111.035429](https://doi.org/10.1103/PhysRevB.111.035429) | 2025 | Metric-corrected tight-binding/Anderson models; A | Metric-dependent renormalization; detailed V convention unverified | Very relevant recent work, but its abstract cannot establish which of the candidate matrices it uses. |
| R31 | Cohan and Weissmann, *Numerical study of systems with diagonal and off-diagonal disorder in a non-orthogonal basis*. [10.1088/0022-3719/12/10/015](https://doi.org/10.1088/0022-3719/12/10/015) | 1979 | Iterative extended Hückel for structurally disordered systems; A | Configuration-dependent diagonal/off-diagonal H and charge redistribution; exact correction unverified | Early direct non-orthogonal disorder precedent, but structural disorder is not the fixed-basis random-potential problem. |

R26 is dated 2021 by the journal issue; OpenAlex returns publication_year=2020, corresponding to its online publication metadata. No later-year inference is made from preprint compilation dates.

### Relevant full-text passages and open versions

These provide the shortest useful reading route. Section and equation identifiers refer to the inspected version; journal pagination can differ.

* **R16:** [DFTB for beginners, open manuscript](https://arxiv.org/pdf/0910.5861), §II.E, Eqs. (36)–(42), for Mulliken populations and their Hamiltonian derivative; **§VI, Eqs. (66)–(67)**, for the external-potential integral and endpoint approximation. Also available from the [authors’ repository](https://jyx.jyu.fi/bitstreams/52379a6c-3df6-4959-ba5a-1b2582bbc7bb/download).
* **R17:** [open manuscript](https://arxiv.org/abs/1004.2812), §II, **Eq. (4)**, explicitly separates the empirical extended-Hückel hopping from the symmetrized Hartree correction; §III treats gates and bias through Poisson electrostatics.
* **R8:** [correct SIESTA manuscript](https://arxiv.org/abs/cond-mat/0111138), §IV, Eq. (16), and **§VI, Grid integrals**. The latter explicitly accumulates weighted orbital products into H. The OpenAlex location cond-mat/0104182 downloaded in the initial extension is a different paper and was excluded from the substantive reading.
* **R9:** [TranSIESTA manuscript](https://arxiv.org/abs/cond-mat/0110650), §III.A on non-orthogonal basis/Green functions and §IV on the effective potential and biased electrodes. This is the real-space DFT branch, not a proof that all NEGF codes use Mulliken potentials.
* **R12:** [open manuscript](https://arxiv.org/abs/cond-mat/0606600) or [published PDF at DTU](https://backend.orbit.dtu.dk/ws/files/4774216/Brand.pdf), **§II.C, Eq. (8)**, Löwdin transformation for Kubo, and **§III.A, Eq. (11)**, diagonal random noise. Eq. (11) uses an unqualified H while the two numerical branches have different basis treatments; it should not be used alone to infer a unique covariant perturbation for LSQUANT.
* **R11:** [open manuscript](https://arxiv.org/abs/cond-mat/0507204), **Eqs. (8)–(12)**, operator transformation and its overlap expansion; discussion immediately following on intra-atomic position elements.
* **R21:** [open manuscript](https://arxiv.org/abs/1102.1920), §II, especially the full/on-site, dual and tensorial representations. Read for why a “local projector” must specify its metric and subspace.
* **R22:** [open manuscript](https://arxiv.org/abs/1404.2043), **§§II–V**, operator reconstruction, dual basis, Hermitian versus oblique projections; §VII for the projected-charge ambiguity.
* **R23:** [open manuscript](https://arxiv.org/abs/1608.05300), §II.A, contrasts the natural mixed representation with the usual generalized matrix equation and their different Hermiticity conditions.
* **R27:** [open manuscript](https://arxiv.org/abs/2309.02558), **§II.A–B, Eqs. (7)–(14)**, pseudo-projectors, operator reconstruction, and Hermitian local partitions with a sum rule. This is corroborating operator algebra, not a scalar-potential prescription.
* **R25:** [open published PDF](https://publikationen.bibliothek.kit.edu/1000118836/71417894), **§II**, energy expansion, Mulliken charges and overlap-weighted Hamiltonian terms; also discusses on-site corrections beyond the elementary approximation.
* **R29:** [open preprint](https://www.cambridge.org/engage/chemrxiv/article-details/6570eb9a29a13c4d47ce6c93), **§2.1.1**, third-centre contributions to the zeroth-order H; **§2.2**, auxiliary AO-product mappings and comparison with monopole/Mulliken approximations. This is the 2023 preprint of the 2024 article.
* **R28:** [open accepted manuscript](https://www.osti.gov/servlets/purl/2251622), **§2.3, Eqs. (25)–(32)**, charge, dipole and quadrupole assignments and integrals; **§2.4, Eq. (33)**, the corresponding Hamiltonian contributions.
* **R18:** [open manuscript](https://arxiv.org/abs/0907.5053), §II.C, **Eq. (19)**, Löwdin construction; §III for joint distributions of site and hopping parameters. Its adapted orbitals are not proposed as a stage-6c requirement.
* **R26:** [open manuscript](https://arxiv.org/abs/1811.07387), numerical examples around **Eq. (114)** and the disorder models section. **R14:** [KPM review](https://arxiv.org/abs/cond-mat/0504627), Chebyshev rescaling and random-vector trace machinery. These are algorithmic context rather than evidence for a particular non-orthogonal site operator.

## 4. What is approximate, disputed, or easy to misread

### The constant-shift test is necessary but does not determine a unique model

Mulliken, Löwdin and exact scalar projection all pass it. There are many other partitions of the identity that also pass. A raw covariant diagonal D fails it, yet remains a mathematically legitimate Hermitian perturbation. The failure disqualifies interpreting all equal diagonal entries as a **uniform scalar energy shift**; it does not invalidate every phenomenological diagonal-matrix disorder calculation.

For illustration take a normalized two-orbital basis
\[
S=\begin{pmatrix}1&s\\s&1\end{pmatrix},\quad |s|<1.
\]
Adding wI to a zero Hamiltonian produces generalized energies w/(1+s) and w/(1−s), rather than two copies of w. Adding wS produces the correct rigid shift. This distinction survives even if the random distribution is chosen to have zero mean.

### Mulliken locality is not positivity or exact spatial locality

For D=diag(w,0),
\[
V^{\rm M}=\begin{pmatrix}w&ws/2\\ws/2&0\end{pmatrix},\qquad
\det V^{\rm M}=-w^2s^2/4.
\]
For w>0 and nonzero s this has a negative direction, also after transformation to an orthonormal basis. With s=0.3 and w=1, the generalized eigenvalues are approximately −0.02414 and 1.02414. Löwdin-diagonal disorder instead has eigenvalues 0 and 1. These are direct algebraic/numerical checks, not a reported benchmark from a cited paper.

Consequently, a Mulliken model is not guaranteed to preserve the pointwise bounds of a physical barrier. For the Löwdin model,
\[
w_{\min}S\preceq V^{\rm L}\preceq w_{\max}S.
\]
The same bound holds for exact projection of a spatial potential bounded between those values. No such general bound holds for the symmetrized model. Hermitian half-sharing and relevant sum rules [R27] do not imply positive local projectors.

There are also three different locality notions: support of matrix entries; spatial support of an operator kernel; and support in a selected orbital subspace. Global inverse overlaps and square roots generally spread orbital components. Conversely, exact scalar-potential integrals vanish for disjoint compact orbital supports, but can be nonzero where Sᵢⱼ is zero by symmetry or cancellation. The underlying multiplication operator V(r) is local; its finite-space reconstruction PVP need not have a local continuum kernel. A uniform electric field can mix different angular orbitals on one atom even when their overlap vanishes. The endpoint approximation misses such dipole/polarization matrix elements. Auxiliary-product and multipole developments address this class of limitations [R20, R28, R29].

### “Diagonal after orthogonalization” specifies an ensemble

For a fixed V, transforming H and V consistently is just a change of representation and preserves all predictions. Defining fresh independent random diagonal entries **after** changing basis specifies a new random operator ensemble. Likewise, the same numerical disorder width W in raw-diagonal, Mulliken and Löwdin models need not produce the same scattering rates or localization lengths. Agreement to first order in small overlap does not make the ensembles identical.

The matrix square roots in the Löwdin prescription are not generally sparse. Decay depends on basis localization, overlap conditioning and the system; a small-overlap expansion or truncation needs an error assessment [R11, R24]. The Löwdin basis is a particular orthonormal basis, not automatically a maximally localized Wannier basis. A different orthonormalization followed by independent diagonal draws generally defines different disorder.

### The meaning of an atom in a multi-zeta SIESTA basis must be explicit

Independent wᵢ on each basis function changes the disorder model when the number of zeta/polarization functions changes. It can also introduce random orbital or spin splittings when only a scalar atom potential was intended. One w_A for every orbital and both spin components of atom A is a clear default. Within an atom, different radial functions need not be orthogonal: even that atom's constant shift contains **w_A S_AA**, not necessarily a diagonal block.

Atomwise scalar D commutes with changes of orbital coordinates confined within each atomic block, making the Mulliken construction covariant under such regroupings. This does not remove dependence on the choice of the atomic subspaces or on basis enlargement. A code should record the random distribution, whether W is a full width or a standard deviation, the atom/orbital assignment and the chosen operator convention. For uniform [−W/2,W/2], the variance is W²/12.

### Impurity chemistry, screening, and basis changes are separate issues

An imposed random potential, a screened charged impurity, a substitutional atom and a vacancy need not have the same matrix structure. A prescribed scalar potential can be projected in a fixed basis, including all induced off-diagonal elements. Changing those matrix elements is consistent with “potential disorder only”: it does not mean the basis or kinetic operator changed. A true substitution or vacancy may require chemical, structural or screening information beyond the present model [R12]. DFTB self-consistency also changes the potential through the density; externally imposed Anderson noise by itself supplies no screening.

The continuum-to-lattice work [R18] usefully illustrates correlated diagonal/off-diagonal parameters, but its disorder-adapted orbitals also change the basis. Historical extended-Hückel disorder work [R31] likewise concerns structurally disordered configurations and their charge redistribution. These studies should not be read as prescriptions for fixed-S potential noise. For the stage-6c question, S remains exactly fixed.

### NEGF and gauge statements need their scope stated

The endpoint average is explicitly used in some semiempirical NEGF models [R17]; NEGF itself imposes no such prescription. TranSIESTA instead obtains local matrix elements from the self-consistent DFT potential [R8, R9]. Global energy-origin invariance requires shifting the Hamiltonian, energy arguments and reservoir chemical potentials consistently. Shifting only a device region is a physical potential step, whose boundary elements must follow the selected potential model. Passing the constant scalar-shift test is weaker than satisfying every electromagnetic gauge/position-operator constraint discussed in optical tight-binding work [R11, R24].

## 5. Consequences for a stage-6c implementation

For the recommended model the covariant update is simply
\[
H'_{ij}=H_{ij}+\tfrac12(w_{A(i)}+w_{A(j)})S_{ij},\qquad S'=S.
\]
An implementation can traverse the overlap entries or use the matrix action
\[
V^{\rm M}x=\tfrac12\{D(Sx)+S(Dx)\}.
\]
The required storage is the union of H and S supports. Entries present in S but absent from the stored H must not be silently discarded. Real atom labels and exact conjugate pairing preserve Hermiticity for complex matrices and spinor calculations.

For KPM, the same physical update must enter whichever representation the existing algorithm uses:
\[
A'=S^{-1}(H+V),\quad\text{or}\quad
h'=S^{-1/2}(H+V)S^{-1/2}.
\]
Adding D directly to the first expression gives a covariant perturbation SD, which is generally non-Hermitian. Adding D directly to the second instead implements the **Löwdin-diagonal model**. Neither is automatically the recommended Mulliken model. Generalized KPM can apply S⁻¹ through solves rather than an explicit inverse [R15]. Likewise, centring and scaling the covariant Hamiltonian for Chebyshev recursion uses (H+V−bS)/a, consistent with scaling the coefficient-action matrix as (A′−bI)/a. Metric-consistent traces, normalizations and current operators remain necessary; this note does not assume LSQUANT's implementation of them.

Appropriate validation checks are: Hermiticity of H′; a uniform w giving H′=H+wS and a rigid generalized spectral shift; the S=I model limit; equality of spectra when the **same** V is transformed; and invariant S across realizations. Recompute or safely bound the KPM spectral range after adding disorder. The Löwdin perturbation is bounded by the range of its wᵢ, whereas a Mulliken perturbation requires an actual bound on S⁻¹/²VS⁻¹/², for example a safe norm estimate involving the conditioning of S. The two-orbital example above shows why merely enlarging the clean range by max|wᵢ| is not always rigorous for Mulliken disorder.

For a physical application, compare at least one small system against a specified smooth-field integral or a Löwdin-disorder calculation, depending on the intended interpretation. This checks the modelling approximation separately from the linear algebra. No code changes or repository operations were performed for this study.

## 6. Search provenance and limits

The first command attempted was `python3 /root/onsite_search.py`, as requested. It failed because pyalex was missing. Pyalex was installed under `/tmp/onsite-python`; the next attempt revealed that `/root/papers` did not exist and its read-only parent prevented creating it in this sandbox. The supplied script was therefore run with its output-directory expression redirected to `/tmp/onsite_nonorth`; the original script was not edited. Its serial download phase was subsequently interrupted and resumed with concurrent fetching, reuse of existing PDFs and shorter network timeouts. The same 24 queries and record construction produced **549 records, 500 indexed abstracts and 200 downloaded PDFs** in `/tmp/onsite_nonorth/corpus.jsonl`; these are retrieval totals, not counts of substantively relevant papers.

The supplied 24 search queries were supplemented using `pyalex.Works` with targeted searches for external scalar potentials in DFTB, non-orthogonal gauge/position operators, local projectors, non-orthogonal KPM, disordered optical-lattice projections and multipole/three-centre methods, plus direct DOI lookups. Backward reference lists and forward citing works were followed for R7, R22, R24 and R30. The extension retrieved 352 distinct works before the final supplementary queries. Retrieved OpenAlex metadata, abstracts, citation relations, download logs and converted open PDFs are held under `/tmp/onsite_nonorth` (`extended_works.jsonl`, `citation_works.jsonl`, `all_extended.jsonl`, `extra_works.jsonl`, and PDF/TXT files).

Abstracts were used for screening; central formulae were checked in the full texts identified above. A retrieved PDF was checked against its title rather than presumed correct from the OpenAlex location. Inaccessible articles are marked A; in particular, this study does not establish the detailed potential convention of the recent R30 paper. Search matches and citation counts were not treated as evidence for an equation. This is a focused literature synthesis, with the default recommendation and the elementary matrix derivations identified as the present analysis, rather than a systematic-review claim of exhaustive coverage.
