"""Symbolic source of the three metric forms (brief sections 3 and 12).

X=S^-1, Y=S^-1/2, Z=S^1/2. Physical matrices O are Hermitian; their
coefficient representation XO is S-Hermitian, not generally Euclidean
Hermitian. The anticommutator is an explicitly labelled approximation for
products and spectral functions. SymPy expressions here state the formulas
implemented numerically in models/nonorthogonal.py and bases.py.
"""
from __future__ import annotations

import sympy as sp

H, S, X, Y, Z, V, A, B, q, dH, dS, dH_cell, dS_cell = sp.symbols(
    "H S X Y Z V A B sigma_z dH dS dH_cell dS_cell", commutative=False)
HERMITIAN = {H, S, X, Y, Z, V, A, B, q, dH, dS, dH_cell, dS_cell}
D = sp.Symbol("D", commutative=False)
HERMITIAN.add(D)


def dagger(expression):
    """Hermitian conjugation: reverse products, conjugate scalars, O^+=O.

    The raw row-centred dipole is excluded: A here is the Hermitian
    midpoint-centred connection returned by siesta.midpoint_connection.
    """
    expression = sp.sympify(expression)
    if expression in HERMITIAN:
        return expression
    if expression.is_commutative:
        return sp.conjugate(expression)
    if expression.is_Add:
        return sp.Add(*(dagger(term) for term in expression.args))
    if expression.is_Mul:
        return sp.Mul(*(dagger(term) for term in reversed(expression.args)))
    if expression.is_Pow and expression.exp.is_Integer:
        return dagger(expression.base)**expression.exp
    return sp.adjoint(expression)


def _reduce_word(word):
    """Reduce adjacent metric factors and commute spin only through the metric."""
    word = list(word)
    rules = {(Y, Y): [X], (Z, Z): [S], (S, X): [], (X, S): [],
             (Y, Z): [], (Z, Y): [], (X, Z): [Y], (Z, X): [Y],
             (S, Y): [Z], (Y, S): [Z]}
    changed = True
    while changed:
        changed = False
        for i in range(len(word)-1):
            pair = tuple(word[i:i+2])
            if pair in rules:
                word[i:i+2] = rules[pair]
                changed = True
                break
            if word[i] == q and word[i+1] in (S, X, Y, Z):
                word[i:i+2] = list(reversed(pair))
                changed = True
                break
    return tuple(word)


def normal_form(expression, *, trace=False):
    """Expand an NC polynomial, reduce YY=X, SX=XS=I and cyclic trace words.

    Z is used solely to express similarity ``f(H_L)=Z f(XH) Y``. The
    rewrite rules also include YZ=ZY=I. No commutation of H, V, A or a
    spectral function with S is assumed. The spin q commutes only with the
    spin-scalar metric, as required by spin_operators.
    """
    expression = sp.expand(expression)
    collected = {}
    for term in sp.Add.make_args(expression):
        coefficient, word = sp.Integer(1), []
        for factor in sp.Mul.make_args(term):
            if factor.is_commutative:
                coefficient *= factor
            elif factor.is_Pow and factor.exp.is_Integer and factor.exp >= 0:
                word.extend([factor.base]*int(factor.exp))
            else:
                word.append(factor)
        reduced = _reduce_word(word)
        if trace and reduced:
            while True:
                rotations = [_reduce_word(reduced[i:]+reduced[:i]) for i in range(len(reduced))]
                best = min(rotations, key=lambda w: (len(w), tuple(map(str, w))))
                if len(best) == len(reduced):
                    reduced = best
                    break
                reduced = best
        collected[reduced] = collected.get(reduced, 0)+coefficient
    return sp.Add(*(sp.simplify(c)*sp.Mul(*word) for word, c in collected.items() if c != 0))


def operator_forms(physical):
    """Return O_c=XO, O_L=YOY, O_sym=(XO+OX)/2 for a physical matrix O."""
    return {"physical": physical, "generalised": X*physical, "lowdin": Y*physical*Y,
            "anticommutator": (X*physical+physical*X)/2}


def hamiltonian_forms():
    """H_c=XH, H_L=YHY, H_sym={X,H}/2; H_sym is not similar to H_c."""
    return operator_forms(H)


def velocity_forms(mode="covariant"):
    """Eq. (4): V=dH-(dS X H+H X dS)/2+i(H X A-A X H).

    bare uses cell derivatives and omits A; berry_connection uses atom
    derivatives and omits A; covariant requires the midpoint dipole A.
    """
    if mode not in ("bare", "berry_connection", "covariant"):
        raise ValueError("unknown velocity mode")
    dh, ds = (dH_cell, dS_cell) if mode == "bare" else (dH, dS)
    physical = dh-(ds*X*H+H*X*ds)/2
    if mode == "covariant":
        physical += sp.I*(H*X*A-A*X*H)
    return operator_forms(physical)


def spin_forms():
    """Sigma_z=S q, Sigma_c=q, Sigma_L=q, Sigma_sym=q when [S,q]=0."""
    return operator_forms(S*q)


def spin_current_forms():
    """J=(V X Sigma+Sigma X V)/2=(V q+q V)/2 for spin-scalar S."""
    return operator_forms((V*q+q*V)/2)


def spin_hall_normalisation():
    """Symbolic source for reconstruct.shc's Pauli-spin Bastin units."""
    bastin = sp.Symbol("KuboBastin", real=True)
    num_cells, cell_area = sp.symbols("N_k A_cell", positive=True)
    return 2*sp.pi*bastin/(num_cells*cell_area)


def hermiticity(forms):
    """Return Euclidean and S-Hermiticity verdicts via symbolic reduction."""
    result = {}
    for name, expression in forms.items():
        result[name] = {
            "hermitian": normal_form(dagger(expression)-expression) == 0,
            "s_hermitian": normal_form(dagger(expression)*S-S*expression) == 0,
        }
    return result


def trace_equivalence():
    """Prove generalised/Lowdin equality for arbitrary catalogue spectral slots.

    H_L=Z (XH) Y and YZ=I imply f(H_L)=Z f(XH) Y by functional calculus
    (polynomial induction, then analytic continuation). Here F_j are
    arbitrary, noncommutative f_j(XH), so this is not a finite-order test.
    Substitute that similarity into the whole word and reduce the cyclic
    polynomial. TD words follow identically; MSD and projectors follow by
    expanding their sums/commutators and applying linearity of the trace.
    """
    result = {"hamiltonian_similarity": normal_form(Z*(X*H)*Y-Y*H*Y) == 0}
    for name, statics, slots in [("dos", 0, 1), ("spectral_operator", 1, 1),
                                  ("greenwood", 2, 2), ("bastin", 2, 2),
                                  ("time_correlation", 2, 3), ("projected_evolution", 3, 3),
                                  ("msd_expanded_term", 2, 3)]:
        operators = sp.symbols(f"O0:{statics}", commutative=False)
        functions = sp.symbols(f"F0:{slots}", commutative=False)
        generalised, lowdin = sp.Integer(1), sp.Integer(1)
        for i in range(max(statics, slots)):
            if i < statics:
                generalised *= X*operators[i]
                lowdin *= Y*operators[i]*Y
            if i < slots:
                generalised *= functions[i]
                lowdin *= Z*functions[i]*Y
        residual = normal_form(lowdin-generalised, trace=True)
        result[name] = {"generalised": generalised, "lowdin_after_similarity": lowdin,
                        "trace_residual": residual, "equal": residual == 0}
    # A single-operator trace is equal. Already with two static operators
    # and no spectral functions the anticommutator is not generally equal.
    single = normal_form(operator_forms(A)["anticommutator"]-X*A, trace=True)
    pair = normal_form(operator_forms(A)["anticommutator"]*operator_forms(B)["anticommutator"]-
                       X*A*X*B, trace=True)
    result["anticommutator"] = {
        "single_operator_equal": single == 0, "pair_residual": pair,
        "generic_products_equal": pair == 0,
        "sufficient_exact_condition": "[S,H]=0 and [S,O_i]=0 for every static slot",
        "classification": "approximation for DOS, Bastin, Greenwood and TD in general",
    }
    return result


def describe():
    """Return printable symbolic forms and verdicts for the one self-check."""
    forms = {"V": velocity_forms(), "Sigma_z": spin_forms(), "J": spin_current_forms()}
    return {name: {"forms": {route: str(normal_form(expr)) for route, expr in item.items()},
                   "hermiticity": hermiticity(item)} for name, item in forms.items()}


def anderson_forms():
    """Derive the symmetric on-site approximation V_ij=(w_i+w_j) S_ij/2.

    Assigning w_j to the ket orbital gives SD; the adjoint assigns w_i to
    the bra and gives DS. Their equal-weight average is the Hermitian
    approximation to a local scalar potential. It is a modelling ansatz,
    not the exact integral of an arbitrary continuum random potential.
    A constant potential must give w*S, since physical identity has matrix
    S. D alone would give w*I and change the generalized spectrum by w*S^-1.
    """
    alpha = sp.Symbol("alpha", real=True)
    trial = alpha*D*S+(1-alpha)*S*D
    # Independent NC words DS and SD have equal real coefficients iff alpha=1/2.
    coefficient = sp.solve(sp.Eq(alpha, 1-alpha), alpha)[0]
    physical = sp.expand(trial.subs(alpha, coefficient))
    w = sp.Symbol("w", real=True)
    forms = operator_forms(physical)
    return dict(forms=forms, hermiticity=hermiticity(forms),
                symmetric_weight=coefficient,
                constant_residual=normal_form(physical.subs(D, w)-w*S),
                orthogonal_residual=sp.expand(physical.subs(S, 1)-D),
                incorrect_diagonal_constant_residual=w-w*S)


def disorder_velocity_forms(mode="covariant"):
    """Insert H+V into the SAME position/velocity ladder as stage 6b.

    The periodic derivative is the derivation delta(M)_ij=i*d_ij*M_ij
    before wrapping hopping translations. D commutes with the orbital
    coordinates, so delta(V)=(D*dS+dS*D)/2. On the infinite covering lattice
    R={position,S}/2+A, and i(H X R-R X H), using the product rule,
    equals delta(H)-{delta(S) X H+H X delta(S)}/2+i(H X A-A X H).
    This is equation (3) with the dipoles, including boundary bonds; a
    commutator with wrapped supercell coordinates would cut those bonds.
    """
    clean = velocity_forms(mode)["physical"]
    ds = dS_cell if mode == "bare" else dS
    potential = anderson_forms()["forms"]["physical"]
    derivative = (D*ds+ds*D)/2
    correction = derivative-(ds*X*potential+potential*X*ds)/2
    if mode == "covariant":
        correction += sp.I*(potential*X*A-A*X*potential)
    substituted = velocity_forms(mode)["physical"].subs(
        {H: H+potential, dH_cell if mode == "bare" else dH:
         (dH_cell if mode == "bare" else dH)+derivative}, simultaneous=True)
    return dict(forms=operator_forms(clean+correction), derivative=derivative,
                correction=correction,
                substitution_residual=normal_form(substituted-clean-correction),
                hermitian=normal_form(dagger(clean+correction)-clean-correction) == 0)


def periodic_position_velocity_proof():
    """Reduce theory.md (3) to (4) plus the midpoint connection.

    P denotes unwrapped orbital coordinates; delta(M)=i[M,P]. Therefore
    PH=HP+i*dH, PS=SP+i*dS, PX=XP-i*X*dS*X. These rules are the product
    rule for the periodic hopping derivative, including inverse overlaps.
    """
    P = sp.Symbol("position", commutative=False)
    R = (P*S+S*P)/2+A
    commutator = sp.I*(H*X*R-R*X*H)
    expanded = normal_form(commutator)
    while True:
        previous = expanded
        expanded = normal_form(expanded.subs(P*H, H*P+sp.I*dH).subs(
            P*S, S*P+sp.I*dS).subs(P*X, X*P-sp.I*X*dS*X))
        if expanded == previous:
            break
    return dict(position=R, commutator=commutator,
                residual=normal_form(expanded-velocity_forms("covariant")["physical"]))
def spectral_projection_forms():
    """Derive stage-6d forms in a complete S-normalised eigenframe.

    C^+ S C=I, C C^+=X, H=S C D C^+ S. W is the diagonal
    0/1 window mask, D_s agrees with D on W and clips its complement.
    O denotes a PHYSICAL matrix; its coefficient action is X O.
    Thus the brief's S P O_c P is P^+ O P, not S P O P.
    No derivative of the modified H is substituted for a projected velocity.
    """
    C, Cd, D, W, Ds, O = sp.symbols("C C_dagger D W D_shift O", commutative=False)
    P = C*W*Cd*S
    Pd = S*C*W*Cd
    hard = S*C*W*D*Cd*S
    shift = S*C*Ds*Cd*S
    projected = Pd*O*P
    return dict(projector=P, projector_adjoint=Pd,
                hard=operator_forms(hard), shift=operator_forms(shift),
                operator=operator_forms(projected),
                shift_update=H+S*C*(Ds-D)*Cd*S,
                identity=Pd*S*P, coefficient_operator=X*O,
                symbols=dict(C=C, Cd=Cd, D=D, W=W, Ds=Ds, O=O))


def spectral_projection_proof():
    """Symbolic identities for arbitrary dimension, not a numeric fixture.

    Reduce NC words with completeness C C^+=X and C^+ S C=I.
    The sole diagonal identity needed for retained eigenpairs is
    D_s W=D W; W D=D W and W^2=W express spectral selection.
    The tester executes these reductions before running numerical checks.
    """
    forms = spectral_projection_forms()
    C, Cd, D, W, Ds, O = (forms["symbols"][n] for n in ("C", "Cd", "D", "W", "Ds", "O"))

    def reduce(expression):
        rules = {(Cd, S, C): [], (C, Cd): [X], (S, X): [],
                 (X, S): [], (W, W): [W], (W, D): [D, W],
                 (Ds, W): [D, W]}
        terms = {}
        for term in sp.Add.make_args(sp.expand(expression)):
            scalar, word = sp.Integer(1), []
            for factor in sp.Mul.make_args(term):
                if factor.is_commutative:
                    scalar *= factor
                elif factor.is_Pow and factor.exp.is_Integer and factor.exp >= 0:
                    word.extend([factor.base]*int(factor.exp))
                else:
                    word.append(factor)
            changed = True
            while changed:
                changed = False
                for i in range(len(word)):
                    for pattern, replacement in rules.items():
                        if tuple(word[i:i+len(pattern)]) == pattern:
                            word[i:i+len(pattern)] = replacement
                            changed = True
                            break
                    if changed:
                        break
            key = tuple(word)
            terms[key] = terms.get(key, 0)+scalar
        return sp.expand(sum(value*sp.Mul(*key) for key, value in terms.items()))

    P, Pd = forms["projector"], forms["projector_adjoint"]
    HP, Hs = forms["hard"]["physical"], forms["shift"]["physical"]
    # Adjoint expressions use C^+=Cd, W^+=W, D^+=D, Ds^+=Ds.
    Oprojected = Pd*O*P
    def projection_dagger(expression):
        return dagger(expression).xreplace({sp.adjoint(C): Cd, sp.adjoint(Cd): C,
                                             sp.adjoint(D): D, sp.adjoint(W): W,
                                             sp.adjoint(Ds): Ds, sp.adjoint(O): O})
    identities = {
        "idempotent": P*P-P,
        "s_self_adjoint": Pd*S-S*P,
        "hard_retained_eigenpairs": HP*C*W-S*C*D*W,
        "shift_retained_eigenpairs": Hs*C*W-S*C*D*W,
        "hard_complement_zero": HP*C*(1-W),
        "physical_operator_equals_coefficient_form": Oprojected-S*P*X*O*P,
        "projected_identity": Pd*S*P-S*P,
        "hard_from_coefficient_hamiltonian": HP-S*P*(C*D*Cd*S)*P,
        "shift_as_update": Hs-(S*C*D*Cd*S+S*C*(Ds-D)*Cd*S),
        "hard_physical_hermitian": HP-projection_dagger(HP),
        "shift_physical_hermitian": Hs-projection_dagger(Hs),
        "projected_physical_hermitian": Oprojected-projection_dagger(Oprojected),
    }
    return {name: reduce(expression) for name, expression in identities.items()}


def orthogonal_projection_forms():
    """Polar reduced-coordinate forms and the Fourier convention of 6d.

    C_W^+ S C_W=I, U unitary => Q=C_W U, Q^+ S Q=I.
    H_W=U^+ diag(E_W) U, O_W=Q^+ O Q. The retained metric is
    identity only in this optional reduced tool; hard/shift retain S.
    The scalar phase and FFT sign below are the source of real_space.py's
    conversion before calling the unchanged origin Model.
    """
    Cw, Cwd, U, Ud, D, O = sp.symbols("C_W C_W_dagger U U_dagger E_W O", commutative=False)
    k, R, tau_i, tau_j = sp.symbols("k R tau_i tau_j", real=True)
    phase = sp.exp(sp.I*k*(R+tau_j-tau_i))
    cell = sp.simplify(sp.exp(sp.I*k*tau_i)*phase*sp.exp(-sp.I*k*tau_j))
    return dict(isometry=Cw*U, overlap=Ud*(Cwd*S*Cw)*U,
                hamiltonian=Ud*D*U, operator=Ud*Cwd*O*Cw*U,
                cell_phase=cell, cell_phase_residual=sp.simplify(cell-sp.exp(sp.I*k*R)),
                inverse_fourier_phase=sp.exp(-sp.I*k*R))
