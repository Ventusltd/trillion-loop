# The second channel

Source: `t02_positron.py`

A verdict that is not trusted unless it is born more than once. This channel
does not re-run the enumeration kernel and does not copy it. It rebuilds the
same two quantities from the formulas, along arithmetic routes that are
algebraically identical and numerically distinct, and it rebuilds the index
decode twice by unrelated methods.

---

## 1. The mathematics

### The space

Six axes, each a finite table of values:

$$\mathcal{X}=\prod_{k=0}^{5}\mathbb{Z}_{r_k},\qquad
\mathcal{T}_k:\mathbb{Z}_{r_k}\to\mathbb{R},\qquad
N=\prod_{k=0}^{5}r_k$$

with $(\mathcal{T}_0,\dots,\mathcal{T}_5)$ the tables for
$(n,\ V,\ \beta,\ \tau,\ I,\ \rho)$ respectively.

### The decode, computed twice

Place values here run the other way: $p_5=1$ and $p_k=\prod_{j>k}r_j$, so axis
five is fastest-varying. Two maps are defined:

$$D_1:\ x\mapsto(d_5,d_4,\dots,d_0)\ \text{by successive division},\qquad
d_k=\Bigl\lfloor\frac{x_k}{1}\Bigr\rfloor\bmod r_k,\ \ x_{k-1}=\Bigl\lfloor\frac{x_k}{r_k}\Bigr\rfloor,\ \ x_5=x$$

$$D_2:\ x\mapsto\bigl(d_k\bigr)_{k=0}^{5},\qquad
d_k=\Bigl\lfloor\frac{x}{p_k}\Bigr\rfloor\bmod r_k .$$

Two conditions are required on every index of the block
$\mathcal{R}=\{x_0,\dots,x_0+w-1\}$:

$$\forall x\in\mathcal{R}:\quad \Bigl\lfloor\frac{x}{\prod_k r_k}\Bigr\rfloor=0
\qquad\text{(no leftover)}$$

$$\forall x\in\mathcal{R}:\quad D_1(x)=D_2(x)
\qquad\text{(the two decodes agree)}$$

together with the containment precondition $x_0+w\le N$, which is checked
here rather than assumed.

### Three arithmetic routes per quantity

With $\tau_0$ the datum and $\varphi$ the design factor, define on
$\mathcal{C}=D(\mathcal{R})$, writing $\Delta\tau=\tau-\tau_0$:

$$\begin{aligned}
\nu_A&=n\,V\,(1+\beta\,\Delta\tau) &
\iota_A&=\varphi\,I\,(1+\rho)\\
\nu_B&=nV+nV\beta\,\Delta\tau &
\iota_B&=\varphi I+\varphi I\rho\\
\nu_C&=\bigl(V+V(\beta\,\Delta\tau)\bigr)\,n &
\iota_C&=\bigl(I(1+\rho)\bigr)\,\varphi
\end{aligned}$$

$\nu_A=\nu_B=\nu_C$ and $\iota_A=\iota_B=\iota_C$ as identities in $\mathbb{R}$:
$B$ applies distributivity, $C$ applies associativity and commutativity of
multiplication. In $\mathbb{F}$, none of the three is required to equal another,
because neither law is exact under rounding.

### The two acceptance conditions

**Magnitude**, with a tolerance. Let $\epsilon$ be the unit roundoff of the
working precision and $\lambda$ the permitted slack in units in the last place.
For each unordered pair $(f,g)$ drawn from $\{\nu_A,\nu_B,\nu_C\}$ and from
$\{\iota_A,\iota_B,\iota_C\}$:

$$M(f,g)\iff \bigl\lvert f-g\bigr\rvert\ \le\ \lambda\,\epsilon\,\max\bigl(\lvert f\rvert,\lvert g\rvert\bigr).$$

**Verdict**, with no tolerance whatever. With the admissible set
$W=[\ell_v^-,\ell_v^+]\times[\ell_i^-,\ell_i^+]$,

$$\omega_X=\bigl[\,(\nu_X,\iota_X)\in W\,\bigr],\qquad X\in\{A,B,C\}$$

$$V\iff \omega_A=\omega_B=\omega_C .$$

### The rule

$$\text{accept}(\mathcal{R})\iff
\underbrace{\forall x:\ D_1=D_2}_{\text{decode}}\ \wedge\
\underbrace{\forall x\ \forall (f,g):\ M(f,g)}_{\text{magnitude}}\ \wedge\
\underbrace{\forall x:\ V}_{\text{verdict}}$$

and the block is atomic: $\neg\text{accept}(\mathcal{R})$ yields **no** partial
result. There is no subset of $\mathcal{R}$ that survives a parting.

### The reported quantities of a block that is accepted

$$\pi=\lvert\{x\in\mathcal{R}:\omega_A(x)\}\rvert,\qquad
\gamma_\nu=\max_{x,(f,g)}\frac{\lvert f-g\rvert}{\max(\lvert f\rvert,\lvert g\rvert)},\qquad
\eta_\nu=\min_{x}\frac{\lvert\nu_A(x)-\ell_v^+\rvert}{\lvert\ell_v^+\rvert}$$

$\gamma$ is the worst relative gap observed between routes; $\eta$ is the closest
any case came to a limit. A block with small $\eta$ is a block that nearly
parted, and the fact is reported even when it did not.

### Why the decode must be checked separately

If $D$ is wrong but shared, then $\nu_A,\nu_B,\nu_C$ are all evaluated at the
same wrong point, and

$$\nu_A\bigl(D'(x)\bigr)=\nu_B\bigl(D'(x)\bigr)=\nu_C\bigl(D'(x)\bigr)$$

holds perfectly for any faulty $D'$. Agreement among routes is invariant under a
shared change of argument, so no number of routes detects a decode fault. Only a
second, differently constructed $D$ can.

---

## 2. The code

| symbol | name in the code |
|---|---|
| $x_0,\ w$ | `base`, `n` — the first index of the block and its length |
| $r_k$ | `pos_radices[k]` |
| $\mathcal{T}_k$ | `pos_tables[k]` |
| $p_k$ | `pos_place[k]` |
| $D_1$ | `pos_dig_a` (successive division) |
| $D_2$ | `pos_dig_b` (place values) |
| $\Delta\tau$ | `pos_dt` |
| $\nu_A,\nu_B,\nu_C$ | `pos_v_ra`, `pos_v_rb`, `pos_v_rc` |
| $\iota_A,\iota_B,\iota_C$ | `pos_i_ra`, `pos_i_rb`, `pos_i_rc` |
| $\epsilon,\ \lambda$ | `pos_eps`, `POS_ULP_SLACK` |
| $M$ | `~pos_bad_mag` |
| $\omega_X$ | `pos_ok_a`, `pos_ok_b`, `pos_ok_c` |
| $\neg V$ | `pos_split` |
| $\pi,\ \gamma,\ \eta$ | `pos_pass_total`, `pos_worst_v_gap`, `pos_closest_v` |

The containment precondition, checked rather than assumed:

```python
pos_space = 1
for pos_r in pos_radices:
    pos_space *= pos_r
if pos_base + pos_n > pos_space:
    return _pos_abort(
        "block [%d, %d) runs past the end of the space (size %d); "
        "a wrapped index decodes to a tuple that is not this block's"
        % (pos_base, pos_base + pos_n, pos_space))
```

The decode, twice, by unrelated methods:

```python
pos_dig_a = [None] * pos_k_axes
pos_rest = pos_idx.copy()
for pos_k in range(pos_k_axes - 1, -1, -1):        # D1: successive divmod
    pos_rad = cp.int64(pos_radices[pos_k])
    pos_dig_a[pos_k] = pos_rest % pos_rad
    pos_rest = pos_rest // pos_rad
if bool((pos_rest != 0).any()):
    return _pos_abort("decode leftover: ...")

for pos_k in range(pos_k_axes):                    # D2: precomputed places
    pos_dig_b = (pos_idx // cp.int64(pos_place[pos_k])) % cp.int64(pos_radices[pos_k])
    if bool((pos_dig_a[pos_k] != pos_dig_b).any()):
        return _pos_abort("decode parting on axis %s ..." % AXES[pos_k])
```

The three routes:

```python
pos_dt   = pos_tmin - POS_T_REF

pos_v_ra = pos_ns * pos_voc * (1.0 + pos_beta * pos_dt)          # factored
pos_v_rb = pos_ns * pos_voc + pos_ns * pos_voc * pos_beta * pos_dt  # distributed
pos_v_rc = (pos_voc + pos_voc * (pos_beta * pos_dt)) * pos_ns    # re-associated

pos_i_ra = POS_I_FACTOR * pos_isc * (1.0 + pos_rear)
pos_i_rb = POS_I_FACTOR * pos_isc + POS_I_FACTOR * pos_isc * pos_rear
pos_i_rc = (pos_isc * (1.0 + pos_rear)) * POS_I_FACTOR
```

Magnitudes with slack; verdicts with none:

```python
pos_scale   = cp.maximum(cp.abs(pos_lhs), cp.abs(pos_rhs))
pos_tol     = POS_ULP_SLACK * pos_eps * cp.maximum(pos_scale, 1e-300)
pos_bad_mag = cp.abs(pos_lhs - pos_rhs) > pos_tol
if bool(pos_bad_mag.any()):
    return _pos_abort("magnitude parting, %s, at index %d: ..." % ...)

pos_ok_a = ((pos_v_ra <= pos_v_max) & (pos_v_ra >= pos_v_min)
            & (pos_i_ra <= pos_i_max) & (pos_i_ra >= pos_i_min))
# ... the same for b and c ...
pos_split = (pos_ok_a != pos_ok_b) | (pos_ok_a != pos_ok_c)
if bool(pos_split.any()):
    return _pos_abort("VERDICT parting at index %d: ...")
```

Note `& | ~` with parentheses throughout: Python's `and` and `or` do not
vectorise, and silently taking the truth value of an array is exactly the kind
of fault this module exists to catch.

### Where the code and the mathematics differ

Two real discrepancies, both about scope rather than arithmetic.

1. **This is not a second channel for the enumeration kernel.** The kernel of
   [document two](02-the-enumeration-kernel.md) has four axes with the *first*
   fastest-varying, and it derives the count in series, the temperature and the
   rear gain from arithmetic progressions. This module has six axes with the
   *last* fastest-varying, and reads all six from tables. The two index the same
   physical quantities but they are different bijections over different spaces,
   so a given index does not name the same case in both. Handing this module a
   range of kernel indices and reading agreement as confirmation would be
   comparing two different questions. The module's own contract does not claim
   otherwise — it takes its radices and tables as arguments — but the pairing it
   provides is a pairing with *itself*, not with the kernel.

2. **The magnitude tolerance is a constant, not a bound.** The permitted slack
   is a fixed number of units in the last place, justified in the source as
   roughly what two roundings of a product plus one of a sum can produce. That
   is a reasonable order of magnitude and it is not a proof: the actual error of
   the distributed route grows with cancellation, and where the correction term
   nearly cancels the leading term, the true relative gap can exceed any fixed
   slack. In that region the module aborts. Aborting is the safe direction, but
   it means the slack sets a boundary of applicability rather than a guarantee.

A structural note rather than a fault: the verdict test uses closed intervals
where the enumeration kernel uses a strict comparison against an upper limit.
A case landing exactly on a limit therefore passes here and fails there.

---

## 3. Plain English

The same sum, written down two ways, is still the same sum. It is not always the
same answer.

A machine holds a number to a fixed number of digits and throws the rest away.
Multiply two numbers and the result is rounded. Add the rounded results and the
total is rounded again. So a calculation written as one multiplication of a
bracketed quantity, and the very same calculation with the bracket multiplied
out, do not round in the same places and can land on either side of a hair. On
paper the two expressions are indistinguishable. In a machine they are two
different experiments.

This channel makes use of that. For each of the two quantities it cares about, it
computes the answer three times: once with the correction gathered in a bracket,
once with the bracket multiplied out, and once with the factors applied in a
different order. All three are the same expression to anyone doing algebra, and
all three round differently.

The three answers are then judged twice over, and the two judgements have
deliberately different strictness.

The first judgement is about the numbers themselves, and it allows a little
slack. A few roundings apart is expected and means nothing. The routes must agree
to within a small allowance, measured relative to how large the numbers are
rather than as a fixed amount, because an allowance that made sense for a small
quantity would be meaningless for a large one.

The second judgement is about the verdict: is this case within its limits or not.
Here there is no slack at all, and that is the deliberate part. Two answers may
differ in their last digit and still mean exactly the same thing, so the
comparison that matters is not between the numbers but between the conclusions
drawn from them. Where the three routes draw different conclusions, the case is
sitting so close to a limit that the arithmetic itself cannot say which side it
is on. That case is not a near miss to be rounded and reported. It is a place
where the question has no answer at the precision available, and treating it as
an answer is how a search invents findings.

When that happens, the whole block is abandoned. Not the offending case, the
block. Nothing from it is kept, nothing is logged and carried on with, and the
failure is announced rather than tallied. The reason is arithmetic of a different
kind: a search running at this scale will meet vast numbers of cases, and a
process that walks past a parting and keeps going produces nonsense in bulk and
at speed. One inconsistency is a curiosity. The same inconsistency repeated
without limit is a report.

There is a second thing this channel checks, and it is the one that is easy to
forget. Every case is found by turning a single counting number into a set of
dial positions. If that translation is wrong, every route computes the answer for
the wrong case, and every route agrees with every other route perfectly, because
they are all being asked about the same wrong thing. Agreement between routes
cannot see a fault they share. So the translation itself is done twice, by two
unrelated methods: once by repeatedly dividing and taking remainders, once by
dividing by precomputed step sizes. If those two disagree by so much as one
position, the block is abandoned in the same way.

Two limitations should be said clearly.

This channel is a second opinion on itself, not on the enumeration described
elsewhere in this repository. It lays its dials out in the opposite order and
reads some of them from tables where the other builds them by counting. A given
case number does not name the same case in both. Running them side by side and
finding agreement would prove nothing, because they would not be answering the
same question. What this channel does establish is that its own answer does not
depend on the order in which its own arithmetic was written.

And the small allowance on the numbers is a judgement, not a theorem. It is set
at about what a few roundings can be expected to produce. Where a correction
almost exactly cancels the quantity it corrects, the true disagreement between
routes can be larger than any fixed allowance, and there the block is abandoned.
Abandoning is the safe direction to fail in. It is still a limit on where this
method applies rather than a guarantee about where it does.
