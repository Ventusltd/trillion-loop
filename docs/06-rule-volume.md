# Rule volume, surface, and the implication lattice

Source: `t12_volume.py`

A rule is not a test. A rule is a region. Once it is treated as a region, two
questions become askable that no count of failures can answer: how much of the
space does it claim, and does any other rule already claim all of it.

---

## 1. The mathematics

### The space

$$\mathcal{X}=\mathcal{B}\times\mathcal{S}\times\mathcal{P}\times\mathcal{T}
=\mathbb{Z}_{m}\times\mathbb{Z}_{N_S}\times\mathbb{Z}_{N_P}\times\mathbb{Z}_{N_T}$$

indexed by bin $i$, series count $\sigma$, rear-gain step $\pi$ and temperature
step $\theta$, with the axis values on arithmetic grids

$$\rho(\pi)=\rho_-+\pi\frac{\rho_+-\rho_-}{N_P-1},\qquad
\tau(\theta)=\tau_-+\theta\frac{\tau_+-\tau_-}{N_T-1},\qquad
n(\sigma)=n_-+\sigma .$$

### Two families of rule

Each rule $r$ is a predicate $\mathcal{X}\to\{0,1\}$ of one of two shapes. With
$\mu_r$ its multiplier, $L_r(i)$ its limit (a constant or a per-bin rating), and
per-bin quantities $\varphi_i,\alpha_i,\beta_i,V_i$:

$$r\in\mathsf{I}:\qquad
F_r(i,\pi,\theta)=\mu_r\Bigl(X_r(i)\bigl(1+\varphi_i\,\rho(\pi)\bigr)\Bigr)\Bigl(1+\tfrac{\alpha_i}{100}\bigl(\tau(\theta)-\tau_0\bigr)\Bigr)>L_r(i)$$

$$r\in\mathsf{V}:\qquad
G_r(i,\sigma,\theta)=\mu_r\Bigl(n(\sigma)\bigl(1+\tfrac{\beta_i}{100}(\tau(\theta)-\tau_0)\bigr)\Bigr)V_i>L_r(i)$$

A rule of family $\mathsf{I}$ does not depend on $\sigma$; a rule of family
$\mathsf{V}$ does not depend on $\pi$. Call the ignored axis the rule's **free
axis** and the axis it varies along its **own axis**.

### Monotonicity, and the suffix structure it forces

Assume, for every bin,

$$\varphi_i\ge0,\qquad 1+\tfrac{\alpha_i}{100}\bigl(\tau(\theta)-\tau_0\bigr)>0,
\qquad 1+\tfrac{\beta_i}{100}\bigl(\tau(\theta)-\tau_0\bigr)>0 .$$

Then at fixed $(i,\theta)$, $F_r$ is nondecreasing in $\pi$ and $G_r$ is
nondecreasing in $\sigma$. A monotone predicate on a totally ordered finite set
is an **up-set**, so there exists a threshold

$$\vartheta_r(i,\theta)\in\{0,1,\dots,n_r\},\qquad
n_r=\begin{cases}N_P,& r\in\mathsf{I}\\ N_S,& r\in\mathsf{V}\end{cases}$$

such that the rule fires exactly on the suffix

$$r(i,\cdot,\theta)^{-1}(1)=\{\,k\ :\ k\ge\vartheta_r(i,\theta)\,\},$$

with $\vartheta_r=n_r$ meaning it never fires at that $(i,\theta)$. The whole
region is therefore described by the integer field $\vartheta_r$ of size
$m\times N_T$, no matter how large $N_P$ and $N_T$ are.

$\vartheta_r$ is obtained by **bisection on the predicate itself**, not by
algebraically inverting it, so no rounding is introduced that brute force would
not also have made.

### Volume

With $M_r$ the cardinality of the free axis,

$$\bigl\lvert r\bigr\rvert=M_r\sum_{i,\theta}\bigl(n_r-\vartheta_r(i,\theta)\bigr),
\qquad \lvert\mathcal{X}\rvert=m\,N_S\,N_P\,N_T .$$

$$r\ \text{is \textbf{dead}}\iff \lvert r\rvert=0\ \ \text{or}\ \ \lvert r\rvert=\lvert\mathcal{X}\rvert .$$

A rule of zero volume never fires anywhere and can never change an answer. A rule
of full volume fires everywhere and separates nothing. Both are invisible to
anyone who only counts failures, and both are worth finding.

### Surface

$$\partial r=\Bigl\{\,x\in r\ :\ \exists\,y\ \text{adjacent to}\ x\ \text{in}\ \mathcal{X},\ y\notin r\,\Bigr\}$$

where adjacency is one step along bin, series, rear gain or temperature, and the
edge of the swept box is **not** counted: that is a boundary of the question, not
of the rule. In terms of $\vartheta$, writing
$\vartheta^{\max}_r(i,\theta)$ for the maximum of $\vartheta_r$ over
$(i,\theta)$ and its four neighbours, clipped to $n_r$:

$$\bigl\lvert\partial r\bigr\rvert=M_r\sum_{i,\theta}\max\Bigl(
\underbrace{\vartheta^{\max}_r-\vartheta_r}_{\text{flips along }i\text{ or }\theta},\ \
\underbrace{\mathbf 1\bigl[0<\vartheta_r<n_r\bigr]}_{\text{flips along its own axis}}\Bigr)$$

$$\text{ratio}(r)=\frac{\lvert\partial r\rvert}{\lvert r\rvert}\in[0,1].$$

Near zero: a fat region with a thin skin, robust to a small change in any input.
Near one: the verdict is settled by the last digit of an input.

### The implication lattice

For every ordered pair,

$$A\Rightarrow B\ \iff\ \lvert A\rvert>0\ \wedge\ \lvert A\rvert-\lvert A\cap B\rvert=0
\ \iff\ \varnothing\neq A\subseteq B .$$

$$A\ \text{strictly inside}\ B\iff A\subseteq B\ \wedge\ B\not\subseteq A$$
$$A\equiv B\iff \lvert A\rvert=\lvert B\rvert=\lvert A\cap B\rvert\neq0$$
$$A\perp B\iff \lvert A\cap B\rvert=0,\ \lvert A\rvert>0,\ \lvert B\rvert>0$$

The joint counts come in closed form from the thresholds. Same family: the
intersection of two suffixes along the same axis is a suffix at the larger
threshold,

$$\lvert A\cap B\rvert=M\sum_{i,\theta}\bigl(n-\max(\vartheta_A,\vartheta_B)\bigr).$$

Different families: at fixed $(i,\theta)$ the two regions are suffixes along
*different* axes, so the intersection is their product,

$$\lvert A\cap B\rvert=\sum_{i,\theta}\bigl(N_P-\vartheta_{\mathsf{I}}(i,\theta)\bigr)\bigl(N_S-\vartheta_{\mathsf{V}}(i,\theta)\bigr).$$

Consequence worth stating: $A$ strictly inside $B$ means $A$ can **never be the
binding constraint**. It refuses nothing that $B$ has not already refused.

### Exactness, and the second implementation

Every count above is the count brute force would return. It is not a sample and
not an estimate: the cost of a trillion cases is the cost of a field of size
$m\times N_T$, because the shape of the region is known rather than explored.

The claim is checked and not asserted. On a box small enough for explicit boolean
arrays, both paths are run:

$$\text{accept}\iff \forall r:\ \lvert r\rvert_{\text{shape}}=\lvert r\rvert_{\text{brute}}
\ \wedge\ \lvert\partial r\rvert_{\text{shape}}=\lvert\partial r\rvert_{\text{brute}}
\ \wedge\ \forall (A,B):\ \lvert A\cap B\rvert_{\text{shape}}=\lvert A\cap B\rvert_{\text{brute}}$$

with **exact** equality, no tolerance, and nothing reported at all if it fails.

### What is not claimed

$$\lvert r\rvert,\ \lvert\partial r\rvert\ \text{are properties of }r\ \text{over the chosen grid.}$$

Halve the step and the surface roughly halves. That is geometry, not physics. And
where $L_r$ is a constant with no source, every statement about $r$ — its volume,
its deadness, its place in the lattice — is a fact about that constant and about
no hardware anywhere.

---

## 2. The code

| symbol | name in the code |
|---|---|
| $m,\ N_S,\ N_P,\ N_T$ | `NB`, `NS`, `NR`, `NT` |
| $\rho(\pi),\ \tau(\theta),\ n(\sigma)$ | `rear`, `temp`, `ser` from `axes()` |
| $\mu_r,\ X_r,\ L_r$ | `r["mult"]`, `r["x"]`, `limit_of(r, b)` |
| $\mathsf{I},\ \mathsf{V}$ | `r["g"] == "i"`, `r["g"] == "v"` |
| $\varphi_i,\alpha_i,\beta_i,V_i$ | `b["phi"]`, `b["alpha"]`, `b["beta"]`, `b["voc"]` |
| $\vartheta_r$ | `starts[r["id"]]`, computed by `_starts()` |
| $n_r,\ M_r$ | `_span(r, NR)`, `_other(r, NR)` |
| $\lvert r\rvert,\ \lvert\partial r\rvert$ | `volume`, `surface` from `volumes()` |
| $\lvert A\cap B\rvert$ | `joint[(a, b)]` from `implications()` |
| the lattice relations | `lattice()` returning `imp, strict, equiv, excl, dead` |
| the second path | `brute()` |
| unsourced limits | `INVENTED`, and `flag()` on every line that touches one |

The threshold, by bisection on the real predicate:

```python
if r["g"] == "i":
    kb   = 1.0 + b["phi"] * rear                     # (NR,)
    kt   = 1.0 + b["alpha"] / 100.0 * (temp - 25.0)  # (NT,)
    base = b[r["x"]] * kb                            # (NR,)
    lo = xp.zeros(NT, dtype=xp.int64)
    hi = xp.full(NT, NR, dtype=xp.int64)
    while bool((lo < hi).any()):
        mid   = (lo + hi) // 2
        val   = r["mult"] * (base[xp.clip(mid, 0, NR - 1)] * kt)
        fires = val > L                              # the predicate itself
        hi = xp.where(fires, mid, hi)
        lo = xp.where(fires, lo, mid + 1)
    s0[bi] = lo
```

Volume and surface from the threshold field alone:

```python
vol = int(other) * int(xp.asnumpy(xp.sum(n - s0)))

m = s0.copy()                                   # neighbours along temperature
if NT > 1:
    m = xp.maximum(m, xp.concatenate([s0[:, :1], s0[:, :-1]], axis=1))
    m = xp.maximum(m, xp.concatenate([s0[:, 1:], s0[:, -1:]], axis=1))
if NB > 1:                                      # and along bin index
    m = xp.maximum(m, xp.concatenate([s0[:1], s0[:-1]], axis=0))
    m = xp.maximum(m, xp.concatenate([s0[1:], s0[-1:]], axis=0))
m = xp.minimum(m, n)

band  = xp.maximum(m - s0, 0)                   # flips along bin/temperature
first = ((s0 > 0) & (s0 < n)).astype(xp.int64)  # the first cell of the suffix
sur   = int(other) * int(xp.asnumpy(xp.sum(xp.maximum(band, first))))
```

The joint counts, the two cases separated:

```python
if a["g"] == b["g"]:                      # same axis: suffix at the larger start
    j = int(_other(a, NR)) * int(xp.asnumpy(xp.sum(n - xp.maximum(sa, sb))))
else:                                     # different axes: a product
    ci = sa if a["g"] == "i" else sb
    cv = sa if a["g"] == "v" else sb
    j = int(xp.asnumpy(xp.sum((NR - ci) * (NS - cv))))
```

The lattice:

```python
for a in ids:
    if   single[a] == 0:     dead.append((a, "zero volume - it never fires anywhere"))
    elif single[a] == total: dead.append((a, "full volume - it fires everywhere, "
                                             "so it separates nothing"))
for a in ids:
    for b in ids:
        if a != b and single[a] and single[a] - joint[(a, b)] == 0:
            imp.append((a, b, single[a], single[b]))
strict = [(a, b, sa, sb) for (a, b, sa, sb) in imp
          if not (single[b] - joint[(a, b)] == 0)]
```

The pair rule, with no tolerance and no partial report:

```python
for r in RULES:
    if v_shape[r["id"]]["volume"]  != v_brute[r["id"]]: bad.append("volume "  + r["id"])
    if v_shape[r["id"]]["surface"] != s_brute[r["id"]]: bad.append("surface " + r["id"])
for k, v in j_shape.items():
    if j_brute[k] != v: bad.append("joint %s & %s" % k)
if bad:
    print("  PAIR RULE BROKEN, nothing reported:")
    return 2
```

### Where the code and the mathematics differ

1. **The monotonicity preconditions are assumed, never checked.** The bisection
   is valid only because the predicate is monotone along the axis being bisected,
   and that requires the rear-gain coefficient to be non-negative and both
   temperature factors to stay positive across the swept range. All three hold
   for the inputs in use. None is asserted in the code. Feed it a bin with a
   negative rear coefficient, or a temperature range wide enough to drive a
   factor through zero, and the bisection returns a threshold for a region that
   is not a suffix — silently, with no error, and with the brute-force check
   passing on a small box where the sign happens not to flip. This is the single
   assumption the exactness claim rests on and it deserves a guard.

2. **The two paths are not fully independent.** Both derive from the same bin
   parsing, the same axis construction and the same expression written in the
   same order of operations. The order is matched deliberately, so that no
   rounding difference is introduced. The price of that choice is that the
   agreement tests the *shape* argument and not the *arithmetic*: a mistake in
   the shared expression is reproduced identically by both paths and confirmed.
   The source says this plainly in its own closing statement; it is repeated
   here because it is the thing most likely to be forgotten.

3. **The check box is small and the reported sweep is not.** Exact agreement is
   established on a box of a few hundred steps per axis and then relied on at a
   scale many orders larger. The shape argument does carry across — that is what
   makes it an argument — but the *verification* does not. What is verified is
   that the two paths agree where both can run.

---

## 3. Plain English

Counting how often a rule refuses something tells you about the cases you
happened to feed it. Asking how much of the space a rule claims tells you about
the rule.

Treat each rule as a region. Somewhere in the space of possible arrangements
there is a set of them that this rule objects to, and everywhere else it is
silent. That set has a size. It has a shape. It has a surface, meaning the part
of it lying right against the boundary where the rule changes its mind. And it
sits in some relation to the regions claimed by every other rule: inside one,
overlapping another, perhaps nowhere near a third.

Two findings fall straight out of that and neither is reachable by counting
failures.

A rule with no volume at all never fires anywhere in the space. It sits in the
list looking like a safety check and it does nothing. It has never changed an
answer and it never will. It passes every test, because it is never invoked.

A rule that claims the entire space is dead in the opposite way. It objects to
everything, which distinguishes nothing, and a check that refuses every case is
as useless as one that refuses none.

The other finding concerns pairs. If one rule's region sits entirely inside
another's, then the inner rule can never be the deciding constraint. Every case
it refuses has already been refused by the outer one. It may be perfectly
correct, carefully reasoned, and completely without effect on any answer the
system ever gives, and nobody looking at pass counts would ever notice.

The obstacle is size. Checking a rule against every arrangement, and against
every neighbour of every arrangement, is far more work than any machine here will
do. It does not have to be done that way, because these rules have a shape, and
the shape can be exploited.

Every rule here compares a quantity that only ever grows as you turn one
particular dial. Turn up the amount the back of the equipment collects and the
current only rises. Add another unit to the chain and the voltage only rises. So
along that dial, once the rule starts objecting it never stops. Its region is not
a scatter of isolated cases; it is everything from some position onward.

That means the entire region can be described by a single whole number for each
combination of the other settings: the position on the dial where objecting
begins. Finding that position is a matter of halving the range and testing,
halving again and testing, a handful of times, and the test performed is the real
rule applied to the real arrangement, not a rearranged formula, so no new
rounding creeps in anywhere. Once those positions are known, the size of the
region, the size of every overlap between two regions, and the size of every
surface all follow by arithmetic, exactly, for any number of cases whatever.
There is no sampling here and no estimating. The answer is the one brute force
would have given, obtained by knowing the shape instead of walking it.

The claim is checked rather than trusted. On a range small enough that every
single arrangement can be written out and examined one at a time, both methods
are run and every number compared. They must match exactly, with no allowance,
and if any one of them does not, nothing is reported at all.

Four cautions, and they belong in the body.

Half the rules here rest on limits that appear in no document. Every statement
about those rules, whether about size, deadness or where they sit relative to
others, is a fact about an invented number and not about any equipment anywhere.
Each line of the report that touches one says so on the line itself, rather than
in a note at the bottom.

A size and a surface are properties of the rule measured against the particular
grid that was swept, at the particular step size chosen. Use a finer step and the
surface shrinks roughly in proportion. That is a fact about the measuring, not
about the thing measured.

The shape shortcut is valid only because the quantity really does rise as the
dial is turned. If it ever fell instead, the halving procedure would return a
confident answer for a region that is not shaped the way it assumes, with no
complaint and no sign that anything was wrong. Nothing in the code checks that
the quantity rises. It happens to, for the inputs in use.

And both methods, the shortcut and the exhaustive one, were written here, from
the same formulas, in the same order. Where they agree, they establish that the
shortcut reasoning about shape is sound. They cannot establish that the formula
they share is right. Agreement between two things built in the same room proves
consistency, never truth.
