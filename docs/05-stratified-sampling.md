# Stratified adaptive sampling

Source: `t04_learn.py`

Steering effort towards the interesting part of a space biases every count taken
from it. This controller steers hard and still reports an unbiased headline,
because the steering is confined to move effort *between* strata and never
*within* one.

---

## 1. The mathematics

### The space, and a partition into strata of equal volume

$$\mathcal{X}=\bigl(\mathbb{Z}_{R}\bigr)^{D},\qquad \lvert\mathcal{X}\rvert=R^{D}$$

Choose $b\mid R$ and put $s=R/b$. The cell map

$$\gamma:\mathcal{X}\to\mathbb{Z}_{C},\qquad
\gamma(x)=\sum_{d=0}^{D-1}\Bigl\lfloor\frac{x_d}{s}\Bigr\rfloor\,b^{\,D-1-d},
\qquad C=b^{D}$$

induces a partition $\mathcal{X}=\bigsqcup_{c} \gamma^{-1}(c)$ with every block of
the same cardinality:

$$\bigl\lvert\gamma^{-1}(c)\bigr\rvert=s^{D}=\frac{\lvert\mathcal{X}\rvert}{C}
\qquad\text{for every }c .$$

Equal volume is the load-bearing property, and it is exactly what $b\mid R$ buys.

### The quantity being estimated

Let $Y:\mathcal{X}\to\{0,1\}$ be the verdict. Define

$$p_c=\frac{1}{s^{D}}\sum_{x\in\gamma^{-1}(c)}Y(x),\qquad
P=\frac{1}{\lvert\mathcal{X}\rvert}\sum_{x\in\mathcal{X}}Y(x).$$

Because the blocks are equal in size,

$$P=\frac{1}{C}\sum_{c=0}^{C-1}p_c . \tag{$\ast$}$$

### The sampling design

At block $t$ the controller chooses an allocation $q^{(t)}\in\Delta^{C-1}$ as a
function of the data observed before block $t$. A draw is then made in two
stages:

$$c\sim q^{(t)},\qquad x\mid c\ \sim\ \mathrm{Unif}\bigl(\gamma^{-1}(c)\bigr).$$

The second stage is the whole argument. Conditional on the cell, the point is
uniform on that cell, so with $n_c$ draws landing in $c$ and $\hat p_c$ their
sample mean,

$$\mathbb{E}\bigl[\hat p_c \,\bigm|\, n_c=n,\ n\ge1\bigr]=p_c
\qquad\text{for every }n,$$

independently of how $q^{(t)}$ was chosen, because $q^{(t)}$ is measurable with
respect to the past and carries no information about *where inside* $c$ the draw
lands. Hence by $(\ast)$

$$\hat P=\frac{1}{C}\sum_{c}\hat p_c
\qquad\text{satisfies}\qquad \mathbb{E}\bigl[\hat P\bigr]=P$$

whenever every $n_c\ge1$. This is the Horvitz–Thompson argument in its simplest
form: the strata have known equal weight $1/C$, so the weights need not be
estimated.

$$\mathrm{Var}\bigl(\hat p_c\bigr)=\frac{p_c(1-p_c)}{n_c},\qquad
\mathrm{SE}\bigl(\hat P\bigr)=\frac{1}{C}\sqrt{\sum_{c}\frac{\hat p_c(1-\hat p_c)}{n_c}} .$$

### The allocation, and why it has a floor

Let $\hat p$ denote the current per-cell rate field and $n$ the per-cell counts.
Define the **surface** field by neighbour disagreement, maximised over axes:

$$\mathrm{surf}_c=\max_{0\le d<D}\ \max\Bigl(\bigl\lvert\hat p_c-\hat p_{c-e_d}\bigr\rvert,\ \bigl\lvert\hat p_c-\hat p_{c+e_d}\bigr\rvert\Bigr)$$

(edges of the grid replicated, so the boundary of the question is not mistaken
for a boundary of the verdict) and the **uncertainty** field

$$\mathrm{unc}_c=\begin{cases}\sqrt{\dfrac{\hat p_c(1-\hat p_c)}{n_c}}, & n_c>0\\[8pt] 1, & n_c=0\end{cases}$$

Then with $\lambda>0$ and a floor share $\phi\in(0,1)$,

$$\mathrm{sc}_c=\mathrm{surf}_c+\lambda\,\mathrm{unc}_c,\qquad
q_c=\frac{\phi}{C}+(1-\phi)\frac{\mathrm{sc}_c}{\sum_{c'}\mathrm{sc}_{c'}} .$$

$$\boxed{\ q_c\ \ge\ \frac{\phi}{C}\ >\ 0\quad\text{for every }c\text{ and every }t\ }$$

so no stratum can be starved, $\sum_t q^{(t)}_c=\infty$, and by Borel–Cantelli
every cell is visited infinitely often almost surely. The sum $(\ast)$ therefore
has no permanent holes.

### The second, unsteered stream

A fixed share $\psi$ of every block is drawn
$x\sim\mathrm{Unif}(\mathcal{X})$ and tallied in counters the controller never
touches, giving

$$\hat P_{U}=\frac{\#\{Y=1\}}{\#\text{drawn}},\qquad
\mathbb{E}\bigl[\hat P_{U}\bigr]=P,\qquad
\mathrm{SE}\bigl(\hat P_{U}\bigr)=\sqrt{\frac{\hat P_U(1-\hat P_U)}{\#\text{drawn}}} .$$

It is slow and wide. Its purpose is to be assumption-free, so that

$$\bigl\lvert\hat P-\hat P_{U}\bigr\rvert\ \le\ 3\sqrt{\mathrm{SE}(\hat P)^2+\mathrm{SE}(\hat P_U)^2}$$

is a falsifiable statement about the stratified estimator rather than a
restatement of it.

### The quantity that must never be quoted

$$\hat P_{A}=\frac{\sum_c \hat t_c^{\,A}}{\sum_c n_c^{\,A}}
\ \xrightarrow{\ \ }\ \frac{\sum_c q_c\,p_c}{\sum_c q_c}\ \neq\ P
\quad\text{in general.}$$

This is the rate of the searchlight, not of the room. It is computed and printed
solely as a diagnostic of how hard the controller is steering.

### What is learned

Per-axis first-order sensitivity, by the law of total variance:

$$\mathcal{S}_d=\frac{\mathrm{Var}_{i}\Bigl(\mathbb{E}\bigl[\hat p\mid \text{axis }d=i\bigr]\Bigr)}{\mathrm{Var}(\hat p)}\ \in[0,1]$$

with $\mathcal{S}_d\approx0$ meaning axis $d$ is inert: conditioning on it does
not move the mean.

Convergence is measured, not asserted:

$$\Delta p^{(t)}=\frac{1}{C}\sum_c\bigl\lvert\hat p_c^{(t)}-\hat p_c^{(t-1)}\bigr\rvert,
\qquad
\Delta w^{(t)}=\bigl\lVert w^{(t)}-w^{(t-1)}\bigr\rVert_1 .$$

### The caveat that survives all of the above

The headline is computed over the **resolved** set
$\mathcal{R}=\{c:n_c>0\}$, not over all of $\mathbb{Z}_C$:

$$\hat P=\frac{1}{\lvert\mathcal{R}\rvert}\sum_{c\in\mathcal{R}}\hat p_c
\qquad\text{estimates}\qquad
\frac{1}{\lvert\mathcal{R}\rvert}\sum_{c\in\mathcal{R}}p_c ,$$

which equals $P$ only when $\mathcal{R}=\mathbb{Z}_C$. Since $\mathcal{R}$ is
random and data-dependent, the unbiasedness claim is conditional on full
resolution, and $\lvert\mathbb{Z}_C\setminus\mathcal{R}\rvert$ is reported
alongside so that the condition can be checked rather than assumed.

---

## 2. The code

| symbol | name in the code |
|---|---|
| $D,\ R,\ b,\ s,\ C$ | `NDIM`, `RES`, `BINS`, `SUB`, `NCELL` |
| $\gamma$ | `cell_of()` |
| $n_c$ (steered / unsteered) | `a_seen` / `u_seen` |
| cell true counts | `a_true` / `u_true` |
| $\hat p_c,\ n_c$ | the pair returned by `rate_field()` |
| $\mathrm{surf}_c$ | `surface()[0]` |
| $\mathrm{unc}_c$ | `unc` inside `next_plan()` |
| $\phi,\ \psi,\ \lambda$ | `FLOOR_SHARE`, `UNIFORM_SHARE`, `UCB_C` |
| $q_c$ | `prob`, cached as `_alloc_p` |
| $\hat P,\ \mathrm{SE}(\hat P)$ | `stratified_rate`, `stratified_stderr` |
| $\hat P_U$ | `uniform_rate` |
| $\hat P_A$ | `adaptive_raw_rate` — tagged `BIASED` |
| $\mathcal{S}_d$ | `sensitivity()` |
| $\Delta p,\ \Delta w$ | the pair returned by `convergence()` |
| $\mathcal{R}$ | `resolved`, counted as `cells_resolved` |

The partition, defined once and asserted:

```python
NDIM, RES, BINS = 6, 100, 10
assert RES % BINS == 0          # equal-volume strata, the load-bearing fact
SUB   = RES // BINS
NCELL = BINS ** NDIM
```

The two-stage draw. The second stage is what keeps the per-cell estimate honest:

```python
if stream == "uniform":
    return rng.integers(0, RES, size=(n, NDIM), dtype=cp.int64)

cdf  = cp.cumsum(self._alloc_p)          # choose a cell from q
u    = rng.random(n, dtype=cp.float64)
cid  = cp.clip(cp.searchsorted(cdf, u).astype(cp.int64), 0, NCELL - 1)
rem  = cid
for d in range(NDIM - 1, -1, -1):        # then UNIFORM inside that cell
    b = rem % BINS
    rem = rem // BINS
    coords[:, d] = b * SUB + rng.integers(0, SUB, size=n, dtype=cp.int64)
return coords
```

The allocation, with the floor written as a separate additive term so that it
cannot be optimised away:

```python
unc   = cp.sqrt(cp.maximum(p * (1.0 - p), 1e-12) / cp.maximum(n, 1.0))
unc   = cp.where(n > 0, unc, 1.0)        # never seen: maximum pull
score = s + UCB_C * unc
tot   = float(score.sum())
prob  = (cp.full(NCELL, 1.0 / NCELL) if tot <= 0
         else (FLOOR_SHARE / NCELL) + (1.0 - FLOOR_SHARE) * (score / tot))
prob  = prob / prob.sum()
```

The three headline numbers, each carrying its own tag:

```python
p_strat = float(p[resolved].mean()) if nres else float("nan")
var     = cp.where(resolved, p * (1 - p) / cp.maximum(n, 1.0), 0.0)
se_strat = float(cp.sqrt(var[resolved].sum())) / max(nres, 1)

p_unif  = ut / us if us else float("nan")
p_adapt = at / as_ if as_ else float("nan")

return {
  "stratified_rate":  p_strat, "stratified_bias": "UNBIASED (equal-volume strata, floor>0)",
  "uniform_rate":     p_unif,  "uniform_bias":    "UNBIASED (never steered)",
  "adaptive_raw_rate":p_adapt, "adaptive_raw_bias":"BIASED - searchlight rate, NOT coverage",
  "cross_check_3sigma_agree": agree,
  "cells_resolved": nres, "cells_unresolved": NCELL - nres,
  ...
}
```

Ingestion is pure counting with no policy in it at all, which is what keeps the
estimator separable from the controller:

```python
def observe(self, br):
    cid  = cell_of(br.coords)
    seen = cp.bincount(cid, minlength=NCELL).astype(cp.float64)
    true = cp.bincount(cid, weights=br.verdict.astype(cp.float64), minlength=NCELL)
    if br.stream == "uniform":
        self.u_seen += seen; self.u_true += true
    else:
        self.a_seen += seen; self.a_true += true
```

### Where the code and the mathematics differ

1. **The headline is an average over resolved cells, and is labelled unbiased
   without qualification.** With a grid of a million cells and a block of a few
   million draws, of which only the floor share is spread flat, a great many
   cells are unvisited for many blocks. The printed tag reads unbiased; the
   property actually held is unbiasedness *for the mean over the cells that have
   been visited*. The count of unvisited cells is reported honestly right beside
   it, so the shortfall is visible, but the tag overstates what has been
   established and a reader who takes the tag at face value will overstate it
   too. The two numbers should be read together or not at all.

2. **The standard error treats the per-cell estimates as independent and the
   allocation as fixed.** Neither holds exactly: the counts are themselves random
   and the allocation at each block depends on the estimates from the previous
   one. The formula is the right first approximation and it is not a confidence
   interval. The unsteered cross-check is the defence against this, which is
   the reason it exists.

3. **The block bookkeeping does not index the space.** Blocks carry a running
   base and a count, in the shape of an enumeration, but every draw is random;
   the base plays no part in choosing points. It survives as a seed ingredient
   and a label. Nothing depends on it, and nothing should be read into it.

---

## 3. Plain English

Point a searchlight at the part of a field where you expect to find something,
count what you find, divide by the number of places you looked, and you have
measured your own expectations. This is the most comfortable mistake in
searching, because the number it produces looks exactly like a measurement and is
usually a flattering one.

The space here is far too large to walk. Almost all of it is flat: the answer is
the same everywhere and confirming that again teaches nothing. What is worth
finding is the surface, the places where the answer changes, and finding it means
deliberately spending most of the effort there. So the steering is not optional.
The question is how to steer and still be able to state an honest figure for the
whole space afterwards.

The answer is to confine the steering. The space is cut into a grid of boxes, all
exactly the same size, and the steering is only ever allowed to choose how many
draws each box receives. Inside a box, the draws are spread evenly and the
controller has no say at all.

That single restriction is what makes the whole thing work. Because the choice of
box carries no information about where inside the box the draw will land, each
box's own answer is a fair estimate of that box's own truth, no matter how many
or how few draws it was given. A box that received a great many draws has a
precise estimate; a box that received very few has a vague one; neither is
skewed. And since the boxes are all the same size, the answer for the whole space
is simply the ordinary average of the boxes' answers, with no weighting needed,
because the weights are known in advance and equal. That average is the figure
worth quoting.

Three further precautions sit underneath it.

Every box keeps a guaranteed minimum share of every round, for ever. The
guarantee is small but it is never zero, so no box can be starved out of the
average entirely, and a box that has been ignored for a long time is eventually
visited. Boxes that have not yet been visited at all are counted and reported as
unresolved rather than quietly filled in with the overall average, which would
be assuming the answer in the places least looked at.

A fixed fraction of every round is drawn completely at random across the whole
space, ignoring everything the controller has learned, and kept in its own
separate tally that the controller never touches. This stream is slow and its
answer is imprecise, but it rests on no assumptions whatever. If the careful
stratified figure and the crude random figure drift apart by more than their
combined uncertainty allows, something is wrong with the reasoning, and the
report says so. A check that can only ever agree is not a check.

And the naive figure, the one obtained by simply dividing everything found by
everything looked at, is computed and printed too, clearly labelled as biased. It
is not a measurement of the space and must never be quoted as one. It is a useful
diagnostic of how hard the controller is currently steering, and nothing else.

What the controller actually learns is two things. First, where the surface is:
for each box it compares its answer with those of its immediate neighbours, and a
box whose neighbours disagree with it is on the surface and is worth returning
to. Second, which of the dials actually matter: if sorting the results by one
particular dial does not move the average at all, that dial is inert, its weight
collapses, and the search stops paying to vary it. Boxes that merely look settled
because they have hardly been sampled are protected by a separate term that pulls
attention towards anything poorly known, so a box is never written off on the
strength of two or three draws.

Whether the search has settled is measured rather than declared. Each round the
controller records how much the picture moved since the last one, and how much
the dial weights moved. Those numbers are reported. Nobody announces convergence.

One caution belongs at the end rather than buried. The headline figure is the
average over the boxes that have actually been visited, and the fairness argument
is an argument about exactly that. When many boxes are still unvisited, the
figure describes the visited part of the space and the rest is an open question.
The number of unvisited boxes is printed right next to the headline for that
reason, and the two are meant to be read together. A reader who takes the
headline alone, on the strength of the word attached to it, will claim more than
has been shown.
