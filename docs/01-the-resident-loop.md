# The resident loop

Source: `trillion_loop.py`

---

## 1. The mathematics

### The space

A **bin** carries four measured quantities of one class of equipment:

$$b=(V,\ \beta,\ I,\ F)\in\mathcal{B}\subset\mathbb{R}^4,\qquad |\mathcal{B}|=m$$

$V$ the open-circuit voltage at the rating datum, $\beta$ its temperature
coefficient, $I$ the short-circuit current, $F$ the overcurrent rating printed
for that bin.

A **gene** carries the three free choices:

$$g=(n,\ \tau,\ \rho)\in\mathcal{G}=[n_-,n_+]\times[\tau_-,\tau_+]\times[\rho_-,\rho_+]$$

$n$ the count in series, $\tau$ the assessment temperature, $\rho$ the rear-side
gain. A **candidate** is a pair:

$$c=(b,g)\in\mathcal{C}=\mathcal{B}\times\mathcal{G}.$$

### Two channels

Fix the datum $\tau_0$ and the design factor $\varphi$. Define two pairs of maps
$\mathcal{C}\to\mathbb{R}$:

$$\nu_A(c)=n\,V\bigl(1+\beta(\tau-\tau_0)\bigr),\qquad
\iota_A(c)=\varphi\bigl(I(1+\rho)\bigr)$$

$$\nu_B(c)=nV+nV\bigl(\beta(\tau-\tau_0)\bigr),\qquad
\iota_B(c)=\varphi I+\varphi I\rho$$

Over $\mathbb{R}$, $\nu_A\equiv\nu_B$ and $\iota_A\equiv\iota_B$. Over the
machine's float set $\mathbb{F}$, with rounding $r:\mathbb{R}\to\mathbb{F}$
applied after each operation, $\nu_B$ performs one extra rounding of a product
and one of a sum, so in general

$$r\circ\nu_A\ \neq\ r\circ\nu_B\quad\text{on }\mathcal{C}.$$

All arithmetic below is in $\mathbb{F}$; the $r$ is suppressed.

### Limits and their provenance

$$\mathcal{L}=\{\ell_v,\ \ell_f,\ \ell_s,\ \ell_m\},\qquad
\chi:\mathcal{L}\to\{0,1\}$$

$\chi(\ell)=1$ if and only if $\ell$ occurs in the evidence base. Write
$\mathcal{S}=\chi^{-1}(1)$ and $\mathcal{U}=\mathcal{L}\setminus\mathcal{S}$.

### The verdict vector

For channel $X\in\{A,B\}$, with $[\,\cdot\,]$ the Iverson bracket,

$$q_X(c)=\bigl(\,[\nu_X>\ell_v],\ [\iota_X>F],\ [\iota_X>\ell_s],\ [\iota_X>\ell_m]\,\bigr)\in\{0,1\}^4 .$$

### The three survival predicates

$$P(c)\iff q_A(c)=q_B(c)$$

$$\sigma(c)\iff[\nu_A>\ell_v]\vee[\iota_A>F]\qquad\text{(a sourced limit fires)}$$

$$u(c)\iff[\iota_A>\ell_s]\vee[\iota_A>\ell_m]\qquad\text{(an unsourced limit fires)}$$

$$\Pi(c)\iff\neg\bigl(u(c)\wedge\neg\sigma(c)\bigr)$$

$$\delta(c)=\min\!\left(\frac{\lvert\nu_A(c)-\ell_v\rvert}{\ell_v},\ \frac{\lvert\iota_A(c)-F\rvert}{F}\right),
\qquad E(c)\iff\delta(c)<\varepsilon$$

$$\boxed{\ \Sigma(c)\iff P(c)\wedge\Pi(c)\wedge E(c)\ }$$

The complement partitions into three disjoint causes of death:

$$\Sigma^{-1}(0)=\neg P\ \sqcup\ \bigl(P\wedge\neg\Pi\bigr)\ \sqcup\ \bigl(P\wedge\Pi\wedge\neg E\bigr).$$

### The generation operator

Let $X_t\subset\mathcal{C}$ with $|X_t|=N$, let
$A_t=\{c\in X_t:\Sigma(c)\}$ and $G_t=\pi_{\mathcal{G}}(A_t)$ its gene
projection. Let $T_N$ repeat a finite multiset to length $N$, let $J$ be the
jitter law on $\mathbb{R}^3$, and let $\Pi_{\mathcal{G}}$ be projection onto the
box $\mathcal{G}$. Then

$$X_{t+1}=\begin{cases}
\bigl(U(\mathcal{B}),\ U(\mathcal{G})\bigr)^{\otimes N}, & A_t=\varnothing\\[6pt]
\Bigl\{\bigl(U(\mathcal{B}),\ \Pi_{\mathcal{G}}(g+j)\bigr)\ :\ g\in T_N(G_t),\ j\sim J\Bigr\}, & A_t\neq\varnothing
\end{cases}$$

where $U(\cdot)$ denotes an independent uniform draw. Note that the bin
coordinate is redrawn uniformly in both branches: only the gene is inherited.

### The write rule and its fixed point

Let $Q(G)\in\mathbb{R}^{5\times3}$ be the componentwise quantiles of $G$ at
$\{0,\tfrac14,\tfrac12,\tfrac34,1\}$, let $R_k$ round to $k$ places, and let $H$
be a cryptographic digest truncated to a fixed width. Define

$$s_t=H\bigl(R_k\,Q(G_t)\bigr).$$

The tape receives a record if and only if $s_t\neq s_{t-1}$, so

$$\#\text{writes}=\bigl\lvert\{\,t\le T\ :\ s_t\neq s_{t-1}\,\}\bigr\rvert .$$

The loop is at a **fixed point of the write rule** whenever $s_t=s_{t-1}$: the
population keeps moving, the drive does not.

### What is and is not established

$$\Sigma(c)\ \Longrightarrow\ c\in\ker\bigl(q_A-q_B\bigr)\quad\text{and nothing more.}$$

No $q_{\text{true}}$ exists in this system, so no statement of the form
$q_A=q_{\text{true}}$ is derivable from $\Sigma$. And since $X_{t+1}$ depends on
$A_t$, the examined multiset is not uniform on $\mathcal{C}$:

$$\frac{1}{\lvert\mathcal{C}\rvert}\sum_{t\le T}\lvert X_t\rvert\ \text{ is \textit{not} a coverage of }\mathcal{C}.$$

---

## 2. The code

| symbol | name in the code |
|---|---|
| $b=(V,\beta,I,F)$ | `voc, beta, isc, fmax` — one row of `pop` |
| $g=(n,\tau,\rho)$ | `ns, tmin, rear` — one row of `gene` |
| $\nu_A,\ \nu_B$ | `v_a`, `v_b` |
| $\iota_A,\ \iota_B$ | `i_a`, `i_b` |
| $\ell_v,\ \ell_s,\ \ell_m$ | `v_max`, `str_a`, `mppt_a` |
| $F$ | `fmax` (per bin, read from the feed) |
| $\chi$ | the boolean second field of each `LIMITS` entry |
| $P,\ \Pi,\ E$ | `!parted`, `!unsourced_only`, `boundary` |
| $\Sigma$ | `survives` |
| $\delta,\ \varepsilon$ | `min(dv, di)`, `edge` |
| $A_t$ | `out_gene[mask]` |
| $s_t$ | `sig` |

The survival predicate, one thread per candidate:

```c
/* CHANNEL A: the factored form. One rounding of the bracket. */
float ku_a = 1.0f + beta * (tmin - 25.0f);
float v_a  = ns * voc * ku_a;
float i_a  = 1.25f * (isc * (1.0f + rear));

/* CHANNEL B: the distributive law. Rounds twice where A rounds once. */
float v_b  = ns * voc + ns * voc * (beta * (tmin - 25.0f));
float i_b  = 1.25f * isc + 1.25f * isc * rear;

/* P: agreement on the VERDICT, not on the digits. */
int parted = ((v_a > v_max) != (v_b > v_max))
           | ((i_a > fmax)  != (i_b > fmax))
           | ((i_a > str_a) != (i_b > str_a))
           | ((i_a > mppt_a)!= (i_b > mppt_a));

/* Pi: a verdict resting only on limits with no source is not knowledge. */
int fires_sourced   = (v_a > v_max) | (i_a > fmax);
int fires_unsourced = (i_a > str_a) | (i_a > mppt_a);
int unsourced_only  = fires_unsourced & !fires_sourced;

/* E: within one step of changing its answer. */
float dv = fabsf(v_a - v_max) / v_max;
float di = fabsf(i_a - fmax)  / fmax;
int   boundary = (dv < edge) | (di < edge);

int survives = (!parted) & (!unsourced_only) & boundary;
alive[i] = survives;
if (survives) { out_gene[i*3+0] = ns;
                out_gene[i*3+1] = tmin;
                out_gene[i*3+2] = rear; }
```

The generation operator and the write rule, on the host, once per generation:

```python
mask    = alive.astype(bool)           # the indicator of Sigma
n_alive = int(mask.sum())
if n_alive == 0:                       # A_t empty: reseed uniformly
    gene[:, 0] = rng.integers(18, 34, size=nc).astype(cp.float32)
    gene[:, 1] = rng.random(nc, dtype=cp.float32) * 40.0 - 30.0
    gene[:, 2] = rng.random(nc, dtype=cp.float32) * 0.35
    continue

surv = out_gene[mask]                  # A_t, compacted on the card
reps = (nc + n_alive - 1) // n_alive
seed = cp.tile(surv, (reps, 1))[:nc]   # T_N(G_t)
jitter = cp.empty_like(seed)           # j ~ J
jitter[:, 0] = rng.integers(-1, 2, size=nc).astype(cp.float32)
jitter[:, 1] = (rng.random(nc, dtype=cp.float32) - 0.5) * 2.0
jitter[:, 2] = (rng.random(nc, dtype=cp.float32) - 0.5) * 0.02
gene  = cp.clip(seed + jitter, lo, hi) # projection onto G
idx   = rng.integers(0, len(B), size=nc)          # the bin is REDRAWN
d_pop = cp.asarray(B)[idx].astype(cp.float32).ravel()

q   = cp.asnumpy(cp.percentile(surv, [0, 25, 50, 75, 100], axis=0))
sig = hashlib.sha256(np.round(q, 3).tobytes()).hexdigest()[:16]
if sig != last_sig:                    # the write rule
    last_sig, writes = sig, writes + 1
    tape.write(json.dumps(record) + "\n")
```

Nothing else crosses from the card to the drive. The dead are never copied
anywhere; they are simply not selected by `mask`.

### Where the code and the mathematics differ

Three honest discrepancies, recorded rather than smoothed over.

1. **A fourth cause of death is announced and never computed.** The module
   header and the banner printed at every run name four ways to die, the fourth
   being a candidate whose neighbourhood has already converged. No such test
   exists in the kernel, no counter is allocated for it, and every candidate
   that is paired, sourced and interior is charged to the interior count. The
   mathematics above lists three causes because three is what the code has.

2. **The survivor's bin is not inherited.** The predicate $\Sigma$ is a property
   of the pair $(b,g)$, but the generation operator carries only $g$ forward and
   redraws $b$ uniformly. A gene that sat on the boundary for one bin is
   therefore re-tested against an unrelated bin, and the claim that the search
   "walks the surface" holds only to the extent that the surface is insensitive
   to the bin coordinate. This is written as $U(\mathcal{B})$ in the generation
   operator above so that the gap is visible rather than assumed away.

3. **One limit is passed and ignored.** The kernel takes $\ell_f$ as an argument
   and never reads it; the fuse comparison uses the per-bin rating $F$ instead.
   The argument is dead. The mathematics omits $\ell_f$ from every predicate for
   that reason, and keeps it in $\mathcal{L}$ only because the reported limit
   table still carries it.

A fourth point is a risk rather than a discrepancy: the two channels are
compiled into one kernel with fast-math enabled, which permits the compiler to
contract and reassociate. If it ever contracts both channels to the same
instruction sequence, they will agree on everything by construction and the
pairing test will silently become free. Nothing in the run would report this.

---

## 3. Plain English

A search is a machine for turning attention into knowledge, and the way to waste
it is to spend attention on things you already know.

This one begins with a population of possible arrangements. Each is a choice of
equipment class together with a few settings: how many units are wired in a
chain, how cold the site is assumed to get, how much extra light the back of the
equipment is assumed to collect. The population is large enough that no person
would ever read it and small enough to sit entirely inside the graphics card's
own memory. That last point is the whole design. The card's memory is very much
faster than the drive underneath it, so a loop that sent every evaluation down to
the drive would crawl along at a small fraction of the card's speed, and would
wear the drive out carrying numbers nobody will ever open. So the population
lives on the card for its whole life, and the drive is written to only when there
is something new to say.

For each arrangement the machine works out two things: the voltage the chain
would reach when it is coldest, and the current it would be asked to carry. It
works each of them out twice, by two different routes to the same answer. One
route keeps the correction gathered inside a bracket and multiplies once. The
other multiplies the bracket out first, so it rounds twice where the first route
rounds once. On paper the two routes are the same expression. Inside a machine
that can only hold so many digits they are not, and where they differ they differ
exactly at the places where the answer is delicate.

The machine then compares not the two numbers but the two verdicts: does this
arrangement go past the limit, yes or no. If the two routes give the same
verdict, the arrangement survives. If they give different verdicts, the pair
destroys itself and nothing at all is recorded about it. Not a warning, not a
note, not a tally in a file next to the word success. An arrangement that cannot
be arrived at twice is not arrived at.

Surviving that comparison is not enough. There are two further tests, and they
are what make this a search for knowledge rather than a search for output.

The first test is about where the limits came from. Some of the limits used here
can be pointed at in a document that somebody wrote while reading an equipment
sheet: the voltage the equipment is rated for, the current its protective device
is rated for. Others cannot. They are numbers that circulate in the trade and
appear in no document at all. An arrangement whose verdict turns only on a limit
with no source is counted and then dropped, because breeding the next generation
from it would fill that generation with the consequences of an invented number.
This is not a caveat added at the end. It is a term in the rule that decides who
lives, sitting beside the arithmetic and exactly as important as the arithmetic.

The second test is about information. An arrangement sitting comfortably inside a
settled region, far from every limit, tells you nothing you did not already know.
It is the expensive kind of nothing, because it looks exactly like work: the card
is busy, the counters climb, and afterwards you know what you knew before. So an
arrangement lives only if it sits within one small step of changing its answer.
Those arrangements trace out the edge of the region where the answer is known,
and the edge is the only place where another generation can teach you anything.

What survives, then, is the set of arrangements that two routes agree about, that
rest on limits somebody can point at, and that are one nudge away from flipping.
Their settings are copied, repeated until they fill the population again, and
nudged slightly, so that the next generation starts from the edge the last one
found rather than scattering across the whole space afresh. If a generation is
wiped out entirely, that is recorded as a fact and the population is drawn again
from nothing; an empty generation is an answer, not a failure.

There is a gap here that should be said out loud. The settings are carried
forward but the choice of equipment class is not: it is drawn again at random
each time. So a set of settings that sat on the edge for one class of equipment
is handed to a different class and asked the question again. That is a weaker
inheritance than the description of the loop implies, and the edge the search
walks is the edge in the settings only.

Finally, the machine takes a compact summary of the survivors, the sort that says
where the middle of the group sits and where its extremes lie, and reduces that
to a short fingerprint. If the fingerprint is the same as last time, nothing is
written. The search is still moving, but it has stopped saying anything new, and
a loop that writes a file every few seconds whether or not anything changed is
not keeping a record, it is keeping a noise. Only when the fingerprint changes
does a line go to the drive.

Two things must be said plainly, and they belong in the body rather than in small
print.

The first is that all of this establishes consistency, not truth. Both routes
were written by the same hands, in the same sitting, out of the same
understanding. If that understanding is wrong, both routes are wrong together,
they agree perfectly, and every arrangement sails through. A second opinion from
yourself is not a second opinion. What the agreement does rule out is narrower
and still worth having: an arrangement whose verdict depends on the order in
which somebody happened to write the arithmetic down.

The second is that the number of generations is not a measure of how much of the
space has been looked at. Later generations return deliberately to the edge that
earlier ones found. Effort is piled up there on purpose, which is the point of
the design, but it means the attention is concentrated in one region rather than
spread evenly, and dividing the number of evaluations by the size of the space
would give a figure that describes the searchlight rather than the room.
