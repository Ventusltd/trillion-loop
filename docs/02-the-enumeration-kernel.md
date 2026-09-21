# The enumeration kernel

Source: `t01_kernel.py`

The loop in [document one](01-the-resident-loop.md) samples. This kernel
*enumerates*: it walks a product space in a fixed order, with no random draws
anywhere, so that a stated range of case indices covers exactly the cases it
claims to and no others.

---

## 1. The mathematics

### The space as a mixed-radix numeral system

Let the space be a product of four finite axes with radices
$R=(r_0,r_1,r_2,r_3)\in\mathbb{N}^4$:

$$\mathcal{X}=\prod_{k=0}^{3}\mathbb{Z}_{r_k},\qquad
N=\lvert\mathcal{X}\rvert=\prod_{k=0}^{3}r_k .$$

Define the place values $p_0=1$ and $p_k=\prod_{j<k}r_j$, and the decoding
bijection

$$D:\ \mathbb{Z}_N\ \longrightarrow\ \mathcal{X},\qquad
D(x)=\bigl(d_0,d_1,d_2,d_3\bigr),\quad
d_k=\Bigl\lfloor \frac{x}{p_k}\Bigr\rfloor \bmod r_k .$$

Digit $0$ is the fastest-varying. $D$ is a bijection; its inverse is
$D^{-1}(d)=\sum_k d_k p_k$.

When every radix below the top is a power of two, $r_k=2^{a_k}$ for $k<3$, the
same $D$ is computed by masking and shifting:

$$d_k=\Bigl(x\gg\textstyle\sum_{j<k}a_j\Bigr)\ \&\ (r_k-1),\qquad k<3 .$$

This is an identity, not an approximation: for $r=2^a$, $y\bmod r=y\ \&\ (r-1)$
and $\lfloor y/r\rfloor=y\gg a$ for $y\ge0$.

### The point associated with an index

Let $\mathcal{B}$ be the bin table, $\mathcal{B}[d_0]=(V,\beta,I)$. The remaining
three digits index arithmetic progressions:

$$n(d_1)=n_0+d_1\,\Delta n,\qquad
\tau(d_2)=\tau_0'+d_2\,\Delta\tau,\qquad
\rho(d_3)=\rho_0+d_3\,\Delta\rho .$$

### The two quantities and the verdict

$$\nu(x)=n\,V\bigl(1+\beta(\tau-\tau_0)\bigr),\qquad
\kappa(x)=\varphi\,I(1+\rho)$$

$$a(x)=[\nu(x)>\ell_v],\qquad b(x)=[\kappa(x)>\ell_i]$$

Headroom, defined only where the case passes:

$$h_\nu(x)=1-\frac{\nu(x)}{\ell_v},\qquad h_\kappa(x)=1-\frac{\kappa(x)}{\ell_i}.$$

### The counters as a measure on a partition

The verdict space is partitioned into four blocks, and the passing block is
further marked by three (not mutually exclusive) headroom predicates. For a
range $\mathcal{R}=\{x_0,\dots,x_0+w-1\}\subset\mathbb{Z}_N$ define the
counting measure $\mu(S)=\lvert S\cap\mathcal{R}\rvert$ and the eight tallies

$$\begin{aligned}
c_0&=\mu(\mathcal{R}) &
c_1&=\mu(\neg a\wedge\neg b) \\
c_2&=\mu(a\wedge\neg b) &
c_3&=\mu(\neg a\wedge b) \\
c_4&=\mu(a\wedge b) &
c_5&=\mu(\neg a\wedge\neg b\wedge[h_\nu<\theta]) \\
c_6&=\mu(\neg a\wedge\neg b\wedge[h_\kappa<\theta]) &
c_7&=\mu(\neg a\wedge\neg b\wedge[h_\nu>\theta_{\text{slack}}])
\end{aligned}$$

with the invariant

$$c_1+c_2+c_3+c_4=c_0=w .$$

### The reduction, and why the order of summation is forced

Each thread accumulates a private vector $c^{(t)}\in\mathbb{N}^8$ over the
indices it owns under the grid-stride partition

$$\mathcal{R}=\bigsqcup_{t} \{\,x_0+t+jS\ :\ j\ge0\,\}\cap\mathcal{R},
\qquad S=\text{(threads per block)}\times\text{(blocks)} ,$$

then block $\mathcal{P}$ reduces to $c^{(\mathcal{P})}=\sum_{t\in\mathcal{P}}c^{(t)}$
in shared memory, and the grid reduces to
$c=\sum_{\mathcal{P}}c^{(\mathcal{P})}$ in global memory. Because
$(\mathbb{N},+)$ is associative and commutative, the result is independent of the
partition and of the completion order; the tally is the *same integer* for every
choice of grid and block shape. This is the property that makes integer counting
safe to parallelise and floating-point summation not.

### The caller's obligation

The kernel evaluates $D$ on $\mathcal{R}$ without a bounds check. The top digit
$d_3$ is never reduced modulo $r_3$ (there is nothing above it to divide by), so

$$x_0+w>N\ \Longrightarrow\ \exists x\in\mathcal{R}:\ d_3(x)\ge r_3 ,$$

and $\rho(d_3)$ leaves its intended interval rather than the enumeration
stopping. $x_0+w\le N$ is a precondition, not a checked invariant.

---

## 2. The code

| symbol | name in the code |
|---|---|
| $r_0,r_1,r_2,r_3$ | `radices[0..3]` |
| $x_0,\ w$ | `base`, `n` — the first index of the range and its length |
| $d_0,d_1,d_2,d_3$ | `d0, d1, d2, d3` |
| $V,\beta,I$ | `voc, bvoc, isc` |
| $n,\tau,\rho$ | `n_s, tmin, rear` |
| $\nu,\kappa$ | `v`, `k_i` |
| $\ell_v,\ell_i$ | `v_max`, `i_max` |
| $h_\nu,h_\kappa$ | `hv`, `hi` |
| $\theta$ | `tight` |
| $c_0\ldots c_7$ | `tally[0..7]`, see `COUNTERS` |
| $c^{(t)}$ | `c_all, c_pass, c_fv, c_fi, c_fb, c_tv, c_ti, c_sv` |
| $c^{(\mathcal{P})}$ | the shared array `s[NC]` |

The decode, both paths:

```c
unsigned int d0, d1, d2, d3;
if (pow2) {                       /* mask and shift: an identity, not a trick */
    unsigned long long u = idx;
    d0 = (unsigned int)(u & (r0 - 1u)); u >>= k0;
    d1 = (unsigned int)(u & (r1 - 1u)); u >>= k1;
    d2 = (unsigned int)(u & (r2 - 1u)); u >>= k2;
    d3 = (unsigned int)u;
} else if (idx <= 0xFFFFFFFFULL) { /* 32-bit divmod: cheaper where it fits */
    unsigned int u = (unsigned int)idx;
    d0 = u % r0; u /= r0;
    d1 = u % r1; u /= r1;
    d2 = u % r2; u /= r2;
    d3 = u;
} else {                           /* 64-bit divmod: always correct */
    unsigned long long u = idx;
    d0 = (unsigned int)(u % r0); u /= r0;
    d1 = (unsigned int)(u % r1); u /= r1;
    d2 = (unsigned int)(u % r2); u /= r2;
    d3 = (unsigned int)u;
}
```

The point, the verdict and the private tally:

```c
const float voc  = bins[d0 * 3u + 0u];
const float bvoc = bins[d0 * 3u + 1u];
const float isc  = bins[d0 * 3u + 2u];

const float n_s  = fmaf((float)d1, ns_s, ns_b);
const float tmin = fmaf((float)d2, t_s,  t_b);
const float rear = fmaf((float)d3, bf_s, bf_b);

const float k_u    = fmaf(bvoc, tmin - 25.0f, 1.0f);
const float v      = n_s * voc * k_u;
const float k_i    = 1.25f * (isc * (1.0f + rear));

const bool ov = (v   > v_max);
const bool oi = (k_i > i_max);

c_all += 1u;
if      (ov & oi) { c_fb += 1u; }
else if (ov)      { c_fv += 1u; }
else if (oi)      { c_fi += 1u; }
else {
    c_pass += 1u;
    const float hv = 1.0f - v   * inv_v;
    const float hi = 1.0f - k_i * inv_i;
    if (hv < tight)  c_tv += 1u;
    if (hi < tight)  c_ti += 1u;
    if (hv > 0.10f)  c_sv += 1u;
}
```

The two-stage reduction, once per thread and once per block, never inside the
loop:

```c
if (c_all)  atomicAdd(&s[0], (unsigned long long)c_all);
/* ... one line per counter ... */
__syncthreads();
if (threadIdx.x < NC) {
    const unsigned long long acc = s[threadIdx.x];
    if (acc) atomicAdd((unsigned long long*)&tally[threadIdx.x], acc);
}
```

The grid is capped, so a single launch covers an arbitrarily large range:

```python
grid = int(min(max_blocks, (n + block - 1) // block))
fn((grid,), (block,), (d_bin, cp.int32(nbin), cp.uint64(base),
                       cp.uint64(n), radices, limits, tally))
```

### Where the code and the mathematics agree, and how that is shown

The module's own self-test recomputes the same enumeration on the host in double
precision and compares the first five counters, allowing a tolerance
proportional to the case count for cases that straddle a limit in single
precision but not in double. It runs the general decode path with radices that
are not powers of two, and it runs a range whose base exceeds the
thirty-two-bit boundary to force the wide decode. All three are checks that the
bijection $D$ is the same bijection on both sides.

The tolerance is the honest part: the two sides are *not* expected to agree
exactly, because the kernel rounds to single precision and the check does not. A
case whose voltage lands within a rounding step of the limit can be counted as
passing on one side and failing on the other. The tolerance says how many such
cases are being tolerated; it does not say which ones they were, and a real
defect small enough to hide under it would not be seen.

---

## 3. Plain English

Every combination this kernel examines is given a number, and the numbering is
the whole trick.

Think of a set of dials. One selects which class of equipment is being
considered. One selects how many units are wired into a chain. One selects the
assumed coldest temperature. One selects how much extra the back of the equipment
is assumed to collect. Each dial has a fixed set of positions. Turning the first
dial through all of its positions, then clicking the second dial on one notch and
going round the first again, and so on, will eventually visit every combination
exactly once, in a fixed order, and never revisit one.

That ordering means a single whole number can stand for a whole combination. Give
the machine a number and it can work out, with no lookup and no memory of what
came before, which position each dial is in. Give it a first number and a count
and it will cover exactly that stretch of combinations, with nothing skipped and
nothing done twice. This is what lets the work be split across many thousands of
workers at once without any of them needing to agree with any other about who is
doing what: each takes its own numbers, and the numbers cannot collide.

Recovering the dial positions from the number is ordinary division with
remainders. When the number of positions on a dial happens to be a power of two,
the same division can be done by chopping the number into bit-sized pieces
instead, which is faster and gives exactly the same answer, not an approximation
of it. The kernel checks once, for the whole run, which of those two cases it is
in, so the choice never varies from worker to worker and none of them has to wait
for the others while a decision is made.

Having recovered the dial positions, the worker computes two quantities: the
voltage the chain would reach at the assumed coldest temperature, and the current
it would be asked to carry. Each is compared with the limit the equipment is
rated for. That gives four possible outcomes: both within their limits, the
voltage over, the current over, or both over. For the combinations that are
within both limits, the worker also asks how much room is left, and notes
separately those with almost none, and those with a great deal, because a design
that only just fits and a design that wastes most of what it paid for are both
worth knowing about and neither shows up in a pass count.

Nothing is written down about any individual combination. The worker keeps a
handful of running counts in its own registers. When it has finished its share it
adds those counts into a small set shared by its immediate neighbours, and once
the neighbourhood has finished, one of them adds the neighbourhood's totals into
the grand totals. The reason for the two stages is contention: if every worker
added directly into the grand totals, they would spend their time queueing behind
one another rather than working.

The counts come out the same regardless of how the work was divided, how many
workers there were, or which finished first. That is a real guarantee and it
comes from the fact that whole numbers can be added in any order and grouped in
any way without changing the sum. It is worth noticing that the same guarantee
would not hold if the workers were adding up fractions, because fractions in a
machine are rounded, and rounding makes the order of addition matter.

There is one duty the kernel places on whoever calls it, and it is not checked.
The top dial has no dial above it, so nothing forces it to stop at its last
position. Ask for a stretch of numbers that runs off the end of the space and the
top dial will simply keep climbing past where it was meant to stop, quietly
producing combinations that were never part of the question. The kernel does not
guard against this, because a guard would cost a comparison on every single
combination, and there are a very great many of them. The caller must not ask.

Finally, the kernel checks itself against a version of the same enumeration
written a completely different way, in a form that carries more digits. They are
not expected to match exactly, and that expectation is itself the point: where a
combination sits almost exactly on a limit, the version carrying fewer digits can
land on one side and the version carrying more on the other. The check allows a
number of such disagreements and says how many. What it cannot do is tell a
rounding disagreement from a real mistake small enough to hide among them.
