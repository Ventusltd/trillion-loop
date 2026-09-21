# Coverage and the chance of missing

Source: `t13_coverage.py`

A count is not a finding. A count of none is a finding only if the search that
produced it had a real chance of finding something. This module attaches that
chance to every count, and refuses to attach it to counts that have not earned
one.

---

## 1. The mathematics

### The setting

A space $\mathcal{X}$ with $N=\lvert\mathcal{X}\rvert$, an unknown offending set
$\mathcal{O}\subseteq\mathcal{X}$ with $k=\lvert\mathcal{O}\rvert$, and a sample
of $n$ cases. The reported count is $\lvert\text{sample}\cap\mathcal{O}\rvert$.

### Miss probability, distinct draws

If the $n$ examined cases are distinct and every $n$-subset is equally likely,
the number of offenders seen is hypergeometric and

$$\Pr\bigl[\text{see }0\ \bigm|\ k\bigr]
=\frac{\binom{N-k}{n}}{\binom{N}{n}}
=\prod_{i=0}^{k-1}\frac{N-n-i}{N-i} .$$

Two special cases are worth having in the hand:

$$k=1:\qquad \Pr[\text{miss}]=\frac{N-n}{N}=1-f,\qquad f=\frac{n}{N}$$

$$n>N-k:\qquad \Pr[\text{miss}]=0 \quad\text{(pigeonhole).}$$

### Miss probability, independent draws

If the $n$ draws are independent and uniform, each misses $\mathcal{O}$ with
probability $1-k/N$, so

$$\Pr\bigl[\text{see }0\bigm|k\bigr]=\Bigl(1-\frac{k}{N}\Bigr)^{n}
=\exp\Bigl(n\log\bigl(1-\tfrac{k}{N}\bigr)\Bigr).$$

### Effective coverage

Distinct cases touch $n$ points. Independent draws revisit, so

$$f_{\text{eff}}(D)=1-\Bigl(1-\frac{1}{N}\Bigr)^{D}$$

is the expected fraction of $\mathcal{X}$ actually touched by $D$ independent
draws. In particular $D=N$ gives $f_{\text{eff}}\to1-e^{-1}$: drawing as many
times as there are cases leaves better than a third of the space untouched.

### Detection floor and the interval for a count of none

$\Pr[\text{see }0\mid k]$ is nonincreasing in $k$, so for a confidence level
$1-\alpha$ the two natural quantities are well defined and found by bisection:

$$k_\alpha=\min\Bigl\{k\ :\ \Pr[\text{see }0\mid k]\le\alpha\Bigr\}
\qquad\text{(the thinnest band this sweep would catch)}$$

$$k^{+}=k_\alpha-1=\max\Bigl\{k\ :\ \Pr[\text{see }0\mid k]>\alpha\Bigr\}
\qquad\text{(the most offenders still consistent with seeing none)}$$

$k^{+}$ is the exact finite-population generalisation of the rule of three; the
familiar $3/n$ bound on a rate is its large-$N$ limit.

### The three verdicts on a count of none

$$\Pr[\text{miss}\mid k{=}1]\ \begin{cases}
\ge\tfrac12 & \text{not evidence of absence}\\
\in[\alpha,\tfrac12) & \text{weak evidence of absence}\\
<\alpha & \text{evidence of absence at the stated width}
\end{cases}$$

### A steered sample has no coverage

If the sample was not uniform on $\mathcal{X}$, then $f$, $\Pr[\text{miss}]$ and
every quantity above are undefined, because they are all functions of a design
that was not the design used. The module returns them as absent rather than as
numbers. Formally the refusal is:

$$\text{uniform}=\text{false}\ \wedge\ \pi\ \text{unknown}
\ \Longrightarrow\ f:=\bot,\quad \Pr[\text{miss}]:=\bot,\quad \text{reportable}:=\text{false}.$$

### With inclusion probabilities, an estimate and not an enumeration

If each examined case $i$ carries a known inclusion probability
$\pi_i\in(0,1]$, the Horvitz–Thompson estimator of the offender total is

$$\hat T_{\mathrm{HT}}=\sum_{i\,\in\,\text{hits}}\frac{1}{\pi_i},
\qquad \mathbb{E}\bigl[\hat T_{\mathrm{HT}}\bigr]=k .$$

It is an estimate carrying sampler variance. It is never an enumeration, and the
word *verified* does not apply to it.

### Bias is contagious

For records $r_1,\dots,r_j$,

$$\text{sample}\bigl(\textstyle\bigoplus_i r_i\bigr)=
\begin{cases}\text{uniform}, & \forall i:\ \text{sample}(r_i)=\text{uniform}\ \wedge\ N_i\ \text{all equal}\\
\text{biased}, & \exists i:\ \text{sample}(r_i)\neq\text{uniform}\end{cases}$$

and combining uniform tallies swept over *different* spaces raises an error
rather than returning a number: they are not the same measurement, and
$\Pr[\text{miss}]$ has no meaning for their sum.

### Discrimination: whether a filter has teeth

Separately from coverage, a filter is characterised by its pass rate on real
cases and on random decoys:

$$\hat r_{\mathrm{real}}=\frac{a}{n_{\mathrm{real}}},\qquad
\hat r_{\mathrm{decoy}}=\frac{b}{n_{\mathrm{decoy}}},\qquad
\text{lift}=\frac{\hat r_{\mathrm{real}}}{\hat r_{\mathrm{decoy}}}$$

$$\text{informative}\iff \hat r_{\mathrm{decoy}}<\tfrac12\,\hat r_{\mathrm{real}}
\ \wedge\ n_{\mathrm{decoy}}\ \text{large enough}$$

$$b=0\ \Longrightarrow\ \hat r_{\mathrm{decoy}}\le \frac{3}{n_{\mathrm{decoy}}}\ \text{at }95\%\ \text{(rule of three)} .$$

A filter with $\hat r_{\mathrm{decoy}}\approx\hat r_{\mathrm{real}}$ carries zero
information about its input. Its output is the expected output of a filter that
accepts everything, whatever its pass count says.

---

## 2. The code

| symbol | name in the code |
|---|---|
| $N,\ n,\ k$ | `space`, `examined`, `k` |
| $\Pr[\text{see }0\mid k]$ | `miss_probability()` |
| $f,\ f_{\text{eff}}$ | `coverage_fraction`, `effective_coverage()` |
| $k_\alpha$ | `detection_floor()` |
| $k^{+}$ | `rule_of_three()` |
| $\bot$ | `None` in the returned record |
| reportable | `rec["reportable"]`, enforced by `require_reportable()` |
| $\pi_i$ | `inclusion_probs` |
| $\hat T_{\mathrm{HT}}$ | `rec["estimate"]` |
| lift, informative | `discrimination()` |
| $\bigoplus$ | `combine()` |

Both draw models, in one function:

```python
def miss_probability(examined, space, k=1, with_replacement=False):
    if k <= 0:      return 1.0
    if k >= space:  return 0.0
    if with_replacement:
        if examined <= 0: return 1.0
        return math.exp(int(examined) * math.log1p(-k / space))
    n = min(int(examined), space)
    if n <= 0:          return 1.0
    if n > space - k:   return 0.0          # pigeonhole
    p = 1.0
    for i in range(k):
        p *= (space - n - i) / (space - i)
        if p == 0.0: break
    return p
```

The two thresholds, by bisection on a monotone function:

```python
def detection_floor(examined, space, alpha=0.05, with_replacement=False):
    lo, hi = 1, int(space)
    if miss_probability(examined, space, hi, with_replacement) > alpha:
        return space
    while lo < hi:
        mid = (lo + hi) // 2
        if miss_probability(examined, space, mid, with_replacement) <= alpha:
            hi = mid
        else:
            lo = mid + 1
    return lo

def rule_of_three(examined, space, with_replacement=False):
    return max(0, detection_floor(examined, space, 0.05, with_replacement) - 1)
```

The refusal, which is the point of the module:

```python
if inclusion_probs is None:
    rec["sample"] = "adaptive-biased"
    rec["sentence"] = (
        "{0}: {1:,} found in {2:,} ADAPTIVELY SAMPLED cases. This is a raw "
        "tally over a steered subset, not a measurement of the {3:,}-case "
        "space. It has no coverage fraction and no miss probability, and it "
        "MUST NOT be quoted as a rate, a percentage, a fraction of the "
        "space, a coverage figure or a clean bill of health. ..."
        .format(label, count, examined, space))
    return rec
```

The gate that stops a steered tally leaving as a fraction:

```python
def require_reportable(record):
    if not record.get("reportable", False):
        raise QualificationError(...)
    return record
```

Contagion, and the refusal to add unlike measurements:

```python
biased = [r for r in records if r.get("sample") != "uniform"]
if biased:
    out = qualify(total, examined, space, uniform=False, label=label)
    out["sentence"] = "MIXED SAMPLE. " + out["sentence"] + " ..."
    return out
spaces = set(int(r["space"]) for r in records)
if len(spaces) != 1:
    raise QualificationError("cannot add uniform tallies swept over different "
                             "spaces ({0}) -- they are not the same "
                             "measurement".format(sorted(spaces)))
```

The self-test reproduces the two failures the module exists to prevent: a count
of none whose chance of having missed a lone offender was better than
nine-hundred-and-ninety-nine in a thousand, and a filter that absorbed every
random decoy handed to it and reported a clean reading.

### Where the code and the mathematics differ

1. **The weighted miss probability is not derived.** With inclusion
   probabilities supplied, the code computes
   $\exp\bigl(\tfrac1N\sum_i\log(1-\pi_i)\bigr)$ and calls it the weighted
   chance of missing a lone offender. That expression is a plausible
   construction — it decreases as the reach of the sample grows, and it tends to
   one as the space grows — but it is not $\Pr[\text{miss}]$ under any sampling
   design stated here. The honest reading is that it is an index of reach, not a
   probability, and the sentence built from it should not use the word chance.
   This is the one place in the module where a number is offered without a
   design behind it.

2. **The estimator is not the Horvitz–Thompson estimator.** The code computes
   $\hat T=\text{count}\times\overline{(1/\pi_i)}$, the count times the *mean*
   design weight, rather than the sum of the weights of the cases that actually
   hit. The two agree only when the hits' weights are typical of the sample. The
   source says so in a comment — the identity of the hits is not available at
   this interface — and the substitution is reasonable there, but the record is
   labelled with the name of an estimator it is not computing.

3. **Coverage of a with-replacement sweep is not $n/N$.** The docstring's summary
   line gives the coverage fraction as examined over space; the code correctly
   uses the effective coverage, which accounts for revisits, whenever the draw
   model is independent. The code is right and the summary line is loose. It
   matters, because the difference is large exactly in the regime this module was
   written for.

4. **The discrimination threshold is a convention.** Requiring the decoy pass
   rate to be below half the real pass rate, with a minimum decoy count, is a
   stated rule of thumb rather than a test with a stated size and power. It is
   reported as a boolean, which reads more decisive than it is.

---

## 3. Plain English

Finding nothing is not the same as there being nothing. Whether it is even
evidence depends entirely on how much of the space was looked at, and that is a
question most reports never ask.

Consider a search that examined a couple of hundred thousand cases out of a space
of several hundred million, and found no disagreements anywhere. Reported as a
clean verdict, this sounds like a result. Worked out properly, a single offending
case sitting somewhere in that space would have been missed almost every time.
The honest sentence is not that the thing was verified. It is that none were
found, and that a lone offender would have gone unseen with near certainty. Those
two sentences describe the same run and they lead to opposite decisions.

This module computes the missing sentence and attaches it to every count.

Its central quantity is the chance of seeing nothing when something is actually
there. That chance depends on three things: how big the space is, how many cases
were examined, and how widespread the thing being hunted would be if it existed.
For a lone offender and a search that examined each case once, the answer is
simple and unforgiving: the chance of missing it is exactly the fraction of the
space that was not looked at. Look at a hundredth of the space and you will miss
a lone offender ninety-nine times in a hundred.

There is a second version for searches that draw cases at random rather than
working through them in order, and the difference is worth understanding. Random
draws revisit. Draw as many times as there are cases in the space and you will
still have left better than a third of it untouched, because the draws pile up on
cases already seen. A count of draws is not a count of cases examined, and
treating one as the other overstates coverage by a wide margin.

From the same arithmetic come two more useful numbers. The first is the thinnest
band of offending cases this search would have been likely to catch. That is a
sensitivity: it says what the search was actually capable of noticing, as opposed
to what it was hoping to notice. The second is the largest number of offenders
that remains consistent with having found none. That is an upper bound, and it is
the number that belongs in any sentence beginning with the word verified.

Depending on where the chance of missing falls, the module issues one of three
verdicts, and it will not be talked out of any of them. If a lone offender would
have been missed more often than not, the result is not evidence of absence at
all. If it would have been missed occasionally, the result is weak evidence. Only
when it would almost certainly have been caught does the result become real
evidence of absence, and even then only at the stated width.

Now the part that took a failure to learn. All of the above assumes every case in
the space had an equal chance of being examined. A search that deliberately
steers towards the interesting part of the space, which is exactly what the
controller described elsewhere in this repository does, breaks that assumption
completely. The count it produces is a real count of real cases, but it is not a
fraction of anything. Both a high count and a count of none are expected
behaviours of a searchlight: the high count because the light was pointed where
the things are, the zero because it was pointed elsewhere.

So for a steered sample, this module does not compute a coverage, does not
compute a chance of missing, and does not return a smaller or more cautious
number in their place. It returns nothing where those numbers would go, marks the
record as not reportable, and writes a sentence saying in plain words that the
figure may not be quoted as a rate, a percentage, a fraction of the space, a
coverage or a clean bill of health. There is a gate that raises an error rather
than let such a record through as a fraction, and the gate exists because a
polite note in a comment would be ignored on a busy night.

The only way back is to supply, for each case examined, how likely it was to have
been examined at all. Given those, a genuine estimate of how many offenders exist
can be made by weighting each finding by how unlikely it was to have been looked
at. That is an estimate carrying sampling uncertainty. It is not an enumeration,
and the word verified does not apply to it.

Bias spreads. Adding a steered tally to a clean one gives a steered total, not a
partly clean one, and the module says so on the combined record. Adding two clean
tallies taken over different spaces raises an error instead of returning a
number, because they are not the same measurement and their sum has no meaning.

There is a last check, aimed at a different failure entirely. A filter that
accepts everything handed to it will report a clean pass rate, and that reading
is worthless. The way to find out is to feed it decoys: cases made up at random,
which it ought to reject. If it accepts the decoys at roughly the rate it accepts
real cases, it is not discriminating, its output is precisely the output expected
of a filter that does nothing, and no count derived from it is evidence of
anything at all, however impressive the count. If it rejects almost all the
decoys while accepting the real cases, it has teeth, and its counts carry
information. If too few decoys were tried to tell either way, that is reported
too, because a filter nobody tested is a filter nobody tested.

Two honest weaknesses to finish on.

When inclusion probabilities are supplied, the figure the module offers as the
weighted chance of missing something is not derived from any stated design. It
behaves sensibly, moving the right way as the search reaches further, but it is
an index of reach rather than a probability, and it should not be read as one.

And this module cannot tell whether the space it was handed is the honest one. A
search can cover the whole of a space that was chosen, deliberately or
accidentally, to contain nothing worth finding, and every number here will come
back immaculate. What the space was, and whether the thing being tested was the
right thing, are the caller's to defend. This module only measures the chance of
missing, given a space and a manner of looking that somebody else stated.
