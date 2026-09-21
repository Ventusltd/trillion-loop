# The provenance gate

Source: `t03_gate.py`

Every other check in this repository compares two things that were written here.
Two things written here can agree perfectly and both be invented. This gate does
not ask whether a verdict is consistent. It asks where the limit came from.

---

## 1. The mathematics

### The evidence base as a relation

Let $K$ be the finite set of canonical limit keys, totally ordered by name, and
let $\mathcal{A}$ be the set of aliases any channel might use. The audit is a
finite multiset of *sites*

$$\mathcal{E}\subset \mathcal{A}\times\mathbb{R}\times\{\mathsf{S},\mathsf{U}\}$$

where $\mathsf{S}$ marks a value present in the evidence base and $\mathsf{U}$ a
value absent from it. Let

$$\kappa:\ \mathcal{A}\longrightarrow K\cup\{\bot\}$$

be the canonicalisation, total on $\mathcal{A}$ and sending anything it does not
recognise to $\bot$. For $k\in K$ define the two value sets

$$S_k=\{\,v\ :\ \exists(a,v,\mathsf{S})\in\mathcal{E},\ \kappa(a)=k\,\},\qquad
U_k=\{\,v\ :\ \exists(a,v,\mathsf{U})\in\mathcal{E},\ \kappa(a)=k\,\}.$$

$S_k$ and $U_k$ need not be disjoint as *names*: a key is **contaminated** when

$$S_k\neq\varnothing\ \wedge\ U_k\neq\varnothing ,$$

which is the case one name carrying two different numbers, one traceable and one
not.

### The second document

Let $\varsigma:K\to\Sigma\cup\{\bot\}$ record what the normative reading says
about each key, with

$$\Sigma=\{\ \mathsf{from\text{-}standard},\ \mathsf{absent},\ \mathsf{not\text{-}a\text{-}standard\text{-}limit},\ \mathsf{recorded\text{-}not\text{-}a\text{-}requirement}\ \}$$

and $\bot$ meaning the documents were never asked about this key.

### Classification

$$\mathrm{ok}(k,v)\ \iff\ \underbrace{v\in S_k}_{\text{present}}\ \wedge\ \underbrace{\varsigma(k)\neq\mathsf{absent}}_{\text{not recorded missing}}$$

Two asymmetries are deliberate and carry the whole design.

$$v\notin S_k\ \Longrightarrow\ \neg\,\mathrm{ok}(k,v)
\qquad\text{(absence of a finding is not a source)}$$

$$\kappa(a)=\bot\ \Longrightarrow\ \neg\,\mathrm{ok}
\qquad\text{(an unknown limit is unsourced, never safe)}$$

Neither follows from $v\in U_k$. The default is refusal, and $U_k$ only supplies
the explanation.

### Encoding a declaration

Index $K$ by its order, $K=\{k_0<k_1<\dots<k_{m-1}\}$, and define
$\beta:K\to\mathbb{N}$, $\beta(k_i)=2^{i}$, together with a distinguished
**poison** bit $\Omega=2^{62}\notin\beta(K)$.

A cycle declares the pairs it rested on, $L=\{(a_1,v_1),\dots,(a_j,v_j)\}$. Its
mask is

$$B(L)=\begin{cases}
-1, & \exists i:\ \kappa(a_i)=\bot\\[4pt]
\displaystyle\bigvee_{i}\beta\bigl(\kappa(a_i)\bigr)\ \vee\
\bigvee_{i\,:\,v_i\neq\bot\,\wedge\,\neg\mathrm{ok}(\kappa(a_i),v_i)}\Omega, & \text{otherwise.}
\end{cases}$$

### The defensible set and the publication rule

$$\Delta=\bigvee_{k\,:\,\mathrm{ok}(k)}\beta(k)\ \in\mathbb{N},\qquad \Omega\not\le\Delta$$

$$\boxed{\ \mathrm{pub}(m)\ \iff\ m\ge 0\ \wedge\ m\neq 0\ \wedge\ \bigl(m\ \wedge\ \neg\Delta\bigr)=0\ }$$

reading $\wedge,\vee,\neg$ as bitwise operations, that is, as meet, join and
complement in the Boolean lattice $2^{K\cup\{\Omega\}}$. Equivalently, with
$\sqsubseteq$ the inclusion order on that lattice,

$$\mathrm{pub}(m)\iff \varnothing\neq m\sqsubseteq\Delta .$$

Two immediate consequences:

**Downward closure.** $m\sqsubseteq m'$ and $\mathrm{pub}(m')$ and $m\neq0$
imply $\mathrm{pub}(m)$. Resting on fewer limits never costs publishability.

**Poison is absorbing.** $\Omega\le m\Rightarrow\neg\mathrm{pub}(m)$, since
$\Omega\not\le\Delta$. A declared value that fails classification cannot be lost
by any subsequent join.

**The empty declaration is refused.** $m=0$ is excluded explicitly: a verdict
resting on no declared limit is not a verdict about the world.

### Refusal as a bottom element

If the audit itself is missing, set $\Delta=0$. Then $\mathrm{pub}(m)$ is false
for every $m$, and the state is announced rather than inferred. The gate does not
assume a limit is sourced because it cannot see the file that would say so.

### Cost

Applied to a device array of masks, the test is

$$\mathbf{1}\bigl[m\neq0\bigr]\ \wedge\ \mathbf{1}\bigl[(m\ \&\ {\sim}\Delta)=0\bigr]$$

— two integer operations per cycle, independent of the number of keys. There is
therefore no performance argument for skipping it.

### What the gate cannot do

$$\mathrm{ok}(k,v)\ \not\Longrightarrow\ v\ \text{is correct},\qquad
\mathrm{ok}(k,v)\ \not\Longrightarrow\ v\ \text{is used correctly}.$$

It separates *present* from *absent*, never *right* from *wrong*.

---

## 2. The code

| symbol | name in the code |
|---|---|
| $\mathcal{E}$ | the `findings` list of the provenance audit |
| $\kappa$ | `canonical()`, backed by the `_ALIAS` table |
| $K$ | `LIMIT_ORDER` |
| $S_k,\ U_k$ | the `sourced` and `unsourced` dicts from `_gather()` |
| contamination | the `CONTAMINATED:` clause in `classify()`'s third return value |
| $\varsigma$ | `_STD[key]["status"]`, built from `_STD_ANCHORS` |
| $\mathrm{ok}$ | the second element of `classify()` |
| $\beta(k)$ | `_BIT[k]` |
| $\Omega$ | `_POISON` |
| $B(L)$ | `bits_for()` |
| $\Delta$ | `defensible_bits()` |
| $\mathrm{pub}$ | `publishable()` / `publishable_mask()` |
| refusal | `_refuse()`, `gate_state() == "REFUSED"` |

The classification, with both asymmetries visible:

```python
d, sourced_v, unsourced_v, sites = _gather(key)
std        = _STD.get(d["std"]) if d["std"] else None
std_absent = bool(std and std["status"] == "ABSENT")

asked = None if value is None else float(value)
if asked is not None:
    verdict = "SOURCED" if asked in sourced_v else "UNSOURCED"
    ev = (sourced_v.get(asked) or unsourced_v.get(asked) or {}).get(
        "evidence",
        "this exact value has no limit site in the audit - absence of "
        "a finding is not a source")
...
ok = (verdict == "SOURCED") and not std_absent
```

An unrecognised name is not an error to be skipped; it is an unsourced limit:

```python
key = canonical(limit_name)
if key is None:
    return (None, False,
            "UNKNOWN LIMIT %r - not one of the %d limits this gate knows. "
            "An unknown limit is treated as UNSOURCED, never as safe."
            % (limit_name, len(LIMIT_ORDER)))
```

The encoding, and the absorbing poison bit:

```python
def bits_for(limits):
    m = 0
    for nm, val in _pairs(limits):
        k = canonical(nm)
        if k is None:
            return -1                      # unknown limit: never publishable
        m |= _BIT[k]
        if val is not None and not classify(k, val)[1]:
            m |= _POISON                   # cannot be lost by a later join
    return m
```

The rule itself, in one line each for the host and the device:

```python
def publishable(verdict_bits) -> bool:
    ...
    if bits < 0:   return False
    if bits == 0:  return False            # rested on nothing: not a verdict
    return (bits & ~defensible_bits()) == 0

def publishable_mask(bits_dev):
    a = _XP.asarray(bits_dev)
    if _REFUSAL is not None:
        return _XP.zeros(a.shape, dtype=bool)
    ok = defensible_bits()
    return (a != 0) & ((a & ~_XP.asarray(ok, dtype=a.dtype)) == 0)
```

Refusal, announced rather than inferred:

```python
if not os.path.isfile(PROVENANCE_JSON):
    _refuse("provenance.json NOT FOUND at: %s\n"
            "The audit of where every limit came from has not been run, or "
            "its output has moved." % PROVENANCE_JSON)
    return False
```

Both evidence files are located through environment variables, so nothing in
this description depends on where any particular machine keeps them.

### Where the code and the mathematics differ

1. **Classification by name alone is not the same function as classification by
   name and value, and the two are both exposed.** Asked for a key without a
   value, the code returns the sourced value *if there is exactly one*, and
   otherwise falls back to a declared default or to nothing. That is a
   convenience, and it is the branch that a caller reaching for the short form
   will hit. The mathematics above defines $\mathrm{ok}$ only on pairs, because
   the pair is the only form that is safe on a contaminated key. Where a key is
   contaminated, the one-argument form reports the traceable value and the
   contamination note together; a caller that reads only the boolean gets the
   answer for a value its cycle may not have used.

2. **Two degrees of severity are collapsed into one.** The normative index
   distinguishes a key the documents record as absent from a key that is merely
   not a limit in those documents at all, and from one recorded explicitly as
   not being a requirement. Only the first of these blocks publication. The
   other two are carried into the explanatory text and have no effect on
   $\mathrm{ok}$. The mathematics reflects the code rather than the intent: the
   status alphabet has four members and exactly one of them is load-bearing.

3. **A missing normative reading disables half the gate without refusing.** If
   the clause-cited document cannot be opened, every absence check is switched
   off, a warning is written to the error stream, and publication continues on
   the audit alone. The audit's absence is a refusal; the normative reading's
   absence is a warning. That asymmetry is deliberate in the source and is worth
   naming, because it means the printed state can read as armed while one of the
   two gates is not running.

---

## 3. Plain English

A number that everybody in a trade uses, that appears in no document anybody can
produce, is not knowledge. It is a habit.

This part of the system does one job and refuses to do any other. It does not ask
whether a calculation is right. It does not ask whether two independent routes
agreed. It asks a single question about every limit that a conclusion leans on:
can somebody point at where this number came from.

To answer it, the gate opens two files, and neither of them was written by a
program. One is an audit: somebody went through the evidence base, took every
number being used as a limit anywhere in the estate, and recorded whether that
number actually appears in an equipment sheet or whether it does not appear
anywhere at all. The other is a reading of the governing documents: somebody read
them and wrote down, clause by clause, what they do and do not say, including the
cases where a number in common use is recorded as being absent from the documents
entirely.

The audit turned up something worth stating. Two of the limits in everyday use
here can be pointed at: the voltage the equipment is rated for, and the current
its protective device is rated for. Both come off the sheets. Two others cannot
be pointed at anywhere. They are simply numbers that circulate. And one limit
name was found carrying two different values at once across its many uses: one of
them traceable to a sheet, the other traceable to nothing. If the gate judged by
name alone it would either publish the untraceable value or bury the traceable
one, so it judges by name and value together, and a name used both ways is
flagged as contaminated, with a demand that the conclusion say which number it
actually used.

Two rules govern the whole thing, and both point the same way.

The first is that absence of a finding is not a source. If the audit has nothing
to say about a particular number, that number is treated as untraceable. It is
not given the benefit of the doubt, because the doubt is the whole point: a
number nobody recorded is exactly the kind of number that gets invented and then
repeated until it feels official.

The second is that an unknown limit is unsourced, never safe. If a calculation
declares that it rested on a limit the gate has never heard of, the gate does not
shrug and let it through. An unrecognised name is treated as the worst case and
the conclusion is held back.

The mechanics are deliberately trivial. Every limit the gate knows is given a
place in a row of switches. A conclusion declares which limits it rested on by
flipping the switches for those limits. The gate keeps a second row showing which
limits are currently traceable. A conclusion may be published if it flipped at
least one switch, and if every switch it flipped is also up in the traceable row.
There is an extra switch reserved for poison, which is flipped when a conclusion
declares a specific value that fails the check, and which is never up in the
traceable row, so a poisoned conclusion can never be rescued by combining it with
anything else. A conclusion that flipped no switches at all is refused too, on
the grounds that something resting on no stated limit is not a statement about the
world.

The cost of this test is two operations, whatever the conclusion is and however
many limits exist. That matters more than it sounds. A gate that was expensive
would be skipped when the machine was busy, which is exactly when it is needed;
a gate that costs nothing has no excuse.

If the audit itself is missing, the gate does not guess. It refuses, loudly,
and everything produced while it is refusing is counted and kept but may not
leave the machine as a claim about anything. Not publishing is a state the system
is willing to be in.

Failing this gate is not a failure of the arithmetic. The numbers are real and
the calculation happened. What is withheld is the right to describe the result as
a fact about equipment in the world, rather than as a fact about the number that
was fed in.

Three things this gate cannot do, and they should be said with the same weight as
what it can.

It can tell a number that is present from a number that is absent. It cannot tell
a right number from a wrong one. A rating copied off the wrong page of the right
sheet passes here without a murmur.

It does not read the governing documents. It trusts one person's reading of them,
written down by hand. A misreading is repeated faithfully, not caught.

And it says nothing about whether a traceable limit was used correctly. A
perfectly sourced protective rating compared against the wrong current gives a
wrong answer that this gate will happily publish.
