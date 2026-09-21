# trillion-loop

A search that runs entirely inside a graphics card's own memory, keeps only the
candidates that survive it, and writes to the drive only when the set of
survivors changes.

Every algorithm here is documented in three layers, always in the same order:

1. **The mathematics.** Set-builder and standard notation. The space, the
   predicates, the survival rule, the fixed point. No code and no prose.
2. **The code.** The Python that implements exactly that mathematics, directly
   below it, short enough that each symbol above has one name below.
3. **Plain English.** No symbols and no numbers. Readable by someone who cannot
   read either layer above it, and still true.

Layer three is the one that matters. It is written last and carries the most
care, because a result nobody outside the room can check is not a result.

## The documents

| | |
|---|---|
| [The resident loop](docs/01-the-resident-loop.md) | the whole cycle: evaluate, annihilate, compact, breed, write only on change |
| [The enumeration kernel](docs/02-the-enumeration-kernel.md) | one launch over a mixed-radix case index, tallied in block-shared counters |
| [The second channel](docs/03-the-second-channel.md) | the same two quantities rebuilt along unrelated arithmetic routes, and a decode checked twice |
| [The provenance gate](docs/04-the-provenance-gate.md) | whether the limit a verdict rests on can be pointed at in a document |
| [Stratified adaptive sampling](docs/05-stratified-sampling.md) | steering effort at the boundary without biasing the headline |
| [Rule volume and the implication lattice](docs/06-rule-volume.md) | how much of the space a rule claims, exactly, by bisection on a monotone predicate |
| [Coverage and the chance of missing](docs/07-coverage.md) | what a count of none is worth, given how much was looked at |

## The four ideas the whole repository rests on

**Annihilation.** Every candidate is computed by two independent routes: one
factored, one expanded by the distributive law, which rounds twice where the
factored form rounds once. Agreement lets the candidate live. Disagreement
destroys it and writes nothing at all: no log line, no warning, no count sitting
in a file beside the word SUCCESS. A number is born as a pair or it is not born.

**Residency.** The card's own memory runs about a hundred times faster than the
drive beneath it. So evaluations live and die in the card's memory, and only
survivors are written. The drive is touched only when the survivor set changes.

**Survival selects for information, not for success.** A candidate lives if two
channels agree on it, if its verdict rests on a limit that traces to a document,
and if it sits within one step of changing its answer. Points deep inside a
settled region repeat what is already known; they are the expensive kind of
nothing, because they look like work.

**Provenance is part of the mathematics.** The survival predicate has a term for
whether a limit appears in any document. A verdict resting on an invented
constant is not knowledge, however cleanly it computes.

## The honest limit

Agreement between two implementations proves consistency, never truth. Both were
built here; both can be wrong in the same way, and a second channel that shares a
mistake with the first will confirm it forever.

A generation count is not a coverage fraction. Later generations deliberately
resample the boundary that earlier ones found, so the evaluations are
concentrated rather than spread, and the total number of evaluations divided by
the size of the space means nothing. The one document that is allowed to turn a
count into a fraction is [Coverage and the chance of
missing](docs/07-coverage.md), and it refuses to do so for a steered sample.

No site, client, project, location, maker or model is named anywhere in this
repository. Limits are described by the role they play - the voltage the
equipment is rated for, the current the protective device is rated for - never by
brand and never by the document they were copied from.
