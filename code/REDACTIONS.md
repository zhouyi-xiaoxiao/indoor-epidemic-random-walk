# Redacted copies of the protocol records

The protocols of the five tests and their post-hoc logs were frozen by SHA-256 hashes of the files on the
author's machine (freeze records `PREREG_FREEZE.json` and similar). The files of this repository are
redacted copies of those records: passages that do not concern the analyses of this article are omitted
(marked `[...]`) or reworded, identifiers are renamed as in the rest of the repository, and the words
for the organisation of the work are replaced as in the other notes (throughout the post-hoc logs, and in
single lines of the protocols where the list below says so). Nothing else was changed. The hashes in the freeze records therefore refer to the unredacted originals, which
are kept by the author and can be provided to the editor and referees; they cannot be checked against
the copies below. The timing evidence for the protocols, which was self-attested in any case, is
weaker for it. The list gives, for each redacted file, the SHA-256 hash of the original, that of the
released copy, the lines of the released copy that differ from the original, and the kinds of change.

## `code/bsc_outbreaks/PROTOCOL_PROPOSAL.md`

* SHA-256 of the original: `61484ecb92b67debab98d6dfb4a4ea90d1b8bd6ec6dff7975eb6f1f53cacfd61` (section 5 of the working notes of the outbreak dataset)
* SHA-256 of this copy: `d2729ff6a39be6c977020be1ff1469376b217cbbf59379a1dc82876098bbfb8e`
* lines changed in this copy: 16, 19, 22, 33, 35-38, 42, 94, 130
* changes:
  * reworded; the copy reads: "This is the per-cell hazard lambda(r)"
  * reworded; the copy reads: "which the spectral decomposition [...] delivers in closed form."
  * reworded; the copy reads: "(the SIR model makes new cases infectious immediately)"
  * reworded; the copy reads: "[...] The spatial model has to beat it on the same data."
  * reworded; the copy reads: "- **M1, same-site model.**"
  * reworded; the copy reads: "fixed at the default values for the mapped scene. Free parameters: ..."
  * reworded; the copy reads: "- **M2, framework with a non-local term.**"
  * reworded; the copy reads: "The propagator and spectral results [...] apply unchanged"
  * reworded; the copy reads: "A statement of agreement with documented outbreaks is supported onl..."
  * passage omitted, marked [...]

## `code/bsc_validation/POSTHOC_LOG.md`

* SHA-256 of the original: `6c7688087fbeb01cab7b95da1c0d878c81462f20c8b6d50c4ffe4b71c169e451`
* SHA-256 of this copy: `e4bc03f70878266e03af4eca479f5401dd5db6f38ee5585fb6d78b2c7241fdcc`
* lines changed in this copy: 11, 15, 17
* changes:
  * reworded; the copy reads: "the exact per-cell rule cannot reach"
  * reworded; the copy reads: "occupancy factor E[1/N] of the per-cell rule"
  * words for the organisation of the work replaced as in the other notes of the release ("re-check" for the second pass, "stretch of work" for an interrupted or resumed working period, the names of the analyses for their labels)

## `code/bsc_validation/PREREGISTRATION.md`

* SHA-256 of the original: `632db9b7b00796a4694d54ef533c3a3a89a3779fe18381a6a456d942b3a3f585`
* SHA-256 of this copy: `eac138ec5decd183ee626113d9a54ffd840546c9cc6c042f6b80b5a26e3b5a59`
* lines changed in this copy: 3, 16, 27-28, 83, 85, 91, 96, 100, 219-220, 251-252, 257-258
* changes:
  * passage omitted, marked [...]
  * passage omitted, marked [...]
  * reworded; the copy reads: "does model M1, or model M2 with the non-local term, predict"
  * reworded; the copy reads: "### 3.2 M1, same-site model"
  * reworded; the copy reads: "(hazard lambda"
  * reworded; the copy reads: "by the cosine-mode expansion [...]."
  * reworded; the copy reads: "per-cell rule, only the index infectious"
  * reworded; the copy reads: "### 3.3 M2, framework with the non-local term"
  * reworded; the copy reads: "the 10-400 quanta/h range [...] read as a 95% range"
  * reworded; the copy reads: "If only M2 passes, the non-local term must be part of the model bef..."
  * reworded; the copy reads: "in the sense of a heterogeneous D)"
  * reworded; the copy reads: "The re-check of the dataset raised binding objections."

## `code/bsc_validation/PREREG_ADDENDUM_A1_power.md`

* SHA-256 of the original: `08a68d93426913e507ce1dd09f37d455b307a4d7f2dabf644a0862ccbafd29a2`
* SHA-256 of this copy: `7fb166ede666b09ad9f5e8e3a359661f8fd59e51bf31aaff28bc7e2b26fb8e3d`
* lines changed in this copy: 62, 72, 74
* changes:
  * reworded; the copy reads: "scene default floor 50 m x 30 m"
  * reworded; the copy reads: "metro scene 22.5 m x 3 m"
  * reworded; the copy reads: "The re-check of the outbreak dataset computed"

## `code/bsc_validation/src/valmod/exact.py`

* SHA-256 of the original: `4d365d1503ff5898237023c73364a1e3aa503833a17c7ba0d04a5cc0625036d6`
* SHA-256 of this copy: `91d1df0f4d223f4d024f33a9f325755bddf8ac16dd24e5bef3815480d3eb83ac`
* lines changed in this copy: 3, 61
* changes:
  * reworded; the copy reads: "Runs the per-cell rule (same site,"
  * identifiers renamed as in the rest of the release (per-cell rule `percell`, two-zone mobility variant `zoneD`, mapping fields `model_scene`, `scene_parameters`)

## `code/bsc_validation/src/valmod/lattice.py`

* SHA-256 of the original: `25369ff046c7cb403422dde3cd8b369282a6fe6d305dcfaa0b05f3869049fc0f`
* SHA-256 of this copy: `ce747f711c9a560885976d66174c726dbd85dffd5c9d0660ba8a38557652fe8f`
* lines changed in this copy: 12
* changes:
  * reworded; the copy reads: "(the cosine-mode expansion)."

## `code/bsc_validation2/a1_more_outbreaks/POSTHOC_LOG.md`

* SHA-256 of the original: `4640073f9256d16d6176742d086b4e09e90e9424484a93e2a42127cd5efa3430`
* SHA-256 of this copy: `23a17d3e8f862f0d345ae246fc1aeaf90a55f6addaea3e8aa658a666fbd8dd7a`
* lines changed in this copy: 1, 12-13, 15
* changes:
  * words for the organisation of the work replaced as in the other notes of the release ("re-check" for the second pass, "stretch of work" for an interrupted or resumed working period, the names of the analyses for their labels)

## `code/bsc_validation2/a1_more_outbreaks/PREREGISTRATION.md`

* SHA-256 of the original: `33f034d580b002057ebc2e13fddc42452971a6e71b7df3a37b4c74c8caed182e`
* SHA-256 of this copy: `687686d17af34fc9400dddf2687537dd506336b67769e656ad07646b34fde217`
* lines changed in this copy: 3, 26, 100, 104, 108
* changes:
  * passage omitted, marked [...]
  * reworded; the copy reads: "Does the non-local airborne term (model M2:"
  * reworded; the copy reads: "- **M2, framework with the non-local term**:"
  * reworded; the copy reads: "- **M1cf, same-site model** (people walk,"
  * reworded; the copy reads: "(not model M1, not part of the decision rule for M2)"

## `code/bsc_validation2/a1_more_outbreaks/src/a1/lattice.py`

* SHA-256 of the original: `25369ff046c7cb403422dde3cd8b369282a6fe6d305dcfaa0b05f3869049fc0f`
* SHA-256 of this copy: `ce747f711c9a560885976d66174c726dbd85dffd5c9d0660ba8a38557652fe8f`
* lines changed in this copy: 12
* changes:
  * reworded; the copy reads: "(the cosine-mode expansion)."

## `code/bsc_validation2/a2_contact_mobility/POSTHOC_LOG.md`

* SHA-256 of the original: `2e8730c01eb5bc395747add53ad0ad0d1446124649c34e760aed72e3a7088a4a`
* SHA-256 of this copy: `b87bbbdd55386121c8d96f834dca8be03929befa181ff13891640fe1684a552a`
* lines changed in this copy: 1, 23, 27, 30, 36, 43, 56
* changes:
  * reworded; the copy reads: "use the harmonic-mean rule. This is exactly the model's"
  * reworded; the copy reads: "outside the model's formulation, recorded as such"
  * reworded; the copy reads: "prefer Z (inside the model's framework)"
  * reworded; the copy reads: "RWstickZ (the model's own D(r) mechanism)"
  * reworded; the copy reads: "inside the model's D(r) formulation"
  * words for the organisation of the work replaced as in the other notes of the release ("re-check" for the second pass, "stretch of work" for an interrupted or resumed working period, the names of the analyses for their labels)

## `code/bsc_validation2/a2_contact_mobility/PREREGISTRATION.md`

* SHA-256 of the original: `8df17115f573dd01e8006f86576e2d31da659888539917fe9acd6e2026a2f66e`
* SHA-256 of this copy: `6e370aced612c359eef308eba55f31a62d30d4e79ee38820ba155a0df030a5d5`
* lines changed in this copy: 9, 20, 66, 68, 70, 72, 91, 129, 138, 165
* changes:
  * reworded; the copy reads: "The model's defining mechanism [...]:"
  * reworded; the copy reads: "* **RW0 — lattice walk, homogeneous.**"
  * reworded; the copy reads: "(discrete-time version of the continuous-time walk)"
  * reworded; the copy reads: "* **RWhet — lattice walk, heterogeneous two-zone.**"
  * reworded; the copy reads: "(harmonic-mean rule for W(r→r'))"
  * reworded; the copy reads: "the test of the lattice-walk mechanism as formulated"
  * reworded; the copy reads: "A0 for the lattice walk as formulated."
  * reworded; the copy reads: "keeps the framework (a Markov walk"
  * reworded; the copy reads: "two purely descriptive scripts were run in an earlier, interrupted ..."
  * reworded; the copy reads: "("matched mean contact rate", as planned for this analysis)."

## `code/bsc_validation2/a2_contact_mobility/src/a2lib.py`

* SHA-256 of the original: `7656161a924e8b8afe2dbd1990abcf4b2e8b89c190b4fdfd93476292dc358f51`
* SHA-256 of this copy: `a9b8dd419695d8b36c9c1143ebd431d0e0a59da9829405bd7ad7616911a7f079`
* lines changed in this copy: 249
* changes:
  * reworded; the copy reads: "(harmonic mean, the interface rule of the model)"

## `code/bsc_validation2/a3_tracer_physics/POSTHOC_LOG.md`

* SHA-256 of the original: `0e7e4bb5531497ee16f0c636771e8b0d1c68ba6f6ad9698564a61e65b8079a88`
* SHA-256 of this copy: `f0a4d50a37271d9356ac0962e1ce77c8baf627826fbe26d75e6e791e272646b2`
* lines changed in this copy: 48, 51, 56, 58
* changes:
  * reworded; the copy reads: "(that stretch of work stopped after writing `src/s07_memo.py`)"
  * words for the organisation of the work replaced as in the other notes of the release ("re-check" for the second pass, "stretch of work" for an interrupted or resumed working period, the names of the analyses for their labels)

## `code/bsc_validation2/a3_tracer_physics/PREREGISTRATION.md`

* SHA-256 of the original: `e721b8c49d4c6dd88480229c1326eaa148f59b825097c83004fcbb93cfad6499`
* SHA-256 of this copy: `5fd0e419171756daafb2ba197310647ae78fad5dc370cb1e4c2a003c3c9afeb7`
* lines changed in this copy: 3, 26, 37-39
* changes:
  * passage omitted, marked [...]
  * reworded; the copy reads: "the non-local airborne term of model M2,"
  * reworded; the copy reads: "(lattice-walk propagator [...] for the quanta"
  * reworded; the copy reads: "This tests the diffusion form of the kernel."

## `code/bsc_validation2/a4_zone_level/POSTHOC_LOG.md`

* SHA-256 of the original: `86dca3798efeabd0b88e0c26bd3aa9faa6c9e8d34ea91f3d126b47ace24d7b68`
* SHA-256 of this copy: `dc08b276820c878c31d6e1f750a4b9ddeeca5183df52c33044d6e43365604380`
* lines changed in this copy: 1, 14-17
* changes:
  * reworded; the copy reads: "The same computation had been made separately before, for the overv..."
  * words for the organisation of the work replaced as in the other notes of the release ("re-check" for the second pass, "stretch of work" for an interrupted or resumed working period, the names of the analyses for their labels)

## `code/bsc_validation2/a4_zone_level/PREREGISTRATION.md`

* SHA-256 of the original: `a434c9500d02ad68ccbb5f4cddd1adf3bae1514def731cd063ce5aff4332bdde`
* SHA-256 of this copy: `08359b1cbfd571cd76c728107235cfb12bd7a36e8c34045471b5e84d067f8c16`
* lines changed in this copy: 26, 29, 35, 79
* changes:
  * reworded; the copy reads: "with same-site, frequency-dependent transmission [...] the force"
  * reworded; the copy reads: "This is the same-site transmission kernel ("M1")"
  * reworded; the copy reads: "| **Z** — zone model,"
  * reworded; the copy reads: "The theory of the model [...] says that"

## `code/bsc_validation2/a4_zone_level/src/01_lyon_mixing.py`

* SHA-256 of the original: `cddb8f476d19bd86efc3c6792f446d83e2f2d13172a981a0be43ecd6c7566fcd`
* SHA-256 of this copy: `c794f65e6025ad45d71b066a5c14da591e67b8578120362884dd832bda68d1e5`
* lines changed in this copy: 30
* changes:
  * reworded; the copy reads: "(duration is the quantity of the model)"
