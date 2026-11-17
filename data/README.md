# Classical data and provenance

These files are preserved observations, not synthetic replacements for the failed
Work-context archive. `source_manifest.json` records all 27 downloaded files,
their exact URLs, byte counts, SHA-256 hashes, and (for CRAN mirrors) Git revisions.
Run `python scripts/fetch_classical_data.py` to restore those exact versions.
`--refresh` deliberately replaces the manifest with new upstream versions.

| Data | Preserved source | What is available |
|---|---|---|
| Corbet butterflies | VGAM `corbet` | Frequencies 1–24; documentation supplies n=9,031 and D=620, including 119 omitted species above frequency 24. |
| Rothamsted moths | sads `moths` | 240 abundances, totaling 15,609; Fisher, Corbet and Williams (1943), Table 3, column I. |
| Preston birds | sads `birds` | 80 abundances, totaling 14,353; Preston (1948), Table IA, from Saunders's census. |
| Shakespeare | vcdExtra `ShakeWords` | Frequencies 1–100; n=884,647 and D=31,534 from the documentation; 846 high-frequency types omitted. Efron and Thisted (1976). |
| BCI trees | vegan `BCI` | 50 one-hectare plots, 225 species, 21,457 trees at least 10 cm DBH; 1982 subset compiled by Condit et al. (2002). Not the complete modern BCI census. |
| Dahl trout | Petersen `data_rodli` | Rodli Tarn two-occasion history frequencies: 109 marked, 177 second catches, 57 recaptures; Ricker (1975) reporting Knut Dahl. |
| Gerking sunfish | fishmethods `Gerking` | 14 catches, recaptures, and net marked increments; known deaths after occasions 11 and 13 are retained in both population and marked-pool updates. |
| Snowshoe hares | Rcapture `hare` | 68 histories over six occasions; documentation cites Cormack (1989), Agresti (1994), and Baillargeon and Rivest (2007). |
| Edwards–Eberhardt cottontails | RMark `edwards.eberhardt` | 76 histories over 18 occasions, derived from Edwards and Eberhardt (1967), as distributed with MARK. |
| McCrea–Morgan collection | textbook companion `schnabel_Tab3.1.xlsx` | 15 rows: two vole, chipmunk, hare, skink, rabbit, squirrel, two deer-mouse, house-mouse, pocket-mouse, wood-mouse, two taxi, and golf-tee series. |

Source `.Rd` documentation and package `DESCRIPTION` files are retained under
`documentation/`. Refer to them for full original citations. The
[McCrea–Morgan companion page](https://www.capturerecapture.co.uk/datasets.html)
requests credit to the original collectors and sources. Our citation is to
McCrea and Morgan (2014), *Analysis of Capture–Recapture Data*, Table 3.1,
and to the original investigators where identified: Edwards and Eberhardt
for rabbits and Carothers for taxis. The spreadsheet itself does not identify
every original collector; no invented attribution is supplied.

`scripts/extract_schnabel.py` reads the original OpenXML cells and preserves
blanks in `processed/schnabel_table_raw.json`. It does not fill in absent
frequencies. The wood-mouse row lists 21 occasions but ends at frequency 18;
its recorded total 334 is incomplete and is not used for a richness estimate.
The hare and rabbit nonzero frequency tables agree exactly with their independent
history files. The cottontail histories give Schnabel 86.106, whereas Chao
(1987, p. 787) quotes 96. This unresolved source/formula discrepancy is retained;
Schumacher–Eschmeyer 96.852 agrees with the quoted rounded 97.

The Gerking table's net increments agree with Krebs's online chapter 2,
Table 2.2, p. 46, including two deaths on June 12 and one on June 14.
Its final zero is a package placeholder for the source table's dash, not
evidence of 16 more deaths. It is excluded from inter-occasion updates.

`processed/frequency_variants.json` preserves additional sourced versions.
Taxicab A has 71 singletons in the workbook versus 142 in Ridout's software
notes and Anan, Böhning and Maruotti (2017). Taxicab B has 142 versus 104 in
C. Yu's 2015 Kent thesis, Table 2.1, p. 13. Golf tees have 104 versus 46 in
Anan et al.'s published Tables 4 and 6 (their Table 5 gives total population
250). Other frequencies match for the alternatives we compare. The source
pages were visually checked. No workbook cells are overwritten; the manuscript
reports both sets of point summaries and explains confidence-set unions for
record uncertainty. Yu's table has an extra golf-tee entry, so its golf row is
not used as the alternative source. These are secondary-version sensitivities,
not a claim to have resolved the original field records.

## Third-party terms

The repository's MIT license applies to original code and writing, **not** to
these third-party files. Package declarations are preserved verbatim:

| Package | Declared license |
|---|---|
| sads, vegan, Rcapture, Petersen | GPL-2 |
| VGAM | GPL-3 |
| vcdExtra, RMark, fishmethods | GPL (>= 2) |

GPL texts are in `licenses/`. The exact package-source revisions and URLs are
recorded in the manifest; the upstream CRAN mirrors retain the corresponding
package sources. No new license is asserted for the McCrea–Morgan workbook:
its public companion page provides the download and attribution request, but
does not give an explicit redistribution license. It is preserved here as
research source material with its original provenance, not relicensed as MIT.

## Serial-number availability

`processed/wartime_comparison.json` preserves the three-row aggregate table
in Miller, Sharma and Yang (2024), PUMP 7, Figure 1, p. 232, attributed there
to Ruggles and Brodie. The page was visually checked. Its descriptive error
comparison is generated in the manuscript; it supplies neither raw serials
nor the historical sample sizes or selection rules.

No authenticated raw wartime serial-number sample was recovered in this audit.
The [Richard Ruggles Papers finding aid](https://hdl.handle.net/10079/fa/mssa.ms.1871)
identifies MS 1871, box 1, folder 2 (1947 article), and folders 4–6 (1943–1944
tank and tire reports). An archival locator is not a raw dataset. The tank
examples elsewhere in the repository are explicitly mathematical examples;
they are not relabeled as historical observations.
