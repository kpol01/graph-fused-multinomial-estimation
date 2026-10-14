# Data descriptions

This folder contains two real-data caches. Chicago includes a prepared analysis
object; India includes source workbooks and district boundaries from which the
preparation function constructs counts and a graph. Simulated observations are
created by `simulations/` and saved under `results/generated/`.

The inventories below were checked against the uploaded files on 2026-10-13.
Original real-data files were retained unchanged during this review.

## Chicago crimes, calendar year 2024

### Files and sources

| File in `chicago_crime_data/` | Contents |
| --- | --- |
| `raw_chicago_2024.rds` | Record data plus API, period, expected-count, maximum-ID, and download metadata |
| `chicago_city_boundary.geojson` | City-boundary snapshot used for filtering and mapping |
| `chicago_2024_grid1000m_min1_all.rds` | Full/training/test counts, graph, geometry, and audits |
| `chicago_2024_grid1000m_min1_all_counts.csv` | One row per retained cell; identifiers, totals, centroid coordinates, and category counts |
| `chicago_2024_grid1000m_min1_all_edges.csv` | Rook edges as one-based count-matrix row indices |

Sources: [Chicago Crimes — 2001 to Present](https://data.cityofchicago.org/d/ijzp-q8t2)
and [city boundary](https://data.cityofchicago.org/d/qqq8-j68g).
The cache records **2026-09-09 16:32:37 UTC** as its download time. Its period
is `date >= 2024-01-01T00:00:00` and `date < 2025-01-01T00:00:00`.
Current live-source counts may differ from this extract.

The raw table has 259,607 rows and five fields:

| Field | Meaning |
| --- | --- |
| `id` | Unique source-record identifier used for pagination and duplicate checks |
| `date` | Source date/time string |
| `primary_type` | Published crime classification |
| `latitude`, `longitude` | Geographic coordinates parsed during preparation |

### Cleaning, grid, and categories

Preparation excludes missing/nonfinite/invalid coordinates and missing
classifications, then applies a bounding-box filter and an exact city-polygon
check. Each remaining record is assigned once to a square grid in **NAD83 /
UTM zone 16N, EPSG:26916**, with **1,000-metre** cell sides. Integer cell
coordinates make assignment deterministic. Occupied cells meeting
`min_count = 1` are retained. Edges are horizontal/vertical rook neighbors.
Grid polygons intersect the city; display clipping is separate from estimation.

| Audit | Uploaded cache |
| --- | ---: |
| Downloaded records | 259,607 |
| Invalid-coordinate exclusions | 2,081 |
| Missing-classification exclusions | 0 |
| Outside-city exclusions | 901 |
| Retained records | 256,625 |
| Zero-count city-intersecting cells omitted from the graph | 64 |
| Retained cells | 652 |
| Undirected edges | 1,176 |
| Components | 2 |
| Isolated retained cells | 0 |

Each source record contributes one count. Homicide records can represent
victims, so these should not be described as distinct incidents for every type.

| Count column | Definition | Total |
| --- | --- | ---: |
| `THEFT` | Source THEFT type | 60,038 |
| `BATTERY` | Source BATTERY type | 45,813 |
| `CRIMINAL DAMAGE` | Source CRIMINAL DAMAGE type | 28,370 |
| `ASSAULT` | Source ASSAULT type | 23,303 |
| `BURGLARY` | Source BURGLARY type | 8,373 |
| `OTHER` | All other retained types; final/reference category | 90,728 |

Use `read.csv(..., check.names = FALSE)` to preserve names containing spaces.

| Cell field | Meaning |
| --- | --- |
| `cell_key` | Unique key formed from integer grid coordinates |
| `grid_x`, `grid_y` | Projected lattice coordinates in cell units |
| `node_id` | One-based prepared-object row index |
| `n_i` | Sum of full-data category counts |
| `x_m`, `y_m` | Projected cell-centroid coordinates, metres |

### Holdout and object fields

Training counts are sampled without replacement within cells, using sequential
hypergeometric draws. The requested fraction is 0.8 and seed is 123. Totals
are determined from cell totals alone, rounded and clamped to retain at least
one observation in each set for non-singletons. Singletons stay in training.

| Split audit | Value |
| --- | ---: |
| Training records | 205,051 |
| Test records | 51,574 |
| Achieved training fraction | 0.7990297 |
| Cells with test observations | 642 |
| Training-only singleton cells | 10 |

`Y_train + Y_test == Y` exactly, with a shared node order and graph. Occupancy
and minimum-count selection use full counts before the split. This tests
prediction within retained cells, rather than unseen locations.

The prepared RDS contains count matrices and totals; graph `g`, adjacency `W`,
incidence `D_G`, and edges; `grid_sf`, `city_boundary`, and `coords`; component
membership/totals; `prep_info`, `split_info`, `download_info`, and `cell_audit`;
and category/configuration/session metadata. Load `R/load_all.R` before using
its package-specific objects. Its stored `output_file` describes the creator's
location: use an explicit new path when saving derived objects.

The main component has 650 cells. The two-cell component contains five full-
data records and four training records, with zero training counts for THEFT,
CRIMINAL DAMAGE, ASSAULT, and BURGLARY. A finite unconstrained logit minimizer
need not exist on the full graph. Raw zeros are preserved and the workflow
reports this condition and iteration-limit status.

`chicago_largest_component()` explicitly retains the 650-cell domain, preserves
its existing holdout assignments, and records the excluded cell keys. It is an
analysis-domain choice, not an automatic graph repair. Both methods must use
the same selected object.

## India mother tongues, Census 2011

### Sources and files

The source is C-16, **Population by Mother Tongue**, Census of India 2011.
Sources and pinned revisions remain those in the supplied preparation script:

- [Census publication](https://censusindia.gov.in/2011Census/Language_MTs.html).
- Workbook mirror: `sunil-dhaka/census-language-analysis`, revision
  `d70c2ae4c18f603e165bf3adb4b1c06566848e4a`, under `datasets/C-16/`.
- [DataMeet district boundaries](https://github.com/datameet/maps/tree/master/Districts/Census_2011),
  revision `b3fbbde595310b397a55d718e0958ce249a4fa1f`.
- Boundary credit retained from the source script: DataMeet, CC BY 2.5 India.

| Files in `india_mother_tongue_cache/` | Contents |
| --- | --- |
| `C16_00.xlsx` | National workbook, state totals, and category catalog |
| `C16_01.xlsx` through `C16_35.xlsx` | State/union-territory workbooks using 2011 Census codes |
| `2011_Dist.shp` | Polygon geometry |
| `2011_Dist.shx` | Geometry index |
| `2011_Dist.dbf` | District attributes and codes |
| `2011_Dist.prj` | Source CRS, WGS 84 |

Keep all four shapefile components together. The source has 641 features:
640 matched district codes and one code `000` feature recorded as unmatched.
Rows are 2011 districts, not present-day districts. State codes likewise follow
2011 geography; for example, code 28 predates the Andhra Pradesh/Telangana split.

### Workbook fields and checks

The parser reads the first sheet as text. It selects `C0116` rows whose
subdistrict code is `0` or `00000`, keeping district/state totals rather than
adding subdistrict observations to their parents.

| Parsed field | Meaning |
| --- | --- |
| `state_code` | Two-digit Census state/union-territory code |
| `district_code` | Three-digit all-India district code; `000` denotes an aggregate |
| `area_name` / `district_name` | Published geographic name |
| `mother_tongue_code` | Six-digit category code |
| `mother_tongue` | Published category label |
| `persons` | Area/category population count |
| `is_language_total` | Code ends in `000` |
| `language_code` | Parent group, formed from the first three digits and `000` |

Checks enforce nonnegative integer counts, persons = male + female,
persons = rural + urban, unique area/category rows, and parent/child
reconciliation. `124000` OTHERS is a standalone total without children.
District totals reconcile with national-workbook state totals. Denominators
sum only top-level language totals; summing all hierarchy rows double-counts.
The parser returns 355 non-total national catalog rows. This describes the
parsed publication, not an assertion that every row is an unaggregated taxon.

### Default categories

The defaults select exact published mother tongues, rather than broader
language totals. OTHER collects unselected categories and is last/reference.

| Category | Mother-tongue code | Population |
| --- | --- | ---: |
| Hindi | `006240` | 322,230,097 |
| Bengali | `002007` | 96,177,835 |
| Marathi | `013071` | 82,801,140 |
| Telugu | `021046` | 80,912,459 |
| Tamil | `020027` | 68,888,839 |
| OTHER | Remaining district population | 559,844,607 |

For example, mother tongue `006240` differs from Hindi language total `006000`.
Custom categories must use distinct leaf codes and unique labels; language
parents cannot be mixed with their children as multinomial categories.

### Geography, graph, and split

District matching uses `censuscode`, the all-India district identifier, and
checks `ST_CEN_CD`. `DT_CEN_CD` is not that national identifier. Geometry is
projected into the original India-centred equal-area CRS:

```text
+proj=laea +lat_0=22 +lon_0=80 +datum=WGS84 +units=m +no_defs
```

Edges use **queen contiguity**, including edge or point contacts. No distance
bridges are added. Full preparation preserves zeros and isolated districts.

| Verified quantity | Full dataset |
| --- | ---: |
| Districts | 640 |
| Default categories | 6 |
| Population | 1,210,854,977 |
| Training population, fraction 0.8 and seed 123 | 968,683,726 |
| Test population | 242,171,251 |
| Undirected edges | 1,580 |
| Components | 13 |
| Isolated districts | 9 |

Training/test counts use the supplied sequential hypergeometric routine,
without replacement within districts. Full counts remain available. This is an
individual-within-district holdout, not a test on unseen districts.

Isolated codes: `175`, `212`, `390`, `587`, `635`, `637`, `638`, `639`, `640`.
`india_largest_component()` retains 625 districts, 1,577 edges, and population
1,192,218,754. It preserves their counts and splits, rebuilds graph matrices,
retains the full boundary, and records all 15 excluded codes. No absent default
category was found in the full graph's component training totals; isolated
vertices still require an explicit domain decision for the MCAR interface.

Preparation returns counts and totals; district attributes/geometries; `W`,
`D_G`, `L`, neighbors, degrees, and components; category catalog, source
metadata, and split audits. No prepared India RDS was present in the original
cache. `applications/india_preparation.R` creates one in `results/generated/`.

## Simulation and summary assets

Simulations generate counts at runtime with totals sampled uniformly from
150 through 300. Grid/kNN truth assigns category one probability 0.08 within
a diagonal X band and 0.70 outside; other categories share the remainder.
SBM truth has block centers and Gaussian logit perturbations. Each generated
analysis saves its settings, truth, counts, and fits.

Supplied `../results/` tables are imported summaries, not raw caches:

| Files | Stored fields |
| --- | --- |
| `mse_{grid,knn,sbm}_n.csv` | n; `gfm`, `mcar`, `mmle`, and corresponding `_sig` fields |
| `mse_{grid,knn,sbm}_D.csv` | Same method fields; D axis |
| `mse_{grid,knn,sbm}_n_D.csv` | Same method fields; n axis; no D column |
| `t_{grid,knn,sbm}_n.csv` | n; method timing values in seconds and `_sig` fields |
| `time_grid_D.csv` | Exported row index, D, `re_ours`, `t_ours`, `re_prop`, `t_prop` |
| `time_grid_n_D.csv` | Exported row index, n, D, and the same error/time fields |

Error-table n values are {49, 100, 196, 289, 484} for grid and {50, 100, 200,
300, 500} for kNN/SBM. Category tables use D = {5, 50, 100, 300, 500}; joint
error tables use n = {49, 144, 256, 361, 484}. All timing-n tables store
{50, 100, 200, 300, 500}, retained literally rather than replaced with grid-error
values. `time_grid_D.csv` uses D = {5, 30, 50, 100, 150}.

Replicate records, complete seed/hardware logs, and all metric-normalization
metadata are absent. Filenames alone do not establish every error unit or
whether `_sig` is SD or SE. The reviewed plots display stored values directly,
with an explicit uncertainty multiplier, and do not invent relative-error
denominators or claim to regenerate the source experiments.

Unrelated course materials and their images were moved from the Chicago cache
to `../archive/coursework/` with their bytes preserved. They are not inputs to
these preparation pipelines. Synthetic workbooks under `../tests/fixtures/`
are software-test fixtures, not Census observations.
