# Records

What the two solvers return for every pair, under the pair's file name. Quantities are defined in [MODEL.md](../MODEL.md); arrays are `[mark 1, mark 2]`, each over the eight channels in the order of `channels`.

| field | content |
|---|---|
| `claimed` | which of the eight channels the pair claims |
| `bie` | boundary-integral solver at its finer setting: `shift` (reading minus placement, nm), `visibility`, `strength`, `imbalance`, the two first-order efficiencies, and `change_with_resolution`, the largest change from the default setting |
| `fem10` | finite elements of order 10, for the claimed channels |
| `fem8` | finite elements of order 8, for all eight channels |
| `agreement_bie_fem10`, `agreement_bie_fem8` | the largest difference between the two solvers |
| `d2` | for a pair within noise: `d^2` on the claimed channels at the pair's noise (indistinguishable when at most 1) |
| `largest_difference_where_equal` | for quantities the pair claims equal: the largest difference between the two marks on the claimed channels |
| `largest_reading_difference_outside_claim` | the largest difference between the two marks' readings in the other channels (nm) |
| `separation`, `separation_between_centroids`, `centroid_offset` | distance between the two placements; the same between the centroids of the two lines' cross-sections; where each centroid sits from the midpoint of its line's feet (nm) |

- `band/`: both marks across 400 to 1100 nm in 50 nm steps (`python check.py pairs/... --band 400 1100 50`). For a mirror pair only mark 1 is listed: its mirror image has the opposite shift and the two orders swapped.
- `spectrum-01.csv`: mark 1 of pair 01 in 1 nm steps, by boundary integrals. `spectrum-01-fem.csv`: the same every 5 nm from 500 to 1000 nm, by finite elements of order 8 (`python check.py pairs/01-r80-all8-identical.txt --band 500 1000 5 --fem`).
- `asymmetries.json`: one asymmetry at a time on the line of pair 01 (`python asymmetries.py`).
