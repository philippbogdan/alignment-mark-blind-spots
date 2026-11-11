# The 51 pairs

One plain-text file per pair: two marks, their placements, the channels that cannot tell them apart and the noise that claim assumes. Definitions are in [MODEL.md](../MODEL.md); both solvers' data for each pair are in [records/](../records/) under the same name.

A **mirror** pair is a mark and its mirror image, which have the same visibility and strength in every channel. For identical data the mark must also shift every blind channel by the same amount: those equal-shift equations were solved by Newton's method on a rigorous solver, and the mirror image is then placed twice the shift away. For a mirror pair within noise the shifts need only agree to within the position noise: a nearly blind mark, scaled towards its symmetric form until the spread fits. Pairs 06, 07, 08, 39 are **general** pairs, two different marks: pair 08 solves for equal readings in all eight channels, pairs 06 and 07 are two nearby marks along the direction the data constrain least, and pair 39 minimises the difference in the data at a fixed distance.

- **cap**: no deformation parameter of either mark exceeds it (nm).
- **blind channels**: filled where the pair claims the channel, in the order 532 TE, 532 TM, 633 TE, 633 TM, 775 TE, 775 TM, 850 TE, 850 TM.
- **data**: *identical*: every position, visibility and strength agrees. *positions identical*: the eight positions agree and the intensities are not compared. *within 0.01 nm*: `d <= 1` with that noise on every position and equal visibility and strength; a second number is the noise on visibility and strength.
- **apart**: distance between the two placements (nm). A lower bound on what those channels can hide: the pairs were built to pass 0.1 nm, not to be as far apart as possible.
- **centroids**: the same distance measured between the centroids of the two lines' cross-sections (nm).
- **imbalance**: the largest `|e+ - e-| / (e+ + e-)` over the blind channels. The two marks of a mirror pair have opposite imbalances.
- **other**: the largest difference between the two marks' readings in those of the eight channels the pair does not claim (nm).
- **median**, **largest**: the median and the largest difference between the two marks' readings at other wavelengths, every 50 nm from 500 to 1000 nm in both polarisations (nm). The largest is usually where that line's first orders are weak.

| pair | cap | blind channels | data | apart | centroids | imbalance | other | median | largest |
|---|---|---|---|---|---|---|---|---|---|
| [01](01-r80-all8-identical.txt) | 80 | ●●●●●●●● | identical | 0.3996 | 0.1676 | 1.3e-03 |  | 0.017 | 0.13 |
| [02](02-r10-all8-noise.txt) | 10 | ●●●●●●●● | within 0.01 nm | 0.0763 | 0.0306 | 1.4e-04 |  | 0.004 | 0.01 |
| [03](03-r20-all8-noise.txt) | 20 | ●●●●●●●● | within 0.01 nm | 0.0923 | 0.0339 | 1.9e-04 |  | 0.005 | 0.02 |
| [04](04-r40-all8-noise.txt) | 40 | ●●●●●●●● | within 0.01 nm | 0.1201 | 0.0343 | 3.6e-04 |  | 0.006 | 0.03 |
| [05](05-r80-all8-noise.txt) | 80 | ●●●●●●●● | within 0.01 nm | 0.5301 | 0.2388 | 1.2e-03 |  | 0.024 | 0.11 |
| [06](06-r10-all8-noise-1e-3.txt) | 10 | ●●●●●●●● | within 0.01 nm, 0.001 | 0.1047 | 0.0445 | 1.2e-02 |  | 0.005 | 0.02 |
| [07](07-r10-all8-noise-1e-2.txt) | 10 | ●●●●●●●● | within 0.01 nm, 0.01 | 0.2012 | 0.1167 | 2.4e-02 |  | 0.011 | 0.05 |
| [08](08-r10-all8-positions.txt) | 10 | ●●●●●●●● | positions identical | 0.4913 | 0.7501 | 4.0e-03 |  | 0.046 | 0.71 |
| [09](09-r20-0134567-identical.txt) | 20 | ●●○●●●●● | identical | 0.3999 | 1.5822 | 1.1e-02 | 0.122 | 0.104 | 0.54 |
| [10](10-r40-0123457-identical.txt) | 40 | ●●●●●●○● | identical | 0.4046 | 0.1207 | 8.6e-03 | 0.126 | 0.094 | 1.10 |
| [11](11-r40-0123567-identical.txt) | 40 | ●●●●○●●● | identical | 0.3999 | 1.5738 | 9.2e-03 | 0.145 | 0.089 | 0.33 |
| [12](12-r10-0123456-noise.txt) | 10 | ●●●●●●●○ | within 0.01 nm | 0.1409 | 0.0428 | 2.7e-04 | 0.033 | 0.008 | 0.04 |
| [13](13-r10-0124567-noise.txt) | 10 | ●●●○●●●● | within 0.01 nm | 0.1221 | 0.0704 | 1.1e-04 | 0.016 | 0.008 | 0.02 |
| [14](14-r10-0234567-noise.txt) | 10 | ●○●●●●●● | within 0.01 nm | 0.1129 | 0.0631 | 8.0e-05 | 0.083 | 0.004 | 0.41 |
| [15](15-r10-1234567-noise.txt) | 10 | ○●●●●●●● | within 0.01 nm | 0.1075 | 0.2912 | 1.6e-03 | 0.052 | 0.016 | 0.13 |
| [16](16-r20-0123457-noise.txt) | 20 | ●●●●●●○● | within 0.01 nm | 0.1087 | 0.0229 | 9.2e-04 | 0.024 | 0.014 | 0.04 |
| [17](17-r20-0123567-noise.txt) | 20 | ●●●●○●●● | within 0.01 nm | 0.1007 | 0.0276 | 2.5e-04 | 0.006 | 0.006 | 0.02 |
| [18](18-r10-123456-identical.txt) | 10 | ○●●●●●●○ | identical | 0.4057 | 0.1191 | 2.1e-03 | 0.176 | 0.059 | 0.31 |
| [19](19-r10-012467-identical.txt) | 10 | ●●●○●○●● | identical | 0.4054 | 0.1140 | 2.3e-03 | 0.079 | 0.025 | 6.60 |
| [20](20-r10-014567-identical.txt) | 10 | ●●○○●●●● | identical | 0.3706 | 0.8568 | 6.6e-03 | 0.160 | 0.076 | 0.81 |
| [21](21-r10-013467-identical.txt) | 10 | ●●○●●○●● | identical | 0.3426 | 0.6797 | 9.0e-03 | 0.151 | 0.091 | 0.85 |
| [22](22-r10-234567-identical.txt) | 10 | ○○●●●●●● | identical | 0.3140 | 0.5404 | 1.8e-03 | 0.309 | 0.043 | 0.69 |
| [23](23-r10-023456-identical.txt) | 10 | ●○●●●●●○ | identical | 0.1200 | 0.0511 | 2.0e-04 | 0.033 | 0.011 | 0.08 |
| [24](24-r10-123457-identical.txt) | 10 | ○●●●●●○● | identical | 0.1200 | 0.6241 | 3.7e-03 | 0.200 | 0.058 | 0.31 |
| [25](25-r10-013456-identical.txt) | 10 | ●●○●●●●○ | identical | 0.1200 | 0.0181 | 8.8e-04 | 0.027 | 0.011 | 0.05 |
| [26](26-r20-012456-identical.txt) | 20 | ●●●○●●●○ | identical | 0.1200 | 0.0778 | 2.4e-04 | 0.022 | 0.005 | 0.03 |
| [27](27-r20-012346-identical.txt) | 20 | ●●●●●○●○ | identical | 0.1200 | 0.1136 | 5.7e-04 | 0.032 | 0.009 | 0.04 |
| [28](28-r20-012457-identical.txt) | 20 | ●●●○●●○● | identical | 0.1200 | 0.0124 | 1.4e-03 | 0.029 | 0.018 | 0.05 |
| [29](29-r20-012345-identical.txt) | 20 | ●●●●●●○○ | identical | 0.1200 | 0.1482 | 1.7e-03 | 0.051 | 0.020 | 0.06 |
| [30](30-r20-012567-identical.txt) | 20 | ●●●○○●●● | identical | 0.1200 | 0.4440 | 2.8e-03 | 0.041 | 0.022 | 0.10 |
| [31](31-r20-012356-identical.txt) | 20 | ●●●●○●●○ | identical | 0.1200 | 0.1380 | 1.6e-03 | 0.034 | 0.020 | 0.05 |
| [32](32-r20-012357-identical.txt) | 20 | ●●●●○●○● | identical | 0.1200 | 0.0271 | 5.0e-03 | 0.182 | 0.156 | 0.27 |
| [33](33-r20-124567-identical.txt) | 20 | ○●●○●●●● | identical | 0.1197 | 0.0795 | 1.1e-04 | 0.022 | 0.006 | 0.03 |
| [34](34-r20-012347-identical.txt) | 20 | ●●●●●○○● | identical | 0.1196 | 0.0912 | 2.6e-03 | 0.056 | 0.035 | 0.10 |
| [35](35-r40-023467-identical.txt) | 40 | ●○●●●○●● | identical | 0.3516 | 0.0852 | 4.5e-03 | 0.257 | 0.024 | 1.76 |
| [36](36-r40-123467-identical.txt) | 40 | ○●●●●○●● | identical | 0.1200 | 0.0850 | 7.6e-05 | 0.014 | 0.006 | 0.03 |
| [37](37-r40-024567-identical.txt) | 40 | ●○●○●●●● | identical | 0.1200 | 0.2673 | 3.7e-03 | 0.144 | 0.035 | 0.13 |
| [38](38-r10-012347-noise.txt) | 10 | ●●●●●○○● | within 0.01 nm | 0.1022 | 0.0214 | 1.4e-03 | 0.015 | 0.013 | 0.40 |
| [39](39-r10-012367-noise-1e-4.txt) | 10 | ●●●●○○●● | within 0.01 nm, 0.0001 | 0.1001 | 0.0032 | 9.1e-04 | 0.017 | 0.009 | 0.21 |
| [40](40-r10-01256-identical.txt) | 10 | ●●●○○●●○ | identical | 0.4004 | 0.1700 | 6.6e-04 | 0.048 | 0.012 | 0.08 |
| [41](41-r10-03567-identical.txt) | 10 | ●○○●○●●● | identical | 0.4000 | 0.5236 | 4.9e-03 | 1.183 | 0.051 | 3.62 |
| [42](42-r10-02367-identical.txt) | 10 | ●○●●○○●● | identical | 0.4000 | 0.1978 | 8.9e-04 | 0.550 | 0.018 | 1.56 |
| [43](43-r10-01245-identical.txt) | 10 | ●●●○●●○○ | identical | 0.4000 | 0.2627 | 4.8e-04 | 0.079 | 0.018 | 0.11 |
| [44](44-r10-01257-identical.txt) | 10 | ●●●○○●○● | identical | 0.3278 | 0.4558 | 2.3e-03 | 0.099 | 0.076 | 0.21 |
| [45](45-r10-01234-identical.txt) | 10 | ●●●●●○○○ | identical | 0.1200 | 0.0538 | 1.1e-04 | 0.030 | 0.008 | 0.04 |
| [46](46-r10-01237-identical.txt) | 10 | ●●●●○○○● | identical | 0.1200 | 0.2315 | 3.0e-03 | 0.076 | 0.057 | 0.20 |
| [47](47-r10-01236-identical.txt) | 10 | ●●●●○○●○ | identical | 0.1200 | 0.0333 | 2.3e-04 | 0.029 | 0.008 | 0.04 |
| [48](48-r10-01357-identical.txt) | 10 | ●●○●○●○● | identical | 0.1200 | 1.0047 | 6.6e-03 | 0.212 | 0.080 | 0.51 |
| [49](49-r10-03457-identical.txt) | 10 | ●○○●●●○● | identical | 0.1200 | 0.1390 | 1.1e-03 | 0.232 | 0.019 | 1.23 |
| [50](50-r10-13567-identical.txt) | 10 | ○●○●○●●● | identical | 0.1200 | 0.5897 | 3.4e-03 | 0.088 | 0.030 | 0.24 |
| [51](51-r20-12367-identical.txt) | 20 | ○●●●○○●● | identical | 0.1200 | 0.0392 | 9.0e-04 | 0.023 | 0.017 | 0.04 |
