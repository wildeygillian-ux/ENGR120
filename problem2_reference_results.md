# Problem 2 reference output (from the included topology interpretation)

These values were generated with the same equations implemented in `solve_problem2_network.m`.

## (a) Turbulence check
- All 22 branches were turbulent (`Re > 4000`), minimum Reynolds number was approximately `5859`.

## (b) Flow rates after convergence (m^3/s)
`Q = [0.068338, 0.028463, 0.019666, 0.011077, 0.031662, 0.031662, 0.039875, 0.031744, 0.005888, -0.000205, 0.008797, 0.008131, 0.017134, 0.008589, 0.005723, 0.005682, 0.001538, 0.007261, 0.042739, 0.025856, 0.025856, 0.004144]`

## (c) Branch head losses after convergence (m)
`Δh = [27.401780, 14.999772, 14.484532, 4.713891, 49.586157, 12.013819, 46.998238, 19.936290, 22.600464, -0.043203, 74.578423, 42.579957, 11.051159, 71.145050, 21.378874, 31.624376, 0.848859, 2.080696, 89.890728, 26.567419, 28.405325, 0.747904]`

## (d) Highest head-loss source-to-outlet path
- Highest-loss directed outlet path was to outlet node 15.
- Branch sequence: `[1, 2, 3, 14, 15, 18]`
- Total head loss along that path: `151.490704 m`

> Note: These numerical values depend on branch connectivity. The connectivity used here is encoded in `problem2_run.m` (`data.from` and `data.to`) based on Figure 2 interpretation.
