# DarkMatterHalos
Calculating half mass radius of dark matter halo simulations.

## Existing Approach
Dark Matter Halos are currently characterized by the [Navarro-Frenk-White Profile](https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile). This profile depends on the change in deinsity gradient which is difficult to calculate precisely for simulations.

## Revised Approach
This project attempts to recharacterize halos using their half-mass-radius. That is, the ellipsoidal radius containing half the particles in the simlation. The first problem is to fit the halo particles to an ellipsoid so the correct radius can be calculated.

An example of a Halo:
![image](https://raw.githubusercontent.com/hazrmard/DarkMatterHalos/master/example_halo_fit.png)

### Attribution:
`ellipsoids.py` and `bgc2.py` taken from [Daniel Sissom's repository](https://github.com/djsissom).
