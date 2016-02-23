# DarkMatterHalos
Calculating half mass radius of dark matter halo simulations.

## Existing Approach
Dark Matter Halos are currently characterized by the [Navarro-Frenk-White Profile](https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile). This profile depends on the change in deinsity gradient which is difficult to calculate precisely for simulations.

## Revised Approach
This project attempts to recharacterize halos using their half-mass-radius. That is, the ellipsoidal radius containing half the particles in the simlation. The first problem is to fit the halo particles to an ellipsoid so the correct radius can be calculated.

An example of a Halo:
![image](https://raw.githubusercontent.com/hazrmard/DarkMatterHalos/master/example_halo_fit.png)
  
## Instructions
### Installation
This package requires `numpy` and python 2.7.x. Clone this repository. In console browse to the repository and run:  
```bash
$ python setup.py install
```
Or you can just import `halos` if you are in the same directory.  
### Usage  
There are several examples present in `testing.py`. One primary interface is the `HalfMassRadius` class. It performs calculations on BGC2 files that contain data for multiple halos.
```python
import halos
H = HalfMassRadius(PATH_TO_BGC2_FILE)     # BGC2 files contain data for multiple halos
H.read_data()                             # Instantiates a Halo object for each halo in file
H.filter(minimum=10)                      # Leave out halos with < 10 particles
H.center_halos()                          # (1) Translate all halos around center points
H.get_covariance_matrices()               # (2)
H.get_eigenvectors()                      # (3)
H.convert_bases()                         # (4)
H.get_radii()                             # (5)
H.get_half_mass_radii()                   # (6) 1-6 are intended to be run in order
```
  
The fundamental object computed on is the `Halo`. The above class `HalfMassRadius` in its function calls simply iterates over `Halo` objects in `HalfMassRadius().halos`. A `Halo` can either be instantiated from an ascii data file or by passing coordinates:
```python
# Generate halo from coordinates
x = [1,2,3,4,5]
y = [5,4,6,3,1]
z = [3,4,6,7,8]
id = 'sample_halo'
center = (0,0,0)
h = halos.helpers.create_halo(id, center, x, y, z)

# OR generate from ascii data file
coords = halos.helpers.read_ascii_pos(FILEPATH)
h = halos.Halo('halo_id', pos=(0,0,0), particles=coords)

# Run computations
h.center_halo()
h.get_covariance_matrix()
h.get_eigenvectors()
h.convert_basis()
h.get_radii()
h.get_half_mass_radius()
```  
  
For higher order fitting:
```python
h.higher_order_fit(ORDER)     # 2,3,4...
```
Finally, the halo can be visualized with an ellipsoidal fit:
```python
# default mode='cleave', also mode='eval'; default transform=False
h.visualize(ellipsoids=True, mode='cleave', transform=True)     
``` 

`cleave` uses absolute maximum projections of particles along each principal axis to draw ellipsoids. `eval` uses halo eigenvalues
to calculate ellipsoid dimensions. `transform=True` rotates inner halo to represent measured distribution of inner particles only,
otherwise the orientation of the inner halo is assumed to be the same as the rest of the particles.  
  
To get stats about the halo on the console:
```python
h.report()
```
  
There are several helper functions in `halos.gendata` module for generating random particle distributions. Please look at source code documentation for more details.  
This project is under active development and function behavior may change significantly in the future.

### Attribution:
`ellipsoids.py` and `bgc2.py` taken from [Daniel Sissom's repository](https://github.com/djsissom).
