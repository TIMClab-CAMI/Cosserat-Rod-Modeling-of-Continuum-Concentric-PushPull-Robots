
## Warning !
This is a temporary version of the code. It is functional yet, it is not documented nor is it well formatted. A code revision will come shortly (within the next few months).


# Cosserat-Rod-Modeling-of-Continuum-Concentric-PushPull-Robots
This code is associated to the paper Tummers et al., “Continuum concentric push-pull robots: a Cosserat rod model,” International Journal of Robotics Research to model continuum concentric push-pull robots based on Cosserat rod theory.

Various approaches and structures emerged recently to design continuum robots. One of the most promising designs regards a new concept of continuum concentric push-pull robots (CPPRs) that have the characteristic of combining several key advantages  of  tendon  actuated,  multi-backbone,  and  concentric  tube  ones  (direct  curvature  actuation,  small  outer/innerdiameter ratio, free lumen, etc.). Geometrically-exact models of such recently introduced robots are yet to be developed togain leverage of their full potential. This article extends beyond usual definitions of Cosserat rod theory in order to take intoaccount this new type of continuum robots, constituted by sliding rods, in a shape of tubes whose cross-sections are neitheruniform nor symmetrical along their entire length. The introduced model is capable of considering versatile design options, external loads, 3D deformations, an arbitrary number of tubes and profiles of the centroid lines, as well as a new actuationmethod consisting of an input rotation. Numerical simulations and experiments on CPPR prototypes validate our model.

![hellicoid_CCPPR](https://github.com/TIMClab-CAMI/Cosserat-Rod-Modeling-of-Continuum-Concentric-PushPull-Robots/assets/127660512/525ff9aa-7f37-4c7b-b9c4-cf2b5b64eba6)

## Structure of the code
* The entry point for the generic code repository is the ‘main.m’ script, which reproduces the results of the associated article (see above).
* The "tools" folder contains various general order tools regarding Lie algebra, Chebyshev grids, Legendre polynomials, spectral integration, quaternion operations, saving, reading, plotting results, etc.
 
## Prerequisites
* MATLAB

## Licence
This project is licensed under the GPL v3.0 License - see the [LICENSE](https://github.com/matthiastummers/Cosserat-Rod-Modeling-of-Tendon-Actuated-Continuum-Robots/blob/main/LICENSE) file for details

## Contributing
Feel free to submit pull requests and use the issue tracker to start a discussion about any bugs you encounter. Please provide a description of your MATLAB version and operating system for any software related bugs.

## Acknowledgements
This work was supported by grants ANR-11-LABX-0004-01, ANR-19-P3IA-0003, ANR-20-CE33-0001, ANR-10-IAHU-02, and ANR-18-CE19-0012


