# bqh-coadjoint

This repository contains the SAGE code used in the paper:

"On the big quantum cohomology of coadjoint varieties" by Nicolas Perrin and Maxim Smirnov.

# Running the code

The code has been tested on SageMath version 9.4.

The easiest way to run the code is by first cloning the whole repository

```
git clone https://github.com/msmirnov18/bqh-coadjoint.git
```
and then going to its subfolder ``code`` and running the scripts from there.
For example, by typing
```
sage e6-p2.sage
```
you will run the computations for the coadjoint variety in type E6.
The first time you run any of the 4 available scripts, the repository
[littlewood-richardson-rule](https://github.com/msmirnov18/littlewood-richardson-rule)
is cloned into the subfolder ``code/littlewood-richardson-rule``.
