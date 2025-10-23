# TIRA: Toolbox for Interval Reachability Analysis

TIRA is a Matlab library gathering several methods to compute an interval over-approximation of the finite-time reachable set for both continous-time and discrete-time systems.

## Documentation

The documentation for TIRA can be found in the [Documentation](https://gitlab.com/pj_meyer/TIRA/tree/master/Documentation) folder. It contains:

* `TIRA_user_manual.pdf`: the user manual of the toolbox detailing the installation instructions, toolbox architecture, an overview of the implemented reachability methods and their requirements, as well as guidelines for three types of use of the toolbox (testing it on pre-defined examples; calling implemented reachability methods on your own system; adding new reachability methods to the toolbox).
* `HSCC2019_paper.pdf`: tool-paper at HSCC19 describing the toolbox architecture and the theory of the 5 main over-approximation methods implemented in version 1 of TIRA (2019-02-19).
* `Springer2021_partial_preprint.pdf`: partial preprint of the book [**Interval Reachability Analysis: Bounding Trajectories of Uncertain Systems with Boxes for Control and Verification**](https://link.springer.com/book/10.1007/978-3-030-65110-7) containing a thourough tutorial-like presentation of all 7 reachability methods implemented in version 2 of TIRA (2021-12-03). The uploaded partial preprint contains the preface, table of content, and Appendix B related to TIRA). If you are interested in this book, feel free to contact the authors directly.

## Changelog

### Version 2 - 2021-12-03

* 2 new reachability methods
	* Quasi-Monte Carlo
	* Monte Carlo
* 1 new submethod for Sampled-data mixed-monotonicity, using second-order sensitivity
* Reworked logger for cleaner console outputs (error, warning, info, runtime)
* Folder `Test_suite` calling all methods and submethods to check that future changes of the toolbox do not break anything
* Tighter method for over-approximation of interval matrix exponential using Horner formulation and Scaling and squaring
* Option to pick or define the ODE solver used in continuous-time methods that need one
* Improved and more thourough user manual in folder [Documentation](https://gitlab.com/pj_meyer/TIRA/tree/master/Documentation)

### Version 1 - 2019-02-19

* Unified framework for interval reachability analysis of continuous-time and discrete-time systems
* 5 main reachability methods:
	* Continuous-time monotonicity
	* Growth-bound and contraction analysis of continuous-time systems
	* Continuous-time mixed-monotonicity
	* Sampled-data mixed-monotonicity, with 2 submethods:
		* Interval arithmetics
		* Sampling and falsification
	* Discrete-time mixed-monotonicity
* Folder `Examples` containing 12 pre-defined example systems

## References

### Citing TIRA

HSCC-2019 tool paper: [**TIRA: Toolbox for Interval Reachability Analysis**](https://arxiv.org/pdf/1902.05204.pdf)
```
@inproceedings{meyer2019hscc,
  title={{TIRA}: Toolbox for Interval Reachability Analysis},
  author={Meyer, Pierre-Jean and Devonport, Alex and Arcak, Murat},
  booktitle={$22^{nd}$ ACM International Conference on Hybrid Systems: Computation and Control},
  pages={224-229},
  year={2019}
}
```

### Reference book

The 2021 book [**Interval Reachability Analysis: Bounding Trajectories of Uncertain Systems with Boxes for Control and Verification**](https://link.springer.com/book/10.1007/978-3-030-65110-7) provides a tutorial-like description of all interval reachability analysis methods implemented in version 2 of TIRA (2021-12-03)
```
@book{meyer2021springer,
  title={Interval Reachability Analysis: Bounding Trajectories of Uncertain Systems with Boxes for Control and Verification},
  author={Meyer, Pierre-Jean and Devonport, Alex and Arcak, Murat},
  year={2021},
  publisher={Springer}
}
```

### Other relevant publications from the authors

IEEE TAC-2019 paper [**Hierarchical decomposition of LTL synthesis problem for nonlinear control systems**](https://arxiv.org/pdf/1712.06014.pdf) for the method on continuous-time mixed-monotonicity
```
@article{meyer2019tac,
  title={Hierarchical decomposition of LTL synthesis problem for nonlinear control systems},
  author={Meyer, Pierre-Jean and Dimarogonas, Dimos V},
  journal={IEEE Transactions on Automatic Control},
  volume={64},
  number={11},
  pages={4676--4683},
  year={2019},
  publisher={IEEE}
}
```

IEEE L-CSS-2018 paper [**Sampled-data reachability analysis using sensitivity and mixed-monotonicity**](https://arxiv.org/pdf/1803.02214.pdf) for the methods on sampled-data mixed-monotonicity (submethods using interval arithmetics, or sampling and falsification) and discrete-time mixed-monotonicity
```
@article{meyer2018lcss,
  title={Sampled-data reachability analysis using sensitivity and mixed-monotonicity},
  author={Meyer, Pierre-Jean and Coogan, Samuel and Arcak, Murat},
  journal={IEEE Control Systems Letters},
  volume={2},
  number={4},
  pages={761--766},
  year={2018}
}
```

IFAC-2020 paper [**Interval Reachability Analysis using Second-Order Sensitivity**](https://arxiv.org/pdf/1911.09775.pdf) for the sampled-data mixed-monotonicity submethod using second-order sensitivity
```
@inproceedings{meyer2020ifacReachability,
  title={Interval Reachability Analysis using Second-Order Sensitivity},
  author={Meyer, Pierre-Jean and Arcak, Murat},
  booktitle={Proceedings of the $21^{st}$ IFAC World Congress},
  pages={1851--1856},
  year={2020},
}
```

L4DC-2020 paper [**Estimating Reachable Sets with Scenario Optimization**](https://par.nsf.gov/servlets/purl/10170424) for the method on quasi-Monte Carlo
```
@inproceedings{devonport2020l4dc,
  title={Estimating reachable sets with scenario optimization},
  author={Devonport, Alex and Arcak, Murat},
  booktitle={Learning for dynamics and control},
  pages={75--84},
  year={2020},
  organization={PMLR}
}
```

## Contributors

* [**Pierre-Jean Meyer**](http://chapal.eu/pierre-jean_meyer/), COSYS-ESTAS, Univ Gustave Eiffel, Lille
	* email @univ-eiffel.fr : pierre-jean.meyer
* **Alex Devonport**, EECS, UC Berkeley
	* email @berkeley.edu : alex_devonport
* **Neelay Junnarkar**, EECS, UC Berkeley
* [**Murat Arcak**](https://people.eecs.berkeley.edu/~arcak/), EECS, UC Berkeley
	* email @berkeley.edu : arcak

## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

