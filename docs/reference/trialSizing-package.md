# trialSizing: Tools for Experimental Design Sizing

Sizes field experiments from uniformity-trial data, following the
relationship between the coefficient of variation and plot size. Checks
a trial for the spatial structure the sizing methods assume
(semivariogram, Moran's I, kriged field map), summarises the coefficient
of variation over every plot shape the grid admits, and estimates the
optimal plot size by the modified maximum curvature method of Meier and
Lessman (1971), by the linear response plateau (LRP) and quadratic
response plateau (QRP) models, and by the closed form of Paranaiba,
Ferreira and Morais (2009), which can be compared side by side. From the
coefficient of variation at the optimum it derives the number of
replications needed to detect a given difference between treatment
means, as in Cargnelutti Filho and others (2014). Every method returns
standardised diagnostic statistics, optional bootstrap uncertainty for
the breakpoint, and publication-style plots.

## See also

Useful links:

- <https://github.com/willyanjnr/trialSizing>

- Report bugs at <https://github.com/willyanjnr/trialSizing/issues>

## Author

**Maintainer**: Willyan Bandeira <bandeira.wjab@gmail.com>
([ORCID](https://orcid.org/0000-0002-9430-3664))

Authors:

- Willyan Bandeira <bandeira.wjab@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9430-3664))

- Leonardo Pradebon <leonardopradebon@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-7827-6312))

- Ivan Carvalho <carvalho.irc@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-7947-4900))

- Murilo Loro <muriloloro@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0241-4226))
