# epidemik

Compartmental Epidemic Models in Python


---

## Table of contents[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#toc)
- [Installation](#installation)
- [Tech Stack](#tech)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

---

## Installation[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#installation)

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install epidemik.

```bash
pip install epidemik
```

<div align="right">[ <a href="#table-of-contents">↑ Back to top ↑</a> ]</div>

---

## Tech Stack[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#tech)


Here's a brief high-level overview of the tech stack the `epidemik` package uses:

- The model is implemented as a directed multigraph using [networkx](https://networkx.org/)
- Ordinary Differential Equations are numerically integrated using [scipy](https://scipy.org/)
- 



<div align="right">[ <a href="#table-of-contents">↑ Back to top ↑</a> ]</div>

---


## Usage[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#usage)

`epidemik` provides three main modules, `EpiModel`, `NetworkEpiModel` and `MetaEpiModel`, usually imported directly from the `epidemik` package using

```python
from epidemik import EpiModel
```

To instanciate a new compartmental model we just need to create a `EpiModel` object and add the relevant transitions:

```python
beta = 0.2
mu = 0.1

SIR = EpiModel()
SIR.add_interaction('S', 'I', 'I', beta)
SIR.add_spontaneous('I', 'R', mu)
```

This fully defines the model. We can get a textual representation of the model using
```python
print(SIR)
```

    Epidemic Model with 3 compartments and 2 transitions:

	S + I = I 0.200000
	I -> R 0.100000

	R0=2.00

or a graphical representation by calling `draw_model`:

```python
SIR.draw_model()
```

<img src="https://raw.githubusercontent.com/DataForScience/epidemik/main/images/SIR.png" />

```python
N = 10_0000
I0 = 10

SIR.integrate(365, S=N-I0, I=I0, R=0)
```

```python
SIR.plot()
```

<img src="https://raw.githubusercontent.com/DataForScience/epidemik/main/images/SIR_results.png" />

<div align="right">[ <a href="#table-of-contents">↑ Back to top ↑</a> ]</div>

---

## Contributing[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#contributing)

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

Join our project and provide assistance by:
* Checking out the list of [open issues](https://github.com/aregtech/areg-sdk/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22) where we need help.
* If you need new features, please open a [new issue](https://github.com/DataForScience/epidemik/issues) or start a [discussion](https://github.com/DataForScience/epidemik/discussions).

Please [![epidemik](https://img.shields.io/github/stars/DataForScience/epidemik.svg?style=social&label=Star%20epidemik)](https://github.com/DataForScience/epidemik/) repository to encourage contributors and helps to grow our community. Contact us for the feedback or new ideas.

<div align="right">[ <a href="#table-of-contents">↑ Back to top ↑</a> ]</div>

---

<div align="right">[ <a href="#table-of-contents">↑ Back to top ↑</a> ]</div>

---

## License[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#license)

`epidemik` is free and open-source software licensed under the [MIT License](https://choosealicense.com/licenses/mit/) [2024]  - Bruno Gonçalves. Please have a look at the [LICENSE.md](LICENSE) for more details.

## Spread The Word[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#spread)

If you want to say thank you and/or support active development of the `epidemik` package:

- Add a GitHub Star to the project!
- Tweet about the project on your Twitter!
	- Tag [@data4sci](https://twitter.com/data4sci) and/or [@bgoncalves](https://twitter.com/bgoncalves)

Thank you so much for your interest in growing the reach of the Well app!

