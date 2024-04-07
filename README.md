# epidemik

Compartmental Epidemic Models in Python

## Installation[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#project-status)

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install epidemik.

```bash
pip install epidemik
```

## Usage[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#project-status)

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


## Contributing[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#project-status)

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License[![](https://raw.githubusercontent.com/DataForScience/epidemik/main/images/pin.svg)](#project-status)

[MIT](https://choosealicense.com/licenses/mit/)