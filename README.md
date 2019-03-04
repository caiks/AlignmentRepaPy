# AlignmentRepaPy

The AlignmentRepaPy repository is a fast Python and C implementation of some of the *practicable inducers* described in the paper *The Theory and Practice of Induction by Alignment* at https://greenlake.co.uk/. The AlignmentRepa repository depends on the [AlignmentPy repository](https://github.com/caiks/AlignmentPy) for the underlying *model* framework. The slower implementations of some of the *practicable inducers* in the Alignment repository can be used to verify the correctness of equivalent faster implementations in AlignmentRepa.

The AlignmentRepaPy repository uses [NumPy](http://www.numpy.org/) high performance arrays. In addition, some compute-intensive array processing is implemented in C using the [NumPy C-API](https://docs.scipy.org/doc/numpy/user/c-info.html). 

## Documentation

The [Python implementation of fast Practicable Inducers](https://greenlake.co.uk/pages/inducer_python_impl_repa) discusses the implementation of the *inducers* using this repository. 

## Installation

The `Alignment` module requires the [Python 3 platform](https://www.python.org/downloads/) to be installed.

For example in Ubuntu,
```sh
sudo apt-get update
sudo apt-get install python3.5
sudo apt install python3-pip
```
Then use the Python installer tool `pip` to install [Sorted Containers](http://www.grantjenks.com/docs/sortedcontainers) and [NumPy and SciPy](https://www.scipy.org/), 
```sh
python3.5 -m pip install --user numpy
python3.5 -m pip install --user scipy
python3.5 -m pip install --user sortedcontainers
```
Then download the zip file or use git to get the repositories -
```sh
cd
git clone https://github.com/caiks/AlignmentPy.git
git clone https://github.com/caiks/AlignmentRepaPy.git
```
Then to build the `AlignmentForeignPy` module that contains the `C` interface,
```sh
cd AlignmentRepaPy
python3 setup.py build
sudo python3 setup.py install
```
In Windows the `AlignmentForeignPy` module must be built explicitly. For example,
```sh
cd /d AlignmentRepaPy
cl.exe -Dinline=__inline /LD /I "." /I "C:\Program Files (x86)\Python\Python37-32\include" /I "C:\Program Files (x86)\Python\Python37-32\Lib\site-packages\numpy\core\include" AlignmentForeignPy.c "C:\Program Files (x86)\Python\Python37-32\libs\python37.lib"
move /Y AlignmentForeignPy.dll  AlignmentForeignPy.pyd
```
Note that if the `AlignmentForeignPy` module fails to build, the `AlignmentRepa` module will still run in pure NumPy Python, albeit more slowly. The success of the `AlignmentForeignPy` module import can be checked at runtime,
```py
AlignmentForeignPy_ok
# True
```

## Usage

Load `AlignmentDevRepa` to import the modules and define various useful abbreviated functions. For example, in Ubuntu,
```sh
export PYTHONPATH=../AlignmentPy:../AlignmentRepaPy
cd ../AlignmentRepaPy
python3
```
```py
from AlignmentRepaPy import *

aarr = systemsHistogramsHistogramRepa
rraa = systemsHistogramRepasHistogram

aa = regdiag(2,2)

aa
# {({(1, 1), (2, 1)}, 1 % 1), ({(1, 2), (2, 2)}, 1 % 1)}

ar = aarr(sys(aa),aa)

ar
# ([1, 2], {(1, 0), (2, 1)}, array([[1., 0.],
#        [0., 1.]]))

rraa(sys(aa),ar)
# {({(1, 1), (2, 1)}, 1 % 1), ({(1, 1), (2, 2)}, 0 % 1), ({(1, 2), (2, 1)}, 0 % 1), ({(1, 2), (2, 2)}, 1 % 1)}
```



