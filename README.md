# AlignmentRepaPy

The AlignmentRepaPy repository is a fast Python and C implementation of some of the *practicable inducers* described in the paper *The Theory and Practice of Induction by Alignment* at https://greenlake.co.uk/. The AlignmentRepa repository depends on the [AlignmentPy repository](https://github.com/caiks/AlignmentPy) for the underlying *model* framework. The slower implementations of some of the *practicable inducers* in the Alignment repository can be used to verify the correctness of equivalent faster implementations in AlignmentRepa.

The AlignmentRepa repository uses high performance arrays. Single-dimensional arrays are defined in the [vector](http://hackage.haskell.org/package/vector) library. See [Numeric Haskell](https://wiki.haskell.org/Numeric_Haskell:_A_Vector_Tutorial). Multi-dimensional shape polymorphic parallel arrays are defined in the [repa](http://hackage.haskell.org/package/repa) library. In addition, some compute-intensive array processing is implemented in C using the [Foreign Function Interface](https://wiki.haskell.org/Foreign_Function_Interface). See also [FFI](http://dev.stephendiehl.com/hask/#ffi) and [Data.Vector.Storable](http://hackage.haskell.org/package/vector-0.12.0.1/docs/Data-Vector-Storable.html).

The *induced models* are made persistent using the JSON format which is implemented in the [aeson](http://hackage.haskell.org/package/aeson) library.

There are a couple of useful libraries that should be installed along with repa and aeson to ensure consistent package versions:

[zlib](http://hackage.haskell.org/package/zlib): Compression and decompression in the gzip and zlib formats

[cassava](http://hackage.haskell.org/package/cassava): A CSV parsing and encoding library

## Documentation

The [Haskell implementation of fast Practicable Inducers](https://greenlake.co.uk/pages/inducer_haskell_impl_repa) discusses the implementation of the *inducers* using this repository. 

## Installation

The `AlignmentRepa` module requires the [Haskell platform](https://www.haskell.org/downloads#platform) to be installed.

For example in Ubuntu,
```
sudo apt-get update
sudo apt-get install haskell-platform
```
Now the libaries not included in Haskell platform must be installed,
```
cabal update
cabal install repa repa-io vector-algorithms zlib cassava aeson aeson-pretty
```
Then download the zip files or use git to get the AlignmentRepa repository and the underlying Alignment repository -
```
cd
git clone https://github.com/caiks/Alignment.git
git clone https://github.com/caiks/AlignmentRepa.git
```

## Usage

Typically we wish to force compilation in ghci in order to have the highest performance. See [Compiling to object code inside GHCi](https://downloads.haskell.org/~ghc/8.4.1/docs/html/users_guide/ghci.html#compiling-to-object-code-inside-ghci).
Load `AlignmentDevRepa` to import the modules and define various useful abbreviated functions,
```sh
cd ../Alignment
rm *.o *.hi

cd ../AlignmentRepa
rm *.o *.hi

gcc -fPIC -c AlignmentForeign.c -o AlignmentForeign.o -O3
ghci -i../Alignment -i../AlignmentRepa ../AlignmentRepa/AlignmentForeign.o
```
```hs
:set -fobject-code
:set +m
:l AlignmentDevRepa

let aa = regdiag 2 2

rp $ aa
"{({(1,1),(2,1)},1 % 1),({(1,2),(2,2)},1 % 1)}"

aa
Histogram (fromList [(State (fromList [(VarInt 1,ValInt 1),(VarInt 2,ValInt 1)]),1 % 1),(State (fromList [(VarInt 1,ValInt 2),(VarInt 2,ValInt 2)]),1 % 1)])

aaar (sys aa) aa
HistogramRepa {histogramRepasVectorVar = [VarInt 1,VarInt 2], histogramRepasMapVarInt = fromList [(VarInt 1,0),(VarInt 2,1)], histogramRepasArray = AUnboxed [2,2] [1.0,0.0,0.0,1.0]}
```
Note that if forcing compilation causes functions to be unresolved, for example,
```hs
rp $ Set.fromList [1,2,3]

<interactive>:9:1: Not in scope: ‘Set.fromList’
```
or 
```hs
rp $ fudEmpty

<interactive>:10:6: Not in scope: ‘fudEmpty’
```
then either (a) import the modules explicitly, for example,
```hs
import qualified Data.Set as Set
import qualified Data.Map as Map
import Alignment

rp $ Set.fromList [1,2,3]
"{1,2,3}"

rp $ fudEmpty
"{}"
```
or (b) interpret module `AlignmentDevRepa` by itself. Exit `ghci` and then delete `AlignmentDevRepa.o`,
```sh
rm AlignmentDevRepa.o

ghci -i../Alignment -i../AlignmentRepa ../AlignmentRepa/AlignmentForeign.o
```
```hs
:set +m
:l AlignmentDevRepa

rp $ Set.fromList [1,2,3]
"{1,2,3}"

rp $ fudEmpty
"{}"
```


