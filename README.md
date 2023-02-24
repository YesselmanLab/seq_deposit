# seq_deposit

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

a set of commands for ensuring that all sequences ordered are properly organized and documented for the Yesselman Lab.

## Install

To install seq_deposit

```shell
python -m pip install git+https://github.com/jyesselm/seq_deposit
```

## How to use
The python package generates a new executable called `seq-deposit` ensure that is exists in your path. 
```shell
$ seq-deposit --help
Usage: seq-deposit [OPTIONS] COMMAND [ARGS]...

  command line interface for seq_deposit script

Options:
  --help  Show this message and exit.

Commands:
  assembly          generate required files for primer assemblies
  opools            generate required files for opools
  set-deposit-path  update deposit path
```

### set-deposit-path
The first step is to set the deposit path. This is the path to the directory where all the sequencing data is stored. Make sure you have 
```shell

## TODO
