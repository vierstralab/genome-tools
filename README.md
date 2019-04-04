# genome-tools

See [this notebook](examples/UsageExamples.ipynb) for example usage scenarios. If people find this useful I will add more content/examples.

Note: Should work with both python versions 2 and 3. 

## Installation
1. Clone git respository
  ```
  [jvierstra@test0 ~]$ module load git
  [jvierstra@test0 ~]$ mkdir -p ~/.local/src && cd ~/.local/src
  [jvierstra@test0 src]$ git clone https://github.com/jvierstra/genome-tools.git
  ```
2. Compile and locally install python package
  ```
  [jvierstra@test0 src]$ cd genome-tools
  [jvierstra@test0 genome-tools]$ python setup.py install --user
  ```
