from distutils.core import setup
setup(name='bioviper',
      version='0.2.0',
      author='Sam Berry',
      author_email="sberry@g.harvard.edu",
      description = "Enhancements to Biopython for working with biological data",
      py_modules=["bioviper", "bioviper.msa", "bioviper.pdb", "bioviper.hmmer_tools",
                  "bioviper.phylo", "bioviper.utils"],
      install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'ete3',
        'matplotlib',
        'tqdm']
      )