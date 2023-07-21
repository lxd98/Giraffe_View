import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="Giraffe_View",
  version="0.0.9.4",
  author="Xudong Liu",
  author_email="xudongliu98@gmail.com",
  description="A small tool help assess and visualize the accuracy of a sequencing dataset, \
  	specifically for Oxford Nanopore Technologies (ONT) long-read sequencing.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/lxd98/Giraffe_View",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
  python_requires = '>=3.7, <=3.10',
  install_requires=[
  'pysam >= 0.17.0',
  'plotnine >= 0.12.1',
  'numpy >= 1.7.0',
  'pandas >= 1.5.0',
  'tqdm == 4.64.0'
  ],
  scripts = ["Giraffe_View/giraffe"]
)
