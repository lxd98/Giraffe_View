import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="Giraffe_View",
  version="0.1.0.1",
  author="Xudong Liu",
  author_email="xudongliu98@gmail.com",
  description="Giraffe_View is specially designed to provide a comprehensive assessment of the accuracy of long-read sequencing datasets obtained from both the PacBio and Nanopore platforms.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/lxd98/Giraffe_View",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
  python_requires = '>=3.8, <=3.11',
  install_requires=[
  'pysam >= 0.17.0',
  'plotnine >= 0.12.1',
  'numpy >= 1.7.0',
  'pandas >= 1.5.0',
  'tqdm == 4.64.0',
  'termcolor >= 2.0.0'
  ],
  scripts = ["Giraffe_View/giraffe"]
)
