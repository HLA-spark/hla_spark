from setuptools import setup, find_packages
def get_readme():
    with open("README.md","rt",encoding="utf-8") as fh:
        return fh.read()

def get_requirements():
    with open("requirements.txt","rt",encoding="utf-8") as fh:
        return [line.strip() for line in fh.readlines()]

def get_version():
    with open("hla_spark/__init__.py", "rt", encoding="utf-8") as fh:
        for line in fh.readlines():
            if line.startswith('__version__'):
                delim='"' if '"' in line else "'"
                return line.split(delim)[1].strip()
    raise RuntimeError("Unable to find version string in HLA_spark/__init__.py")

def get_license():
    with open("LICENSE","rt",encoding="utf-8") as fh:
        return fh.read()

setup(
    name='hla_spark',
    version=get_version(),
    description='This is a package for calling coding SNPs in human leukocyte antigen (HLA) genes from next-generation sequencing (NGS) data',
    author='zhoulin',
    author_email='270194248@qq.com',
    url='https://github.com/zhoulin8908/HLA-spark',
    packages=find_packages(),
    package_dir={"hla_spark":"hla_spark"},
    package_data={
      'hla_spark': ['data/*','data/hla_alignment_matrices/*','data/hla_allele_references/*','data/hla_allele_references/bed/*','config/config_HLA_spark.ini'],
    },
    include_package_data=True,
    license=get_license(),
    install_requires=['numpy','pandas','scipy'],
    entry_points={
        'console_scripts':[
            'hla_spark=hla_spark.hla_spark:hla_spark'
        ],
    },
    requires=get_requirements()

)