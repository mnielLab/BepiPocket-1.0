from setuptools import setup, find_packages

setup(
    name='bepipocket',
    version='0.1',
    packages=find_packages(),  
    include_package_data=True,  # Include non-code files specified in MANIFEST.in
    package_data={
        'bepipocket': [
            'models/hmm_antibody_identification/*.hmm' # include hmm models for antibody identification
            'models/bepipocket/models/discotope3_models/*' # include DiscoTope-3 models for DiscoPocket Method

        ]
    },
    install_requires=[
        # List your package's dependencies here
    ],
    author='Joakim Clifford',
    author_email='cliffordjoakim@gmail.com',
    description='BepiPocket: Predict antibody-antigen structures using Chai-1 guided by antibody-epitope restraints predicted by BepiPred-3.0 or DiscoTope-3.0.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/mnielLab/BepiPocket-1.0.git',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
