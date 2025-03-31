from setuptools import setup, find_packages

# setup(
#     name="ssdna_tool",
#     version="0.1",
#     description="Package to provide modules to facilitate ssDNA virus annotation.",
#     author="Pavan R. & Tisza M.",
#     author_email="pavan.4@osu.edu",
#     url="https://github.com/ricrocha82/ssdna_tool",
#     packages=["ssdna_tool"],
#     zip_safe=False,
# )

setup(
    name="cressent",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "cressent=cressent_core.cli:cli",
        ],
    },
    install_requires=[
        "biopython",
        "pandas",
        "matplotlib",
        "seaborn",
        "numpy",
        "gffutils",
        "click",
        "viennarna",
    ],
    python_requires=">=3.6",
    description="A comprehensive toolkit for ssDNA virus analysis",
    author="Pavan R. & Tisza M.",
    author_email="pavan.4@osu.edu",
    url="https://github.com/ricrocha82/cressent",
)